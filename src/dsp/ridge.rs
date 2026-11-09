//! Ridge detection: link spectral-contrast local maxima into moving-tone tracks.
//!
//! A flyover produces a near-monotonic *downward* frequency sweep (approach ->
//! closest approach -> recession). Candidate ridges must therefore be smooth,
//! bounded in slope and curvature, and survive a robust Doppler-model fit;
//! unconstrained linking produces zigzags across the whole band that no
//! constant-velocity source could emit.

use crate::model::{Ridge, RidgeSample, Spectrogram};

/// Max allowed frequency change per second during linking (Hz/s).
///
/// Physical bound: a jet at 250 m/s and 500 m slant produces a steady tone near
/// 1 kHz with `|df/dt| ~ f0*v^2/(c*d) ~ 360 Hz/s`, so anything far beyond this
/// is noise-driven linking rather than a source signature.
const MAX_SLOPE_HZ_S: f32 = 400.0;
/// Max change in per-frame step between consecutive links (Hz) — enforces
/// curvature smoothness so the track cannot zigzag.
const MAX_CURVATURE_HZ: f32 = 20.0;
/// Max upward frequency movement per frame (Hz). A single flyover's observed
/// frequency only decreases, so allowing significant upward steps lets the DP
/// wander through unrelated loud peaks.
const MAX_UP_STEP_HZ: f32 = 2.0;
/// Max frames a ridge may bridge without a detected peak.
const MAX_GAP_FRAMES: usize = 2;
/// Minimum ridge duration (s). Flyover sweeps span many seconds; short tracks
/// are almost always noise flicker.
const MIN_DURATION_S: f64 = 15.0;
/// Minimum net downward frequency change (Hz).
const MIN_NET_DROP_HZ: f32 = 4.0;
/// Minimum net drop as a fraction of the tone frequency (a real flyover sweeps).
const MIN_REL_DROP: f32 = 0.02;
/// Evidence a track must earn per frame to keep extending, as a multiple of the
/// seed contrast. Without a cost the DP joins unrelated loud peaks into one
/// arbitrary mega-path; with a cost this high a track has to sustain a genuine
/// narrowband tone to survive.
const CONTINUATION_COST_FACTOR: f32 = 4.0;
/// Maximum number of ridges kept after de-duplication.
const MAX_PATHS: usize = 400;
/// Extraction passes (each re-solves the DP on the residual peaks).
const MAX_EXTRACTIONS: usize = 24;
/// Evidence retained by a peak already used by an earlier track (soft masking).
const USED_GAIN_FACTOR: f32 = 0.2;
/// Minimum fraction of samples that must move downward.
const MIN_DOWNWARD_FRACTION: f32 = 0.7;
/// Minimum mean spectral contrast for a ridge to be considered a real tone.
const MIN_MEAN_CONTRAST: f32 = 2.0;

#[derive(Clone, Copy)]
struct Peak {
    freq_hz: f32,
    contrast: f32,
}

/// Extract ridges from a contrast spectrogram.
pub fn detect_ridges(spec: &Spectrogram, seed_contrast: f32) -> Vec<Ridge> {
    let bench = std::env::var_os("AVALON_TIMING").is_some();
    let t0 = std::time::Instant::now();
    let n_frames = spec.time_axis_s.len();
    let n_freq = spec.freq_axis_hz.len();
    if n_frames == 0 || n_freq == 0 {
        return Vec::new();
    }
    let seed = seed_contrast.max(MIN_MEAN_CONTRAST);

    // 1. Per-frame local maxima above the seed contrast, then non-maximum
    // suppression at the true resolution. Without NMS, a zero-padded FFT
    // presents several correlated bins across one Hann mainlobe (and its
    // sidelobes) as separate "tonals", producing dozens of parallel ridges.
    let bin_hz = spec
        .freq_axis_hz
        .get(1)
        .map(|f| *f - spec.freq_axis_hz[0])
        .filter(|d| *d > 0.0)
        .unwrap_or(1.0);
    let min_sep_hz = 0.75 * spec.resolution_hz.max(bin_hz);
    let local_window = ((min_sep_hz / bin_hz / 2.0).round() as usize).clamp(1, 64);

    // Every frame is independent here (read-only contrast column, per-frame peak
    // list), so this is the analysis hot loop and it runs on the whole pool.
    use rayon::prelude::*;
    let peaks: Vec<Vec<Peak>> = (0..n_frames)
        .into_par_iter()
        .map(|t| {
            let mut out: Vec<Peak> = Vec::new();
            let mut candidates: Vec<Peak> = Vec::new();
            for k in 0..n_freq {
                let c = spec.contrast[k][t];
                if c < seed {
                    continue;
                }
                let mut is_max = true;
                for d in 1..=local_window {
                    if k >= d && spec.contrast[k - d][t] > c {
                        is_max = false;
                        break;
                    }
                    if k + d < n_freq && spec.contrast[k + d][t] > c {
                        is_max = false;
                        break;
                    }
                }
                if is_max {
                    candidates.push(Peak {
                        freq_hz: spec.freq_axis_hz[k],
                        contrast: c,
                    });
                }
            }
            // Greedy NMS: strongest first, drop anything within `min_sep_hz`.
            candidates.sort_by(|a, b| {
                b.contrast
                    .partial_cmp(&a.contrast)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            for cand in candidates.drain(..) {
                if out
                    .iter()
                    .all(|p| (p.freq_hz - cand.freq_hz).abs() >= min_sep_hz)
                {
                    out.push(cand);
                }
            }
            out.sort_by(|a, b| {
                a.freq_hz
                    .partial_cmp(&b.freq_hz)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            out
        })
        .collect();

    let hop_s = if n_frames > 1 {
        spec.time_axis_s[1] - spec.time_axis_s[0]
    } else {
        0.5
    };
    let max_step_per_frame = MAX_SLOPE_HZ_S * hop_s as f32;
    let min_len = (MIN_DURATION_S / hop_s).ceil().max(2.0) as usize;
    if bench {
        let total: usize = peaks.iter().map(|p| p.len()).sum();
        eprintln!(
            "[ridge] peaks: {total} in {:.2}s",
            t0.elapsed().as_secs_f64()
        );
    }
    let t_dp = std::time::Instant::now();

    // 2. Repeated Viterbi extraction over the residual peak set.
    //
    // Greedy "nearest peak" linking lets a loud noise peak steal the path from a
    // faint but consistent flyover tone. A Viterbi pass maximises accumulated
    // evidence along a whole track instead. After each extraction the winning
    // track's peaks are masked and the DP is re-solved, because a single pass
    // only ever yields prefixes of one trunk — the second curve lives in the
    // residual, not in a variant of the first.
    let mut used: Vec<Vec<bool>> = peaks.iter().map(|p| vec![false; p.len()]).collect();
    let mut ridges: Vec<Ridge> = Vec::new();
    let min_path_gain = 2.0 * seed;

    for _ in 0..MAX_EXTRACTIONS {
        let Some((acc, prev, plen)) = viterbi(
            &peaks,
            n_frames,
            max_step_per_frame,
            seed,
            &used,
        ) else {
            break;
        };

        // Best endpoint that already terminates a long enough path.
        let mut best: Option<(f32, usize, usize)> = None;
        for t in 0..n_frames {
            for i in 0..peaks[t].len() {
                if used[t][i] || plen[t][i] < min_len {
                    continue;
                }
                if best.map_or(true, |(bs, _, _)| acc[t][i] > bs) {
                    best = Some((acc[t][i], t, i));
                }
            }
        }
        let Some((score, end_t, end_i)) = best else { break };
        if score < min_path_gain {
            break;
        }

        // Backtrack. Masking the whole track is only safe when it really is a
        // flyover curve: a junk path that merely crosses a faint real curve
        // would otherwise consume the peaks the real curve needs. Always mask
        // the endpoint so the next pass cannot return the identical path.
        let mut chain: Vec<(usize, usize)> = Vec::new();
        let mut cur = (end_t, end_i);
        loop {
            chain.push(cur);
            match prev[cur.0][cur.1] {
                Some(p) => cur = p,
                None => break,
            }
        }
        chain.reverse();
        for (t, i) in &chain {
            used[*t][*i] = true;
        }
        if chain.len() < min_len {
            continue;
        }

        let samples: Vec<RidgeSample> = chain
            .iter()
            .map(|(t, i)| RidgeSample {
                t_s: spec.time_axis_s[*t],
                freq_hz: peaks[*t][*i].freq_hz,
                weight: peaks[*t][*i].contrast,
            })
            .collect();

        let Some(mut ridge) = build_ridge(samples) else {
            continue;
        };
        // Provisional fit decides whether this track is a real flyover.
        ridge.fit = crate::dsp::doppler::fit_doppler(&ridge);
        ridge.plausible = crate::dsp::doppler::is_plausible(&ridge);
        if bench {
            let first = &chain[0];
            let last = &chain[chain.len() - 1];
            eprintln!(
                "[ridge]   pass len {} score {:.1} {:.1}->{:.1} Hz {}",
                chain.len(),
                score,
                peaks[first.0][first.1].freq_hz,
                peaks[last.0][last.1].freq_hz,
                if ridge.plausible { "PLAUSIBLE" } else { "weak" },
            );
        }
        ridges.push(ridge);
    }

    ridges.sort_by(|a, b| a.start_s.partial_cmp(&b.start_s).unwrap_or(std::cmp::Ordering::Equal));
    if bench {
        eprintln!(
            "[ridge] extract: {:.2}s ({} paths) dp {:.2}s",
            t_dp.elapsed().as_secs_f64(),
            ridges.len(),
            t_dp.elapsed().as_secs_f64()
        );
    }
    // Follow each plausible curve to its natural ends. A Viterbi path stops
    // where the tone fades relative to the per-frame cost, which is typically
    // mid-sweep; the fitted Doppler model predicts where the tone should be, so
    // the track can be extended along that prediction into the faint tail.
    // Fit and extend alternate because each extension constrains (v, d) better,
    // which in turn predicts the flat recession tail more accurately.
    for r in ridges.iter_mut().filter(|r| r.plausible) {
        let original = r.clone();
        for _ in 0..3 {
            let before = r.samples.len();
            r.fit = crate::dsp::doppler::fit_doppler(r);
            extend_along_model(r, &peaks, &spec.time_axis_s);
            if r.samples.len() == before {
                break;
            }
        }
        r.fit = crate::dsp::doppler::fit_doppler(r);
        r.plausible = crate::dsp::doppler::is_plausible(r);
        if std::env::var_os("AVALON_TIMING").is_some() {
            eprintln!(
                "[ridge] extend {:?}->{:?} samples, rmse {:?} -> {:?}, inliers {:?} -> {:?}",
                original.samples.len(),
                r.samples.len(),
                original.fit.map(|f| f.rmse_hz),
                r.fit.map(|f| f.rmse_hz),
                original.fit.map(|f| f.inlier_fraction),
                r.fit.map(|f| f.inlier_fraction),
            );
        }

        // Keep the extended track as long as it still constitutes a plausible
        // flyover; otherwise fall back to the fragment that was verified.
        if !r.plausible {
            *r = original;
        }
    }

    let mut ridges = dedupe(ridges, 3.0 * spec.resolution_hz.max(0.5));
    ridges.truncate(MAX_PATHS);
    for (i, r) in ridges.iter_mut().enumerate() {
        r.id = i;
    }
    ridges
}

/// Extend a ridge past where the Viterbi path stopped.
///
/// The extension follows the *observed* tone from sample to sample (a flyover
/// frequency decreases monotonically, so the next peak is just below/above the
/// last one); the fitted model is only a tie-breaker. Pure model following fails
/// here: a fragment pinning (v, d) predicts the wrong recession asymptote, so
/// the true flatter tail drifts out of tolerance and the track stops early.
fn extend_along_model(ridge: &mut Ridge, peaks: &[Vec<Peak>], times: &[f64]) {
    if peaks.is_empty() || times.len() < 2 {
        return;
    }
    let fit = ridge.fit;
    let hop = (times[1] - times[0]).abs().max(1e-3);
    // Frames are indexed from the *centre* time of frame 0, not from zero: a
    // sample's timestamp already includes the half-window offset, so dividing it
    // by the hop would place every sample several frames late (and rewriting the
    // timestamps from those keys shifts the whole track, cumulatively per pass).
    let frame_of = |t: f64| -> isize {
        (((t - times[0]) / hop).round() as isize).clamp(0, times.len() as isize - 1)
    };
    let mut start_idx = frame_of(ridge.start_s) as usize;
    let mut end_idx = frame_of(ridge.end_s) as usize;

    // Existing samples keyed by frame so extension merges rather than duplicates.
    let mut by_frame: std::collections::BTreeMap<usize, (f32, f32)> = std::collections::BTreeMap::new();
    for s in &ridge.samples {
        let idx = frame_of(s.t_s) as usize;
        by_frame.insert(idx, (s.freq_hz, s.weight));
    }

    let grow = |dir: isize,
                from: usize,
                last_f: f32,
                by_frame: &mut std::collections::BTreeMap<usize, (f32, f32)>|
     -> usize {
        let mut idx = from as isize;
        let mut cur_f = last_f;
        let mut misses = 0;
        let mut last_good = from;
        loop {
            idx += dir;
            if idx < 0 || idx as usize >= times.len() {
                break;
            }
            // Local search window around the last observed frequency. Forward
            // in time the tone descends, so only a little upward jitter is
            // allowed; backwards it ascends.
            let gate = (0.03 * cur_f).max(4.0);
            let (lo, hi) = if dir > 0 {
                (cur_f - gate, cur_f + 2.0)
            } else {
                (cur_f - 2.0, cur_f + gate)
            };
            let pred = fit.map(|f| crate::dsp::doppler::evaluate_fit(&f, times[idx as usize]) as f32);
            let mut best: Option<(f32, f32, f32)> = None; // (cost, freq, contrast)
            for p in &peaks[idx as usize] {
                if p.freq_hz < lo || p.freq_hz > hi {
                    continue;
                }
                // Prefer the model when it agrees; otherwise stay continuous.
                let cost = match pred {
                    Some(pf) => (p.freq_hz - pf).abs().min((p.freq_hz - cur_f).abs()),
                    None => (p.freq_hz - cur_f).abs(),
                };
                if best.map_or(true, |(c, _, _)| cost < c) {
                    best = Some((cost, p.freq_hz, p.contrast));
                }
            }
            match best {
                Some((_, freq, contrast)) => {
                    by_frame.insert(idx as usize, (freq, contrast));
                    cur_f = freq;
                    misses = 0;
                    last_good = idx as usize;
                }
                None => {
                    misses += 1;
                    if misses > 2 {
                        break;
                    }
                }
            }
        }
        last_good
    };

    let first_f = ridge.samples.first().map(|s| s.freq_hz).unwrap_or(0.0);
    let last_f = ridge.samples.last().map(|s| s.freq_hz).unwrap_or(0.0);
    start_idx = grow(-1, start_idx, first_f, &mut by_frame);
    end_idx = grow(1, end_idx, last_f, &mut by_frame);
    let _ = (start_idx, end_idx);
    if by_frame.len() <= ridge.samples.len() {
        return;
    }

    let samples: Vec<RidgeSample> = by_frame
        .into_iter()
        .map(|(idx, (freq, weight))| RidgeSample {
            t_s: times[idx],
            freq_hz: freq,
            weight,
        })
        .collect();
    if let Some(ext) = build_ridge(samples) {
        *ridge = ext;
    }
}

/// One Viterbi pass over the sparse peak set, ignoring masked peaks.
///
/// Returns accumulated evidence, back-pointers and path length per node.
fn viterbi(
    peaks: &[Vec<Peak>],
    n_frames: usize,
    max_step_per_frame: f32,
    seed: f32,
    used: &[Vec<bool>],
) -> Option<(
    Vec<Vec<f32>>,
    Vec<Vec<Option<(usize, usize)>>>,
    Vec<Vec<usize>>,
)> {
    let gain = |p: &Peak| (p.contrast - seed).max(0.0);
    let cost = std::env::var("AVALON_COST")
        .ok()
        .and_then(|v| v.parse::<f32>().ok())
        .unwrap_or(CONTINUATION_COST_FACTOR)
        * seed;
    // Peaks already used by an earlier track keep a fraction of their evidence
    // instead of being removed: later passes prefer new regions, but a real
    // curve that merely crosses a junk path can still be recovered.
    let gain_at = |t: usize, i: usize| {
        let g = gain(&peaks[t][i]);
        if used[t][i] {
            g * USED_GAIN_FACTOR
        } else {
            g
        }
    };
    let mut acc: Vec<Vec<f32>> = peaks
        .iter()
        .enumerate()
        .map(|(t, ps)| (0..ps.len()).map(|i| gain_at(t, i)).collect())
        .collect();
    let mut prev: Vec<Vec<Option<(usize, usize)>>> =
        peaks.iter().map(|ps| vec![None; ps.len()]).collect();
    let mut plen: Vec<Vec<usize>> = peaks.iter().map(|ps| vec![1usize; ps.len()]).collect();

    for t in 1..n_frames {
        for i in 0..peaks[t].len() {
            let q = peaks[t][i];
            let own = gain_at(t, i);
            let mut best_chain = 0.0f32;
            let mut best_prev: Option<(usize, usize)> = None;
            let mut best_len = 1usize;
            for gap in 1..=MAX_GAP_FRAMES {
                if t < gap {
                    break;
                }
                let pt = t - gap;
                if peaks[pt].is_empty() {
                    continue;
                }
                let bound = max_step_per_frame * gap as f32;
                let lo = q.freq_hz - bound;
                let hi = q.freq_hz + bound;
                // Peaks are sorted by frequency: only scan the reachable window.
                let start = peaks[pt].partition_point(|p| p.freq_hz < lo);
                let end = peaks[pt].partition_point(|p| p.freq_hz <= hi);
                for j in start..end {
                    let step = (q.freq_hz - peaks[pt][j].freq_hz) / gap as f32;
                    // Monotonic descent: the observed tone never moves up
                    // appreciably during a single pass.
                    if step > MAX_UP_STEP_HZ {
                        continue;
                    }
                    // Curvature: the per-frame step must stay consistent.
                    if let Some((ppt, pj)) = prev[pt][j] {
                        let dt = (pt - ppt).max(1) as f32;
                        let prev_step = (peaks[pt][j].freq_hz - peaks[ppt][pj].freq_hz) / dt;
                        if (step - prev_step).abs() > MAX_CURVATURE_HZ {
                            continue;
                        }
                    }
                    // Every extra frame must be paid for: a track has to earn more
                    // than CONTINUATION_COST per frame to keep extending.
                    let cand = acc[pt][j] - cost * gap as f32 - 0.01 * step.abs();
                    if cand > best_chain + 1e-6 {
                        best_chain = cand;
                        best_prev = Some((pt, j));
                        best_len = plen[pt][j] + 1;
                    }
                }
            }
            acc[t][i] = own + best_chain;
            prev[t][i] = best_prev;
            plen[t][i] = best_len;
        }
    }
    Some((acc, prev, plen))
}

/// Merge tracks that describe the same physical curve, for callers that detect
/// in separate blocks and need cross-block de-duplication.
pub fn dedupe_tracks(ridges: Vec<Ridge>, tol_hz: f32) -> Vec<Ridge> {
    dedupe(ridges, tol_hz)
}

/// Merge tracks that describe the same physical curve.
fn dedupe(mut ridges: Vec<Ridge>, tol_hz: f32) -> Vec<Ridge> {
    // Keep the longest / strongest representative of each curve.
    ridges.sort_by(|a, b| {
        b.samples
            .len()
            .cmp(&a.samples.len())
            .then_with(|| {
                b.mean_contrast
                    .partial_cmp(&a.mean_contrast)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    let mut kept: Vec<Ridge> = Vec::with_capacity(ridges.len());
    'candidate: for r in ridges {
        for k in &kept {
            let (t0, t1) = (r.start_s.max(k.start_s), r.end_s.min(k.end_s));
            if t1 <= t0 {
                continue;
            }
            // Compare the two tracks where they overlap: paths that describe the
            // same curve agree in frequency there, even if they end at different
            // times (the DP re-solves from a masked peak set).
            if let (Some(fr), Some(fk)) = (freq_at(&r, t0), freq_at(k, t0)) {
                if (fr - fk).abs() < tol_hz {
                    continue 'candidate;
                }
            }
        }
        kept.push(r);
    }
    kept
}

/// Frequency of the track nearest to time `t`.
fn freq_at(r: &Ridge, t: f64) -> Option<f32> {
    r.samples
        .iter()
        .min_by(|a, b| {
            (a.t_s - t)
                .abs()
                .partial_cmp(&(b.t_s - t).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|s| s.freq_hz)
}

fn build_ridge(samples: Vec<RidgeSample>) -> Option<Ridge> {
    if samples.len() < 3 {
        return None;
    }
    let mut min_f = f32::MAX;
    let mut max_f = f32::MIN;
    let mut sum_c = 0.0f32;
    for s in &samples {
        min_f = min_f.min(s.freq_hz);
        max_f = max_f.max(s.freq_hz);
        sum_c += s.weight;
    }
    let net_drop = max_f - min_f;
    if net_drop < MIN_NET_DROP_HZ || net_drop < MIN_REL_DROP * max_f {
        return None;
    }
    let mean_contrast = sum_c / samples.len() as f32;
    if mean_contrast < MIN_MEAN_CONTRAST {
        return None;
    }
    let ts: Vec<f64> = samples.iter().map(|s| s.t_s).collect();
    let fs: Vec<f32> = samples.iter().map(|s| s.freq_hz).collect();
    let slope = linear_slope(&ts, &fs);
    if slope >= 0.0 {
        // Reject upward / stationary sweeps: a flyover is a downward sweep.
        return None;
    }
    let mut down = 0usize;
    for w in samples.windows(2) {
        // Flat steps count as compliant: with fine bins a slow sweep often
        // repeats the same peak frequency between frames.
        if w[1].freq_hz <= w[0].freq_hz + 0.05 {
            down += 1;
        }
    }
    let downward_fraction = down as f32 / (samples.len() - 1).max(1) as f32;
    if downward_fraction < MIN_DOWNWARD_FRACTION {
        return None;
    }
    let start_s = samples.first().unwrap().t_s;
    let end_s = samples.last().unwrap().t_s;
    if (end_s - start_s) < MIN_DURATION_S {
        return None;
    }
    Some(Ridge {
        id: 0, // assigned later
        samples,
        start_s,
        end_s,
        mean_contrast,
        min_freq_hz: min_f,
        max_freq_hz: max_f,
        net_drop_hz: net_drop,
        slope_hz_s: slope,
        downward_fraction,
        fit: None,
        plausible: false,
    })
}

fn linear_slope(ts: &[f64], fs: &[f32]) -> f32 {
    let n = ts.len() as f64;
    let (sx, sy) = ts
        .iter()
        .zip(fs)
        .fold((0.0f64, 0.0f64), |(a, b), (&t, &f)| (a + t, b + f as f64));
    let mx = sx / n;
    let my = sy / n;
    let (mut sxx, mut sxy) = (0.0f64, 0.0f64);
    for (&t, &f) in ts.iter().zip(fs) {
        sxx += (t - mx) * (t - mx);
        sxy += (t - mx) * (f as f64 - my);
    }
    if sxx <= 1e-12 {
        0.0
    } else {
        (sxy / sxx) as f32
    }
}
