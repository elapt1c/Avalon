//! Where does each extracted ridge sit relative to the bright spectral band?
//!
//! Recomputes a raw-magnitude STFT over each event's neighbourhood and, for
//! every ridge sample, compares the sample frequency with the local magnitude
//! maximum. Prints the offset distribution so a systematic bias is visible.
//!
//! Usage: bandcheck [clip]

use std::path::Path;

use avalon::model::{AnalysisSettings, Metadata};
use avalon::{decode, dsp, pipeline};
use rustfft::num_complex::Complex;
use rustfft::FftPlanner;

const CLIP: &str = "/home/tpg/Downloads/Unknown location - raw audio - colsto (128p).mp4";
const SR: u32 = 4000;

fn main() -> anyhow::Result<()> {
    let clip = std::env::args().nth(1).unwrap_or_else(|| CLIP.into());
    let settings = AnalysisSettings { target_sample_rate: SR, ..Default::default() };
    let meta = Metadata::default();
    println!("analysing");
    let a = pipeline::run(Path::new(&clip), &settings, &meta, None)?;
    let (audio, _) = decode::decode_mono_at(Path::new(&clip), SR)?;
    println!(
        "{} ridges, {} events, {:.1}s audio",
        a.ridges.len(),
        a.events.len(),
        audio.duration_s
    );

    for e in &a.events {
        if e.id == 0 {
            continue;
        }
        let Some(fit) = e.doppler else { continue };
        if !fit.f0_hz.is_finite() || fit.f0_hz <= 0.0 {
            continue;
        }
        let Some(rep) = a.ridges.iter().find(|r| r.id == e.representative_ridge_id) else { continue };
        println!(
            "\nevent {} reps ridge {} plausible {} f0 {:.1} Hz v {:.0} m/s d {:.0} m rmse {:.2} inl {:.2}",
            e.id,
            rep.id,
            rep.plausible,
            fit.f0_hz,
            fit.speed_mps,
            fit.slant_distance_m,
            fit.rmse_hz,
            fit.inlier_fraction
        );
        println!(
            "  ridge {:.1}->{:.1}s f {:.1}->{:.1} Hz drop {:.1} mean_contrast {:.2} samples {}",
            rep.start_s,
            rep.end_s,
            rep.max_freq_hz,
            rep.min_freq_hz,
            rep.net_drop_hz,
            rep.mean_contrast,
            rep.samples.len()
        );

        // Raw magnitude over a generous slice around the event.
        let slice_pad = 30.0;
        let a0 = (((rep.start_s - slice_pad).max(0.0)) * SR as f64) as usize;
        let b0 = (((rep.end_s + slice_pad) * SR as f64) as usize).min(audio.samples.len());
        if b0 <= a0 + SR as usize {
            println!("  slice too short");
            continue;
        }
        let lo = (0.55 * fit.f0_hz) as f32;
        let hi = (1.65 * fit.f0_hz) as f32;
        let (mag, t_axis, f_axis) =
            raw_stft(&audio.samples[a0..b0], SR, 3.0, 0.25, lo, hi);
        if mag.is_empty() {
            println!("  no stft");
            continue;
        }
        let nf = f_axis.len();
        let dt = t_axis[1] - t_axis[0];
        let t_off = a0 as f64 / SR as f64;

        // Per-bin median baseline, exactly as the library computes contrast.
        let baselines: Vec<f32> = (0..nf)
            .map(|k| {
                let mut v: Vec<f32> = mag[k].clone();
                v.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
                v[v.len() / 2].max(1e-6)
            })
            .collect();

        let mut offs: Vec<f32> = Vec::new();
        let mut abs_offs: Vec<f32> = Vec::new();
        let mut miss = 0usize;
        for s in rep.samples.iter().step_by(4) {
            let t = s.t_s - t_off;
            if t < 0.0 || t > (t_axis.len() as f64 - 1.0) * dt {
                continue;
            }
            let col = nearest_t(&t_axis, t);
            // Global brightest bin in the analysed band.
            let mut bk = 0usize;
            for k in 1..nf {
                if mag[k][col] > mag[bk][col] {
                    bk = k;
                }
            }
            // Brightest bin within +-6 Hz of the ridge sample.
            let gate = 6.0 / (f_axis[1] - f_axis[0]).max(1e-3);
            let mut lk: Option<usize> = None;
            let k0 = ((s.freq_hz - lo) / (f_axis[1] - f_axis[0])).round();
            if k0.is_finite() && k0 >= 0.0 {
                let start = (k0 as isize - gate as isize).max(0) as usize;
                let end = ((k0 as isize + gate as isize) as usize).min(nf - 1);
                for k in start..=end {
                    if lk.map_or(true, |b| mag[k][col] > mag[b][col]) {
                        lk = Some(k);
                    }
                }
            }
            offs.push(s.freq_hz - f_axis[bk]);
            match lk {
                Some(k) => abs_offs.push((s.freq_hz - f_axis[k]).abs()),
                None => miss += 1,
            }
        }
        if offs.is_empty() {
            println!("  no overlapping samples");
            continue;
        }
        offs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        abs_offs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        let q = |v: &[f32], p: f64| v[((v.len() - 1) as f64 * p) as usize];
        println!(
            "  offset vs global band peak: median {:+.1} Hz  p10 {:+.1}  p90 {:+.1}  |x| med {:.1}",
            q(&offs, 0.5),
            q(&offs, 0.1),
            q(&offs, 0.9),
            q(&abs_offs, 0.5)
        );
        println!("  near-miss samples {miss}");

        // Contrast/magnitude breakdown at a few ridge samples: shows whether the
        // ridge rides a genuine amplitude peak or a quiet-bin contrast artifact.
        let n_show = rep.samples.len();
        for frac in [0.05, 0.35, 0.65, 0.95] {
            let s = &rep.samples[((n_show - 1) as f64 * frac) as usize];
            let t = s.t_s - t_off;
            if t < 0.0 {
                continue;
            }
            let col = nearest_t(&t_axis, t);
            let fmax = (0..nf).map(|k| mag[k][col]).fold(0.0f32, f32::max).max(1e-9);
            let mut rows: Vec<(f32, f32, f32)> = (0..nf)
                .filter_map(|k| {
                    let c = mag[k][col] / baselines[k];
                    if c > 1.5 {
                        Some((f_axis[k], c, mag[k][col] / fmax))
                    } else {
                        None
                    }
                })
                .collect();
            rows.sort_by(|x, y| y.1.partial_cmp(&x.1).unwrap_or(std::cmp::Ordering::Equal));
            rows.truncate(6);
            println!(
                "  t={:.1}s ridge {:.1} Hz (contrast {:.1}, rel {:.2}) -> {:?}",
                s.t_s,
                s.freq_hz,
                mag[nearest(&f_axis, s.freq_hz)][col] / baselines[nearest(&f_axis, s.freq_hz)],
                mag[nearest(&f_axis, s.freq_hz)][col] / fmax,
                rows.iter()
                    .map(|(f, c, r)| format!("{f:.1}Hz c{c:.1} a{r:.2}"))
                    .collect::<Vec<_>>()
            );
        }

        // First/last sample detail.
        for s in [rep.samples.first(), rep.samples.last()].into_iter().flatten() {
            let t = s.t_s - t_off;
            if t < 0.0 {
                continue;
            }
            let col = nearest_t(&t_axis, t);
            let mut tops: Vec<(f32, f32)> = (0..nf).map(|k| (f_axis[k], mag[k][col])).collect();
            tops.sort_by(|x, y| y.1.partial_cmp(&x.1).unwrap_or(std::cmp::Ordering::Equal));
            tops.truncate(4);
            println!(
                "  t={:.1}s ridge {:.1} Hz; top bins {:?}",
                s.t_s,
                s.freq_hz,
                tops.iter()
                    .map(|(f, m)| format!("{f:.1}Hz/{m:.0}"))
                    .collect::<Vec<_>>()
            );
        }
    }
    Ok(())
}

/// Raw-magnitude STFT (no per-bin normalisation), for band comparison only.
fn raw_stft(
    samples: &[f32],
    sr: u32,
    window_s: f32,
    hop_s: f32,
    lo: f32,
    hi: f32,
) -> (Vec<Vec<f32>>, Vec<f64>, Vec<f32>) {
    let win = (window_s * sr as f32).round() as usize;
    if win == 0 || samples.len() < win {
        return (Vec::new(), Vec::new(), Vec::new());
    }
    let hop = (hop_s * sr as f32).round().max(1.0) as usize;
    let fft_size = (win * 4).next_power_of_two();
    let n_frames = (samples.len() - win) / hop + 1;
    let window: Vec<f32> = (0..win)
        .map(|i| 0.5 - 0.5 * (2.0 * std::f32::consts::PI * i as f32 / (win as f32 - 1.0)).cos())
        .collect();
    let bin_hz = sr as f32 / fft_size as f32;
    let k0 = (lo / bin_hz).floor().max(1.0) as usize;
    let k1 = ((hi / bin_hz).ceil() as usize).min(fft_size / 2 - 1);
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(fft_size);
    let mut buf = vec![Complex::new(0.0f32, 0.0f32); fft_size];
    let mut mag = vec![vec![0.0f32; n_frames]; k1 - k0 + 1];
    for f in 0..n_frames {
        let start = f * hop;
        buf.fill(Complex::new(0.0, 0.0));
        for (i, w) in window.iter().enumerate() {
            buf[i] = Complex::new(samples[start + i] * w, 0.0);
        }
        fft.process(&mut buf);
        for k in k0..=k1 {
            mag[k - k0][f] = buf[k].norm();
        }
    }
    let t: Vec<f64> = (0..n_frames)
        .map(|f| (f * hop + win / 2) as f64 / sr as f64)
        .collect();
    let fr: Vec<f32> = (k0..=k1).map(|k| k as f32 * bin_hz).collect();
    (mag, t, fr)
}

fn nearest(axis: &[f32], f: f32) -> usize {
    axis.iter()
        .enumerate()
        .min_by(|a, b| {
            (a.1 - f)
                .abs()
                .partial_cmp(&(b.1 - f).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0)
}

fn nearest_t(axis: &[f64], t: f64) -> usize {
    match axis.binary_search_by(|v| v.partial_cmp(&t).unwrap_or(std::cmp::Ordering::Equal)) {
        Ok(i) => i,
        Err(i) => {
            if i == 0 {
                0
            } else if i >= axis.len() {
                axis.len() - 1
            } else if (axis[i] - t).abs() < (t - axis[i - 1]).abs() {
                i
            } else {
                i - 1
            }
        }
    }
}
