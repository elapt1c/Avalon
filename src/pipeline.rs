//! Analysis pipeline orchestration: decode -> resample -> spectrogram ->
//! ridge detection -> harmonic grouping -> Doppler parameterization.

use std::path::Path;
use std::sync::mpsc::Sender;

use anyhow::Result;
use serde::Serialize;

use crate::decode;
use crate::dsp::{self, doppler, harmonics, ridge, spectrogram};
use crate::model::{AnalysisResult, AnalysisSettings, Metadata};

#[derive(Debug, Clone)]
pub enum Progress {
    Stage { name: &'static str, fraction: f32 },
    Log(String),
    Done(AnalysisResult),
    Error(String),
}

pub fn run(
    path: &Path,
    settings: &AnalysisSettings,
    meta: &Metadata,
    tx: Option<Sender<Progress>>,
) -> Result<AnalysisResult> {
    let send = |p: Progress| {
        if let Some(tx) = &tx {
            let _ = tx.send(p);
        }
    };
    let t_start = std::time::Instant::now();
    let mut t_mark = t_start;
    let lap = |label: &str, send: &dyn Fn(Progress), t_mark: &mut std::time::Instant| {
        let secs = t_mark.elapsed().as_secs_f64();
        send(Progress::Log(format!("timing {label}: {secs:.2}s")));
        crate::log::info(
            "pipeline",
            &format!(
                "{label}: {secs:.2}s, rss {:.0} MB",
                crate::log::rss_bytes() as f64 / 1e6
            ),
        );
        *t_mark = std::time::Instant::now();
    };

    // 1. Decode.
    send(Progress::Stage { name: "Decoding media", fraction: 0.0 });
    // Decode straight to the analysis rate: the native-rate buffer for a long
    // recording is the single largest allocation in the whole pipeline.
    let (audio, backend) = decode::decode_mono_at(path, settings.target_sample_rate)?;
    send(Progress::Log(format!(
        "decoded via {backend}: {} Hz, {:.1} s, {} samples",
        audio.sample_rate,
        audio.duration_s,
        audio.samples.len()
    )));
    lap("decode", &send, &mut t_mark);

    // 2. Resample to the analysis rate.
    send(Progress::Stage { name: "Resampling", fraction: 0.05 });
    let samples = if audio.sample_rate == settings.target_sample_rate {
        audio.samples
    } else {
        dsp::resample_mono(&audio.samples, audio.sample_rate, settings.target_sample_rate)
    };
    let sr = settings.target_sample_rate;
    lap("resample", &send, &mut t_mark);

    // Frequency band: metadata hint overrides settings.
    let (min_f, max_f) = meta
        .freq_band_hz
        .unwrap_or((settings.min_freq_hz, settings.max_freq_hz));
    let (min_f, max_f) = (min_f.min(max_f), min_f.max(max_f));

    // 3./4. Spectrogram + ridge detection.
    //
    // Long recordings are analysed in overlapping blocks. A one-hour file
    // computed as a single spectrogram dilutes a faint flyover twice over: the
    // per-bin baselines are dominated by the whole hour, and one global DP pass
    // never reaches a curve that a 5-minute window resolves trivially.
    send(Progress::Stage { name: "Computing spectrogram", fraction: 0.15 });
    let params = spectrogram::StftParams {
        window_s: settings.window_s,
        hop_s: settings.hop_s,
        min_freq_hz: min_f,
        max_freq_hz: max_f,
        fft_pad: settings.fft_pad,
    };
    let duration = samples.len() as f64 / sr as f64;
    let (spec, ridges) = if duration <= BLOCK_S {
        let spec = spectrogram::compute_spectrogram(&samples, sr, &params);
        let ridges = ridge::detect_ridges(&spec, settings.seed_contrast);
        (spec, ridges)
    } else {
        analyse_blocks(&samples, sr, &params, settings.seed_contrast, &send)
    };
    let mut ridges = ridges;
    send(Progress::Log(format!(
        "spectrogram: {} bins x {} frames",
        spec.freq_axis_hz.len(),
        spec.time_axis_s.len()
    )));
    lap("spectrogram", &send, &mut t_mark);
    send(Progress::Stage { name: "Detecting ridges", fraction: 0.55 });
    send(Progress::Log(format!("detected {} ridges", ridges.len())));
    lap("ridge tracking", &send, &mut t_mark);

    // 5. Per-ridge Doppler fit + plausibility gate.
    send(Progress::Stage { name: "Fitting Doppler curves", fraction: 0.68 });
    for r in ridges.iter_mut() {
        r.fit = doppler::fit_doppler(r);
        r.plausible = doppler::is_plausible(r);
    }
    let n_plausible = ridges.iter().filter(|r| r.plausible).count();
    send(Progress::Log(format!(
        "{n_plausible}/{} ridges fit a plausible flyover curve",
        ridges.len()
    )));
    lap("doppler fit", &send, &mut t_mark);

    // 6. Harmonic grouping.
    send(Progress::Stage { name: "Grouping harmonics", fraction: 0.85 });
    let mut events = harmonics::group_harmonics(&ridges);
    for ev in events.iter_mut() {
        ev.doppler = ridges
            .iter()
            .find(|r| r.id == ev.representative_ridge_id)
            .and_then(|r| r.fit);
    }
    // Order events so plausible flyovers come first; the weak ones remain
    // selectable in the UI but are excluded from search by default.
    events.sort_by_key(|ev| {
        let plausible = ridges
            .iter()
            .find(|r| r.id == ev.representative_ridge_id)
            .map(|r| r.plausible)
            .unwrap_or(false);
        (!plausible, ev.start_s as i64)
    });
    for (i, ev) in events.iter_mut().enumerate() {
        ev.id = i;
    }
    send(Progress::Log(format!(
        "grouped into {} events ({n_plausible} plausible)",
        events.len()
    )));
    lap("harmonics", &send, &mut t_mark);
    send(Progress::Log(format!(
        "timing total: {:.2}s",
        t_start.elapsed().as_secs_f64()
    )));

    send(Progress::Stage { name: "Done", fraction: 1.0 });

    Ok(AnalysisResult {
        source_path: path.display().to_string(),
        sample_rate: audio.sample_rate,
        duration_s: audio.duration_s,
        backend,
        spectrogram: spec,
        ridges,
        events,
        settings: settings.clone(),
    })
}

/// Analysis block length for long recordings (s). Shortened automatically when
/// RAM is tight; see `plan_block_seconds`.
const BLOCK_S: f64 = 300.0;
/// Cap on the display spectrogram (bins x frames). The stitched full-resolution
/// matrix for an hour at 1600 Hz was ~1.4 GB; the app only ever renders
/// 2048x1024, so we pool while stitching instead of holding it all.
const DISPLAY_MAX_BINS: usize = 1024;
const DISPLAY_MAX_FRAMES: usize = 2048;

/// How many blocks to process at once: use the machine's cores, but cap it so a
/// 60 s minimum block length is still respected on huge machines.
fn block_workers() -> usize {
    std::env::var("AVALON_BLOCK_WORKERS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
                .clamp(1, 4)
        })
}

/// Pick a block length whose spectrogram plus FFT scratch fits the memory
/// budget. Cost per block ≈ bins x frames x 4 bytes, doubled for scratch.
fn plan_block_seconds(sr: u32, params: &spectrogram::StftParams) -> f64 {
    let win = (params.window_s * sr as f32).max(64.0) as usize;
    let fft = (win * params.fft_pad.max(1) as usize).next_power_of_two();
    let bin_hz = sr as f32 / fft as f32;
    let bins = (((params.max_freq_hz - params.min_freq_hz) / bin_hz) as f64).max(16.0);
    let frames_per_s = (1.0 / params.hop_s.max(0.01)) as f64;
    let bytes_per_block_s = bins * frames_per_s * 4.0 * 2.0; // matrix + scratch
    let budget = crate::log::budget_bytes() as f64;
    let fit = budget / bytes_per_block_s.max(1.0);
    fit.clamp(45.0, BLOCK_S)
}

/// Analyse a long recording block by block, returning a stitched *display*
/// spectrogram and every ridge found, with absolute times.
fn analyse_blocks(
    samples: &[f32],
    sr: u32,
    params: &spectrogram::StftParams,
    seed: f32,
    send: &dyn Fn(Progress),
) -> (crate::model::Spectrogram, Vec<crate::model::Ridge>) {
    let duration = samples.len() as f64 / sr as f64;
    // Run several blocks at once instead of one giant block serially: the RAM
    // budget is what limits block size, so splitting it across N workers keeps
    // peak memory the same while using N cores. This is the analysis stage's
    // real parallelism — the per-block FFT pass alone barely scaled.
    let workers = block_workers();
    let block_s = (plan_block_seconds(sr, params) / workers as f64).max(60.0);
    // Overlap stays at half a block. Measured on the reference clip: dropping it
    // to 20 % saves 15 s of the 56 s analysis but loses 9 of 37 ridges and 3 of
    // 17 events, i.e. it buys speed by missing weak flyovers. Keep the accuracy.
    let overlap_s = (block_s / 2.0).max(10.0);
    let step = (block_s - overlap_s).max(1.0);
    crate::log::info(
        "pipeline",
        &format!(
            "block plan: {block_s:.0}s blocks, {overlap_s:.0}s overlap, {:.0} MB budget",
            crate::log::budget_bytes() as f64 / 1e6
        ),
    );
    // Block ranges first, then process them on `workers` threads. Blocks are
    // independent (they only overlap for ridge continuity), so this is the
    // parallelism that matters: the per-block FFT pass alone barely scaled.
    let ranges: Vec<(usize, f64, f64)> = {
        let mut v = Vec::new();
        let mut start_s = 0.0f64;
        while start_s < duration {
            let end_s = (start_s + block_s).min(duration);
            let a = (start_s * sr as f64) as usize;
            let b = ((end_s * sr as f64) as usize).min(samples.len());
            if b <= a {
                break;
            }
            v.push((v.len() + 1, start_s, end_s));
            if end_s >= duration {
                break;
            }
            start_s = end_s - overlap_s;
        }
        v
    };
    let total_blocks = ranges.len().max(1);
    crate::log::info(
        "pipeline",
        &format!(
            "block plan: {total_blocks} blocks of {block_s:.0}s ({overlap_s:.0}s overlap) on {workers} worker(s)"
        ),
    );

    struct BlockOut {
        n: usize,
        start_s: f64,
        ridges: Vec<crate::model::Ridge>,
        grid: Option<crate::model::Spectrogram>,
    }

    let done = std::sync::atomic::AtomicUsize::new(0);
    let mut outs: Vec<BlockOut> = {
        use rayon::prelude::*;
        ranges
            .par_iter()
            .map(|&(n, start_s, end_s)| {
                let a = (start_s * sr as f64) as usize;
                let b = ((end_s * sr as f64) as usize).min(samples.len());
                let t0 = std::time::Instant::now();
                let spec = spectrogram::compute_spectrogram(&samples[a..b], sr, params);
                let t_spec = t0.elapsed();
                if spec.time_axis_s.is_empty() {
                    return BlockOut { n, start_s, ridges: Vec::new(), grid: None };
                }
                let t1 = std::time::Instant::now();
                let mut ridges = ridge::detect_ridges(&spec, seed);
                let t_ridge = t1.elapsed();
                for r in ridges.iter_mut() {
                    r.start_s += start_s;
                    r.end_s += start_s;
                    for smp in r.samples.iter_mut() {
                        smp.t_s += start_s;
                    }
                }
                // Each block pools into its own display grid; the merge below is
                // a max-fold, so order does not matter.
                let t2 = std::time::Instant::now();
                let mut grid = make_display_grid(params, duration);
                {
                    let rows = grid.contrast.len();
                    let cols = grid.time_axis_s.len();
                    let n_freq = spec.freq_axis_hz.len().min(spec.contrast.len());
                    for (fi, &t) in spec.time_axis_s.iter().enumerate() {
                        let abs_t = start_s + t;
                        let col = ((abs_t / duration) * (cols - 1) as f64).round() as usize;
                        let col = col.min(cols - 1);
                        for k in 0..n_freq {
                            let row = (k * rows / n_freq).min(rows - 1);
                            let v = spec.contrast[k][fi];
                            if v > grid.contrast[row][col] {
                                grid.contrast[row][col] = v;
                            }
                        }
                    }
                }
                // Progress is reported by the merge loop below: `send` is not
                // Sync, so it cannot be called from the worker threads.
                done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if std::env::var("AVALON_TIMING").is_ok() {
                    crate::log::info(
                        "pipeline",
                        &format!(
                            "block {n}: spectrogram {:.1}s, ridges {:.1}s, display-pool {:.1}s",
                            t_spec.as_secs_f64(),
                            t_ridge.as_secs_f64(),
                            t2.elapsed().as_secs_f64()
                        ),
                    );
                }
                BlockOut { n, start_s, ridges, grid: Some(grid) }
            })
            .collect()
    };
    outs.sort_by_key(|o| o.n);

    let mut ridges_all: Vec<crate::model::Ridge> = Vec::new();
    let mut merged: Option<crate::model::Spectrogram> = None;
    let mut merged_blocks = 0usize;
    for out in outs {
        merged_blocks += 1;
        send(Progress::Stage {
            name: "Analysing blocks",
            fraction: 0.55 + 0.4 * (merged_blocks as f32 / total_blocks as f32),
        });
        send(Progress::Log(format!(
            "block {}: {:.0}-{:.0}s, {} ridges",
            out.n,
            out.start_s,
            out.start_s + block_s,
            out.ridges.len()
        )));
        ridges_all.append(&mut { out.ridges });
        if let Some(g) = out.grid {
            match merged.as_mut() {
                Some(m) => {
                    for (row, rows_g) in g.contrast.iter().enumerate() {
                        let dst = &mut m.contrast[row];
                        for (col, v) in rows_g.iter().enumerate() {
                            if *v > dst[col] {
                                dst[col] = *v;
                            }
                        }
                    }
                }
                None => merged = Some(g),
            }
        }
    }

    let spec = match merged {
        Some(grid) => grid,
        None => pool_display(spectrogram::compute_spectrogram(samples, sr, params), duration),
    };
    let tol = 3.0 * spec.resolution_hz.max(0.5);
    let mut ridges = ridge::dedupe_tracks(ridges_all, tol);
    for (i, r) in ridges.iter_mut().enumerate() {
        r.id = i;
    }
    crate::log::info(
        "pipeline",
        &format!(
            "{} ridges, display {} bins x {} frames, rss {:.0} MB",
            ridges.len(),
            spec.freq_axis_hz.len(),
            spec.time_axis_s.len(),
            crate::log::rss_bytes() as f64 / 1e6
        ),
    );
    (spec, ridges)
}

/// Fixed-shape display grid: every block writes into the same rows/columns, so
/// the result is always rectangular.
fn make_display_grid(params: &spectrogram::StftParams, duration: f64) -> crate::model::Spectrogram {
    let rows = DISPLAY_MAX_BINS;
    let cols = DISPLAY_MAX_FRAMES;
    let f_min = params.min_freq_hz.max(5.0);
    let f_max = params.max_freq_hz.max(f_min + 1.0);
    let freq_axis: Vec<f32> = (0..rows)
        .map(|r| f_min + (f_max - f_min) * (r as f32 + 0.5) / rows as f32)
        .collect();
    let time_axis: Vec<f64> = (0..cols)
        .map(|c| duration * (c as f64 + 0.5) / cols as f64)
        .collect();
    crate::model::Spectrogram {
        time_axis_s: time_axis,
        freq_axis_hz: freq_axis,
        contrast: vec![vec![0.0f32; cols]; rows],
        resolution_hz: 0.0,
    }
}

/// Max-pool a spectrogram down to the display grid so long recordings do not
/// keep a gigabyte-scale matrix alive just to draw it.
fn pool_display(spec: crate::model::Spectrogram, _duration: f64) -> crate::model::Spectrogram {
    let nf = spec.freq_axis_hz.len();
    let nt = spec.time_axis_s.len();
    if nf <= DISPLAY_MAX_BINS && nt <= DISPLAY_MAX_FRAMES {
        return spec;
    }
    let fs = (nf + DISPLAY_MAX_BINS - 1) / DISPLAY_MAX_BINS;
    let ts = (nt + DISPLAY_MAX_FRAMES - 1) / DISPLAY_MAX_FRAMES;
    let rows = (nf + fs - 1) / fs;
    let cols = (nt + ts - 1) / ts;
    let mut contrast = vec![vec![0.0f32; cols]; rows];
    let mut freq_axis = Vec::with_capacity(rows);
    let mut time_axis = Vec::with_capacity(cols);
    for (r, row) in contrast.iter_mut().enumerate() {
        let k = (r * fs).min(nf - 1);
        freq_axis.push(spec.freq_axis_hz[k]);
        for (c, out) in row.iter_mut().enumerate() {
            let t0 = c * ts;
            let t1 = (t0 + ts).min(nt);
            let mut peak = 0.0f32;
            for kk in k..(k + fs).min(nf) {
                for tt in t0..t1 {
                    peak = peak.max(spec.contrast[kk][tt]);
                }
            }
            *out = peak;
        }
    }
    for c in 0..cols {
        time_axis.push(spec.time_axis_s[(c * ts).min(nt - 1)]);
    }
    crate::model::Spectrogram {
        time_axis_s: time_axis,
        freq_axis_hz: freq_axis,
        contrast,
        resolution_hz: spec.resolution_hz,
    }
}

// ---------------------------------------------------------------------------
// JSON export (summary, no full spectrogram).
// ---------------------------------------------------------------------------

#[derive(Debug, Serialize)]
pub struct SummaryReport<'a> {
    pub source: &'a str,
    pub backend_note: &'a str,
    pub sample_rate: u32,
    pub duration_s: f64,
    pub settings: &'a AnalysisSettings,
    pub metadata: &'a Metadata,
    pub ridge_count: usize,
    pub event_count: usize,
    pub ridges: &'a [crate::model::Ridge],
    pub events: &'a [crate::model::DetectedEvent],
}

pub fn save_summary(result: &AnalysisResult, meta: &Metadata, out_path: &Path) -> Result<()> {
    let report = SummaryReport {
        source: &result.source_path,
        backend_note: result.backend,
        sample_rate: result.sample_rate,
        duration_s: result.duration_s,
        settings: &result.settings,
        metadata: meta,
        ridge_count: result.ridges.len(),
        event_count: result.events.len(),
        ridges: &result.ridges,
        events: &result.events,
    };
    let json = serde_json::to_string_pretty(&report)?;
    std::fs::write(out_path, json)?;
    Ok(())
}
