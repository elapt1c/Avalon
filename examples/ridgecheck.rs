//! Isolate ridge detection: run the library STFT + detect_ridges on one clean
//! slice around an event, and report where the extracted tracks actually sit
//! relative to the strongest component of each frame.
//!
//! Usage: ridgecheck [start_s] [end_s]

use std::path::Path;

use avalon::dsp::{ridge, spectrogram};
use avalon::{decode, pipeline};
use avalon::model::{AnalysisSettings, Metadata};

const CLIP: &str = "/home/tpg/Downloads/Unknown location - raw audio - colsto (128p).mp4";
const SR: u32 = 4000;

fn main() -> anyhow::Result<()> {
    let t0: f64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(2060.0);
    let t1: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(2260.0);
    let clip = CLIP;

    let (audio, _) = decode::decode_mono_at(Path::new(clip), SR)?;
    let a0 = (t0 * SR as f64) as usize;
    let b0 = ((t1 * SR as f64) as usize).min(audio.samples.len());
    let slice = &audio.samples[a0..b0];
    println!("slice {t0}..{t1}s, {} samples", slice.len());

    let settings = AnalysisSettings { target_sample_rate: SR, ..Default::default() };
    let params = spectrogram::StftParams {
        window_s: settings.window_s,
        hop_s: settings.hop_s,
        min_freq_hz: settings.min_freq_hz,
        max_freq_hz: settings.max_freq_hz,
        fft_pad: settings.fft_pad,
    };
    println!(
        "stft window {}s hop {}s band {}-{} pad {}",
        params.window_s, params.hop_s, params.min_freq_hz, params.max_freq_hz, params.fft_pad
    );
    let spec = spectrogram::compute_spectrogram(slice, SR, &params);
    println!(
        "{} bins x {} frames, resolution {:.2} Hz",
        spec.freq_axis_hz.len(),
        spec.time_axis_s.len(),
        spec.resolution_hz
    );

    let ridges = ridge::detect_ridges(&spec, settings.seed_contrast);
    println!("{} ridges", ridges.len());
    let dt = spec.time_axis_s[1] - spec.time_axis_s[0];
    let fmin = spec.freq_axis_hz[0];
    let df = spec.freq_axis_hz[1] - spec.freq_axis_hz[0];
    for r in ridges.iter().take(8) {
        // For a few samples: contrast at the ridge vs the frame's strongest bin.
        let mut rows = Vec::new();
        for frac in [0.1, 0.5, 0.9] {
            let s = &r.samples[((r.samples.len() - 1) as f64 * frac) as usize];
            let col = nearest_t(&spec.time_axis_s, s.t_s);
            let k = ((s.freq_hz - fmin) / df).round() as usize;
            let k = k.min(spec.freq_axis_hz.len() - 1);
            let (bk, bc) = (0..spec.freq_axis_hz.len())
                .map(|i| (i, spec.contrast[i][col]))
                .fold((0usize, 0.0f32), |acc, x| if x.1 > acc.1 { x } else { acc });
            rows.push(format!(
                "t{:6.1} f{:7.1} c{:5.1} | peak {:7.1} c{:5.1}",
                s.t_s + t0,
                s.freq_hz,
                spec.contrast[k][col],
                spec.freq_axis_hz[bk],
                bc
            ));
        }
        println!(
            "  ridge {:?} plausible {} f {:.1}->{:.1} n{} mc{:.1} slope {:.1} | {}",
            r.id,
            r.plausible,
            r.max_freq_hz,
            r.min_freq_hz,
            r.samples.len(),
            r.mean_contrast,
            r.slope_hz_s,
            rows.join("  ")
        );
        if let Some(f) = r.fit {
            println!(
                "    fit f0 {:.1} v {:.0} d {:.0} rmse {:.2} inl {:.2}",
                f.f0_hz, f.speed_mps, f.slant_distance_m, f.rmse_hz, f.inlier_fraction
            );
        }
    }

    // What does the pipeline produce over the whole clip for comparison?
    if std::env::var_os("FULL").is_some() {
        let meta = Metadata::default();
        let a = pipeline::run(Path::new(clip), &settings, &meta, None)?;
        println!("pipeline ridges {}", a.ridges.len());
        for r in a.ridges.iter().filter(|r| r.end_s > t0 && r.start_s < t1) {
            println!(
                "  p ridge {:?} plausible {} {:.1}->{:.1}s f {:.1}->{:.1} n{} mc{:.1}",
                r.id,
                r.plausible,
                r.start_s,
                r.end_s,
                r.max_freq_hz,
                r.min_freq_hz,
                r.samples.len(),
                r.mean_contrast
            );
        }
    }
    Ok(())
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
