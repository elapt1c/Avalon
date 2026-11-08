//! Headless pipeline inspection: run analysis on a media file and print the
//! detected ridges/events with their Doppler fit quality.
//!
//! Usage: inspect <file> [--seed C] [--min HZ] [--max HZ] [--window S]
//!                       [--hop S] [--pad N] [--rate HZ]

use std::path::Path;

use avalon::model::{AnalysisSettings, Metadata};
use avalon::pipeline;

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() {
        eprintln!("usage: inspect <file> [--seed C] [--min HZ] [--max HZ] [--window S] [--hop S] [--pad N] [--rate HZ]");
        std::process::exit(2);
    }
    let path = args.remove(0);
    let mut settings = AnalysisSettings::default();
    let mut i = 0;
    while i + 1 < args.len() {
        let v = args[i + 1].parse::<f32>().unwrap_or(f32::NAN);
        match args[i].as_str() {
            "--seed" => settings.seed_contrast = v,
            "--min" => settings.min_freq_hz = v,
            "--max" => settings.max_freq_hz = v,
            "--window" => settings.window_s = v,
            "--hop" => settings.hop_s = v,
            "--pad" => settings.fft_pad = v as u32,
            "--rate" => settings.target_sample_rate = v as u32,
            other => {
                eprintln!("unknown flag {other}");
                std::process::exit(2);
            }
        }
        i += 2;
    }

    let (tx, rx) = std::sync::mpsc::channel();
    let result = pipeline::run(Path::new(&path), &settings, &Metadata::default(), Some(tx))?;
    while let Ok(p) = rx.try_recv() {
        if let avalon::pipeline::Progress::Log(msg) = p {
            println!("  · {msg}");
        }
    }

    println!("source: {}", result.source_path);
    println!(
        "decoded: {} Hz, {:.1}s via {}",
        result.sample_rate, result.duration_s, result.backend
    );
    println!(
        "spectrogram: {} bins x {} frames",
        result.spectrogram.freq_axis_hz.len(),
        result.spectrogram.time_axis_s.len()
    );
    println!("ridges: {}", result.ridges.len());
    let n_plausible = result.ridges.iter().filter(|r| r.plausible).count();
    println!("plausible flyover ridges: {n_plausible}");
    println!("events: {}", result.events.len());
    println!();

    println!("all ridge fits (pre-filter):");
    for r in &result.ridges {
        let lin = avalon::dsp::doppler::linear_rmse_hz(r).unwrap_or(f32::NAN);
        match r.fit {
            Some(f) => println!(
                "  ridge {:>2}: {:6.1}..{:6.1}s  {:7.1}->{:7.1} Hz (drop {:6.1})  \
                 f0 {:7.1}  v {:5.1}  d {:7.0}  rmse {:6.2}  lin {:6.2}  \
                 inliers {:4.0}%  down {:.2}  contrast {:.1}  {} {}",
                r.id,
                r.start_s,
                r.end_s,
                r.max_freq_hz,
                r.min_freq_hz,
                r.net_drop_hz,
                f.f0_hz,
                f.speed_mps,
                f.slant_distance_m,
                f.rmse_hz,
                lin,
                f.inlier_fraction * 100.0,
                r.downward_fraction,
                r.mean_contrast,
                if f.complete_sweep { "complete" } else { "partial" },
                if r.plausible { "PLAUSIBLE" } else { "weak" },
            ),
            None => println!("  ridge {:>2}: no fit", r.id),
        }
    }
    println!();

    for ev in &result.events {
        let rep = result
            .ridges
            .iter()
            .find(|r| r.id == ev.representative_ridge_id);
        match (ev.doppler, rep) {
            (Some(f), Some(r)) => println!(
                "event {:>2}: t {:6.1}..{:6.1}s  ridges {}  f0 {:7.1} Hz  v {:5.1} m/s  slant {:7.0} m  \
                 drop {:7.1} Hz  rmse {:5.2} Hz  inliers {:4.0}%  down {:.2}  {}",
                ev.id,
                ev.start_s,
                ev.end_s,
                ev.ridge_ids.len(),
                f.f0_hz,
                f.speed_mps,
                f.slant_distance_m,
                r.net_drop_hz,
                f.rmse_hz,
                f.inlier_fraction * 100.0,
                r.downward_fraction,
                if f.complete_sweep { "complete" } else { "partial" },
            ),
            _ => println!(
                "event {:>2}: t {:6.1}..{:6.1}s  ridges {}  (no fit)",
                ev.id,
                ev.start_s,
                ev.end_s,
                ev.ridge_ids.len()
            ),
        }
    }
    Ok(())
}
