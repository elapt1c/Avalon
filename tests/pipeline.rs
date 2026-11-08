//! End-to-end pipeline smoke test: synthesize a constant-velocity aircraft
//! Doppler sweep, write it as a WAV, and confirm the pipeline detects a ridge,
//! groups an event, and parameterizes a plausible Doppler fit.

use std::f64::consts::PI;
use std::path::PathBuf;

use avalon::model::{AnalysisSettings, Metadata};
use avalon::pipeline;

fn doppler_freq(t: f64, f0: f64, v: f64, d: f64, tcpa: f64) -> f64 {
    let c = 343.0;
    let x = v * (t - tcpa);
    let r = (d * d + x * x).sqrt();
    let ur = v * x / r;
    f0 * c / (c + ur)
}

fn synth_sweep(sr: u32, dur_s: f64, f0: f64, v: f64, d: f64, tcpa: f64) -> Vec<f32> {
    let n = (sr as f64 * dur_s) as usize;
    let mut phase = 0.0f64;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / sr as f64;
        let f = doppler_freq(t, f0, v, d, tcpa);
        out.push((phase.sin() * 0.5) as f32);
        phase += 2.0 * PI * f / sr as f64;
    }
    out
}

fn write_wav(path: &PathBuf, sr: u32, samples: &[f32]) {
    let mut data = Vec::with_capacity(samples.len() * 2);
    for s in samples {
        let v = (s.clamp(-1.0, 1.0) * 32767.0) as i16;
        data.extend_from_slice(&v.to_le_bytes());
    }
    let mut bytes = Vec::with_capacity(44 + data.len());
    bytes.extend_from_slice(b"RIFF");
    bytes.extend_from_slice(&((36 + data.len()) as u32).to_le_bytes());
    bytes.extend_from_slice(b"WAVE");
    bytes.extend_from_slice(b"fmt ");
    bytes.extend_from_slice(&16u32.to_le_bytes());
    bytes.extend_from_slice(&1u16.to_le_bytes()); // PCM
    bytes.extend_from_slice(&1u16.to_le_bytes()); // mono
    bytes.extend_from_slice(&sr.to_le_bytes());
    bytes.extend_from_slice(&(sr * 2).to_le_bytes()); // byte rate
    bytes.extend_from_slice(&2u16.to_le_bytes()); // block align
    bytes.extend_from_slice(&16u16.to_le_bytes()); // bits per sample
    bytes.extend_from_slice(b"data");
    bytes.extend_from_slice(&(data.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&data);
    std::fs::write(path, bytes).unwrap();
}

#[test]
fn detects_and_parameterizes_synthetic_flyover() {
    let sr = 44_100;
    let (f0, v, d, tcpa) = (100.0, 100.0, 1000.0, 20.0);
    let samples = synth_sweep(sr, 40.0, f0, v, d, tcpa);

    let mut path = std::env::temp_dir();
    path.push(format!("avalon_test_{}.wav", std::process::id()));
    write_wav(&path, sr, &samples);

    let result = pipeline::run(
        &path,
        &AnalysisSettings::default(),
        &Metadata::default(),
        None,
    )
    .expect("pipeline should run");

    let _ = std::fs::remove_file(&path);

    assert!(!result.ridges.is_empty(), "expected at least one ridge");
    assert!(!result.events.is_empty(), "expected at least one event");

    let fitted: Vec<_> = result.events.iter().filter_map(|e| e.doppler).collect();
    assert!(!fitted.is_empty(), "expected at least one Doppler fit");

    // Any plausible fit near the true 100 Hz tone satisfies the check (robust to
    // event ordering).
    let good = fitted.iter().any(|fit| {
        (70.0..=140.0).contains(&fit.f0_hz)
            && fit.speed_mps > 0.0
            && fit.speed_mps < 400.0
            && fit.slant_distance_m > 0.0
    });
    assert!(
        good,
        "expected a plausible flyover fit near 100 Hz, got {:?}",
        fitted.iter().map(|f| (f.f0_hz, f.speed_mps)).collect::<Vec<_>>()
    );

    // The plausibility gate must accept a real flyover, not just avoid false
    // positives: at least one ridge has to be classified as a flyover curve.
    let plausible = result.ridges.iter().filter(|r| r.plausible).count();
    assert!(
        plausible > 0,
        "plausibility gate rejected every synthetic flyover ridge ({} ridges)",
        result.ridges.len()
    );
}
