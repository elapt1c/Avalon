//! Measurement calibration: synthesise a flyover from a *known* track, run the
//! real spectrogram -> ridge -> Doppler-fit path on it, and compare the
//! recovered geometry with the truth.
//!
//! The scorer is already validated by synthtest; this checks the half that
//! produces the curve, which is where a bias shows up as a wrong distance.
//!
//! Usage: fitcheck [radius_km]

use std::sync::atomic::AtomicBool;

use avalon::dsp::{doppler, ridge, spectrogram};
use avalon::model::{AnalysisSettings, Bounds};
use avalon::score;

const TRUTH: (f64, f64) = (40.87, -73.52);
const SR: u32 = 4000;
const F0: f64 = 110.0;

fn main() -> anyhow::Result<()> {
    let radius: f64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(15.0);
    let region = Bounds {
        west: TRUTH.1 - 0.9555,
        south: TRUTH.0 - 0.7226,
        east: TRUTH.1 + 0.9555,
        north: TRUTH.0 + 0.7226,
    };
    let window = (
        avalon::dateutil::civil_to_unix(2026, 10, 31, 11, 0, 0),
        avalon::dateutil::civil_to_unix(2026, 10, 31, 23, 0, 0),
    );
    let cancel = AtomicBool::new(false);
    let mut prog = |_p: avalon::flights::FlightProgress| {};
    let fl = avalon::flights::fetch_flights(window.0, window.1, &[region], "staging-0", &mut prog, &cancel)?;
    println!("{} flights in window", fl.len());

    let obs = avalon::search::ecef(TRUTH.0, TRUTH.1, 0.0);
    // One flight per distance band, so a distance-dependent bias is visible.
    let mut picks: Vec<(f64, String, score::Track, f64)> = Vec::new();
    for f in &fl {
        let mut pts: Vec<(f64, [f64; 3])> = Vec::new();
        let mut last = f64::NEG_INFINITY;
        for p in &f.points {
            if p.t_s - last < 2.0 {
                continue;
            }
            last = p.t_s;
            pts.push((p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m)));
        }
        let Some(tr) = score::Track::from_points(&pts) else { continue };
        let Some((d, t_cpa)) = tr.cpa(obs) else { continue };
        if d > radius * 1000.0 {
            continue;
        }
        let v = tr.vel(t_cpa).map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()).unwrap_or(0.0);
        if v < 40.0 {
            continue;
        }
        // Keep the closest-to-truth flight in each 2 km shell.
        let shell = (d / 2000.0).floor() as i64;
        if let Some(slot) = picks.iter_mut().find(|p| (p.0 / 2000.0).floor() as i64 == shell) {
            if d < slot.3 {
                *slot = (d, f.icao.clone(), tr, d);
            }
        } else {
            picks.push((d, f.icao.clone(), tr, d));
        }
    }
    picks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    println!("{} calibration flights", picks.len());

    let settings = AnalysisSettings { target_sample_rate: SR, ..Default::default() };
    let params = spectrogram::StftParams {
        window_s: settings.window_s,
        hop_s: settings.hop_s,
        min_freq_hz: settings.min_freq_hz,
        max_freq_hz: settings.max_freq_hz,
        fft_pad: settings.fft_pad,
    };

    println!("\n  truth(d km, v m/s)   | recovered (f0, v, d, rmse, inl)   | err d");
    for (d_true, icao, track, _) in &picks {
        let (_, t_cpa) = track.cpa(obs).unwrap();
        let dur = 240.0f64;
        let t0 = t_cpa - dur / 2.0;
        let n = (SR as f64 * dur) as usize;
        let mut audio = vec![0.0f32; n];
        let mut phases = [0.0f64; 3];
        for (i, s) in audio.iter_mut().enumerate() {
            let t = t0 + i as f64 / SR as f64;
            let Some(ratio) = score::doppler_ratio(track, obs, t, score::Atmos::default()) else {
                continue;
            };
            let p = track.at(t).unwrap();
            let r = ((p[0] - obs[0]).powi(2) + (p[1] - obs[1]).powi(2) + (p[2] - obs[2]).powi(2)).sqrt();
            let amp = 1500.0 / r.max(200.0);
            let mut v = 0.0f64;
            for (h, ph) in [1.0f64, 2.0, 3.0].iter().zip(phases.iter_mut()) {
                v += amp * (0.6 / h) * ph.sin();
                *ph += 2.0 * std::f64::consts::PI * F0 * h * ratio / SR as f64;
            }
            *s = (v * 0.35).clamp(-1.0, 1.0) as f32;
        }
        let spec = spectrogram::compute_spectrogram(&audio, SR, &params);
        let ridges = ridge::detect_ridges(&spec, settings.seed_contrast);
        let best = ridges
            .iter()
            .filter(|r| r.plausible)
            .filter_map(|r| r.fit.map(|f| (f, r)))
            .min_by(|a, b| a.0.rmse_hz.partial_cmp(&b.0.rmse_hz).unwrap_or(std::cmp::Ordering::Equal));
        let v_true = track
            .vel(t_cpa)
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .unwrap_or(0.0);
        match best {
            Some((f, _)) => println!(
                "  {:5.1} km {:5.0} m/s {icao} | f0 {:6.1} v {:5.0} d {:5.1} km rmse {:4.2} inl {:3.2} | {:+.1} km ({:+.0}%)",
                d_true / 1000.0,
                v_true,
                f.f0_hz,
                f.speed_mps,
                f.slant_distance_m / 1000.0,
                f.rmse_hz,
                f.inlier_fraction,
                f.slant_distance_m / 1000.0 - (*d_true as f32) / 1000.0,
                100.0 * (f.slant_distance_m / (*d_true as f32).max(1.0) - 1.0),
            ),
            None => {
                println!(
                    "  {:5.1} km {:5.0} m/s {icao} | no plausible ridge ({} ridges)",
                    d_true / 1000.0,
                    v_true,
                    ridges.len()
                );
                for r in ridges.iter().take(4) {
                    println!(
                        "      ridge f {:.0}->{:.0} n{} mc{:.1} fit {:?} why {}",
                        r.max_freq_hz,
                        r.min_freq_hz,
                        r.samples.len(),
                        r.mean_contrast,
                        r.fit.map(|f| format!(
                            "f0 {:.0} v {:.0} d {:.1}km rmse {:.2} inl {:.2}",
                            f.f0_hz,
                            f.speed_mps,
                            f.slant_distance_m / 1000.0,
                            f.rmse_hz,
                            f.inlier_fraction
                        )),
                        why(r)
                    );
                }
            }
        }
        let _ = doppler::is_plausible;
    }
    Ok(())
}

/// Which plausibility test rejects this ridge.
fn why(r: &avalon::model::Ridge) -> String {
    let Some(f) = r.fit else {
        return "no fit".into();
    };
    if !(40.0..=300.0).contains(&f.speed_mps) {
        return format!("speed {:.0}", f.speed_mps);
    }
    if f.slant_distance_m < 200.0 {
        return "slant".into();
    }
    if f.f0_hz < r.min_freq_hz * 0.96 || f.f0_hz > r.max_freq_hz * 1.04 {
        return format!("f0 {:.0} outside {:.0}..{:.0}", f.f0_hz, r.min_freq_hz, r.max_freq_hz);
    }
    if f.inlier_fraction < 0.8 {
        return format!("inliers {:.2}", f.inlier_fraction);
    }
    if f.rmse_hz > (0.025 * f.f0_hz).max(2.5) {
        return format!("rmse {:.2}", f.rmse_hz);
    }
    match doppler::linear_rmse_hz(r) {
        Some(l) if f.rmse_hz <= 0.75 * l => "ok".into(),
        Some(l) => format!("doppler {:.2} vs linear {:.2}", f.rmse_hz, l),
        None => "linear".into(),
    }
}
