//! Is the recording's timebase slightly off? Sweep a playback-speed multiplier k
//! and see which value makes the two published aircraft fit the two measured
//! ridges simultaneously, at the known receiver position.
//!
//! Usage: speedfit [analysis.json] [lat,lon] [icao0,icao1]
use std::sync::atomic::AtomicBool;

use avalon::flights;
use avalon::model::{DetectedEvent, Flight as ArchiveFlight, Ridge};
use avalon::score::{self, Atmos, Track};
use avalon::search::ecef;

#[derive(serde::Deserialize)]
struct Cached {
    events: Vec<DetectedEvent>,
    ridges: Vec<Ridge>,
}

fn build(f: &ArchiveFlight, dt: f64, center: (f64, f64), radius_m: f64) -> Option<Track> {
    let near: Vec<&avalon::model::FlightPoint> = f
        .points
        .iter()
        .filter(|p| {
            let dlat = (p.lat - center.0) * 111_320.0;
            let dlon = (p.lon - center.1) * 111_320.0 * center.0.to_radians().cos();
            (dlat * dlat + dlon * dlon).sqrt() <= radius_m
        })
        .collect();
    if near.len() < 4 {
        return None;
    }
    let (t_first, t_last) = (near[0].t_s, near[near.len() - 1].t_s);
    let mut pos = Vec::new();
    let mut t = t_first;
    while t <= t_last {
        let i = near.partition_point(|p| p.t_s < t).min(near.len() - 1);
        let b = near[i];
        let a = near[i.saturating_sub(1)];
        let span = (b.t_s - a.t_s).max(1e-6);
        let w = ((t - a.t_s) / span).clamp(0.0, 1.0);
        pos.push(ecef(
            a.lat + (b.lat - a.lat) * w,
            a.lon + (b.lon - a.lon) * w,
            a.alt_m + (b.alt_m - a.alt_m) * w,
        ));
        t += dt;
    }
    if pos.len() < 4 {
        return None;
    }
    let origin = pos[0];
    let rel: Vec<[f32; 3]> = pos
        .iter()
        .map(|p| {
            [
                (p[0] - origin[0]) as f32,
                (p[1] - origin[1]) as f32,
                (p[2] - origin[2]) as f32,
            ]
        })
        .collect();
    Some(Track { t0: t_first, dt, origin, pos: rel, bbox: [0.0; 4] })
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let cache = args.get(1).cloned().unwrap_or_else(|| "/root/projects/analysis-default.json".into());
    let truth_s = args.get(2).cloned().unwrap_or_else(|| "40.876466,-73.524951".into());
    let icaos: Vec<String> = args
        .get(3)
        .cloned()
        .unwrap_or_else(|| "a1dfeb,ace391".into())
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    let (tlat, tlon): (f64, f64) = {
        let (a, b) = truth_s.split_once(',').expect("lat,lon");
        (a.trim().parse()?, b.trim().parse()?)
    };

    let cached: Cached = serde_json::from_slice(&std::fs::read(&cache)?)?;
    let region = vec![avalon::model::Bounds {
        west: tlon - 1.5,
        south: tlat - 1.5,
        east: tlon + 1.5,
        north: tlat + 1.5,
    }];
    let mut prog = |_p: flights::FlightProgress| {};
    let cancel = AtomicBool::new(false);
    let day = flights::fetch_flights(1788890400, 1788910200, &region, "staging-0", &mut prog, &cancel)?;

    let obs = ecef(tlat, tlon, 2.0);
    let base0 = 1788898640.0;

    // Gather the two (aircraft, ridge) legs once.
    let mut legs: Vec<(String, Track, Vec<(f64, f32)>, f64)> = Vec::new();
    for (k, icao) in icaos.iter().enumerate() {
        let Some(f) = day.iter().find(|f| f.icao.eq_ignore_ascii_case(icao)) else {
            println!("{icao}: not in archive");
            continue;
        };
        let Some(track) = build(f, 2.0, (tlat, tlon), 90_000.0) else {
            println!("{icao}: no passage");
            continue;
        };
        let Some(ev) = cached.events.get(k) else { continue };
        let Some(ridge) = cached.ridges.iter().find(|r| r.id == ev.representative_ridge_id) else {
            continue;
        };
        let curve: Vec<(f64, f32)> = ridge.samples.iter().map(|s| (s.t_s, s.freq_hz)).collect();
        let f0 = ev.doppler.as_ref().map(|d| d.f0_hz as f64).unwrap_or(120.0);
        legs.push((icao.clone(), track, curve, f0));
    }
    if legs.len() < 2 {
        anyhow::bail!("need two legs, got {}", legs.len());
    }

    // For each speed multiplier: fit each leg independently (free f0, free shift)
    // with the audio time axis stretched by k, and report the joint residual.
    println!("{:>8} {:>10} {:>10} {:>10} {:>12}", "k", "res0 Hz", "res1 Hz", "joint Hz", "start split s");
    let mut best = (f64::INFINITY, 1.0);
    let mut k = 0.96f64;
    while k <= 1.0601 {
        let mut out = Vec::new();
        let mut starts = Vec::new();
        for (_icao, track, curve, f0_audio) in &legs {
            // Stretch the audio timeline (and, optionally, the frequencies).
            let scaled: Vec<(f64, f32)> = curve
                .iter()
                .map(|(t, f)| (t * k, if std::env::var("SPEED_PITCH").is_ok() { f * k as f32 } else { *f }))
                .collect();
            let mut best_leg = (f64::INFINITY, 0.0f64, 0.0f64);
            for f0 in score::f0_grid_rel(*f0_audio, 13, 0.18) {
                for sh in score::shift_grid(400.0, 81) {
                    if let Some(v) = score::curve_rmse_offset(
                        &scaled,
                        track,
                        obs,
                        f0,
                        base0,
                        sh,
                        0.0,
                        Atmos::default(),
                        60.0,
                    ) {
                        if v < best_leg.0 {
                            best_leg = (v, f0, sh);
                        }
                    }
                }
            }
            out.push(best_leg);
            starts.push(base0 + best_leg.2);
        }
        let joint = out.iter().map(|v| v.0).sum::<f64>() / out.len() as f64;
        let split = (starts[0] - starts[1]).abs();
        println!(
            "{k:>8.3} {:>10.2} {:>10.2} {joint:>10.2} {split:>12.1}",
            out[0].0, out[1].0
        );
        // Prefer a low joint residual whose two fits also agree on the start.
        let score_k = joint + split * 0.05;
        if score_k < best.0 {
            best = (score_k, k);
        }
        k += 0.005;
    }
    println!("\nbest k = {:.3} (audio timeline stretched by {:.1}%)", best.1, (best.1 - 1.0) * 100.0);
    Ok(())
}
