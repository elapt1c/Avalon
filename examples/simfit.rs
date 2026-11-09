//! Decisive test: at the KNOWN receiver position, how well can the two published
//! aircraft explain the two measured Doppler ridges?
//!
//! If this residual is comparable to the one the search finds at a wrong place,
//! the model/assumption is fine and the search/scoring is at fault. If it is far
//! worse, we are comparing against the wrong signal and no amount of tuning helps.
//!
//! Usage: simfit [analysis.json] [lat,lon] [icao0,icao1]

use std::path::PathBuf;
use std::sync::atomic::AtomicBool;

use avalon::flights;
use avalon::model::Flight as ArchiveFlight;
use avalon::model::{DetectedEvent, Ridge};
use avalon::score::{self, Atmos, Track};
use avalon::search::ecef;

#[derive(serde::Deserialize)]
struct Cached {
    events: Vec<DetectedEvent>,
    ridges: Vec<Ridge>,
}

fn build(f: &ArchiveFlight, dt: f64, center: (f64, f64), radius_m: f64) -> Option<Track> {
    // Sample the passage that comes near the receiver, the way the search does.
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
        // Linear interpolation between the bracketing samples.
        let i = near.partition_point(|p| p.t_s < t).min(near.len() - 1);
        let b = near[i];
        let a = near[i.saturating_sub(1)];
        let span = (b.t_s - a.t_s).max(1e-6);
        let w = ((t - a.t_s) / span).clamp(0.0, 1.0);
        let lat = a.lat + (b.lat - a.lat) * w;
        let lon = a.lon + (b.lon - a.lon) * w;
        let alt = a.alt_m + (b.alt_m - a.alt_m) * w;
        pos.push(ecef(lat, lon, alt));
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
    Some(Track {
        t0: t_first,
        dt,
        origin,
        pos: rel,
        bbox: [0.0, 0.0, 0.0, 0.0],
    })
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let cache = PathBuf::from(
        args.get(1)
            .cloned()
            .unwrap_or_else(|| "/root/projects/analysis-default.json".into()),
    );
    let truth_s = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "40.876466,-73.524951".into());
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
    let day = flights::fetch_flights(
        1788890400,
        1788910200,
        &region,
        "staging-0",
        &mut prog,
        &cancel,
    )?;
    println!("{} flights loaded", day.len());

    let obs = ecef(tlat, tlon, 2.0);
    for (k, icao) in icaos.iter().enumerate() {
        let Some(f) = day.iter().find(|f| f.icao.eq_ignore_ascii_case(icao)) else {
            println!("{icao}: not in archive");
            continue;
        };
        let Some(track) = build(f, 2.0, (tlat, tlon), 90_000.0) else {
            println!("{icao}: no passage near the receiver");
            continue;
        };
        let Some(ev) = cached.events.get(k) else {
            println!("{icao}: no matching event index");
            continue;
        };
        let Some(ridge) = cached
            .ridges
            .iter()
            .find(|r| r.id == ev.representative_ridge_id)
        else {
            println!("{icao}: no ridge");
            continue;
        };
        let curve: Vec<(f64, f32)> = ridge.samples.iter().map(|s| (s.t_s, s.freq_hz)).collect();
        let spec = ev.doppler.as_ref().map(|d| (d.f0_hz, d.t_cpa_s)).unwrap_or((120.0, 0.0));

        // Sweep the blade-pass frequency and the shared clock, at the TRUE place.
        let mut best = (f64::INFINITY, 0.0f64);
        for f0 in score::f0_grid_rel(spec.0 as f64, 41, 0.25) {
            for sh in score::shift_grid(60.0, 61) {
                if let Some(v) = score::curve_rmse_offset(
                    &curve,
                    &track,
                    obs,
                    f0,
                    1788898640.0, // close to the true recording start; the shift grid covers the rest
                    sh,
                    0.0,
                    Atmos::default(),
                    60.0,
                ) {
                    if v < best.0 {
                        best = (v, f0);
                    }
                }
            }
        }
        println!(
            "{icao}: <{:.2} Hz> best residual at the true location, f0 {:.1} Hz (event {k}, {} samples)",
            best.0,
            best.1,
            curve.len()
        );
    }
    Ok(())
}
