//! Where in a region does the *known* flight pair best explain the measured
//! Doppler curves? If the answer is not near the true receiver, the pairing or
//! the acoustic measurement is wrong — no scoring tweak can fix that.
//!
//! Usage: truthscore [W,S,E,N] [ICAO_event0,ICAO_event1]

use std::sync::atomic::AtomicBool;

use avalon::model::{Bounds, DetectedEvent, Ridge};
use avalon::{flights, score};

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut boxs = (-74.48f64, 40.15f64, -72.56f64, 41.59f64);
    if let Some(v) = args.first() {
        let p: Vec<f64> = v.split(',').filter_map(|x| x.parse().ok()).collect();
        if p.len() == 4 {
            boxs = (p[0], p[1], p[2], p[3]);
        }
    }
    let known: Vec<String> = args
        .get(1)
        .map(|s| s.split(',').map(|x| x.trim().to_lowercase()).collect())
        .unwrap_or_else(|| vec!["a67ba8".into(), "a659a2".into()]);

    let win = (
        avalon::dateutil::civil_to_unix(2026, 10, 30, 11, 0, 0),
        avalon::dateutil::civil_to_unix(2026, 10, 31, 23, 0, 0),
    );
    let region = Bounds { west: -74.48, south: 40.15, east: -72.56, north: 41.59 };
    let v: serde_json::Value =
        serde_json::from_slice(&std::fs::read("/root/projects/analysis-default.json")?)?;
    let events: Vec<DetectedEvent> = serde_json::from_value(v["events"].clone())?;
    let ridges: Vec<Ridge> = serde_json::from_value(v["ridges"].clone())?;

    let cancel = AtomicBool::new(false);
    let mut prog = |_p: flights::FlightProgress| {};
    let fl = flights::fetch_flights(win.0, win.1, &[region], "staging-0", &mut prog, &cancel)?;

    // Measured curves, built like search.rs does (high-contrast core only).
    let mut curves: Vec<Vec<(f64, f32)>> = Vec::new();
    for e in events.iter().take(2) {
        let r = ridges
            .iter()
            .find(|r| r.id == e.representative_ridge_id)
            .ok_or_else(|| anyhow::anyhow!("ridge"))?;
        let mut ws: Vec<f32> = r.samples.iter().map(|s| s.weight).collect();
        ws.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let med = ws[ws.len() / 2].max(1.0);
        let core: Vec<(f64, f32)> = r
            .samples
            .iter()
            .filter(|s| s.weight >= med)
            .map(|s| (s.t_s, s.freq_hz))
            .collect();
        let step = (core.len() / 120).max(1);
        curves.push(core.into_iter().step_by(step).collect());
    }

    let track_of = |icao: &str| -> Option<score::Track> {
        let f = fl.iter().find(|f| f.icao.to_lowercase() == icao)?;
        let pts: Vec<(f64, [f64; 3])> = f
            .points
            .iter()
            .map(|p| (p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m)))
            .collect();
        score::Track::from_points(&pts)
    };
    let mut pair = Vec::new();
    for (k, icao) in known.iter().enumerate() {
        let e = events.get(k).ok_or_else(|| anyhow::anyhow!("event"))?;
        let fit = e.doppler.ok_or_else(|| anyhow::anyhow!("fit"))?;
        let tr = track_of(icao).ok_or_else(|| anyhow::anyhow!("{icao} not in archive"))?;
        pair.push((tr, fit));
    }

    // Coarse scan: one shared clock offset per point, derived from where the
    // first flight's CPA falls; f0 free over a modest grid.
    let step_deg = 2.0 / 111.32;
    let mut hits: Vec<(f64, f64, f64, f64)> = Vec::new(); // (rmse, lat, lon, shift)
    let mut lat = boxs.1;
    while lat <= boxs.3 {
        let mut lon = boxs.0;
        while lon <= boxs.2 {
            let obs = avalon::search::ecef(lat, lon, 0.0);
            if let Some((_, t0)) = pair[0].0.cpa(obs) {
                let base = t0 - pair[0].1.t_cpa_s;
                let shifts = score::shift_grid(10.0, 5);
                let mut best = (f64::INFINITY, 0.0);
                for &sh in &shifts {
                    let mut sum = 0.0;
                    let mut ok = true;
                    for (k, (tr, fit)) in pair.iter().enumerate() {
                        let mut e = f64::INFINITY;
                        for &f0 in &score::f0_grid(fit.f0_hz as f64, 5) {
                            if let Some(v) = score::curve_rmse_offset(
                                &curves[k], tr, obs, f0, base, sh, 0.0,
                                score::Atmos::default(), 400.0,
                            ) {
                                if v < e {
                                    e = v;
                                }
                            }
                        }
                        if !e.is_finite() {
                            ok = false;
                            break;
                        }
                        sum += e;
                    }
                    if ok && sum / 2.0 < best.0 {
                        best = (sum / 2.0, sh);
                    }
                }
                if best.0.is_finite() {
                    hits.push((best.0, lat, lon, best.1));
                }
            }
            lon += step_deg;
        }
        lat += step_deg * 0.75; // ~2 km
    }
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    println!("{} grid points, best 10 for the KNOWN pair:", hits.len());
    for (r, la, lo, sh) in hits.iter().take(10) {
        let d = haversine((*la, *lo), (40.87, -73.52)) / 1000.0;
        println!("  {la:.3},{lo:.3} rmse {r:.2} Hz shift {sh:+.0}s  ({d:.0} km from truth)");
    }
    Ok(())
}

fn haversine(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6_371_000.0f64;
    let (la1, la2) = (a.0.to_radians(), b.0.to_radians());
    let h = ((b.0 - a.0).to_radians() * 0.5).sin().powi(2)
        + la1.cos() * la2.cos() * ((b.1 - a.1).to_radians() * 0.5).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}
