//! How close does the BEST-fitting location for the KNOWN flight pair get to the
//! true receiver? This isolates localization error from candidate-selection
//! error: whatever the search's #1 pick is, this tells us where the correct
//! flights actually put the receiver.
//!
//! Usage: closestspike [W,S,E,N] [ICAO0,ICAO1] [coarse_km]
//! default box = CONUS, default pair = a67ba8,a659a2, default coarse 40 km

use std::sync::atomic::AtomicBool;

use avalon::model::{Bounds, DetectedEvent, Ridge};
use avalon::{flights, score};

fn haversine(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6_371_000.0f64;
    let (la1, la2) = (a.0.to_radians(), b.0.to_radians());
    let h = ((b.0 - a.0).to_radians() * 0.5).sin().powi(2)
        + la1.cos() * la2.cos() * ((b.1 - a.1).to_radians() * 0.5).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (w, s, e, n) = match args.first() {
        Some(v) => {
            let p: Vec<f64> = v.split(',').filter_map(|x| x.parse().ok()).collect();
            if p.len() == 4 { (p[0], p[1], p[2], p[3]) } else { (-124.8, 24.5, -66.9, 49.4) }
        }
        None => (-124.8, 24.5, -66.9, 49.4),
    };
    let known: Vec<String> = args.get(1).map(|s| s.split(',').map(|x| x.trim().to_lowercase()).collect())
        .unwrap_or_else(|| vec!["a67ba8".into(), "a659a2".into()]);
    let coarse_km: f64 = args.get(2).and_then(|v| v.parse().ok()).unwrap_or(40.0);

    let truth = (40.87f64, -73.52f64);

    let win = (
        avalon::dateutil::civil_to_unix(2026, 10, 30, 11, 0, 0),
        avalon::dateutil::civil_to_unix(2026, 10, 31, 23, 0, 0),
    );
    let region = Bounds { west: w, south: s, east: e, north: n };
    let v: serde_json::Value =
        serde_json::from_slice(&std::fs::read("/root/projects/analysis-default.json")?)?;
    let events: Vec<DetectedEvent> = serde_json::from_value(v["events"].clone())?;
    let ridges: Vec<Ridge> = serde_json::from_value(v["ridges"].clone())?;

    let cancel = AtomicBool::new(false);
    let mut prog = |_p: flights::FlightProgress| {};
    eprintln!("loading archive for {}..{} over the box", win.0, win.1);
    let fl = flights::fetch_flights(win.0, win.1, &[region], "staging-0", &mut prog, &cancel)?;
    eprintln!("{} flights in box", fl.len());

    let mut curves: Vec<Vec<(f64, f32)>> = Vec::new();
    for ev in events.iter().take(2) {
        let r = ridges.iter().find(|r| r.id == ev.representative_ridge_id)
            .ok_or_else(|| anyhow::anyhow!("ridge"))?;
        let mut ws: Vec<f32> = r.samples.iter().map(|x| x.weight).collect();
        ws.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let med = ws[ws.len() / 2].max(1.0);
        let core: Vec<(f64, f32)> = r.samples.iter().filter(|x| x.weight >= med)
            .map(|x| (x.t_s, x.freq_hz)).collect();
        let step = (core.len() / 120).max(1);
        curves.push(core.into_iter().step_by(step).collect());
    }

    let track_of = |icao: &str| -> Option<score::Track> {
        let f = fl.iter().find(|f| f.icao.to_lowercase() == icao)?;
        let pts: Vec<(f64, [f64; 3])> = f.points.iter()
            .map(|p| (p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m))).collect();
        score::Track::from_points(&pts)
    };
    let mut pair = Vec::new();
    for (k, icao) in known.iter().enumerate() {
        let ev = events.get(k).ok_or_else(|| anyhow::anyhow!("event"))?;
        let fit = ev.doppler.ok_or_else(|| anyhow::anyhow!("fit"))?;
        let tr = track_of(icao).ok_or_else(|| anyhow::anyhow!("{icao} not in archive"))?;
        eprintln!("  {} track points: {}", icao, fl.iter().find(|f| f.icao.to_lowercase() == *icao).map(|f| f.points.len()).unwrap_or(0));
        pair.push((tr, fit));
    }

    // Score one point: best over the shared-clock shift grid.
    let eval = |lat: f64, lon: f64| -> Option<(f64, f64)> {
        let obs = avalon::search::ecef(lat, lon, 0.0);
        let (_, t0) = pair[0].0.cpa(obs)?;
        let base = t0 - pair[0].1.t_cpa_s;
        let mut best = (f64::INFINITY, 0.0);
        for &sh in &score::shift_grid(10.0, 5) {
            let mut sum = 0.0;
            let mut ok = true;
            for (k, (tr, fit)) in pair.iter().enumerate() {
                let mut err = f64::INFINITY;
                for &f0 in &score::f0_grid(fit.f0_hz as f64, 5) {
                    if let Some(val) = score::curve_rmse_offset(
                        &curves[k], tr, obs, f0, base, sh, 0.0,
                        score::Atmos::default(), 400.0) {
                        if val < err { err = val; }
                    }
                }
                if !err.is_finite() { ok = false; break; }
                sum += err;
            }
            if ok && sum / 2.0 < best.0 { best = (sum / 2.0, sh); }
        }
        if best.0.is_finite() { Some(best) } else { None }
    };

    let step_deg = coarse_km / 111.32;
    let mut hits: Vec<(f64, f64, f64)> = Vec::new(); // rmse, lat, lon
    let mut lat = s;
    let mut done = 0;
    while lat <= n {
        let mut lon = w;
        while lon <= e {
            if let Some((r, _)) = eval(lat, lon) { hits.push((r, lat, lon)); }
            lon += step_deg;
            done += 1;
        }
        lat += step_deg * 0.77;
    }
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    eprintln!("{done} coarse points, best rmse {:.2} Hz", hits[0].0);

    // refine the 8 best coarse cells at 1/8 the step
    let mut refined: Vec<(f64, f64, f64)> = Vec::new();
    for &(_, la, lo) in hits.iter().take(8) {
        let mut i = -1.0;
        while i <= 1.0 {
            let mut j = -1.0;
            while j <= 1.0 {
                let (cla, clo) = (la + i * step_deg * 0.5, lo + j * step_deg * 0.5);
                if let Some((r, _)) = eval(cla, clo) { refined.push((r, cla, clo)); }
                j += 0.25;
            }
            i += 0.25;
        }
    }
    refined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    println!("\nknown pair {} -> best-fit locations (truth {:.2},{:.2}):", known.join(","), truth.0, truth.1);
    println!("{:<6} {:<10} {:<10} {:>8} {:>10}", "rank", "lat", "lon", "rmse Hz", "km from truth");
    for (k, (r, la, lo)) in refined.iter().take(10).enumerate() {
        println!("{:<6} {:<10.4} {:<10.4} {:>8.2} {:>10.0}", k + 1, la, lo, r, haversine((*la, *lo), truth) / 1000.0);
    }
    let best = refined[0];
    println!("\nCLOSEST HIGH-CONFIDENCE SPIKE: {:.0} km from truth, rmse {:.2} Hz", haversine((best.1, best.2), truth) / 1000.0, best.0);
    Ok(())
}
