//! Does a finer clock-shift search alone fix the location?
//!
//! The production scorer evaluates 5 clock shifts (±10 s, 5 s steps) and 5 f0
//! values per candidate. This reproduces that scoring exactly (per-candidate
//! clock anchor from the candidate's own CPA time) but swaps in a dense shift
//! sweep and a finer f0 grid, then refines the winners locally. If the winner
//! jumps to the truth, the 118 km answer was pure grid discretisation.
//!
//! Usage: fit4 [W,S,E,N] [ICAO0,ICAO1] [coarse_km]

use std::sync::atomic::AtomicBool;

use avalon::model::{Bounds, DetectedEvent, Ridge};
use avalon::{flights, score};
use rayon::prelude::*;

fn haversine_km(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6371.0f64;
    let (la1, la2) = (a.0.to_radians(), b.0.to_radians());
    let h = ((b.0 - a.0).to_radians() * 0.5).sin().powi(2)
        + la1.cos() * la2.cos() * ((b.1 - a.1).to_radians() * 0.5).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

struct P {
    curves: Vec<Vec<(f64, f32)>>,
    tracks: Vec<score::Track>,
    f0: [f64; 2],
    t_cpa: [f64; 2],
    atm: score::Atmos,
}

impl P {
    /// Production-style score: anchor follows the candidate location.
    fn score(&self, lat: f64, lon: f64, shifts: &[f64], f0_steps: usize) -> Option<(f64, f64)> {
        let obs = avalon::search::ecef(lat, lon, 0.0);
        let (_, t0) = self.tracks[0].cpa(obs)?;
        let base = t0 - self.t_cpa[0];
        let mut best = (f64::INFINITY, 0.0);
        for &sh in shifts {
            let mut sum = 0.0;
            let mut ok = true;
            for k in 0..2 {
                let mut err = f64::INFINITY;
                for f0 in score::f0_grid(self.f0[k], f0_steps) {
                    if let Some(v) = score::curve_rmse_offset(
                        &self.curves[k], &self.tracks[k], obs, f0, base, sh, 0.0,
                        self.atm, 400.0)
                    {
                        if v < err { err = v; }
                    }
                }
                if !err.is_finite() { ok = false; break; }
                sum += err;
            }
            if ok && sum / 2.0 < best.0 { best = (sum / 2.0, sh); }
        }
        if best.0.is_finite() { Some(best) } else { None }
    }
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (w, s, e, n) = match args.first() {
        Some(v) if !v.is_empty() => {
            let p: Vec<f64> = v.split(',').filter_map(|x| x.parse().ok()).collect();
            if p.len() == 4 { (p[0], p[1], p[2], p[3]) } else { (-124.8, 24.5, -66.9, 49.4) }
        }
        _ => (-124.8, 24.5, -66.9, 49.4),
    };
    let known: Vec<String> = args.get(1)
        .map(|x| x.split(',').map(|y| y.trim().to_lowercase()).collect())
        .unwrap_or_else(|| vec!["a67ba8".into(), "a659a2".into()]);
    let coarse_km: f64 = args.get(2).and_then(|v| v.parse().ok()).unwrap_or(60.0);
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
    let mut prog = |_q: flights::FlightProgress| {};
    let fl = flights::fetch_flights(win.0, win.1, &[region], "staging-0", &mut prog, &cancel)?;

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
    let mut tracks = Vec::new();
    let mut f0s = [0.0; 2];
    let mut tcpa = [0.0; 2];
    for (k, icao) in known.iter().enumerate() {
        let ev = events.get(k).ok_or_else(|| anyhow::anyhow!("event"))?;
        let fit = ev.doppler.ok_or_else(|| anyhow::anyhow!("fit"))?;
        f0s[k] = fit.f0_hz as f64;
        tcpa[k] = fit.t_cpa_s;
        tracks.push(track_of(icao).ok_or_else(|| anyhow::anyhow!("{icao} missing"))?);
    }
    let probe = P { curves, tracks, f0: f0s, t_cpa: tcpa, atm: score::Atmos::default() };

    let dense: Vec<f64> = score::shift_grid(60.0, 61);   // 1 s steps
    let coarse = score::shift_grid(10.0, 5);             // production grid

    println!("\n=== production grid (5 shifts) vs dense (1 s) at the truth ===");
    for (label, grid, f0n) in [("production", &coarse, 5usize), ("dense", &dense, 9)] {
        let r = probe.score(truth.0, truth.1, grid, f0n);
        println!("  {:<11} rmse {:>7.2} Hz shift {:+.1} s", label,
                 r.map(|x| x.0).unwrap_or(f64::NAN), r.map(|x| x.1).unwrap_or(f64::NAN));
    }

    // AOS_LOCAL=lat,lon,km : dense scoring in a window around a seed, skipping CONUS
    if let Ok(spec) = std::env::var("AOS_LOCAL") {
        let p: Vec<f64> = spec.split(',').filter_map(|x| x.parse().ok()).collect();
        let (clat, clon, span) = (p[0], p[1], p[2]);
        let st = 15.0 / 111.32;
        let mut cells: Vec<(f64, f64)> = Vec::new();
        let mut la = clat - span / 111.0;
        while la <= clat + span / 111.0 {
            let mut lo = clon - span / 111.0 * 1.3;
            while lo <= clon + span / 111.0 * 1.3 {
                cells.push((la, lo));
                lo += st;
            }
            la += st;
        }
        eprintln!("{} local cells", cells.len());
        let mut sc: Vec<(f64, f64, f64, f64)> = cells.par_iter()
            .filter_map(|&(la, lo)| probe.score(la, lo, &dense, 9).map(|(r, sh)| (r, la, lo, sh)))
            .collect();
        sc.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        println!("\n=== LOCAL SWEEP around {:.2},{:.2} +/-{:.0} km (15 km cells, dense shift) ===", clat, clon, span);
        let mut keep: Vec<(f64, f64, f64, f64)> = Vec::new();
        for c in sc.into_iter() {
            if keep.iter().any(|k| haversine_km((k.1, k.2), (c.1, c.2)) < 30.0) { continue; }
            keep.push(c);
            if keep.len() >= 5 { break; }
        }
        for (i, (r, la, lo, sh)) in keep.iter().enumerate() {
            println!("  #{} {:.3},{:.3} rmse {:.2} Hz shift {:+.0}s  {:.0} km from truth",
                     i + 1, la, lo, r, sh, haversine_km((*la, *lo), truth));
        }
        let (br, bla, blo, bsh) = keep[0];
        let fs = 3.0 / 111.32;
        let mut loc: Vec<(f64, f64, f64)> = Vec::new();
        let mut la = bla - 0.25;
        while la <= bla + 0.25 {
            let mut lo = blo - 0.32;
            while lo <= blo + 0.32 {
                if let Some((r, _)) = probe.score(la, lo, &[bsh - 3.0, bsh - 1.5, bsh, bsh + 1.5, bsh + 3.0], 11) {
                    loc.push((r, la, lo));
                }
                lo += fs;
            }
            la += fs;
        }
        loc.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        println!("\n  refined best: {:.4},{:.4} rmse {:.2} Hz ({:.1} km from truth)",
                 loc[0].1, loc[0].2, loc[0].0, haversine_km((loc[0].1, loc[0].2), truth));
        return Ok(());
    }

    let step = coarse_km / 111.32;
    let mut cells: Vec<(f64, f64)> = Vec::new();
    let mut lat = s;
    while lat <= n {
        let mut lon = w;
        while lon <= e {
            cells.push((lat, lon));
            lon += step;
        }
        lat += step * 0.77;
    }
    eprintln!("{} coarse cells", cells.len());

    // pass 1: production grid (cheap) to shortlist
    let mut scored: Vec<(f64, f64, f64)> = cells.par_iter()
        .filter_map(|&(la, lo)| probe.score(la, lo, &coarse, 5).map(|(r, _)| (r, la, lo)))
        .collect();
    scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut shortlist: Vec<(f64, f64)> = Vec::new();
    for &(_, la, lo) in scored.iter() {
        if shortlist.iter().any(|p| haversine_km(*p, (la, lo)) < 60.0) { continue; }
        shortlist.push((la, lo));
        if shortlist.len() >= 25 { break; }
    }
    println!("\n=== pass 1 (production grid, {} cells): best {:.1} Hz at {:.0} km from truth ===",
             scored.len(), scored[0].0, haversine_km((scored[0].1, scored[0].2), truth));

    // pass 2: dense shift + finer f0 on the shortlist
    let mut refined: Vec<(f64, f64, f64, f64)> = shortlist.par_iter()
        .filter_map(|&(la, lo)| probe.score(la, lo, &dense, 9).map(|(r, sh)| (r, la, lo, sh)))
        .collect();
    refined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    println!("\n=== pass 2 (dense shift) ===");
    for (i, (r, la, lo, sh)) in refined.iter().take(6).enumerate() {
        println!("  #{} {:.3},{:.3} rmse {:.2} Hz shift {:+.0}s  {:.0} km from truth",
                 i + 1, la, lo, r, sh, haversine_km((*la, *lo), truth));
    }

    // pass 3: local 3 km refinement around the new best
    let (_, bla, blo, bsh) = refined[0];
    let fs = 3.0 / 111.32;
    let mut local: Vec<(f64, f64, f64)> = Vec::new();
    let mut la = bla - 0.5;
    while la <= bla + 0.5 {
        let mut lo = blo - 0.6;
        while lo <= blo + 0.6 {
            if let Some((r, _)) = probe.score(la, lo, &[bsh - 2.0, bsh, bsh + 2.0], 9) {
                local.push((r, la, lo));
            }
            lo += fs;
        }
        la += fs;
    }
    local.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    println!("\n=== pass 3 (local 3 km refine) ===");
    for (i, (r, la, lo)) in local.iter().take(4).enumerate() {
        println!("  #{} {:.4},{:.4} rmse {:.2} Hz  {:.1} km from truth",
                 i + 1, la, lo, r, haversine_km((*la, *lo), truth));
    }
    Ok(())
}
