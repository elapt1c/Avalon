//! Is the ground truth even findable? Check the two-pass structure at the true
//! receiver location against the cached flight data.
//!
//! Usage: paircheck <clip> [--region W,S,E,N] [--window A,B] [--truth LAT,LON]

use std::path::Path;
use std::sync::atomic::AtomicBool;

use avalon::model::{AnalysisSettings, Bounds, Metadata, UtcWindow};
use avalon::{flights, pipeline, score};

const DEFAULT_CLIP: &str =
    "/home/tpg/Downloads/Unknown location - raw audio - colsto (128p).mp4";
const TRUTH: (f64, f64) = (40.87, -73.52);

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let clip = args.first().cloned().unwrap_or_else(|| DEFAULT_CLIP.to_string());
    let mut region = Bounds {
        west: TRUTH.1 - 0.9555,
        south: TRUTH.0 - 0.7226,
        east: TRUTH.1 + 0.9555,
        north: TRUTH.0 + 0.7226,
    };
    // Must match the GUI defaults exactly or the parsed cache key differs.
    let mut window = (
        avalon::dateutil::civil_to_unix(2026, 10, 30, 10, 0, 0),
        avalon::dateutil::civil_to_unix(2026, 11, 1, 12, 0, 0),
    );
    let mut truth = TRUTH;
    let mut i = 1;
    while i + 1 < args.len() {
        match args[i].as_str() {
            "--region" => {
                let v: Vec<f64> = args[i + 1].split(',').filter_map(|s| s.parse().ok()).collect();
                if v.len() == 4 {
                    region = Bounds { west: v[0], south: v[1], east: v[2], north: v[3] };
                }
            }
            "--window" => {
                let v: Vec<i64> = args[i + 1].split(',').filter_map(|s| s.parse().ok()).collect();
                if v.len() == 2 {
                    window = (v[0], v[1]);
                }
            }
            "--truth" => {
                let v: Vec<f64> = args[i + 1].split(',').filter_map(|s| s.parse().ok()).collect();
                if v.len() == 2 {
                    truth = (v[0], v[1]);
                }
            }
            _ => {}
        }
        i += 2;
    }

    // 1. Analyse the clip for its two flyovers.
    println!("analysing {clip}");
    let settings = AnalysisSettings { target_sample_rate: 4000, ..Default::default() };
    let meta = Metadata { region: vec![region], time_window: Some(UtcWindow { start_unix_s: window.0, end_unix_s: window.1 }), ..Default::default() };
    let a = pipeline::run(Path::new(&clip), &settings, &meta, None)?;
    let mut specs: Vec<(f64, f64, f64)> = Vec::new(); // (cpa_s, speed, slant)
    for e in &a.events {
        if e.id == 0 {
            continue;
        }
        let rep = a.ridges.iter().find(|r| r.id == e.representative_ridge_id);
        if rep.map(|r| r.plausible).unwrap_or(false) {
            if let Some(f) = e.doppler {
                specs.push((f.t_cpa_s, f.speed_mps as f64, f.slant_distance_m as f64));
            }
        }
    }
    println!("plausible events: {specs:?}");

    // 2. Flights for the window (cached).
    let cancel = AtomicBool::new(false);
    let mut prog = |_p: flights::FlightProgress| {};
    let fl = flights::fetch_flights(window.0, window.1, &[region], "staging-0", &mut prog, &cancel)?;
    println!("{} flights", fl.len());

    // 3. At the truth point, list every pass and look for the two-pass structure.
    let obs = avalon::search::ecef(truth.0, truth.1, 0.0);
    let mut passes: Vec<(f64, f64, f64, String)> = Vec::new(); // (recv, dist, speed, icao)
    for f in &fl {
        // Subsample to 10 s: a full day at 1 s is ~86k points per flight.
        let mut pts: Vec<(f64, [f64; 3])> = Vec::new();
        let mut last = f64::NEG_INFINITY;
        for p in &f.points {
            if p.t_s - last < 10.0 {
                continue;
            }
            last = p.t_s;
            pts.push((p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m)));
            if pts.len() > 20_000 {
                break;
            }
        }
        let Some(tr) = score::Track::from_points(&pts) else { continue };
        let Some((d, t)) = tr.cpa(obs) else { continue };
        if d > 25_000.0 {
            continue;
        }
        let v = tr.vel(t).map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()).unwrap_or(0.0);
        passes.push((t + d / score::C0, d, v, f.icao.clone()));
    }
    passes.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    println!("\npasses within 25 km of truth: {}", passes.len());
    for p in passes.iter().take(12) {
        println!("  t={:.0} d={:.1} km v={:.0} m/s {}", p.0, p.1 / 1000.0, p.2, p.3);
    }

    if specs.len() < 2 {
        println!("\nneed two plausible events to test pairing");
        return Ok(());
    }
    // Match the GUI's autorun: event 0 is skipped, so pair the last two.
    let (e0, e1) = if specs.len() >= 3 {
        (specs[1], specs[2])
    } else {
        (specs[0], specs[1])
    };
    let dt_audio = (e1.0 - e0.0).abs();
    println!(
        "\naudio Δt = {dt_audio:.0}s · speeds {:.0}/{:.0} m/s · slants {:.0}/{:.0} m",
        e0.1, e1.1, e0.2, e1.2
    );
    println!("(pairs below are filtered to BOTH speeds within 0.75-1.35x)");
    let mut hits = 0;
    let mut pair_match: Option<(String, String, f64)> = None;
    for (i, a) in passes.iter().enumerate() {
        for b in passes.iter().skip(i + 1) {
            let dt = (b.0 - a.0).abs();
            let speed_ok = |v: f64, target: f64| v >= 0.75 * target && v <= 1.35 * target;
            let paired = (speed_ok(a.2, e0.1) && speed_ok(b.2, e1.1))
                || (speed_ok(a.2, e1.1) && speed_ok(b.2, e0.1));
            if (dt - dt_audio).abs() <= 10.0 && paired {
                hits += 1;
                if pair_match.is_none() {
                    pair_match = Some((a.3.clone(), b.3.clone(), a.0 - e0.0));
                }
                println!(
                    "  PAIR Δt={dt:.0}s (err {:+.0}s) {}+{} d={:.1}/{:.1} km v={:.0}/{:.0} m/s",
                    dt - dt_audio,
                    a.3,
                    b.3,
                    a.1 / 1000.0,
                    b.1 / 1000.0,
                    a.2,
                    b.2
                );
            }
        }
    }
    println!("\n{hits} pair(s) at truth match the audio Δt within 10 s");

    // Residual comparison: the same curve-agreement metric the search ranks on.
    let mut curves: Vec<Vec<(f64, f32)>> = Vec::new();
    let mut fits: Vec<(f64, f64, f64, f64)> = Vec::new(); // (cpa, speed, slant, f0)
    for e in &a.events {
        if e.id == 0 {
            continue;
        }
        let rep = a.ridges.iter().find(|r| r.id == e.representative_ridge_id);
        if rep.map(|r| r.plausible).unwrap_or(false) {
            if let (Some(f), Some(r)) = (e.doppler, rep) {
                fits.push((
                    f.t_cpa_s,
                    f.speed_mps as f64,
                    f.slant_distance_m as f64,
                    f.f0_hz as f64,
                ));
                curves.push(r.samples.iter().map(|s| (s.t_s, s.freq_hz)).collect());
            }
        }
    }
    if curves.len() < 2 {
        return Ok(());
    }
    let (ci, cj) = if curves.len() >= 3 { (1, 2) } else { (0, 1) };

    let mut best_at = |lat: f64, lon: f64, label: &str| {
        let o = avalon::search::ecef(lat, lon, 0.0);
        let mut ranked: Vec<(f64, f64, String)> = Vec::new(); // (total rmse, shift, icao)
        for f in &fl {
            let mut pts: Vec<(f64, [f64; 3])> = Vec::new();
            let mut last = f64::NEG_INFINITY;
            for p in &f.points {
                if p.t_s - last < 5.0 { continue; }
                last = p.t_s;
                pts.push((p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m)));
            }
            let Some(tr) = score::Track::from_points(&pts) else { continue };
            let Some((d, t_cpa)) = tr.cpa(o) else { continue };
            if d > 30_000.0 { continue; }
            let base = t_cpa - fits[ci].0;
            let mut best0 = f64::INFINITY;
            for &f0 in &score::f0_grid(fits[ci].3.max(20.0), 9) {
                for &sh in &score::shift_grid(10.0, 9) {
                    for &df in &[-3.0f64, -1.5, 0.0, 1.5, 3.0] {
                        if let Some(v) = score::curve_rmse_offset(
                            &curves[ci], &tr, o, f0, base, sh, df, score::Atmos::default(), 25.0,
                        ) {
                            if v < best0 { best0 = v; }
                        }
                    }
                }
            }
            if best0.is_finite() {
                ranked.push((best0, d, f.icao.clone()));
            }
        }
        ranked.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        println!("\n{label} ({lat:.4},{lon:.4}) best single-event residuals:");
        for r in ranked.iter().take(5) {
            println!("  rmse {:.2} Hz  d={:.1} km  {}", r.0, r.1 / 1000.0, r.2);
        }
    };
    best_at(truth.0, truth.1, "TRUTH");
    best_at(41.1334, -73.6557, "BEST-CELL");

    // Joint pair residuals at the truth: the two events must be explained by two
    // flights that share one recording start, which is what score_cell gates on.
    {
        let o = avalon::search::ecef(truth.0, truth.1, 0.0);
        let track_of = |f: &avalon::model::Flight| {
            let mut pts: Vec<(f64, [f64; 3])> = Vec::new();
            let mut last = f64::NEG_INFINITY;
            for p in &f.points {
                if p.t_s - last < 5.0 {
                    continue;
                }
                last = p.t_s;
                pts.push((p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m)));
            }
            score::Track::from_points(&pts)
        };
        // Per event: best residual per flight, with its own geometry.
        let mut per: Vec<Vec<(f64, f64, String, f64, f64)>> = Vec::new(); // rmse, d, icao, cpa, base
        for k in 0..2 {
            let mut rows = Vec::new();
            for f in &fl {
                let Some(tr) = track_of(f) else { continue };
                let Some((d, t_cpa)) = tr.cpa(o) else { continue };
                if d > 30_000.0 {
                    continue;
                }
                let base = t_cpa - fits[k].0;
                if base < window.0 as f64 - 120.0 || base > window.1 as f64 + 120.0 {
                    continue;
                }
                let mut best = f64::INFINITY;
                for &f0 in &score::f0_grid(fits[k].3.max(20.0), 9) {
                    for &sh in &score::shift_grid(10.0, 9) {
                        for &df in &[-3.0f64, 0.0, 3.0] {
                            if let Some(v) = score::curve_rmse_offset(
                                &curves[k], &tr, o, f0, base, sh, df,
                                score::Atmos::default(), 25.0,
                            ) {
                                if v < best { best = v; }
                            }
                        }
                    }
                }
                if best.is_finite() {
                    rows.push((best, d, f.icao.clone(), t_cpa, base));
                }
            }
            rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            per.push(rows);
        }
        println!("\nper-event top residuals at truth (event0 / event1):");
        for (k, rows) in per.iter().enumerate() {
            for r in rows.iter().take(4) {
                println!(
                    "  e{k} rmse {:.2} Hz d={:.1} km v? base {:+.1}s {}",
                    r.0,
                    r.1 / 1000.0,
                    r.4,
                    r.2
                );
            }
        }
        // Pairs sharing a recording start (within the search's tolerance).
        let mut pairs: Vec<(f64, usize, usize)> = Vec::new();
        for (i, a) in per[0].iter().enumerate().take(60) {
            for (j, b) in per[1].iter().enumerate().take(60) {
                if a.2 == b.2 {
                    continue;
                }
                if (a.4 - b.4).abs() > 12.0 {
                    continue;
                }
                pairs.push((0.5 * (a.0 + b.0), i, j));
            }
        }
        pairs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
        println!("\njoint pair residuals at truth (shared base, {} pairs):", pairs.len());
        for (mean, i, j) in pairs.iter().take(6) {
            let (a, b) = (&per[0][*i], &per[1][*j]);
            println!(
                "  mean {:.2} Hz | e0 {:.2} {} d{:.1}km | e1 {:.2} {} d{:.1}km | Δbase {:+.1}s",
                mean,
                a.0,
                a.2,
                a.1 / 1000.0,
                b.0,
                b.2,
                b.1 / 1000.0,
                b.4 - a.4
            );
        }
    }

    // Curve overlay: measured ridge vs the simulated curve from the pair that
    // matches the audio's timing and speeds at the true location.
    if let Some((i0, i1, base0)) = pair_match {
        let o = avalon::search::ecef(truth.0, truth.1, 0.0);
        for (k, &(ti, t_abs)) in [(ci, fits[ci].0), (cj, fits[cj].0)].iter().enumerate() {
            let _ = t_abs;
            let icao = if k == 0 { &i0 } else { &i1 };
            let Some(f) = fl.iter().find(|f| &f.icao == icao) else { continue };
            let mut pts: Vec<(f64, [f64; 3])> = Vec::new();
            let mut last = f64::NEG_INFINITY;
            for p in &f.points {
                if p.t_s - last < 5.0 { continue; }
                last = p.t_s;
                pts.push((p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m)));
            }
            let Some(tr) = score::Track::from_points(&pts) else { continue };
            let idx = if k == 0 { ci } else { cj };
            let base = base0;
            println!("\nEVENT {idx} measured vs simulated ({icao}, base {base:.0}):");
            println!("     t_rel    measured   simulated   delta");
            let mut tau = 0.0;
            let n = curves[idx].len();
            for j in (0..n).step_by((n / 10).max(1)) {
                let (t_rel, f_meas) = curves[idx][j];
                let model_t = base + t_rel;
                if j == 0 {
                    tau = model_t - 2000.0 / score::C0;
                }
                let sim = score::doppler_ratio_warm(&tr, o, model_t, &mut tau, score::Atmos::default())
                    .map(|r| fits[idx].3 * r);
                match sim {
                    Some(s) => println!(
                        "  {:8.1}  {:8.1}   {:8.1}   {:+.1}",
                        t_rel,
                        f_meas,
                        s,
                        s - f_meas as f64
                    ),
                    None => println!("  {:8.1}  {:8.1}        --", t_rel, f_meas),
                }
            }
        }
    }
    Ok(())
}
