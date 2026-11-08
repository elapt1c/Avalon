//! At the true receiver, how well can the known pair explain each measured
//! curve if the clock offset is free? Distinguishes "wrong clock" from "wrong
//! shape", which is the difference between a fixable bug and a bad trace.
//!
//! Usage: truthpair [a67ba8] [a659a2]

use std::path::Path;
use std::sync::atomic::AtomicBool;

use avalon::model::{AnalysisSettings, Bounds, Metadata, UtcWindow};
use avalon::{flights, pipeline, score};

const CLIP: &str = "/home/tpg/Downloads/Unknown location - raw audio - colsto (128p).mp4";
const TRUTH: (f64, f64) = (40.87, -73.52);

fn main() -> anyhow::Result<()> {
    let a = std::env::args().nth(1).unwrap_or_else(|| "a67ba8".into());
    let b = std::env::args().nth(2).unwrap_or_else(|| "a659a2".into());
    let region = Bounds {
        west: TRUTH.1 - 0.9555,
        south: TRUTH.0 - 0.7226,
        east: TRUTH.1 + 0.9555,
        north: TRUTH.0 + 0.7226,
    };
    // Three days: the pair may sit just outside the search's one-day window, and
    // that is exactly what this check is for.
    let window = (
        avalon::dateutil::civil_to_unix(2026, 10, 30, 10, 0, 0),
        avalon::dateutil::civil_to_unix(2026, 11, 1, 12, 0, 0),
    );
    let search_window = (
        avalon::dateutil::civil_to_unix(2026, 10, 31, 11, 0, 0),
        avalon::dateutil::civil_to_unix(2026, 10, 31, 23, 0, 0),
    );
    let settings = AnalysisSettings { target_sample_rate: 4000, ..Default::default() };
    let meta = Metadata {
        region: vec![region],
        time_window: Some(UtcWindow { start_unix_s: search_window.0, end_unix_s: search_window.1 }),
        ..Default::default()
    };
    println!("analysing");
    let an = pipeline::run(Path::new(CLIP), &settings, &meta, None)?;

    let mut events = Vec::new();
    for e in &an.events {
        if e.id == 0 {
            continue;
        }
        let Some(rep) = an.ridges.iter().find(|r| r.id == e.representative_ridge_id) else { continue };
        if !rep.plausible {
            continue;
        }
        events.push((e.clone(), rep.clone()));
    }
    println!("{} selected events", events.len());
    for (e, r) in &events {
        println!(
            "  event {} span {:.0}-{:.0}s samples {} fit {:?}",
            e.id,
            e.start_s,
            e.end_s,
            r.samples.len(),
            e.doppler.map(|f| (f.f0_hz, f.speed_mps as i32, (f.slant_distance_m / 1000.0 * 10.0) as i32 as f32 / 10.0, f.rmse_hz))
        );
    }

    let cancel = AtomicBool::new(false);
    let mut prog = |_p: flights::FlightProgress| {};
    let fl = flights::fetch_flights(window.0, window.1, &[region], "staging-0", &mut prog, &cancel)?;
    let obs = avalon::search::ecef(TRUTH.0, TRUTH.1, 0.0);
    let track_of = |icao: &str| {
        let f = fl.iter().find(|f| f.icao == icao)?;
        let mut pts: Vec<(f64, [f64; 3])> = Vec::new();
        let mut last = f64::NEG_INFINITY;
        for p in &f.points {
            if p.t_s - last < 1.0 {
                continue;
            }
            last = p.t_s;
            pts.push((p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m)));
        }
        score::Track::from_points(&pts)
    };
    for (icao, ei) in [(&a, 0usize), (&b, 1usize)] {
        let Some((e, r)) = events.get(ei) else { continue };
        let Some(fit) = e.doppler else { continue };
        let Some(tr) = track_of(icao) else {
            println!("{icao}: not in archive");
            continue;
        };
        let (d, t_cpa) = tr.cpa(obs).unwrap();
        let v = tr.vel(t_cpa).map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()).unwrap_or(0.0);
        // Absolute time of the recording's t=0 if this pass is the event's pass.
        let base = t_cpa - fit.t_cpa_s;
        let (dy, dm, dd, hh, mi, ss) = avalon::dateutil::unix_to_civil(base as i64);
        let (cy, cm, cd, ch, cmi, _) = avalon::dateutil::unix_to_civil(t_cpa as i64);
        println!(
            "\n{icao} event {}: d {:.1} km v {:.0} m/s",
            e.id,
            d / 1000.0,
            v
        );
        println!(
            "  pass at {cy}-{cm:02}-{cd:02} {ch:02}:{cmi:02}Z; implied recording start {dy}-{dm:02}-{dd:02} {hh:02}:{mi:02}:{ss:02}Z"
        );
        println!(
            "  search window {}-{}-{}..{}:{}Z  -> start {}the window",
            "2026",
            "10",
            "31",
            "23",
            "00",
            if base < window.0 as f64 || base > window.1 as f64 { "OUTSIDE " } else { "inside " }
        );
        let curve: Vec<(f64, f32)> = r.samples.iter().map(|s| (s.t_s, s.freq_hz)).collect();
        // Fine clock sweep, wide: is the mismatch a clock or a shape?
        let mut rows: Vec<(f64, f64, f64)> = Vec::new(); // shift, rmse, f0
        for k in 0..=240 {
            let sh = -120.0 + k as f64;
            let mut best = f64::INFINITY;
            let mut bf0 = 0.0;
            for &f0 in &score::f0_grid(fit.f0_hz as f64, 17) {
                if let Some(v) = score::curve_rmse(
                    &curve, &tr, obs, f0, base, sh, score::Atmos::default(), 40.0,
                ) {
                    if v < best {
                        best = v;
                        bf0 = f0;
                    }
                }
            }
            if best.is_finite() {
                rows.push((sh, best, bf0));
            }
        }
        rows.sort_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal));
        println!("  best clock shifts:");
        for r in rows.iter().take(6) {
            println!("    shift {:+.0}s rmse {:.2} Hz f0 {:.0}", r.0, r.1, r.2);
        }
        if let Some(best) = rows.first() {
            let (sh, _, f0) = *best;
            println!("  curve at shift {:+.0}s, f0 {:.0}:", sh, f0);
            println!("       t_rel  measured  simulated  delta");
            let mut tau = 0.0;
            let n = curve.len();
            for j in (0..n).step_by((n / 12).max(1)) {
                let (t_rel, f_meas) = curve[j];
                let model_t = base + t_rel - sh;
                if j == 0 {
                    tau = model_t - 2000.0 / score::C0;
                }
                let sim = score::doppler_ratio_warm(&tr, obs, model_t, &mut tau, score::Atmos::default())
                    .map(|r| f0 * r);
                match sim {
                    Some(s) => println!("    {:8.1} {:8.1}  {:8.1}  {:+.1}", t_rel, f_meas, s, s - f_meas as f64),
                    None => println!("    {:8.1} {:8.1}       --", t_rel, f_meas),
                }
            }
        }
    }
    Ok(())
}
