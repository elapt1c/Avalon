//! Which flights does the per-cell shortlist keep at one ground point?
//!
//! The locate stage ranks every cell by the best *flight pair* that explains
//! both measured Doppler curves. If the true pair never survives the shortlist
//! at the true receiver, no amount of grid refinement can find it. This prints
//! the shortlist exactly as `score_cell` builds it.
//!
//! Usage: pairdiag [lat,lon] [--pair ICAO1,ICAO2]

use std::sync::atomic::AtomicBool;

use avalon::model::{AnalysisSettings, Bounds, DetectedEvent, Metadata, Ridge, UtcWindow};
use avalon::{dateutil, flights, pipeline, score};

const CLIP: &str = "/home/tpg/Downloads/Unknown location - raw audio - colsto (128p).mp4";
const TRUTH: (f64, f64) = (40.87, -73.52);
/// Flights the project's own truth tools name as the pair behind the footage.
const KNOWN_PAIR: (&str, &str) = ("a67ba8", "a659a2");

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut here = TRUTH;
    let mut pair = KNOWN_PAIR;
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--pair" {
            i += 1;
            if let Some(p) = args.get(i) {
                let (a, b) = p.split_once(',').unwrap_or((p.as_str(), ""));
                pair = (
                    Box::leak(a.to_string().into_boxed_str()),
                    Box::leak(b.to_string().into_boxed_str()),
                );
            }
        } else if let Some((a, b)) = args[i].split_once(',') {
            here = (a.parse()?, b.parse()?);
        }
        i += 1;
    }
    // Must match the box the bench harness uses, or the parsed-flight cache
    // (keyed on region + window) is missed and the archive is re-parsed.
    let region = Bounds { west: -74.48, south: 40.15, east: -72.56, north: 41.59 };
    let win = (
        dateutil::civil_to_unix(2026, 10, 30, 11, 0, 0),
        dateutil::civil_to_unix(2026, 10, 31, 23, 0, 0),
    );
    let cache = "/root/projects/analysis-default.json";
    let (events, ridges): (Vec<DetectedEvent>, Vec<Ridge>) = match std::fs::read(cache) {
        Ok(b) => {
            let v: serde_json::Value = serde_json::from_slice(&b)?;
            (
                serde_json::from_value(v["events"].clone())?,
                serde_json::from_value(v["ridges"].clone())?,
            )
        }
        Err(_) => {
            let settings = AnalysisSettings::default();
            let meta = Metadata {
                region: vec![region],
                time_window: Some(UtcWindow { start_unix_s: win.0, end_unix_s: win.1 }),
                ..Default::default()
            };
            let a = pipeline::run(std::path::Path::new(CLIP), &settings, &meta, None)?;
            let mut ev = Vec::new();
            let mut rg = Vec::new();
            for e in &a.events {
                if e.id == 0 {
                    continue;
                }
                let Some(rep) = a.ridges.iter().find(|r| r.id == e.representative_ridge_id) else {
                    continue;
                };
                if !rep.plausible {
                    continue;
                }
                ev.push(e.clone());
                rg.push(rep.clone());
            }
            (ev, rg)
        }
    };
    let _ = ridges;

    let cancel = AtomicBool::new(false);
    let mut prog = |_p: flights::FlightProgress| {};
    let fl = flights::fetch_flights(win.0, win.1, &[region], "staging-0", &mut prog, &cancel)?;
    println!("{} flights in region/window", fl.len());

    let obs = avalon::search::ecef(here.0, here.1, 0.0);
    let truth_pair = (pair.0.to_lowercase(), pair.1.to_lowercase());
    for (k, e) in events.iter().enumerate().take(2) {
        let Some(fit) = e.doppler else { continue };
        println!(
            "\nevent {} (index {k}): f0 {:.1} Hz v {:.0} m/s slant {:.2} km cpa {:.1}s",
            e.id,
            fit.f0_hz,
            fit.speed_mps,
            fit.slant_distance_m / 1000.0,
            fit.t_cpa_s
        );
        let mut rows: Vec<(f64, f64, f64, f64, String, bool)> = Vec::new();
        for f in &fl {
            let mut pts: Vec<(f64, [f64; 3])> = Vec::new();
            for p in &f.points {
                pts.push((p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m)));
            }
            let Some(tr) = score::Track::from_points(&pts) else { continue };
            let Some((d, t_cpa)) = tr.cpa(obs) else { continue };
            if d > 40_000.0 {
                continue;
            }
            let gs = tr.vel(t_cpa).map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()).unwrap_or(0.0);
            if gs < 0.5 * fit.speed_mps as f64 || gs > 2.0 * fit.speed_mps as f64 {
                continue;
            }
            if d < 0.1 * fit.slant_distance_m as f64 || d > 8.0 * fit.slant_distance_m as f64 {
                continue;
            }
            let implied = t_cpa - fit.t_cpa_s;
            if implied < win.0 as f64 - 120.0 || implied > win.1 as f64 + 120.0 {
                continue;
            }
            let icao = f.icao.to_lowercase();
            let known = icao == truth_pair.0 || icao == truth_pair.1;
            rows.push((d, gs, implied, (gs - fit.speed_mps as f64).abs(), icao, known));
        }
        println!("  {} flights survive the geometric gates", rows.len());
        let known_rank = {
            let mut r = rows.clone();
            r.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));
            r.iter().position(|x| x.5)
        };
        println!("  speed-ranked: known pair entry rank {:?} (shortlist keeps 40)", known_rank);
        rows.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));
        for (d, gs, implied, _, icao, known) in rows.iter().take(12) {
            let (cy, cm, cd, ch, cmi, cs) = dateutil::unix_to_civil(*implied as i64);
            println!(
                "    {icao}{} d {:6.2} km v {:3.0} m/s implied start {cy}-{cm:02}-{cd:02} {ch:02}:{cmi:02}:{cs:02}Z",
                if *known { " <= KNOWN" } else { "        " },
                d / 1000.0,
                gs
            );
        }
        for (d, gs, implied, _, icao, known) in rows.iter().filter(|r| r.5) {
            let (cy, cm, cd, ch, cmi, cs) = dateutil::unix_to_civil(*implied as i64);
            println!(
                "  known {icao} d {:.2} km v {:.0} m/s implied start {cy}-{cm:02}-{cd:02} {ch:02}:{cmi:02}:{cs:02}Z rank {}",
                d / 1000.0,
                gs,
                rows.iter().position(|x| x.4 == *icao).unwrap_or(usize::MAX)
            );
            let _ = known;
        }
    }
    Ok(())
}
