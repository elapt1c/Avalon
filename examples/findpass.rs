//! Which aircraft actually passed near a known point at a known time?
//!
//! Downloads (and caches) the archive for the window in AVALON_WIN, then reports
//! every flight whose closest approach to the reference point is within
//! `--max-km`, with the slant distance and the UTC time of that approach.
//!
//! Usage: findpass [W,S,E,N] [max_km]
//!   AVALON_WIN=unix_start,unix_end
//!   AVALON_TRUTH=lat,lon

use std::sync::atomic::AtomicBool;

use avalon::model::Bounds;
use avalon::{flights, search};

fn haversine_km(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6371.0f64;
    let (la1, la2) = (a.0.to_radians(), b.0.to_radians());
    let h = ((b.0 - a.0).to_radians() * 0.5).sin().powi(2)
        + la1.cos() * la2.cos() * ((b.1 - a.1).to_radians() * 0.5).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

fn iso(t: f64) -> String {
    let secs = t as i64;
    let days = secs.div_euclid(86_400);
    let sod = secs.rem_euclid(86_400);
    // civil_from_days (Howard Hinnant)
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!("{:04}-{:02}-{:02} {:02}:{:02}:{:02}Z", y, m, d, sod / 3600, (sod % 3600) / 60, sod % 60)
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (w, s, e, n) = match args.first() {
        Some(v) if !v.is_empty() => {
            let p: Vec<f64> = v.split(',').filter_map(|x| x.parse().ok()).collect();
            if p.len() == 4 { (p[0], p[1], p[2], p[3]) } else { (-74.6, 39.8, -71.8, 42.0) }
        }
        _ => (-74.6, 39.8, -71.8, 42.0),
    };
    let max_km: f64 = args.get(1).and_then(|v| v.parse().ok()).unwrap_or(25.0);

    let win_spec = std::env::var("AVALON_WIN").unwrap_or_default();
    let win: Vec<f64> = win_spec.split(',').filter_map(|x| x.trim().parse().ok()).collect();
    let (t0, t1) = if win.len() == 2 { (win[0] as i64, win[1] as i64) } else { (0, 0) };
    anyhow::ensure!(t0 > 0 && t1 > t0, "AVALON_WIN=unix_start,unix_end required");

    let truth_spec = std::env::var("AVALON_TRUTH").unwrap_or_else(|_| "40.876466,-73.524951".into());
    let tp: Vec<f64> = truth_spec.split(',').filter_map(|x| x.parse().ok()).collect();
    let truth = (tp[0], tp[1]);

    println!("window {} .. {}", iso(t0 as f64), iso(t1 as f64));
    println!("reference point {:.6},{:.6}, max slant {:.0} km", truth.0, truth.1, max_km);

    let region = Bounds { west: w, south: s, east: e, north: n };
    let cancel = AtomicBool::new(false);
    let mut prog = |_p: flights::FlightProgress| {};
    let fl = flights::fetch_flights(t0, t1, &[region], "staging-0", &mut prog, &cancel)?;
    println!("{} flights in box\n", fl.len());

    let obs = search::ecef(truth.0, truth.1, 0.0);
    let mut rows: Vec<(f64, f64, String, usize, f64, f64)> = Vec::new(); // dist, t, icao, pts, alt, cpa_t
    for f in fl.iter() {
        let mut best = (f64::INFINITY, 0.0f64, 0.0f64);
        for p in f.points.iter() {
            let d = haversine_km(truth, (p.lat, p.lon));
            if d < best.0 {
                best = (d, p.t_s, p.alt_m);
            }
        }
        if best.0 <= max_km {
            rows.push((best.0, best.1, f.icao.clone(), f.points.len(), best.2, best.1));
        }
    }
    rows.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    println!("{:<10} {:<22} {:>8} {:>7} {:>8}", "icao", "closest-approach (UTC)", "km", "alt m", "pts");
    for (d, t, icao, pts, alt, _) in rows.iter() {
        println!("{:<10} {:<22} {:>8.2} {:>7.0} {:>8}", icao, iso(*t), d, alt, pts);
    }
    let _ = obs;
    Ok(())
}
