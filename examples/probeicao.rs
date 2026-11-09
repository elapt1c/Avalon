//! What actually flew near the receiver at the measured event times?
//!
//! Usage: probeicao [truth_lat,truth_lon] [event_unix=+off_s,...]
//!
//! Prints, for each candidate event, every flight that came within 25 km of the
//! receiver while airborne, sorted by closest approach — the physical ground
//! truth the scoring has to find.
use avalon::flights;
use avalon::model::Bounds;

fn haversine(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6_371_000.0f64;
    let (la1, la2) = (a.0.to_radians(), b.0.to_radians());
    let dla = (b.0 - a.0).to_radians();
    let dlo = (b.1 - a.1).to_radians();
    let h = (dla * 0.5).sin().powi(2) + la1.cos() * la2.cos() * (dlo * 0.5).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

fn main() -> anyhow::Result<()> {
    let truth = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "40.876466,-73.524951".into());
    let (tlat, tlon): (f64, f64) = {
        let (a, b) = truth.split_once(',').expect("lat,lon");
        (a.trim().parse()?, b.trim().parse()?)
    };
    // Events are given as unix seconds of closest approach.
    let events: Vec<f64> = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "1788898180,1788899110".into())
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    let region = vec![Bounds { west: tlon - 1.5, south: tlat - 1.5, east: tlon + 1.5, north: tlat + 1.5 }];
    let t0 = 1788890400.0;
    let t1 = 1788910200.0;
    println!("loading archive …");
    let mut prog = |p: flights::FlightProgress| eprintln!("  {} traces, {} kept", p.traces_seen, p.flights_kept);
    let cancel = std::sync::atomic::AtomicBool::new(false);
    let days = flights::fetch_flights(t0 as i64, t1 as i64, &region, "staging-0", &mut prog, &cancel)?;
    println!("{} flights loaded", days.len());

    // Named aircraft: closest approach anywhere in the window.
    if let Ok(list) = std::env::var("AVALON_DEBUG_ICAO") {
        for name in list.split(',').map(|x| x.trim()).filter(|x| !x.is_empty()) {
            let mut best = (f64::INFINITY, 0.0, 0.0, false);
            let mut seen = false;
            for f in &days {
                if !f.icao.eq_ignore_ascii_case(name) {
                    continue;
                }
                seen = true;
                for p in &f.points {
                    let d = haversine((p.lat, p.lon), (tlat, tlon));
                    if d < best.0 {
                        best = (d, p.t_s, p.alt_m, p.alt_m > 50.0);
                    }
                }
            }
            if seen {
                println!(
                    "{name}: closest {:.1} km at t={:.0} ({:+.0}s from window start), alt {:.0} m",
                    best.0 / 1000.0,
                    best.1,
                    best.1 - t0,
                    best.2
                );
            } else {
                println!("{name}: not in archive");
            }
        }
    }

    for (i, ev) in events.iter().enumerate() {
        let mut near: Vec<(f64, String, f64, f64)> = Vec::new();
        for f in &days {
            let mut cpa = (f64::INFINITY, 0.0);
            for p in &f.points {
                if (p.t_s - ev).abs() > 90.0 || p.alt_m < 50.0 {
                    continue;
                }
                let d = haversine((p.lat, p.lon), (tlat, tlon));
                if d < cpa.0 {
                    cpa = (d, p.t_s);
                }
            }
            if cpa.0 < 25_000.0 {
                let alt = f
                    .points
                    .iter()
                    .filter(|p| (p.t_s - ev).abs() <= 20.0)
                    .map(|p| p.alt_m)
                    .fold(0.0f64, f64::max);
                near.push((cpa.0, f.icao.clone(), cpa.1, alt));
            }
        }
        near.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        println!("\n--- event {i} at t={ev:.0}: {} flights within 25 km", near.len());
        for (d, icao, t, alt) in near.iter().take(15) {
            println!("  {icao:8} {:.1} km  t={:.0} ({:+.0}s)  alt {:.0} m", d / 1000.0, t, t - ev, alt);
        }
    }
    Ok(())
}
