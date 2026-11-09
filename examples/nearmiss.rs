//! Which aircraft actually passed near the receiver, and which pairs of them are
//! ~920 s apart? Enumerates every flight whose closest approach falls inside a
//! window, sorted by CPA time, so the pair matching the measured inter-event gap
//! can be picked by eye instead of assumed.
//!
//! Usage: nearmiss [lat,lon] [radius_km] [unix_from] [unix_to]
use std::sync::atomic::AtomicBool;

use avalon::flights;
use avalon::model::Bounds;
use avalon::search::ecef;

fn haversine(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6_371_000.0f64;
    let (la1, la2) = (a.0.to_radians(), b.0.to_radians());
    let dla = (b.0 - a.0).to_radians();
    let dlo = (b.1 - a.1).to_radians();
    let h = (dla * 0.5).sin().powi(2) + la1.cos() * la2.cos() * (dlo * 0.5).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let truth = args.get(1).cloned().unwrap_or_else(|| "40.876466,-73.524951".into());
    let radius_km: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(15.0);
    let t_from: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1788900000.0);
    let t_to: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(1788902600.0);
    let gap_target: f64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(921.85);

    let (tlat, tlon): (f64, f64) = {
        let (a, b) = truth.split_once(',').expect("lat,lon");
        (a.trim().parse()?, b.trim().parse()?)
    };
    let region = vec![Bounds { west: tlon - 1.5, south: tlat - 1.5, east: tlon + 1.5, north: tlat + 1.5 }];
    let mut prog = |_p: flights::FlightProgress| {};
    let cancel = AtomicBool::new(false);
    let day = flights::fetch_flights(1788890400, 1788910200, &region, "staging-0", &mut prog, &cancel)?;

    // Closest approach of every flight to the receiver, restricted to the window.
    let mut hits: Vec<(f64, String, f64, f64)> = Vec::new(); // (t_cpa, icao, dist_km, alt_m)
    for f in &day {
        let mut best = (f64::INFINITY, 0.0, 0.0);
        for p in &f.points {
            if p.alt_m < 50.0 {
                continue;
            }
            let d = haversine((p.lat, p.lon), (tlat, tlon));
            if d < best.0 {
                best = (d, p.t_s, p.alt_m);
            }
        }
        if best.0.is_finite()
            && best.0 / 1000.0 <= radius_km
            && best.1 >= t_from
            && best.1 <= t_to
        {
            hits.push((best.1, f.icao.clone(), best.0 / 1000.0, best.2));
        }
    }
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    println!("{} flights with CPA < {radius_km:.0} km in the window:", hits.len());
    println!("{:>12} {:>8} {:>8} {:>9} {:>10}", "cpa time", "icao", "dist km", "alt m", "gap to next");
    for (i, (t, icao, d, alt)) in hits.iter().enumerate() {
        let gap = hits
            .get(i + 1)
            .map(|n| n.0 - t)
            .map(|g| format!("{g:+.0}"))
            .unwrap_or_else(|| "—".into());
        println!("{t:>12.0} {icao:>8} {d:>8.2} {alt:>9.0} {gap:>10}");
    }

    // Every pair whose separation matches the measured inter-event gap.
    println!("\npairs separated by {gap_target:.0} s +/- 5 s:");
    let mut found = 0;
    for i in 0..hits.len() {
        for j in i + 1..hits.len() {
            let gap = hits[j].0 - hits[i].0;
            if (gap - gap_target).abs() <= 5.0 {
                found += 1;
                println!(
                    "  {} ({:.1} km) + {} ({:.1} km)  gap {:.0} s  at {:.0} / {:.0}",
                    hits[i].1, hits[i].2, hits[j].1, hits[j].2, gap, hits[i].0, hits[j].0
                );
            }
        }
    }
    if found == 0 {
        println!("  none");
    }
    Ok(())
}
