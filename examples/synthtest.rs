//! Known-answer test: synthesise audio from a real ADS-B track, then check
//! whether the scorer ranks that track's true receiver first.
//!
//! This removes the two unknowns of the real clip (is the date right, are the
//! aircraft in the archive) and measures the locator itself.
//!
//! Usage: synthtest [--truth LAT,LON] [--radius-km N]

use std::sync::atomic::AtomicBool;

use avalon::model::Bounds;
use avalon::{flights, score};

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let mut truth = (40.87, -73.52);
    let mut radius_km = 30.0;
    let mut i = 1;
    while i + 1 < args.len() {
        match args[i].as_str() {
            "--truth" => {
                let v: Vec<f64> = args[i + 1].split(',').filter_map(|s| s.parse().ok()).collect();
                if v.len() == 2 {
                    truth = (v[0], v[1]);
                }
            }
            "--radius-km" => {
                if let Ok(v) = args[i + 1].parse() {
                    radius_km = v;
                }
            }
            _ => {}
        }
        i += 2;
    }

    let region = Bounds {
        west: truth.1 - 0.9555,
        south: truth.0 - 0.7226,
        east: truth.1 + 0.9555,
        north: truth.0 + 0.7226,
    };
    let win = (
        avalon::dateutil::civil_to_unix(2026, 10, 30, 10, 0, 0),
        avalon::dateutil::civil_to_unix(2026, 11, 1, 12, 0, 0),
    );

    let cancel = AtomicBool::new(false);
    let mut prog = |_p: flights::FlightProgress| {};
    let fl = flights::fetch_flights(win.0, win.1, &[region], "staging-0", &mut prog, &cancel)?;
    println!("{} flights", fl.len());

    let obs = avalon::search::ecef(truth.0, truth.1, 0.0);

    // Pick the flight that passes closest to the receiver.
    let mut best: Option<(f64, String, score::Track)> = None;
    for f in &fl {
        let mut pts: Vec<(f64, [f64; 3])> = Vec::new();
        let mut last = f64::NEG_INFINITY;
        for p in &f.points {
            if p.t_s - last < 5.0 {
                continue;
            }
            last = p.t_s;
            pts.push((p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m)));
        }
        let Some(tr) = score::Track::from_points(&pts) else { continue };
        let Some((d, _t)) = tr.cpa(obs) else { continue };
        if d > 8000.0 {
            continue;
        }
        if best.as_ref().map(|(bd, _, _)| d < *bd).unwrap_or(true) {
            best = Some((d, f.icao.clone(), tr));
        }
    }
    let Some((d_true, icao, track)) = best else {
        anyhow::bail!("no flight passes within 8 km of the truth point");
    };
    let (_, t_cpa) = track.cpa(obs).unwrap();
    println!("using {icao}: closest approach {:.2} km at t={t_cpa:.0}", d_true / 1000.0);

    // Synthesise: 3 harmonics, amplitude 1/R, f0 = 110 Hz, 150 s window.
    let sr = 4000u32;
    let dur = 150.0f64;
    let t0 = t_cpa - dur / 2.0;
    let f0 = 110.0;
    let harmonics = [1.0f64, 2.0, 3.0];
    let mut audio = vec![0.0f32; (sr as f64 * dur) as usize];
    let mut phases = [0.0f64; 3];
    for (i, s) in audio.iter_mut().enumerate() {
        let t = t0 + i as f64 / sr as f64;
        let Some(ratio) = score::doppler_ratio(&track, obs, t, score::Atmos::default()) else {
            continue;
        };
        let p = track.at(t).unwrap();
        let r = ((p[0] - obs[0]).powi(2) + (p[1] - obs[1]).powi(2) + (p[2] - obs[2]).powi(2)).sqrt();
        let amp = (1500.0 / r.max(200.0)) as f64;
        let mut v = 0.0f64;
        for (h, ph) in harmonics.iter().zip(phases.iter_mut()) {
            v += amp * (0.6 / h) * ph.sin();
            *ph += 2.0 * std::f64::consts::PI * f0 * h * ratio / sr as f64;
        }
        *s = (v * 0.35).clamp(-1.0, 1.0) as f32;
    }

    // Brightness image of the synthetic clip.
    let spec = score::EventSpec {
        f0: f0 as f32,
        harmonics: harmonics.iter().map(|h| *h as f32).collect(),
        cpa_s: t_cpa - t0,
        half_width_s: 60.0,
    };
    let band = spec.band();
    let b = score::Brightness::build(&audio, sr, 0.0, band, 1.0, 0.25);
    println!("brightness {} frames, band {:?}", b.frames, band);

    // Score a grid around the receiver; report the rank of the true cell.
    let step_deg = 2.0 / 111.0; // ~2 km
    let n = (radius_km / 2.0) as i32;
    let mut cells: Vec<(f32, f64, f64, f64)> = Vec::new(); // (score, lat, lon, km from truth)
    for dy in -n..=n {
        for dx in -n..=n {
            let lat = truth.0 + dy as f64 * step_deg;
            let lon = truth.1 + dx as f64 * step_deg * 1.33;
            let o = avalon::search::ecef(lat, lon, 0.0);
            let mut best_s = f32::NEG_INFINITY;
            for &f in &score::f0_grid(f0, 9) {
                for &sh in &score::shift_grid(4.0, 9) {
                    let s = score::score_event(
                        &b,
                        &track,
                        o,
                        &spec,
                        f,
                        // Absolute time of frame 0 — not 0.0, or every model
                        // lookup falls outside the track and scores zero.
                        t0,
                        sh,
                        score::Atmos::default(),
                        2,
                    ) as f32;
                    if s > best_s {
                        best_s = s;
                    }
                }
            }
            let dlat = (lat - truth.0) * 111.32;
            let dlon = (lon - truth.1) * 111.32 * truth.0.to_radians().cos();
            cells.push((best_s, lat, lon, (dlat * dlat + dlon * dlon).sqrt()));
        }
    }
    let mut sorted = cells.clone();
    sorted.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let rank = sorted
        .iter()
        .position(|c| c.3 < 1.0)
        .map(|i| i + 1)
        .unwrap_or(0);
    println!(
        "\ncells {} · best {:.2} at {:.3},{:.3} ({:.1} km from truth)",
        sorted.len(),
        sorted[0].0,
        sorted[0].1,
        sorted[0].2,
        sorted[0].3
    );
    println!("true cell rank: {rank} of {}", sorted.len());
    for c in sorted.iter().take(6) {
        println!("  {:.2}  {:.4},{:.4}  {:.1} km", c.0, c.1, c.2, c.3);
    }
    Ok(())
}
