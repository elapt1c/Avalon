//! Bounded joint refit of the known pair.
//!
//! The previous attempts let f0/frequency-offset run free and produced absurd
//! nuisance values (f0 +202 %, offsets ±270 Hz) that "explain" any location —
//! degeneracy, not accuracy. This version keeps the nuisance parameters inside
//! physically defensible bounds and asks the two questions that matter:
//!
//!   1. floor at the TRUE receiver with bounded, well-seeded nuisance params
//!   2. does a joint (position + nuisance) fit started from a random cell across
//!      CONUS converge back to the truth? — i.e. is refinement enough to hit 10 km
//!
//! Usage: fit3 [W,S,E,N] [ICAO0,ICAO1] [seed_step_km]

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

type X = [f64; 7];

const F0_HALF: f64 = 0.03;
const OFF_HALF: f64 = 3.0;
const SHIFT_MAX: f64 = 90.0;

struct Probe {
    curves: Vec<Vec<(f64, f32)>>,
    tracks: Vec<score::Track>,
    f0: [f64; 2],
    anchor: f64,
    atm: score::Atmos,
    max_err: f64,
}

impl Probe {
    fn rmse(&self, x: &X) -> f64 {
        let obs = avalon::search::ecef(x[0], x[1], 0.0);
        let alpha = [x[3], x[4]];
        let off = [x[5], x[6]];
        let mut sum = 0.0;
        let mut n = 0usize;
        for k in 0..2 {
            let tr = &self.tracks[k];
            let f0 = self.f0[k] * (1.0 + alpha[k]);
            let mut tau = self.anchor + self.curves[k][0].0 - x[2] - 2000.0 / self.atm.c;
            for &(t_rel, f_meas) in &self.curves[k] {
                let model_t = self.anchor + t_rel - x[2];
                match score::doppler_ratio_warm(tr, obs, model_t, &mut tau, self.atm) {
                    Some(r) => {
                        let e = f0 * r + off[k] - f_meas as f64;
                        let a = e.abs();
                        sum += if a > self.max_err { self.max_err * self.max_err } else { e * e };
                    }
                    None => {
                        tau = model_t - 2000.0 / self.atm.c;
                        sum += self.max_err * self.max_err;
                    }
                }
                n += 1;
            }
        }
        (sum / n.max(1) as f64).sqrt()
    }

    /// RMSE plus soft walls that keep the nuisance parameters physical.
    fn cost(&self, x: &X) -> f64 {
        let mut c = self.rmse(x);
        let mut pen = 0.0;
        for a in [x[3], x[4]] {
            let ov = (a.abs() - F0_HALF).max(0.0) / F0_HALF;
            pen += ov * ov;
        }
        for o in [x[5], x[6]] {
            let ov = (o.abs() - OFF_HALF).max(0.0) / OFF_HALF;
            pen += ov * ov;
        }
        let ov = (x[2].abs() - SHIFT_MAX).max(0.0) / 20.0;
        pen += ov * ov;
        c + 200.0 * pen
    }

    fn nm(&self, seed: X, pos_step_deg: f64) -> (X, f64) {
        let n = 8;
        let mut simplex: Vec<(X, f64)> = vec![(seed, self.cost(&seed))];
        for i in 0..7 {
            let mut p = seed;
            p[i] += match i {
                0 | 1 => pos_step_deg,
                2 => 3.0,
                3 | 4 => 0.01,
                _ => 0.5,
            };
            simplex.push((p, self.cost(&p)));
        }
        for _ in 0..900 {
            simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            let best = simplex[0].0;
            let worst = simplex[n - 1].0;
            let mut cen = [0.0f64; 7];
            for s in simplex.iter().take(n - 1) {
                for i in 0..7 {
                    cen[i] += s.0[i] / (n - 1) as f64;
                }
            }
            let refl = std::array::from_fn(|i| cen[i] + (cen[i] - worst[i]));
            let fr = self.cost(&refl);
            if fr < simplex[0].1 {
                let exp = std::array::from_fn(|i| cen[i] + 2.0 * (cen[i] - worst[i]));
                let fe = self.cost(&exp);
                simplex[n - 1] = if fe < fr { (exp, fe) } else { (refl, fr) };
            } else if fr < simplex[n - 2].1 {
                simplex[n - 1] = (refl, fr);
            } else {
                let con = std::array::from_fn(|i| cen[i] + 0.5 * (worst[i] - cen[i]));
                let fc = self.cost(&con);
                if fc < simplex[n - 1].1 {
                    simplex[n - 1] = (con, fc);
                } else {
                    for j in 1..n {
                        let p = std::array::from_fn(|i| best[i] + 0.5 * (simplex[j].0[i] - best[i]));
                        simplex[j] = (p, self.cost(&p));
                    }
                }
            }
            let spread = simplex.iter().map(|s| (s.1 - simplex[0].1).abs()).fold(0.0, f64::max);
            if spread < 1e-11 {
                break;
            }
        }
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        simplex[0]
    }

    /// Floor at a fixed location, cheapest seed: coarse shift x f0 scan then NM.
    fn floor_at(&self, lat: f64, lon: f64) -> (X, f64) {
        let mut best: Option<(X, f64)> = None;
        for sh in [-60.0, -40.0, -20.0, 0.0, 20.0, 40.0, 60.0] {
            for a in [-0.02, 0.0, 0.02] {
                let s: X = [lat, lon, sh, a, a, 0.0, 0.0];
                let (x, c) = self.nm(s, 0.0);
                if best.as_ref().map_or(true, |b| c < b.1) {
                    best = Some((x, c));
                }
            }
        }
        best.unwrap()
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
        .map(|s| s.split(',').map(|x| x.trim().to_lowercase()).collect())
        .unwrap_or_else(|| vec!["a67ba8".into(), "a659a2".into()]);
    let seed_km: f64 = args.get(2).and_then(|v| v.parse().ok()).unwrap_or(400.0);
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
    let fl = flights::fetch_flights(win.0, win.1, &[region], "staging-0", &mut prog, &cancel)?;

    let decim: usize = std::env::var("AOS_FIT_DECIM").ok().and_then(|x| x.parse().ok()).unwrap_or(200);
    let mut curves: Vec<Vec<(f64, f32)>> = Vec::new();
    for ev in events.iter().take(2) {
        let r = ridges.iter().find(|r| r.id == ev.representative_ridge_id)
            .ok_or_else(|| anyhow::anyhow!("ridge"))?;
        let step = (r.samples.len() / decim).max(1);
        curves.push(r.samples.iter().step_by(step).map(|x| (x.t_s, x.freq_hz)).collect());
    }

    let track_of = |icao: &str| -> Option<score::Track> {
        let f = fl.iter().find(|f| f.icao.to_lowercase() == icao)?;
        let pts: Vec<(f64, [f64; 3])> = f.points.iter()
            .map(|p| (p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m))).collect();
        score::Track::from_points(&pts)
    };
    let mut tracks = Vec::new();
    let mut f0s = [0.0f64; 2];
    for (k, icao) in known.iter().enumerate() {
        let ev = events.get(k).ok_or_else(|| anyhow::anyhow!("event"))?;
        f0s[k] = ev.doppler.ok_or_else(|| anyhow::anyhow!("fit"))?.f0_hz as f64;
        tracks.push(track_of(icao).ok_or_else(|| anyhow::anyhow!("{icao} missing"))?);
    }
    let obs_truth = avalon::search::ecef(truth.0, truth.1, 0.0);
    let (_, t0) = tracks[0].cpa(obs_truth).ok_or_else(|| anyhow::anyhow!("cpa"))?;
    let anchor = t0 - events[0].doppler.ok_or_else(|| anyhow::anyhow!("fit"))?.t_cpa_s;

    let probe = Probe { curves, tracks, f0: f0s, anchor, atm: score::Atmos::default(), max_err: 400.0 };

    let (xt, ct) = probe.floor_at(truth.0, truth.1);
    println!("\n=== FLOOR AT TRUTH ({:.2},{:.2}) ===", truth.0, truth.1);
    println!("rmse {:.3} Hz | shift {:+.2} s | f0 {:+.2}% / {:+.2}% | off {:+.2} / {:+.2} Hz",
             probe.rmse(&xt), xt[2], xt[3] * 100.0, xt[4] * 100.0, xt[5], xt[6]);
    let (xf, cf) = probe.nm([truth.0, truth.1, xt[2], xt[3], xt[4], xt[5], xt[6]], 0.03);
    println!("joint fit from truth -> {:.3},{:.3} ({:.1} km) rmse {:.3} Hz",
             xf[0], xf[1], haversine_km((xf[0], xf[1]), truth), probe.rmse(&xf));

    // explicit seeds (semicolon separated "lat,lon") when provided
    if let Ok(list) = std::env::var("AOS_SEEDS") {
        println!("\n=== REFINE FROM EXPLICIT SEEDS ===");
        for part in list.split(';').filter(|x| !x.trim().is_empty()) {
            let p: Vec<f64> = part.split(',').filter_map(|x| x.trim().parse().ok()).collect();
            if p.len() != 2 { continue; }
            let seed: X = [p[0], p[1], 0.0, 0.0, 0.0, 0.0, 0.0];
            let (x, _) = probe.nm(seed, 0.05);
            println!("  seed {:.2},{:.2} ({:.0} km) -> {:.3},{:.3} rmse {:.2} Hz  {:.1} km from truth",
                     p[0], p[1], haversine_km((p[0], p[1]), truth),
                     x[0], x[1], probe.rmse(&x), haversine_km((x[0], x[1]), truth));
        }
    }

    // global multi-start
    let d = seed_km / 111.32;
    let mut seeds: Vec<X> = Vec::new();
    let mut lat = s;
    while lat <= n {
        let mut lon = w;
        while lon <= e {
            for sh in [-30.0, 0.0, 30.0] {
                seeds.push([lat, lon, sh, 0.0, 0.0, 0.0, 0.0]);
            }
            lon += d;
        }
        lat += d * 0.77;
    }
    if std::env::var("AOS_NOSEED").is_ok() {
        return Ok(());
    }
    eprintln!("{} joint starts", seeds.len());
    let mut fits: Vec<([f64; 7], f64)> = seeds.par_iter().map(|s| probe.nm(*s, d * 0.5)).collect();
    fits.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    println!("\n=== JOINT FITS FROM {} CONUS SEEDS ===", seeds.len());
    let mut kept: Vec<([f64; 7], f64)> = Vec::new();
    for f in fits.into_iter() {
        if kept.iter().any(|k| haversine_km((k.0[0], k.0[1]), (f.0[0], f.0[1])) < 40.0) { continue; }
        kept.push(f);
        if kept.len() >= 6 { break; }
    }
    for (i, (x, c)) in kept.iter().enumerate() {
        println!("  #{} {:.3},{:.3} rmse {:.2} Hz  {:.0} km from truth  (shift {:+.1}s f0 {:+.1}%/{:+.1}%)",
                 i + 1, x[0], x[1], probe.rmse(x), haversine_km((x[0], x[1]), truth),
                 x[2], x[3] * 100.0, x[4] * 100.0);
    }
    let best = &kept[0];
    println!("\nBEST JOINT: {:.0} km from truth (rmse {:.2} Hz); floor at truth {:.2} Hz",
             haversine_km((best.0[0], best.0[1]), truth), probe.rmse(&best.0), probe.rmse(&xt));
    Ok(())
}
