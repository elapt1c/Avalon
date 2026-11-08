//! Continuous refinement of the known-pair fit.
//!
//! `closestspike` scores on a discrete grid (5 clock shifts x 5 f0 values), which
//! can bias the answer by tens of km. This fits the same two curves with a
//! Nelder-Mead search over [lat, lon, shift, f0_1, f0_2, off_1, off_2] and asks
//! the question that decides the viability of the method: at the TRUE receiver,
//! how low can the residual go, and is that optimum competitive with the best
//! location anywhere in CONUS?
//!
//! Usage: fitpair [W,S,E,N] [ICAO0,ICAO1]

use std::sync::atomic::AtomicBool;

use avalon::model::{Bounds, DetectedEvent, Ridge};
use avalon::{flights, score};

fn haversine_km(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6371.0f64;
    let (la1, la2) = (a.0.to_radians(), b.0.to_radians());
    let h = ((b.0 - a.0).to_radians() * 0.5).sin().powi(2)
        + la1.cos() * la2.cos() * ((b.1 - a.1).to_radians() * 0.5).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

#[derive(Clone)]
struct Probe {
    curves: Vec<Vec<(f64, f32)>>,
    tracks: Vec<score::Track>,
    f0: [f64; 2],
    t_shift_anchor: f64,
    atm: score::Atmos,
    max_err: f64,
}

impl Probe {
    /// Residual RMSE (Hz) for receiver at lat/lon with nuisance params.
    fn cost(&self, lat: f64, lon: f64, shift: f64, a1: f64, a2: f64, o1: f64, o2: f64) -> f64 {
        let obs = avalon::search::ecef(lat, lon, 0.0);
        let alpha = [a1, a2];
        let off = [o1, o2];
        let mut sum = 0.0;
        let mut n = 0usize;
        for k in 0..2 {
            let tr = &self.tracks[k];
            let f0 = self.f0[k] * alpha[k].exp();
            let mut tau = self.t_shift_anchor + self.curves[k][0].0 - shift - 2000.0 / self.atm.c;
            for &(t_rel, f_meas) in &self.curves[k] {
                let model_t = self.t_shift_anchor + t_rel - shift;
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

    /// Cost with lat/lon free, nuisance params optimized (starts at `seed`).
    fn fit_from(&self, seed: [f64; 7]) -> ([f64; 7], f64) {
        nelder_mead(seed, |x| self.cost(x[0], x[1], x[2], x[3], x[4], x[5], x[6]), 1.0)
    }

    /// Cost at a fixed location, nuisance-only optimization.
    fn fit_nuisance(&self, lat: f64, lon: f64, seed: [f64; 5]) -> f64 {
        let (x, c) = nelder_mead(
            [lat, lon, seed[0], seed[1], seed[2], seed[3], seed[4]],
            |x| self.cost(lat, lon, x[2], x[3], x[4], x[5], x[6]),
            0.5,
        );
        let _ = x;
        c
    }
}

fn nelder_mead<F: Fn(&[f64; 7]) -> f64>(x0: [f64; 7], f: F, scale: f64) -> ([f64; 7], f64) {
    let n = 7;
    let mut simplex: Vec<([f64; 7], f64)> = Vec::with_capacity(n + 1);
    simplex.push((x0, f(&x0)));
    for i in 0..n {
        let mut p = x0;
        let d = if i < 2 { scale / 111.0 } else if i == 2 { 2.0 * scale } else { 0.02 * scale };
        p[i] += d;
        let v = f(&p);
        simplex.push((p, v));
    }
    for _ in 0..600 {
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let best = simplex[0].0;
        let worst = simplex[n].0;
        let centroid: [f64; 7] = {
            let mut c = [0.0; 7];
            for s in simplex.iter().take(n) {
                for i in 0..7 {
                    c[i] += s.0[i] / n as f64;
                }
            }
            c
        };
        let refl: [f64; 7] = std::array::from_fn(|i| centroid[i] + (centroid[i] - worst[i]));
        let fr = f(&refl);
        if fr < simplex[0].1 {
            let exp: [f64; 7] = std::array::from_fn(|i| centroid[i] + 2.0 * (centroid[i] - worst[i]));
            let fe = f(&exp);
            simplex[n] = if fe < fr { (exp, fe) } else { (refl, fr) };
        } else if fr < simplex[n - 1].1 {
            simplex[n] = (refl, fr);
        } else {
            let con: [f64; 7] = std::array::from_fn(|i| centroid[i] + 0.5 * (worst[i] - centroid[i]));
            let fc = f(&con);
            if fc < simplex[n].1 {
                simplex[n] = (con, fc);
            } else {
                for j in 1..=n {
                    let p: [f64; 7] = std::array::from_fn(|i| best[i] + 0.5 * (simplex[j].0[i] - best[i]));
                    let v = f(&p);
                    simplex[j] = (p, v);
                }
            }
        }
        let spread: f64 = simplex
            .iter()
            .map(|s| (s.1 - simplex[0].1).abs())
            .fold(0.0, f64::max);
        if spread < 1e-9 {
            break;
        }
    }
    simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    simplex[0]
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
    eprintln!("loading archive…");
    let fl = flights::fetch_flights(win.0, win.1, &[region], "staging-0", &mut prog, &cancel)?;
    eprintln!("{} flights in box", fl.len());

    let use_core = std::env::var("AOS_FIT_ALL").is_err();
    let decim: usize = std::env::var("AOS_FIT_DECIM").ok().and_then(|s| s.parse().ok()).unwrap_or(240);
    let mut curves: Vec<Vec<(f64, f32)>> = Vec::new();
    for ev in events.iter().take(2) {
        let r = ridges.iter().find(|r| r.id == ev.representative_ridge_id)
            .ok_or_else(|| anyhow::anyhow!("ridge"))?;
        let mut ws: Vec<f32> = r.samples.iter().map(|x| x.weight).collect();
        ws.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let med = ws[ws.len() / 2].max(1.0);
        let src: Vec<(f64, f32)> = if use_core {
            r.samples.iter().filter(|x| x.weight >= med).map(|x| (x.t_s, x.freq_hz)).collect()
        } else {
            r.samples.iter().map(|x| (x.t_s, x.freq_hz)).collect()
        };
        eprintln!("ridge {} : {} samples (using {}, weight filter {})",
                  r.id, r.samples.len(), src.len(), if use_core { ">= median" } else { "none" });
        let step = (src.len() / decim).max(1);
        curves.push(src.into_iter().step_by(step).collect());
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
        let fit = ev.doppler.ok_or_else(|| anyhow::anyhow!("fit"))?;
        f0s[k] = fit.f0_hz as f64;
        tracks.push(track_of(icao).ok_or_else(|| anyhow::anyhow!("{icao} missing"))?);
    }

    // Shared clock anchor: flight 1's CPA maps recording time -> world time.
    let obs_truth = avalon::search::ecef(truth.0, truth.1, 0.0);
    let (_, t0) = tracks[0].cpa(obs_truth).ok_or_else(|| anyhow::anyhow!("cpa"))?;
    let anchor = t0 - events[0].doppler.ok_or_else(|| anyhow::anyhow!("fit"))?.t_cpa_s;

    let probe = Probe {
        curves,
        tracks,
        f0: f0s,
        t_shift_anchor: anchor,
        atm: score::Atmos::default(),
        max_err: 400.0,
    };

    // 1) How low does the residual go AT TRUTH?
    let seed_nuisance = [0.0, 0.0, 0.0, 0.0, 0.0];
    let c_truth = probe.fit_nuisance(truth.0, truth.1, seed_nuisance);
    println!("\nresidual floor AT TRUTH ({:.2},{:.2}): {:.3} Hz", truth.0, truth.1, c_truth);

    // 2) Multi-start over CONUS: best location the curves allow.
    let mut starts: Vec<[f64; 7]> = Vec::new();
    let mut lat = s;
    while lat <= n {
        let mut lon = w;
        while lon <= e {
            starts.push([lat, lon, 0.0, 0.0, 0.0, 0.0, 0.0]);
            lon += 6.0;
        }
        lat += 5.0;
    }
    eprintln!("{} starts", starts.len());
    use rayon::prelude::*;
    let mut fits: Vec<([f64; 7], f64)> = starts.par_iter().map(|s| probe.fit_from(*s)).collect();
    fits.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    println!("\nbest locations found (continuous nuisance fit):");
    println!("{:<5} {:<10} {:<10} {:>9} {:>10} {:>8} {:>8}", "rank", "lat", "lon", "rmse Hz", "km", "shift s", "f0 dev%");
    let mut shown: Vec<(f64, f64)> = Vec::new();
    let mut rank = 0;
    for (x, c) in fits.iter() {
        if shown.iter().any(|p| haversine_km(*p, (x[0], x[1])) < 25.0) {
            continue;
        }
        shown.push((x[0], x[1]));
        rank += 1;
        println!("{:<5} {:<10.4} {:<10.4} {:>9.3} {:>10.0} {:>8.1} {:>8.1}",
                 rank, x[0], x[1], c, haversine_km((x[0], x[1]), truth), x[2],
                 (x[3].exp() - 1.0) * 100.0);
        if rank >= 8 {
            break;
        }
    }
    let (bx, bc) = fits[0];
    println!("\nBEST ANYWHERE: {:.0} km from truth, rmse {:.3} Hz (at truth: {:.3} Hz)",
             haversine_km((bx[0], bx[1]), truth), bc, c_truth);
    println!("verdict: {}", if bc + 0.02 < c_truth {
        "curves prefer a WRONG location — acoustic/model bias, not geometry"
    } else {
        "truth is the global optimum of these curves"
    });
    Ok(())
}
