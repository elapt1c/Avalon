//! Continuous refit of the known-pair geometry, then a surface map.
//!
//!  1. residual floor AT the true receiver with continuous nuisance parameters
//!  2. is that the global minimum across the box (nuisance held at the truth fit)
//!  3. how far the receiver can move before the residual degrades by 1/3/5 Hz
//!
//! Usage: fit2 [W,S,E,N] [ICAO0,ICAO1] [scan_km] [fine_km]

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

type X = [f64; 7];

struct Probe {
    curves: Vec<Vec<(f64, f32)>>,
    tracks: Vec<score::Track>,
    f0: [f64; 2],
    anchor: f64,
    atm: score::Atmos,
    max_err: f64,
}

impl Probe {
    #[inline]
    fn cost(&self, x: &X) -> f64 {
        self.cost_at(x[0], x[1], x[2], x[3], x[4], x[5], x[6])
    }

    fn cost_at(&self, lat: f64, lon: f64, shift: f64, a0: f64, a1: f64, o0: f64, o1: f64) -> f64 {
        let obs = avalon::search::ecef(lat, lon, 0.0);
        let alpha = [a0, a1];
        let off = [o0, o1];
        let mut sum = 0.0;
        let mut n = 0usize;
        for k in 0..2 {
            let tr = &self.tracks[k];
            let f0 = self.f0[k] * alpha[k].exp();
            let mut tau = self.anchor + self.curves[k][0].0 - shift - 2000.0 / self.atm.c;
            for &(t_rel, f_meas) in &self.curves[k] {
                let model_t = self.anchor + t_rel - shift;
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

    fn fit_nuisance(&self, lat: f64, lon: f64, seed: [f64; 5]) -> (X, f64) {
        let s: X = [lat, lon, seed[0], seed[1], seed[2], seed[3], seed[4]];
        let n = 6; // 1 vertex + 5 nuisance dimensions
        let mut simplex: Vec<(X, f64)> = vec![(s, self.cost(&s))];
        for i in 2..7 {
            let mut p = s;
            p[i] += match i {
                2 => 4.0,
                3 | 4 => 0.02,
                _ => 0.5,
            };
            simplex.push((p, self.cost(&p)));
        }
        for _ in 0..800 {
            simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            let best = simplex[0].0;
            let worst = simplex[n - 1].0;
            let mut centroid = [0.0f64; 7];
            for s in simplex.iter().take(n - 1) {
                for i in 0..7 {
                    centroid[i] += s.0[i] / (n - 1) as f64;
                }
            }
            let fix = |p: X| -> X { [lat, lon, p[2], p[3], p[4], p[5], p[6]] };
            let refl = fix(std::array::from_fn(|i| centroid[i] + (centroid[i] - worst[i])));
            let fr = self.cost(&refl);
            if fr < simplex[0].1 {
                let exp = fix(std::array::from_fn(|i| centroid[i] + 2.0 * (centroid[i] - worst[i])));
                let fe = self.cost(&exp);
                simplex[n - 1] = if fe < fr { (exp, fe) } else { (refl, fr) };
            } else if fr < simplex[n - 2].1 {
                simplex[n - 1] = (refl, fr);
            } else {
                let con = fix(std::array::from_fn(|i| centroid[i] + 0.5 * (worst[i] - centroid[i])));
                let fc = self.cost(&con);
                if fc < simplex[n - 1].1 {
                    simplex[n - 1] = (con, fc);
                } else {
                    for j in 1..n {
                        let p = fix(std::array::from_fn(|i| best[i] + 0.5 * (simplex[j].0[i] - best[i])));
                        let v = self.cost(&p);
                        simplex[j] = (p, v);
                    }
                }
            }
            let spread = simplex.iter().map(|s| (s.1 - simplex[0].1).abs()).fold(0.0, f64::max);
            if spread < 1e-10 {
                break;
            }
        }
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        simplex[0]
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
    let scan_km: f64 = args.get(2).and_then(|v| v.parse().ok()).unwrap_or(25.0);
    let fine_km: f64 = args.get(3).and_then(|v| v.parse().ok()).unwrap_or(140.0);
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

    let use_core = std::env::var("AOS_FIT_ALL").is_err();
    let decim: usize = std::env::var("AOS_FIT_DECIM").ok().and_then(|x| x.parse().ok()).unwrap_or(240);
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
        let step = (src.len() / decim).max(1);
        let curve: Vec<(f64, f32)> = src.into_iter().step_by(step).collect();
        eprintln!("ridge {}: {} pts from {} samples{}", r.id, curve.len(), r.samples.len(),
                  if use_core { " (weight>=median)" } else { " (all)" });
        curves.push(curve);
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
    let obs_truth = avalon::search::ecef(truth.0, truth.1, 0.0);
    let (_, t0) = tracks[0].cpa(obs_truth).ok_or_else(|| anyhow::anyhow!("cpa"))?;
    let anchor = t0 - events[0].doppler.ok_or_else(|| anyhow::anyhow!("fit"))?.t_cpa_s;

    let probe = Probe { curves, tracks, f0: f0s, anchor, atm: score::Atmos::default(), max_err: 400.0 };

    let (xt, c_truth) = probe.fit_nuisance(truth.0, truth.1, [0.0; 5]);
    println!("\n=== AT TRUTH ({:.2},{:.2}) ===", truth.0, truth.1);
    println!("rmse {:.3} Hz | shift {:+.2} s | f0 dev {:+.2}% / {:+.2}% | offset {:+.2} / {:+.2} Hz",
             c_truth, xt[2], (xt[3].exp() - 1.0) * 100.0, (xt[4].exp() - 1.0) * 100.0, xt[5], xt[6]);
    let nuis = [xt[2], xt[3], xt[4], xt[5], xt[6]];

    let step = scan_km / 111.32;
    let mut cells: Vec<(f64, f64, f64)> = Vec::new();
    let mut lat = s;
    while lat <= n {
        let mut lon = w;
        while lon <= e {
            cells.push((probe.cost_at(lat, lon, nuis[0], nuis[1], nuis[2], nuis[3], nuis[4]), lat, lon));
            lon += step;
        }
        lat += step * 0.77;
    }
    cells.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    println!("\n=== GLOBAL SCAN ({} cells @ {:.0} km, nuisance fixed) ===", cells.len(), scan_km);
    let mut kept: Vec<(f64, f64, f64)> = Vec::new();
    for &(c, la, lo) in cells.iter() {
        if kept.iter().any(|k| haversine_km((k.1, k.2), (la, lo)) < 60.0) { continue; }
        kept.push((c, la, lo));
        if kept.len() >= 6 { break; }
    }
    for (i, (c, la, lo)) in kept.iter().enumerate() {
        println!("  #{} {:.3},{:.3}  rmse {:.2} Hz  {:.0} km from truth", i + 1, la, lo, c,
                 haversine_km((*la, *lo), truth));
    }

    let fs = 3.0 / 111.32;
    let mut map: Vec<(f64, f64, f64)> = Vec::new();
    let mut la = truth.0 - fine_km / 111.0;
    while la <= truth.0 + fine_km / 111.0 {
        let mut lo = truth.1 - fine_km / 111.0 * 1.3;
        while lo <= truth.1 + fine_km / 111.0 * 1.3 {
            map.push((probe.cost_at(la, lo, nuis[0], nuis[1], nuis[2], nuis[3], nuis[4]), la, lo));
            lo += fs;
        }
        la += fs;
    }
    println!("\n=== UNCERTAINTY MAP (±{:.0} km, {:.0} km grid, {} cells) ===", fine_km, fs * 111.0, map.len());
    for (thr, label) in [(1.0, "+1 Hz"), (3.0, "+3 Hz"), (5.0, "+5 Hz")] {
        let sel: Vec<&(f64, f64, f64)> = map.iter().filter(|c| c.0 <= c_truth + thr).collect();
        let maxd = sel.iter().map(|c| haversine_km((c.1, c.2), truth)).fold(0.0, f64::max);
        println!("  rmse <= floor{}: {:>4} cells, extends {:.0} km", label, sel.len(), maxd);
    }
    map.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let (bc, bla, blo) = map[0];
    println!("  local best: {:.3},{:.3} rmse {:.3} Hz ({:.1} km from truth)", bla, blo, bc,
             haversine_km((bla, blo), truth));
    Ok(())
}
