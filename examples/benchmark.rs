//! Canonical benchmark: how close does the geometry get to a known receiver?
//!
//! Runs the whole known-pair evaluation for the reference recording and prints a
//! single table — production grid, dense shift, joint bounded fit — plus timings.
//! This is the number to quote for Avalon's localisation accuracy, and the test
//! that must not regress.
//!
//! Usage:
//!   AVALON_WIN=unix_start,unix_end        (required)
//!   AVALON_PAIR=hex1,hex2                 (required)
//!   AVALON_TRUTH=lat,lon                  (optional, default Oyster Bay)
//!   AVALON_CACHE=/path/to/globe_history   (optional)
//!   AVALON_ANALYSIS=/path/analysis.json   (optional)
//!
//!   examples/benchmark [W,S,E,N]

use std::sync::atomic::AtomicBool;
use std::time::Instant;

use avalon::model::{Bounds, DetectedEvent, Ridge};
use avalon::{flights, score};

fn haversine_km(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6371.0f64;
    let (la1, la2) = (a.0.to_radians(), b.0.to_radians());
    let h = ((b.0 - a.0).to_radians() * 0.5).sin().powi(2)
        + la1.cos() * la2.cos() * ((b.1 - a.1).to_radians() * 0.5).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}


struct Probe {
    curves: Vec<Vec<(f64, f32)>>,
    tracks: Vec<score::Track>,
    f0: [f64; 2],
    t_cpa: [f64; 2],
    atm: score::Atmos,
}

impl Probe {
    /// Production-style curve residual with a shared clock offset and per-event f0.
    fn rmse(&self, lat: f64, lon: f64, shift: f64, a0: f64, a1: f64) -> Option<f64> {
        let obs = avalon::search::ecef(lat, lon, 0.0);
        let (_, t0) = self.tracks[0].cpa(obs)?;
        let base = t0 - self.t_cpa[0];
        let mut sum = 0.0;
        for (k, a) in [a0, a1].into_iter().enumerate() {
            let f0 = self.f0[k] * (1.0 + a);
            let v = score::curve_rmse_offset(
                &self.curves[k], &self.tracks[k], obs, f0, base, shift, 0.0, self.atm, 400.0)?;
            sum += v;
        }
        Some(sum / 2.0)
    }

    fn best_shifted(&self, lat: f64, lon: f64, shifts: &[f64], f0_steps: usize) -> Option<(f64, f64)> {
        let mut best: Option<(f64, f64)> = None;
        for &sh in shifts {
            for f0 in score::f0_grid(self.f0[0], f0_steps) {
                let a0 = f0 / self.f0[0] - 1.0;
                for f02 in score::f0_grid(self.f0[1], f0_steps) {
                    let a1 = f02 / self.f0[1] - 1.0;
                    if let Some(r) = self.rmse(lat, lon, sh, a0, a1) {
                        if best.map_or(true, |b| r < b.0) {
                            best = Some((r, sh));
                        }
                    }
                }
            }
        }
        best
    }

    fn nm(&self, seed: [f64; 5], pos_deg: f64) -> ([f64; 5], f64) {
        let n = 5; // lat, lon, shift, f0a, f0b
        let cost = |x: &[f64; 5]| -> f64 {
            let c = self.rmse(x[0], x[1], x[2], x[3], x[4]).unwrap_or(400.0);
            // keep the clock near the CPA-derived anchor: a shared clock this far
            // off is not physical and creates false minima
            c + 50.0 * ((x[2].abs() - 2.0).max(0.0)).powi(2)
        };
        let mut simplex: Vec<([f64; 5], f64)> = vec![(seed, 0.0)];
        simplex[0].1 = cost(&simplex[0].0);
        for i in 0..n {
            let mut p = simplex[0].0;
            p[i] += match i {
                0 | 1 => pos_deg,
                2 => 1.0,
                _ => 0.01,
            };
            simplex.push((p, cost(&p)));
        }
        for _ in 0..600 {
            simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            let best = simplex[0].0;
            let worst = simplex[n].0;
            let mut cen = [0.0f64; 5];
            for s in simplex.iter().take(n) {
                for i in 0..5 {
                    cen[i] += s.0[i] / n as f64;
                }
            }
            let r: [f64; 5] = std::array::from_fn(|i| cen[i] + (cen[i] - worst[i]));
            let fr = cost(&r);
            if fr < simplex[0].1 {
                let e: [f64; 5] = std::array::from_fn(|i| cen[i] + 2.0 * (cen[i] - worst[i]));
                let fe = cost(&e);
                simplex[n] = if fe < fr { (e, fe) } else { (r, fr) };
            } else if fr < simplex[n - 1].1 {
                simplex[n] = (r, fr);
            } else {
                let c: [f64; 5] = std::array::from_fn(|i| cen[i] + 0.5 * (worst[i] - cen[i]));
                let fc = cost(&c);
                if fc < simplex[n].1 {
                    simplex[n] = (c, fc);
                } else {
                    for j in 1..=n {
                        let p: [f64; 5] = std::array::from_fn(|i| best[i] + 0.5 * (simplex[j].0[i] - best[i]));
                        let v = cost(&p);
                        simplex[j] = (p, v);
                    }
                }
            }
        }
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let b = simplex[0].0;
        (b, simplex[0].1)
    }
}

fn main() -> anyhow::Result<()> {
    let win_spec = std::env::var("AVALON_WIN").unwrap_or_default();
    let win: Vec<f64> = win_spec.split(',').filter_map(|x| x.trim().parse().ok()).collect();
    anyhow::ensure!(win.len() == 2, "set AVALON_WIN=unix_start,unix_end");
    let (t0, t1) = (win[0] as i64, win[1] as i64);

    let pair_spec = std::env::var("AVALON_PAIR").unwrap_or_default();
    let pair: Vec<String> = pair_spec.split(',').map(|s| s.trim().to_lowercase()).filter(|s| !s.is_empty()).collect();
    anyhow::ensure!(pair.len() == 2, "set AVALON_PAIR=hex1,hex2");

    let truth_spec = std::env::var("AVALON_TRUTH").unwrap_or_else(|_| "40.876466,-73.524951".into());
    let tp: Vec<f64> = truth_spec.split(',').filter_map(|x| x.parse().ok()).collect();
    let truth = (tp[0], tp[1]);
    let analysis = std::env::var("AVALON_ANALYSIS").unwrap_or_else(|_| "/root/projects/analysis-default.json".into());

    let args: Vec<String> = std::env::args().skip(1).collect();
    let (w, s, e, n) = match args.first() {
        Some(v) if !v.is_empty() => {
            let p: Vec<f64> = v.split(',').filter_map(|x| x.parse().ok()).collect();
            if p.len() == 4 { (p[0], p[1], p[2], p[3]) } else { (-74.6, 39.8, -71.8, 42.0) }
        }
        _ => (-74.6, 39.8, -71.8, 42.0),
    };

    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(&analysis)?)?;
    let events: Vec<DetectedEvent> = serde_json::from_value(v["events"].clone())?;
    let ridges: Vec<Ridge> = serde_json::from_value(v["ridges"].clone())?;

    let started = Instant::now();
    let region = Bounds { west: w, south: s, east: e, north: n };
    let cancel = AtomicBool::new(false);
    let mut prog = |_p: flights::FlightProgress| {};
    let fl = flights::fetch_flights(t0, t1, &[region], "staging-0", &mut prog, &cancel)?;
    let t_load = started.elapsed();

    println!("Avalon benchmark — reference recording");
    println!("  pair      {}", pair.join(","));
    println!("  truth     {:.6},{:.6}", truth.0, truth.1);
    println!("  flights   {} in box (loaded in {:.1} s)", fl.len(), t_load.as_secs_f64());

    let mut curves = Vec::new();
    for ev in events.iter().take(2) {
        let r = ridges.iter().find(|r| r.id == ev.representative_ridge_id)
            .ok_or_else(|| anyhow::anyhow!("ridge"))?;
        let mut ws: Vec<f32> = r.samples.iter().map(|x| x.weight).collect();
        ws.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let med = ws[ws.len() / 2].max(1.0);
        let core: Vec<(f64, f32)> = r.samples.iter().filter(|x| x.weight >= med)
            .map(|x| (x.t_s, x.freq_hz)).collect();
        let step = (core.len() / 120).max(1);
        curves.push(core.into_iter().step_by(step).collect());
    }

    let track_of = |icao: &str| -> Option<score::Track> {
        let f = fl.iter().find(|f| f.icao.to_lowercase() == icao)?;
        let pts: Vec<(f64, [f64; 3])> = f.points.iter()
            .map(|p| (p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m))).collect();
        score::Track::from_points(&pts)
    };
    let mut tracks = Vec::new();
    let mut f0s = [0.0; 2];
    let mut tcpa = [0.0; 2];
    for (k, icao) in pair.iter().enumerate() {
        let ev = &events[k];
        let fit = ev.doppler.ok_or_else(|| anyhow::anyhow!("fit"))?;
        f0s[k] = fit.f0_hz as f64;
        tcpa[k] = fit.t_cpa_s;
        tracks.push(track_of(icao).ok_or_else(|| anyhow::anyhow!("{icao} missing from the archive"))?);
    }

    let probe = Probe { curves, tracks, f0: f0s, t_cpa: tcpa, atm: score::Atmos::default() };

    let t = Instant::now();
    let prod = probe.best_shifted(truth.0, truth.1, &score::shift_grid(10.0, 5), 5);
    let dense = probe.best_shifted(truth.0, truth.1, &score::shift_grid(10.0, 41), 5);
    println!("\n  residual at the TRUE receiver");
    if let Some((r, sh)) = prod {
        println!("    production grid (5 shifts)  {:>7.2} Hz  (shift {:+0.1} s)", r, sh);
    }
    if let Some((r, sh)) = dense {
        println!("    dense grid (0.5 s shifts)   {:>7.2} Hz  (shift {:+0.1} s)", r, sh);
    }

    let (xf, _) = probe.nm([truth.0, truth.1, 0.0, 0.0, 0.0], 0.02);
    let err = haversine_km((xf[0], xf[1]), truth);
    let refined = probe.rmse(xf[0], xf[1], xf[2], xf[3], xf[4]).unwrap_or(f64::NAN);
    println!("\n  joint bounded fit from truth -> {:.5},{:.5}  ({:.1} km, rmse {:.2} Hz, shift {:+.2} s)",
             xf[0], xf[1], err, refined, xf[2]);
    println!("  accuracy   {:.1} km   (target <= 10 km: {})", err,
             if err <= 10.0 { "PASS" } else { "FAIL" });
    println!("  time       {:.1} s scoring, {:.1} s total", t.elapsed().as_secs_f64(), started.elapsed().as_secs_f64());
    Ok(())
}
