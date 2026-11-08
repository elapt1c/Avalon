//! Constant-velocity Doppler parameterization of a detected ridge.
//!
//! Physical model (stationary receiver, straight constant-velocity pass,
//! retarded sound-arrival time neglected since v << c):
//!
//!   f(t) = f0 * c / (c + u_r(t))
//!   u_r(t) = v * x / r,   x = v (t - t_cpa),   r = sqrt(d^2 + x^2)
//!
//! The fit estimates rest frequency `f0`, speed `v`, closest-approach slant
//! distance `d`, and closest-approach time `t_cpa`. This is the same first-stage
//! "bound speed and slant distance" fit DOPPLERGEIST derives from a ridge.

use crate::dsp::median_f32;
use crate::model::{DopplerFit, Ridge, SOUND_SPEED};

/// Minimum number of ridge samples required to attempt a fit.
const MIN_SAMPLES: usize = 8;
/// Max samples used for fitting (keeps the multi-start search cheap).
const MAX_FIT_SAMPLES: usize = 120;
/// Inlier tolerance floor (Hz) for the robust refit.
const INLIER_TOL_HZ: f64 = 2.0;

/// Fit a constant-velocity Doppler model, robustly.
///
/// Multi-start Nelder-Mead over a deterministic grid of plausible
/// `(f0, v, d)` seeds, then iteratively reweighted refits that discard gross
/// outliers. A ridge that is not a real sigmoid will converge to a poor fit and
/// report a low inlier fraction, letting the pipeline reject it.
pub fn fit_doppler(ridge: &Ridge) -> Option<DopplerFit> {
    let samples = &ridge.samples;
    if samples.len() < MIN_SAMPLES {
        return None;
    }

    // Subsample uniformly to keep the search bounded.
    let stride = (samples.len() / MAX_FIT_SAMPLES).max(1);
    let sub: Vec<&crate::model::RidgeSample> = samples.iter().step_by(stride).collect();
    if sub.len() < MIN_SAMPLES {
        return None;
    }
    let ts: Vec<f64> = sub.iter().map(|s| s.t_s).collect();
    let fs: Vec<f64> = sub.iter().map(|s| s.freq_hz as f64).collect();

    // Asymptotic extremes.
    let mut freqs: Vec<f32> = sub.iter().map(|s| s.freq_hz).collect();
    freqs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let k = (freqs.len() as f32 * 0.15).ceil() as usize;
    let low_tail = &freqs[..k.max(1)];
    let high_tail = &freqs[freqs.len() - k.max(1)..];
    let f_min = median_f32(low_tail) as f64;
    let f_max = median_f32(high_tail) as f64;
    if f_max - f_min < 3.0 {
        return None;
    }

    let f0_mid = 2.0 * f_max * f_min / (f_max + f_min);
    let v_mid = (SOUND_SPEED * (f_max - f_min) / (f_max + f_min)).clamp(20.0, 300.0);

    // Deterministic multi-start grid. `f0` is constrained to the observed band:
    // the emitted tone is always inside the sweep, and letting it drift produces
    // degenerate fits for one-sided sweeps.
    let f0_lo = f_min * 0.97;
    let f0_hi = f_max * 1.03;
    let mut best: Option<([f64; 4], f64)> = None; // (params, sse)
    for f0_scale in [0.92, 1.0, 1.08] {
        for v_scale in [0.7, 1.0, 1.4] {
            for d in [150.0, 500.0, 1500.0, 5000.0] {
                let f0 = (f0_mid * f0_scale).clamp(f0_lo, f0_hi);
                let v = (v_mid * v_scale).clamp(20.0, 300.0);
                let t_cpa = crossing_time(&ts, &fs, f0);
                let p0 = [f0, v, d, t_cpa];
                let p = nelder_mead(
                    |p| residual_bounded(&ts, &fs, p, f0_lo, f0_hi),
                    p0,
                    1e-9,
                    220,
                );
                let sse = residual_bounded(&ts, &fs, &p, f0_lo, f0_hi);
                if best.map_or(true, |(_, b)| sse < b) {
                    best = Some((p, sse));
                }
            }
        }
    }
    let (mut p, _) = best?;

    // Iteratively reweighted refit: drop points far from the current curve so a
    // few spurious samples cannot drag the parameters.
    for _ in 0..4 {
        let tol = inlier_tolerance(&ts, &fs, &p, f0_lo, f0_hi);
        let (ts_in, fs_in) = inliers(&ts, &fs, &p, tol, f0_lo, f0_hi);
        if ts_in.len() < MIN_SAMPLES {
            break;
        }
        p = nelder_mead(
            |q| residual_bounded(&ts_in, &fs_in, q, f0_lo, f0_hi),
            p,
            1e-9,
            260,
        );
    }

    let tol = inlier_tolerance(&ts, &fs, &p, f0_lo, f0_hi);
    let n_in = inliers(&ts, &fs, &p, tol, f0_lo, f0_hi).0.len();
    let inlier_fraction = n_in as f32 / ts.len() as f32;
    let rmse = (residual_bounded(&ts, &fs, &p, f0_lo, f0_hi) / fs.len() as f64).sqrt();
    let complete_sweep = (f_max - p[0]).abs() > 0.02 * p[0] && (p[0] - f_min).abs() > 0.02 * p[0];

    Some(DopplerFit {
        f0_hz: p[0] as f32,
        speed_mps: p[1] as f32,
        slant_distance_m: p[2] as f32,
        t_cpa_s: p[3],
        rmse_hz: rmse as f32,
        inlier_fraction,
        n_samples: samples.len(),
        complete_sweep,
    })
}

/// Robust inlier band: max(floor, ~2.5 sigma via median absolute deviation).
fn inlier_tolerance(ts: &[f64], fs: &[f64], p: &[f64; 4], f0_lo: f64, f0_hi: f64) -> f64 {
    let res: Vec<f32> = ts
        .iter()
        .zip(fs)
        .map(|(&t, &f)| (model_bounded(t, p, f0_lo, f0_hi) - f).abs() as f32)
        .collect();
    let med = median_f32(&res) as f64;
    let dev: Vec<f32> = res.iter().map(|r| (*r as f64 - med).abs() as f32).collect();
    let mad = median_f32(&dev) as f64;
    (2.5 * 1.4826 * mad).max(INLIER_TOL_HZ)
}

fn inliers<'a>(
    ts: &'a [f64],
    fs: &'a [f64],
    p: &[f64; 4],
    tol: f64,
    f0_lo: f64,
    f0_hi: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut a = Vec::with_capacity(ts.len());
    let mut b = Vec::with_capacity(fs.len());
    for (&t, &f) in ts.iter().zip(fs) {
        if (model_bounded(t, p, f0_lo, f0_hi) - f).abs() <= tol {
            a.push(t);
            b.push(f);
        }
    }
    (a, b)
}

/// Evaluate the fitted Doppler curve at time `t` (for plotting overlays).
pub fn evaluate_fit(fit: &DopplerFit, t: f64) -> f64 {
    model(
        t,
        &[
            fit.f0_hz as f64,
            fit.speed_mps as f64,
            fit.slant_distance_m as f64,
            fit.t_cpa_s,
        ],
    )
}

/// Physical model with parameters clamped to plausible ranges.
fn model(t: f64, p: &[f64]) -> f64 {
    model_bounded(t, p, 5.0, 20000.0)
}

/// Physical model with the rest frequency constrained to the observed band.
///
/// A flyover's emitted tone always lies inside the observed sweep, so allowing
/// the fit to wander outside it (which the optimiser will happily do for a
/// one-sided sweep) produces degenerate solutions that then fail plausibility.
fn model_bounded(t: f64, p: &[f64], f0_lo: f64, f0_hi: f64) -> f64 {
    let f0 = p[0].clamp(f0_lo, f0_hi);
    let v = p[1].clamp(10.0, 300.0);
    let d = p[2].clamp(5.0, 200000.0);
    let tcpa = p[3];
    let x = v * (t - tcpa);
    let r = (d * d + x * x).sqrt();
    let ur = v * x / r;
    f0 * SOUND_SPEED / (SOUND_SPEED + ur)
}

fn residual_bounded(ts: &[f64], fs: &[f64], p: &[f64], f0_lo: f64, f0_hi: f64) -> f64 {
    ts.iter()
        .zip(fs)
        .map(|(&t, &f)| {
            let e = model_bounded(t, p, f0_lo, f0_hi) - f;
            e * e
        })
        .sum()
}

/// Interpolate the time where the ridge crosses `f0`.
fn crossing_time(ts: &[f64], fs: &[f64], f0: f64) -> f64 {
    for i in 0..fs.len().saturating_sub(1) {
        let (f1, f2) = (fs[i], fs[i + 1]);
        if (f1 - f0) * (f2 - f0) <= 0.0 {
            let t = if (f2 - f1).abs() > 1e-9 {
                ts[i] + (ts[i + 1] - ts[i]) * (f0 - f1) / (f2 - f1)
            } else {
                (ts[i] + ts[i + 1]) / 2.0
            };
            return t;
        }
    }
    (ts.first().unwrap() + ts.last().unwrap()) / 2.0
}

/// Nelder-Mead simplex refinement (derivative-free, scale-tolerant).
fn nelder_mead<F: Fn(&[f64]) -> f64>(f: F, x0: [f64; 4], tol: f64, max_iter: usize) -> [f64; 4] {
    let n = 4;
    let mut simplex: Vec<[f64; 4]> = Vec::with_capacity(n + 1);
    simplex.push(x0);
    for i in 0..n {
        let mut p = x0;
        p[i] += (x0[i].abs() * 0.1).max(1e-3);
        simplex.push(p);
    }
    let mut vals: Vec<f64> = simplex.iter().map(|p| f(p)).collect();

    for _ in 0..max_iter {
        let mut idx: Vec<usize> = (0..=n).collect();
        idx.sort_by(|&a, &b| vals[a].partial_cmp(&vals[b]).unwrap_or(std::cmp::Ordering::Equal));
        let best = idx[0];
        let worst = idx[n];

        // Centroid of all but the worst.
        let mut centroid = [0.0f64; 4];
        for &i in &idx[..n] {
            for k in 0..n {
                centroid[k] += simplex[i][k];
            }
        }
        for k in 0..n {
            centroid[k] /= n as f64;
        }

        // Reflect.
        let mut refl = [0.0f64; 4];
        for k in 0..n {
            refl[k] = centroid[k] + 1.0 * (centroid[k] - simplex[worst][k]);
        }
        let rval = f(&refl);

        if rval < vals[best] {
            // Expand.
            let mut exp = [0.0f64; 4];
            for k in 0..n {
                exp[k] = centroid[k] + 2.0 * (refl[k] - centroid[k]);
            }
            let eval = f(&exp);
            if eval < rval {
                simplex[worst] = exp;
                vals[worst] = eval;
            } else {
                simplex[worst] = refl;
                vals[worst] = rval;
            }
        } else if rval < vals[idx[n - 1]] {
            simplex[worst] = refl;
            vals[worst] = rval;
        } else {
            // Contract.
            let mut con = [0.0f64; 4];
            for k in 0..n {
                con[k] = centroid[k] + 0.5 * (simplex[worst][k] - centroid[k]);
            }
            let cval = f(&con);
            if cval < vals[worst] {
                simplex[worst] = con;
                vals[worst] = cval;
            } else {
                // Shrink toward best.
                for &i in &idx[1..] {
                    for k in 0..n {
                        simplex[i][k] = simplex[best][k] + 0.5 * (simplex[i][k] - simplex[best][k]);
                    }
                    vals[i] = f(&simplex[i]);
                }
            }
        }

        let spread = vals
            .iter()
            .map(|v| (v - vals[best]).abs())
            .fold(0.0f64, f64::max);
        if spread < tol {
            break;
        }
    }

    let mut best_val = f64::INFINITY;
    let mut best_p = simplex[0];
    for (i, p) in simplex.iter().enumerate() {
        if vals[i] < best_val {
            best_val = vals[i];
            best_p = *p;
        }
    }
    best_p
}

/// Keep only ridges whose Doppler fit is physically plausible: the tone must
/// map to a realistic aircraft speed, lie inside the observed sweep, show a
/// resolved sweep (both approach and recession), be explained by the model for
/// most of its samples — and fit clearly better than a straight line, which is
/// what a slow noise drift looks like.
pub fn is_plausible(r: &Ridge) -> bool {
    let Some(fit) = r.fit else {
        return false;
    };
    // Aircraft speeds. Below this the fit is trading speed against distance on a
    // partial sweep to explain slow drift, which is not a flyover.
    if !(40.0..=300.0).contains(&fit.speed_mps) {
        return false;
    }
    // A closest approach within a couple of hundred metres is not an aircraft
    // pass; the fit reaches such solutions by trading speed against distance on
    // a partial sweep.
    if fit.slant_distance_m < 200.0 {
        return false;
    }
    // Rest frequency must lie inside the observed sweep. The tolerance matches
    // the fit's own clamp, so a solution sitting on the boundary is not
    // rejected for rounding.
    if fit.f0_hz < r.min_freq_hz * 0.96 || fit.f0_hz > r.max_freq_hz * 1.04 {
        return false;
    }
    // The model must explain nearly all samples, not just the smooth middle.
    if fit.inlier_fraction < 0.8 {
        return false;
    }
    // Residual must be a small fraction of the tone's own frequency.
    if fit.rmse_hz > (0.025 * fit.f0_hz).max(2.5) {
        return false;
    }
    // An (v, d) solution is only identifiable when the tone is seen on both
    // sides of the closest approach. Partial sweeps (a flyover recorded during
    // approach or recession only) are still real detections — keep them and
    // report the partial flag as uncertainty rather than discarding evidence.
    // A Doppler sweep is an S-curve: it must beat a straight-line fit by a
    // clear margin. Constant-slope drifts (noise bands, hum wander) do not.
    if let Some(lin_rmse) = linear_rmse_hz(r) {
        if fit.rmse_hz > 0.75 * lin_rmse {
            return false;
        }
    } else {
        return false;
    }
    true
}

/// RMSE of a least-squares straight line through the ridge (Hz).
pub fn linear_rmse_hz(r: &Ridge) -> Option<f32> {
    let n = r.samples.len();
    if n < 4 {
        return None;
    }
    let nf = n as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    for s in &r.samples {
        sx += s.t_s;
        sy += s.freq_hz as f64;
    }
    let mx = sx / nf;
    let my = sy / nf;
    let (mut sxx, mut sxy) = (0.0, 0.0);
    for s in &r.samples {
        sxx += (s.t_s - mx) * (s.t_s - mx);
        sxy += (s.t_s - mx) * (s.freq_hz as f64 - my);
    }
    if sxx <= 1e-12 {
        return None;
    }
    let slope = sxy / sxx;
    let mut sse = 0.0;
    for s in &r.samples {
        let pred = my + slope * (s.t_s - mx);
        let e = pred - s.freq_hz as f64;
        sse += e * e;
    }
    Some((sse / nf).sqrt() as f32)
}
