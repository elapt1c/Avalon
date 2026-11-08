//! Render a PNG of the spectrogram with the measured ridge and the *simulated*
//! Doppler curve at a given candidate cell.
//!
//! Usage: overlay <out.png> [--cell LAT,LON] [--event N]
//!
//! Measured = green, simulated = cyan. If cyan does not sit on the bright
//! spectral line, the model or the trace is wrong, and you can see which.

use std::path::Path;
use std::sync::atomic::AtomicBool;

use avalon::model::{AnalysisSettings, Bounds, Metadata, UtcWindow};
use avalon::{dateutil, dsp, flights, pipeline, score};

const CLIP: &str = "/home/tpg/Downloads/Unknown location - raw audio - colsto (128p).mp4";
const TRUTH: (f64, f64) = (40.87, -73.52);

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let out = args.first().cloned().unwrap_or_else(|| "/tmp/overlay.png".into());
    let mut cell = TRUTH;
    let mut event_idx: Option<usize> = None;
    let mut i = 1;
    while i + 1 < args.len() {
        match args[i].as_str() {
            "--cell" => {
                let v: Vec<f64> = args[i + 1].split(',').filter_map(|s| s.parse().ok()).collect();
                if v.len() == 2 {
                    cell = (v[0], v[1]);
                }
            }
            "--event" => event_idx = args[i + 1].parse().ok(),
            _ => {}
        }
        i += 2;
    }

    let region = Bounds {
        west: TRUTH.1 - 0.9555,
        south: TRUTH.0 - 0.7226,
        east: TRUTH.1 + 0.9555,
        north: TRUTH.0 + 0.7226,
    };
    let win = (
        dateutil::civil_to_unix(2026, 10, 31, 11, 0, 0),
        dateutil::civil_to_unix(2026, 10, 31, 23, 0, 0),
    );

    // 1. Analysis (cached by the OS page cache; ~2 min).
    let settings = AnalysisSettings { target_sample_rate: 4000, ..Default::default() };
    let meta = Metadata {
        region: vec![region],
        time_window: Some(UtcWindow { start_unix_s: win.0, end_unix_s: win.1 }),
        ..Default::default()
    };
    println!("analysing clip");
    let a = pipeline::run(Path::new(CLIP), &settings, &meta, None)?;

    // 2. Flights (parsed cache).
    let cancel = AtomicBool::new(false);
    let mut prog = |_p: flights::FlightProgress| {};
    let fl = flights::fetch_flights(win.0, win.1, &[region], "staging-0", &mut prog, &cancel)?;
    println!("{} flights", fl.len());

    let obs = avalon::search::ecef(cell.0, cell.1, 0.0);
    let mut tracks: Vec<(String, score::Track, f64, f64)> = Vec::new(); // (icao, track, cpa, speed)
    for f in &fl {
        let mut pts: Vec<(f64, [f64; 3])> = Vec::new();
        let mut last = f64::NEG_INFINITY;
        for p in &f.points {
            if p.t_s - last < 2.0 {
                continue;
            }
            last = p.t_s;
            pts.push((p.t_s, avalon::search::ecef(p.lat, p.lon, p.alt_m)));
        }
        let Some(tr) = score::Track::from_points(&pts) else { continue };
        let Some((d, t)) = tr.cpa(obs) else { continue };
        if d > 40_000.0 {
            continue;
        }
        let v = tr
            .vel(t)
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .unwrap_or(0.0);
        tracks.push((f.icao.clone(), tr, d, v));
    }
    println!("{} tracks within 40 km", tracks.len());

    // 3. Render each selected event.
    let mut done = 0usize;
    for ev in &a.events {
        if let Some(k) = event_idx {
            if ev.id != k {
                continue;
            }
        }
        let Some(rep) = a.ridges.iter().find(|r| r.id == ev.representative_ridge_id) else {
            continue;
        };
        let Some(fit) = ev.doppler else { continue };
        // Same high-contrast core the search uses.
        let med = {
            let mut w: Vec<f32> = rep.samples.iter().map(|s| s.weight).collect();
            w.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            w[w.len() / 2].max(1.0)
        };
        let core: Vec<&avalon::model::RidgeSample> =
            rep.samples.iter().filter(|s| s.weight >= med).collect();
        let curve: Vec<(f64, f32)> = core
            .iter()
            .step_by((core.len() / 120).max(1))
            .map(|s| (s.t_s, s.freq_hz))
            .collect();

        // Best flight + parameters at this cell for this event.
        let mut best: Option<(f64, usize, f64, f64, f64)> = None; // rmse, idx, f0, base, shift
        for (ti, (_icao, tr, _d, _v)) in tracks.iter().enumerate() {
            let Some((_, t_cpa)) = tr.cpa(obs) else { continue };
            let base = t_cpa - fit.t_cpa_s;
            for &f0 in &score::f0_grid(fit.f0_hz as f64, 25) {
                for &sh in &score::shift_grid(10.0, 17) {
                    for &df in &[-3.0f64, 0.0, 3.0] {
                        if let Some(v) = score::curve_rmse_offset(
                            &curve, tr, obs, f0, base, sh, df, score::Atmos::default(), 25.0,
                        ) {
                            if best.map_or(true, |(b, _, _, _, _)| v < b) {
                                best = Some((v, ti, f0, base, sh));
                            }
                        }
                    }
                }
            }
        }
        let Some((rmse, ti, f0, base, shift)) = best else {
            println!("event {}: no candidate flight", ev.id);
            continue;
        };
        println!(
            "event {}: best {} rmse {:.2} Hz at d={:.1} km",
            ev.id, tracks[ti].0, rmse, tracks[ti].2 / 1000.0
        );

        // Spectrogram of the event window.
        let (audio, _sr) = avalon::decode::decode_mono_at(Path::new(CLIP), 4000)?;
        let sr = 4000u32;
        let half = ((ev.end_s - ev.start_s) * 0.5).clamp(45.0, 180.0);
        let a0 = (((fit.t_cpa_s - half) * sr as f64).max(0.0)) as usize;
        let b0 = (((fit.t_cpa_s + half) * sr as f64) as usize).min(audio.samples.len());
        let slice = &audio.samples[a0..b0];
        let (lo, hi) = (0.6 * fit.f0_hz * 0.9, 1.6 * fit.f0_hz);
        let spec = dsp::spectrogram::compute_spectrogram(
            slice,
            sr,
            &dsp::spectrogram::StftParams {
                window_s: 1.0,
                hop_s: 0.25,
                min_freq_hz: lo.max(15.0),
                max_freq_hz: hi,
                fft_pad: 4,
            },
        );
        let (nf, nt) = (spec.freq_axis_hz.len(), spec.time_axis_s.len());
        if nf == 0 || nt == 0 {
            continue;
        }
        let t_off = a0 as f64 / sr as f64;

        // Colourise with a robust ceiling.
        let mut vals: Vec<f32> = Vec::with_capacity(nf * nt);
        for row in &spec.contrast {
            vals.extend_from_slice(row);
        }
        let idx = (((vals.len() as f64) * 0.995) as usize).min(vals.len() - 1);
        let ceil = *vals
            .select_nth_unstable_by(idx, |x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal))
            .1;
        let ceil = ceil.max(3.0);
        let mut img = image::RgbImage::new(nt as u32, nf as u32);
        for (k, row) in spec.contrast.iter().enumerate() {
            let y = (nf - 1 - k) as u32;
            for (f, &v) in row.iter().enumerate() {
                let t = ((v - 1.0) / (ceil - 1.0)).clamp(0.0, 1.0);
                img.put_pixel(f as u32, y, image::Rgb(colormap(t)));
            }
        }

        // Map (time, freq) -> pixel.
        let to_px = |t_rel: f64, f: f64| -> Option<(i32, i32)> {
            let t = t_rel - t_off;
            if t < 0.0 || t > (nt as f64 - 1.0) * 0.25 {
                return None;
            }
            let x = (t / 0.25) as i32;
            let k = ((f - lo.max(15.0) as f64) / (hi - lo.max(15.0)) as f64 * (nf as f64 - 1.0)) as i32;
            if k < 0 || k >= nf as i32 {
                return None;
            }
            Some((x, nf as i32 - 1 - k))
        };

        // Measured ridge in green.
        for w in curve.windows(2) {
            if let (Some(p0), Some(p1)) = (to_px(w[0].0, w[0].1 as f64), to_px(w[1].0, w[1].1 as f64)) {
                draw_line(&mut img, p0, p1, [80, 255, 120]);
            }
        }
        // Simulated curve in cyan.
        let tr = &tracks[ti].1;
        let mut tau = base + curve[0].0 - shift - 2000.0 / score::C0;
        let mut prev: Option<(i32, i32)> = None;
        for &(t_rel, _) in &curve {
            let model_t = base + t_rel - shift;
            if let Some(r) = score::doppler_ratio_warm(tr, obs, model_t, &mut tau, score::Atmos::default()) {
                let f = f0 * r;
                if let Some(p) = to_px(t_rel, f) {
                    if let Some(q) = prev {
                        draw_line(&mut img, q, p, [80, 220, 255]);
                    }
                    prev = Some(p);
                }
            }
        }

        let path = if a.events.len() > 1 {
            format!("{}.{}", out.trim_end_matches(".png"), ev.id).replace("..png", ".png")
        } else {
            out.clone()
        };
        let path = if path.ends_with(".png") { path } else { format!("{path}.png") };
        img.save(&path)?;
        println!("wrote {path}");
        done += 1;
    }
    println!("{done} image(s)");
    Ok(())
}

fn draw_line(img: &mut image::RgbImage, a: (i32, i32), b: (i32, i32), c: [u8; 3]) {
    let (mut x0, mut y0) = a;
    let (x1, y1) = b;
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    loop {
        if x0 >= 0 && y0 >= 0 && (x0 as u32) < img.width() && (y0 as u32) < img.height() {
            img.put_pixel(x0 as u32, y0 as u32, image::Rgb(c));
        }
        if x0 == x1 && y0 == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x0 += sx;
        }
        if e2 <= dx {
            err += dx;
            y0 += sy;
        }
    }
}

fn colormap(t: f32) -> [u8; 3] {
    let stops: [(f32, [u8; 3]); 5] = [
        (0.0, [8, 8, 12]),
        (0.25, [60, 18, 92]),
        (0.5, [150, 30, 90]),
        (0.75, [230, 110, 40]),
        (1.0, [250, 230, 140]),
    ];
    for i in 0..stops.len() - 1 {
        let (t0, c0) = stops[i];
        let (t1, c1) = stops[i + 1];
        if t <= t1 {
            let f = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
            let mut o = [0u8; 3];
            for k in 0..3 {
                o[k] = (c0[k] as f32 + (c1[k] as f32 - c0[k] as f32) * f) as u8;
            }
            return o;
        }
    }
    [250, 230, 140]
}
