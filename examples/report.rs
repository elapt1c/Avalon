//! One-image diagnostic report: spectrogram with measured and simulated curves,
//! plus the locate heatmap with the ground-truth marker.
//!
//! Usage: report [out.png]

use std::path::Path;
use std::sync::atomic::AtomicBool;

use avalon::model::{
    AnalysisSettings, Bounds, Metadata, SearchRequest, UtcWindow,
};
use avalon::search;
use avalon::{dateutil, decode, dsp, flights, pipeline};

const CLIP: &str = "/home/tpg/Downloads/Unknown location - raw audio - colsto (128p).mp4";
const TRUTH: (f64, f64) = (40.87, -73.52);
const SR: u32 = 4000;

fn main() -> anyhow::Result<()> {
    let out = std::env::args().nth(1).unwrap_or_else(|| "/tmp/report.png".into());

    let region = Bounds {
        west: TRUTH.1 - 0.9555,
        south: TRUTH.0 - 0.7226,
        east: TRUTH.1 + 0.9555,
        north: TRUTH.0 + 0.7226,
    };
    // Same default the GUI now uses: the archive's only Δt/base-consistent pair
    // for this footage is on 2026-10-30, so cover both candidate days.
    let win = (
        dateutil::civil_to_unix(2026, 10, 30, 11, 0, 0),
        dateutil::civil_to_unix(2026, 10, 31, 23, 0, 0),
    );

    // --- analysis -------------------------------------------------------
    let settings = AnalysisSettings { target_sample_rate: SR, ..Default::default() };
    let meta = Metadata {
        region: vec![region],
        time_window: Some(UtcWindow { start_unix_s: win.0, end_unix_s: win.1 }),
        ..Default::default()
    };
    println!("analysing");
    let a = pipeline::run(Path::new(CLIP), &settings, &meta, None)?;

    // Same selection the autorun makes: plausible events, skipping event 0.
    let mut events = Vec::new();
    let mut ridges = Vec::new();
    for e in &a.events {
        if e.id == 0 {
            continue;
        }
        let Some(rep) = a.ridges.iter().find(|r| r.id == e.representative_ridge_id) else { continue };
        if !rep.plausible {
            continue;
        }
        events.push(e.clone());
        ridges.push(rep.clone());
    }
    println!("{} selected event(s)", events.len());
    if events.is_empty() {
        anyhow::bail!("no plausible events");
    }

    // --- locate ---------------------------------------------------------
    let req = SearchRequest {
        region: vec![region],
        time_window: UtcWindow { start_unix_s: win.0, end_unix_s: win.1 },
        event_ids: events.iter().map(|e| e.id).collect(),
        source_path: CLIP.to_string(),
    };
    let cancel = AtomicBool::new(false);
    let mut prog = |_s: search::Snapshot| {};
    let mut log = |m: String| println!("  {m}");
    println!("locating");
    let res = search::run_search(&req, &events, &ridges, "staging-0", &mut prog, &mut log, &cancel)?;
    let survivors = res.cells.iter().filter(|c| c.score_hz > 0.0).count();
    if let Some(b) = res.best.first() {
        let d = haversine((b.lat, b.lon), TRUTH) / 1000.0;
        println!(
            "best {:.4},{:.4} score {:.2} rmse {:.2} Hz -> {:.1} km from truth ({survivors} survivors)",
            b.lat, b.lon, b.score_hz, b.rmse_hz, d
        );
    }
    println!("top cells:");
    for c in res.cells.iter().take(10) {
        println!(
            "  {:.4},{:.4} score {:.3} rmse {:.2} Hz {:>11?} {:.1} km from truth",
            c.lat,
            c.lon,
            c.score_hz,
            c.rmse_hz,
            c.pair,
            haversine((c.lat, c.lon), TRUTH) / 1000.0
        );
    }
    // Did the two-flight pairing produce anything at all? Single-event matches
    // alone are not evidence of a location.
    let pairs: Vec<&avalon::model::HeatmapCell> =
        res.cells.iter().filter(|c| c.pair.contains('+')).collect();
    println!("{} cell(s) explained by a flight pair", pairs.len());
    for c in pairs.iter().take(8) {
        println!(
            "  {:.4},{:.4} score {:.3} rmse {:.2} Hz {:>11?} {:.1} km from truth",
            c.lat,
            c.lon,
            c.score_hz,
            c.rmse_hz,
            c.pair,
            haversine((c.lat, c.lon), TRUTH) / 1000.0
        );
    }
    // Nearest scored cell to the truth: is the right answer even in the list?
    let near = res
        .cells
        .iter()
        .min_by(|a, b| {
            haversine((a.lat, a.lon), TRUTH)
                .partial_cmp(&haversine((b.lat, b.lon), TRUTH))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|c| (c.clone(), haversine((c.lat, c.lon), TRUTH) / 1000.0));
    if let Some((c, d)) = near {
        println!(
            "truth cell {:.4},{:.4} ({:.1} km off) score {:.3} rmse {:.2} Hz {}",
            c.lat, c.lon, d, c.score_hz, c.rmse_hz, c.pair
        );
    }

    // --- canvas ---------------------------------------------------------
    let (w, h) = (1400u32, 1000u32);
    let mut img = image::RgbImage::from_pixel(w, h, image::Rgb([18, 18, 22]));
    let spec_h = 600u32;

    // Spectrogram panel for the first selected event.
    let ev = &events[0];
    let rep = &ridges[0];
    if let Some(fit) = ev.doppler {
        let (audio, _) = decode::decode_mono_at(Path::new(CLIP), SR)?;
        let half = ((ev.end_s - ev.start_s) * 0.5).clamp(45.0, 180.0);
        let a0 = (((fit.t_cpa_s - half) * SR as f64).max(0.0)) as usize;
        let b0 = (((fit.t_cpa_s + half) * SR as f64) as usize).min(audio.samples.len());
        let lo = (0.6 * fit.f0_hz * 0.9).max(15.0);
        let hi = 1.6 * fit.f0_hz;
        let spec = dsp::spectrogram::compute_spectrogram(
            &audio.samples[a0..b0],
            SR,
            &dsp::spectrogram::StftParams {
                window_s: 1.0,
                hop_s: 0.25,
                min_freq_hz: lo,
                max_freq_hz: hi,
                fft_pad: 4,
            },
        );
        let (nf, nt) = (spec.freq_axis_hz.len(), spec.time_axis_s.len());
        if nf > 0 && nt > 0 {
            let mut vals: Vec<f32> = Vec::new();
            for row in &spec.contrast {
                vals.extend_from_slice(row);
            }
            let idx = (((vals.len() as f64) * 0.995) as usize).min(vals.len() - 1);
            let ceil = *vals
                .select_nth_unstable_by(idx, |x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal))
                .1;
            let ceil = ceil.max(3.0);
            for y in 0..spec_h {
                let k = ((y as f32 / spec_h as f32) * (nf as f32 - 1.0)) as usize;
                let k = nf - 1 - k.min(nf - 1);
                for x in 0..w {
                    let f = ((x as f32 / w as f32) * (nt as f32 - 1.0)) as usize;
                    let v = spec.contrast[k][f.min(nt - 1)];
                    let t = ((v - 1.0) / (ceil - 1.0)).clamp(0.0, 1.0);
                    img.put_pixel(x, y, image::Rgb(colormap(t)));
                }
            }
            let t_off = a0 as f64 / SR as f64 + spec.time_axis_s[0];
            let t_last = spec.time_axis_s[nt - 1];
            let f_lo = spec.freq_axis_hz[0] as f64;
            let f_hi = spec.freq_axis_hz[nf - 1] as f64;
            let to_px = |t_rel: f64, f: f64| -> Option<(i32, i32)> {
                let t = t_rel - t_off;
                if t < 0.0 || t > t_last {
                    return None;
                }
                let x = (t / t_last.max(1e-6) * (w as f64 - 1.0)) as i32;
                let k = (f - f_lo) / (f_hi - f_lo).max(1e-9);
                let y = spec_h as i32 - 1 - (k * (nf as f64 - 1.0)) as i32;
                if x < 0 || x >= w as i32 || y < 0 || y >= spec_h as i32 {
                    return None;
                }
                Some((x, y))
            };
            // Measured ridge (green), high-contrast core only.
            let mut ws: Vec<f32> = rep.samples.iter().map(|s| s.weight).collect();
            ws.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
            let med = ws[ws.len() / 2].max(1.0);
            let core: Vec<&avalon::model::RidgeSample> =
                rep.samples.iter().filter(|s| s.weight >= med).collect();
            for win2 in core.windows(2) {
                if let (Some(p), Some(q)) = (to_px(win2[0].t_s, win2[0].freq_hz as f64), to_px(win2[1].t_s, win2[1].freq_hz as f64)) {
                    line(&mut img, p, q, [80, 255, 120]);
                }
            }
            // Simulated curve (cyan) from the search result.
            if let Some(curve) = res.model_curves.first() {
                let mut prev: Option<(i32, i32)> = None;
                for &(t, f) in curve {
                    if let Some(p) = to_px(t, f as f64) {
                        if let Some(q) = prev {
                            line(&mut img, q, p, [80, 220, 255]);
                        }
                        prev = Some(p);
                    }
                }
            }
        }
    }

    // Heatmap panel: bounds -> lower half.
    let map_y0 = spec_h + 20;
    let map_h = h - map_y0 - 20;
    let b = res.bounds.unwrap_or(region);
    // The searched area is green, not neutral: green is the low end of the
    // confidence ramp, so a green field with coloured spots is the legend the
    // user reads at a glance, and purple never disappears against a dark canvas.
    for y in map_y0..(map_y0 + map_h).min(h) {
        for x in 0..w {
            img.put_pixel(x, y, image::Rgb([28, 96, 56]));
        }
    }
    let maxs = res.cells.iter().map(|c| c.score_hz).fold(0.0f32, f32::max).max(1e-6);
    // Survivors are a handful of cells in a 60k-cell grid, so a one-cell mark is
    // invisible. Splat a Gaussian per survivor (radius ~3 km) to get a heat field
    // with an opaque core, matching the in-app rasteriser.
    let km_per_px = {
        let mid = 0.5 * (b.north + b.south);
        let km_w = (b.east - b.west).abs() * 111.3 * mid.to_radians().cos();
        (km_w / w as f64).max(1e-6)
    };
    let sigma = (3.0 / km_per_px / 2.5) as f32;
    let radius = (3.0 * sigma) as i32;
    let mut field = vec![0.0f32; (w * h) as usize];
    for c in &res.cells {
        if c.score_hz <= 0.0 {
            continue;
        }
        let fx = ((c.lon - b.west) / (b.east - b.west)) as f32;
        let fy = ((b.north - c.lat) / (b.north - b.south)) as f32;
        let cx = (fx * (w as f32 - 1.0)) as i32;
        let cy = map_y0 as i32 + (fy * (map_h as f32 - 1.0)) as i32;
        let t = (c.score_hz / maxs).clamp(0.0, 1.0);
        let weight = 0.35 + 0.65 * t;
        for dy in -radius..=radius {
            for dx in -radius..=radius {
                let (x, y) = (cx + dx, cy + dy);
                if x < 0 || y < map_y0 as i32 || x >= w as i32 || y >= h as i32 {
                    continue;
                }
                let d2 = (dx * dx + dy * dy) as f32;
                field[(y as u32 * w + x as u32) as usize] +=
                    weight * (-d2 / (2.0 * sigma * sigma)).exp();
            }
        }
    }
    for y in map_y0..(map_y0 + map_h).min(h) {
        for x in 0..w {
            let v = field[(y * w + x) as usize].min(1.0);
            if v <= 0.004 {
                continue;
            }
            let a = 0.47 + 0.53 * v.sqrt();
            let c = heat(v);
            let p = img.get_pixel_mut(x, y);
            for k in 0..3 {
                p[k] = (p[k] as f32 * (1.0 - a) + c[k] as f32 * a) as u8;
            }
        }
    }
    // Truth cross (white) and best cell (magenta).
    let mark = |img: &mut image::RgbImage, lat: f64, lon: f64, col: [u8; 3], r: i32| {
        let fx = ((lon - b.west) / (b.east - b.west)) as f32;
        let fy = ((b.north - lat) / (b.north - b.south)) as f32;
        let x = (fx * (w as f32 - 1.0)) as i32;
        let y = map_y0 as i32 + (fy * (map_h as f32 - 1.0)) as i32;
        for d in -r..=r {
            line(img, (x + d, y), (x + d, y), col);
            line(img, (x, y + d), (x, y + d), col);
        }
    };
    mark(&mut img, TRUTH.0, TRUTH.1, [255, 255, 255], 10);
    if let Some(bcell) = res.best.first() {
        mark(&mut img, bcell.lat, bcell.lon, [255, 80, 255], 7);
    }

    img.save(&out)?;
    println!("wrote {out}");
    let _ = std::process::Command::new("xdg-open").arg(&out).spawn();
    Ok(())
}

fn line(img: &mut image::RgbImage, a: (i32, i32), b: (i32, i32), c: [u8; 3]) {
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

fn disc(img: &mut image::RgbImage, cx: i32, cy: i32, r: i32, c: [u8; 3]) {
    for dy in -r..=r {
        for dx in -r..=r {
            if dx * dx + dy * dy > r * r {
                continue;
            }
            let (x, y) = (cx + dx, cy + dy);
            if x >= 0 && y >= 0 && (x as u32) < img.width() && (y as u32) < img.height() {
                img.put_pixel(x as u32, y as u32, image::Rgb(c));
            }
        }
    }
}

fn heat(t: f32) -> [u8; 3] {
    let stops: [(f32, [u8; 3]); 5] = [
        (0.0, [70, 200, 110]),
        (0.35, [240, 170, 50]),
        (0.65, [225, 70, 50]),
        (0.85, [170, 60, 160]),
        (1.0, [110, 40, 190]),
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
    [110, 40, 190]
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

fn haversine(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6_371_000.0;
    let (la1, lo1) = (a.0.to_radians(), a.1.to_radians());
    let (la2, lo2) = (b.0.to_radians(), b.1.to_radians());
    let h = ((la2 - la1) / 2.0).sin().powi(2)
        + la1.cos() * la2.cos() * ((lo2 - lo1) / 2.0).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}
