//! Lightweight OpenStreetMap slippy map for egui: pan/zoom + rectangle region
//! selection. Tiles are fetched on a background thread; the UI thread only
//! uploads finished images, so panning and zooming never block on the network.

use std::collections::{HashMap, HashSet};
use std::sync::mpsc::{ channel, Receiver, Sender };

use egui::{Color32, Pos2, Rect, TextureHandle};

use crate::model::{Bounds, HeatmapCell};

/// Max tile requests in flight at once.
const MAX_IN_FLIGHT: usize = 24;
/// Max cached tile textures before a full clear.
const MAX_CACHE: usize = 1024;
/// Number of tiles at a given zoom (longitude axis).
const TILE: f64 = 256.0;
/// Scroll pixels per zoom level. One wheel notch delivers ~80 px through
/// egui's smoothed delta (measured), so this is "one notch = one level".
const NOTCH_PX: f64 = 80.0;
/// Heatmap raster is downsampled to at most this many texels per side.
const HEAT_MAX_TEXELS: usize = 1024;

type TileKey = (u32, u32, u32);
type TileResult = (TileKey, Vec<u8>, (usize, usize));

#[derive(Clone, Copy, PartialEq)]
pub enum Mode {
    Pan,
    Select,
}

/// A rasterised candidate field drawn under the selection rectangles.
struct Heat {
    bounds: Bounds,
    tex: TextureHandle,
    /// Survivors as a point cloud: (lat, lon, confidence). Drawn as fixed-size
    /// screen-space dots, so zooming in separates them instead of smearing the
    /// field — the QGIS read: solid at a distance, individual points up close.
    points: Vec<[f32; 3]>,
    survivors: usize,
    /// Best (lowest) residual among survivors, Hz — the legend's headline.
    best_rmse: f32,
}

pub struct MapState {
    // Centre in normalised Web-Mercator coords (0..1). Storing pixels meant
    // every zoom step had to rescale the centre, and a single bad step threw
    // the view thousands of km off.
    cx: f64,
    cy: f64,
    zoom: f32,
    pub mode: Mode,
    tile_cache: HashMap<TileKey, TextureHandle>,
    /// Requests handed to the loader and not yet answered.
    pending: HashSet<TileKey>,
    /// Tiles the server refused; never retried, so we don't hammer it.
    failed: HashSet<TileKey>,
    tile_tx: Sender<TileKey>,
    tile_rx: Receiver<TileResult>,
    rectangles: Vec<Bounds>,
    heat: Option<Heat>,
    /// Optional ground-truth point drawn as a crosshair (accuracy check).
    marker: Option<(f64, f64)>,
    /// Animated camera target: (cx, cy, zoom) and the interpolation state.
    cam_from: Option<(f64, f64, f32)>,
    cam_to: Option<(f64, f64, f32)>,
    cam_t: f32,
    /// Live peaks (lat, lon, score) drawn as HUD markers while solving.
    pub peaks: Vec<[f32; 3]>,
    /// Final answer, highlighted when the run finishes.
    pub winner: Option<(f64, f64)>,
    // interaction
    zoom_accum: f64,
    last_pointer: Option<Pos2>,
    draw_start: Option<Pos2>,
    draw_end: Option<Pos2>,
}

impl MapState {
    pub fn new(lat: f64, lon: f64, zoom: f32) -> Self {
        let (req_tx, req_rx) = channel::<TileKey>();
        let (res_tx, res_rx) = channel::<TileResult>();
        // One loader thread per map; it lives as long as the app does.
        std::thread::Builder::new()
            .name("avalon-tiles".into())
            .spawn(move || {
                while let Ok(key) = req_rx.recv() {
                    if let Some((bytes, dims)) = fetch_tile(key.0, key.1, key.2) {
                        if res_tx.send((key, bytes, dims)).is_err() {
                            break;
                        }
                    } else {
                        // Signal failure with an empty payload.
                        if res_tx.send((key, Vec::new(), (0, 0))).is_err() {
                            break;
                        }
                    }
                }
            })
            .ok();
        Self {
            cx: lon_to_world_x(lon),
            cy: lat_to_world_y(lat),
            zoom,
            mode: Mode::Pan,
            tile_cache: HashMap::new(),
            pending: HashSet::new(),
            failed: HashSet::new(),
            tile_tx: req_tx,
            tile_rx: res_rx,
            rectangles: Vec::new(),
            heat: None,
            marker: None,
            cam_from: None,
            cam_to: None,
            cam_t: 1.0,
            peaks: Vec::new(),
            winner: None,
            zoom_accum: 0.0,
            last_pointer: None,
            draw_start: None,
            draw_end: None,
        }
    }

    /// Screen pixels per unit of normalised world space.
    fn scale(&self) -> f64 {
        TILE * 2f64.powf(self.zoom as f64)
    }

    pub fn rectangles(&self) -> &[Bounds] {
        &self.rectangles
    }

    pub fn set_rectangles(&mut self, r: Vec<Bounds>) {
        self.rectangles = r;
    }

    pub fn clear_rectangles(&mut self) {
        self.rectangles.clear();
    }

    /// Rasterise the candidate grid once. Scoring spans orders of magnitude, so
    /// the colour scale is clipped to the 2nd–98th percentile of the cells that
    /// actually matched; otherwise one outlier flattens the whole map.
    ///
    /// Confidence, not error, drives the colour: a cell whose fit is near the
    /// best of the field is the one worth looking at, so it gets the hot end.
    pub fn set_heatmap(&mut self, ctx: &egui::Context, cells: &[HeatmapCell], step_deg: f64) {
        if cells.is_empty() {
            self.heat = None;
            return;
        }
        let (mut w, mut s, mut e, mut n) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for c in cells {
            w = w.min(c.lon);
            e = e.max(c.lon);
            s = s.min(c.lat);
            n = n.max(c.lat);
        }
        let step = step_deg.max(1e-4);
        let cols = (((e - w) / step).round() as usize + 1).clamp(1, 4096);
        let rows = (((n - s) / step).round() as usize + 1).clamp(1, 4096);

        // Split the grid: survivors carry a score, everything else is either an
        // elimination or unscored. Only survivors define the colour scale.
        let mut elim = vec![false; cols * rows];
        let mut vals = vec![0.0f32; cols * rows];
        for c in cells {
            let cx = (((c.lon - w) / step).round() as isize).clamp(0, cols as isize - 1) as usize;
            let cy = (((c.lat - s) / step).round() as isize).clamp(0, rows as isize - 1) as usize;
            let idx = cy * cols + cx;
            if c.score_hz > 0.0 {
                vals[idx] = vals[idx].max(c.score_hz);
            } else {
                elim[idx] = true;
            }
        }

        let mut sorted: Vec<f32> = vals.iter().copied().filter(|v| *v > 0.0).collect();
        if sorted.is_empty() {
            self.heat = None;
            return;
        }
        let survivors = sorted.len();
        let best_rmse = cells
            .iter()
            .map(|c| c.rmse_hz)
            .filter(|r| r.is_finite() && *r > 0.0)
            .fold(f32::INFINITY, f32::min);
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let lo = sorted[((sorted.len() as f64 * 0.02) as usize).min(sorted.len() - 1)];
        let hi = sorted[((sorted.len() as f64 * 0.98) as usize).min(sorted.len() - 1)];
        let span = (hi - lo).max(1e-6);

        // Survivors are a handful of cells in a 60k-cell grid, so painting one
        // texel each is invisible at the zoom a searched 2-degree box is viewed
        // at. Splat a Gaussian per survivor: the result reads as a heat field
        // with an opaque core, and overlapping survivors reinforce each other.
        const SIGMA: f32 = 1.8;
        const RADIUS: isize = 5;
        let two_sigma_sq = 2.0 * SIGMA * SIGMA;
        let mut field = vec![0.0f32; cols * rows];
        for cy in 0..rows {
            for cx in 0..cols {
                let v = vals[cy * cols + cx];
                if v <= 0.0 {
                    continue;
                }
                // score_hz is a composite where higher is better, so the
                // normalised value *is* the confidence.
                let conf = ((v - lo) / span).clamp(0.0, 1.0);
                let weight = 0.30 + 0.70 * conf;
                for dy in -RADIUS..=RADIUS {
                    let yy = cy as isize + dy;
                    if yy < 0 || yy >= rows as isize {
                        continue;
                    }
                    for dx in -RADIUS..=RADIUS {
                        let xx = cx as isize + dx;
                        if xx < 0 || xx >= cols as isize {
                            continue;
                        }
                        let d2 = (dx * dx + dy * dy) as f32;
                        field[yy as usize * cols + xx as usize] +=
                            weight * (-d2 / two_sigma_sq).exp();
                    }
                }
            }
        }

        // Max-pool down to a texture size every GPU accepts, keeping peaks.
        let ds = ((cols.max(rows) as f64) / HEAT_MAX_TEXELS as f64).ceil().max(1.0) as usize;
        let out_cols = cols.div_ceil(ds);
        let out_rows = rows.div_ceil(ds);
        // Scale to the field's own maximum: overlapping splats saturate the core
        // instantly, which collapsed the whole ramp into purple.
        let field_max = field.iter().cloned().fold(0.0f32, f32::max).max(1e-6);
        let mut pixels = Vec::with_capacity(out_cols * out_rows);
        for oy in 0..out_rows {
            for ox in 0..out_cols {
                let mut best = 0.0f32;
                let mut any_elim = false;
                for dy in 0..ds {
                    let y = oy * ds + dy;
                    if y >= rows {
                        break;
                    }
                    for dx in 0..ds {
                        let x = ox * ds + dx;
                        if x >= cols {
                            break;
                        }
                        let i = y * cols + x;
                        best = best.max(field[i]);
                        any_elim |= elim[i];
                    }
                }
                if best > 0.004 {
                    // Gamma < 1 spreads the mid-range so the ramp is legible.
                    let v = (best / field_max).clamp(0.0, 1.0).powf(0.55);
                    let c = heat_color(v);
                    let a = (150.0 + 105.0 * v.powf(0.6)) as u8;
                    pixels.push(Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), a));
                } else if any_elim {
                    // Eliminated: a faint wash so the searched box still reads as
                    // "covered", without drowning the survivors.
                    pixels.push(Color32::from_rgba_unmultiplied(58, 168, 96, 34));
                } else {
                    pixels.push(Color32::TRANSPARENT);
                }
            }
        }

        // Point cloud for the dot rendering (full resolution, not pooled).
        let mut points: Vec<[f32; 3]> = Vec::with_capacity(survivors);
        for cy in 0..rows {
            for cx in 0..cols {
                let v = vals[cy * cols + cx];
                if v <= 0.0 {
                    continue;
                }
                let conf = ((v - lo) / span).clamp(0.0, 1.0);
                let lat = s + cy as f64 * step;
                let lon = w + cx as f64 * step;
                points.push([lat as f32, lon as f32, conf]);
            }
        }

        let image = egui::ColorImage::new([out_cols, out_rows], pixels);
        let handle = ctx.load_texture("search_heatmap", image, egui::TextureOptions::LINEAR);
        // Store the covered box, expanded by half a cell so cells sit centred.
        self.heat = Some(Heat {
            bounds: Bounds {
                west: w - step / 2.0,
                south: s - step / 2.0,
                east: e + step / 2.0,
                north: n + step / 2.0,
            },
            tex: handle,
            points,
            survivors,
            best_rmse,
        });
    }

    pub fn set_marker(&mut self, p: Option<(f64, f64)>) {
        self.marker = p;
    }

    /// Current (fractional) zoom level, for the header readout.
    pub fn zoom(&self) -> f32 {
        self.zoom
    }

    /// Smoothly fly to a region: centre and zoom are eased over ~0.9 s.
    pub fn fly_to_bounds(&mut self, b: Bounds) {
        if !b.valid() {
            return;
        }
        // Fit in *Mercator* space: the geometric midpoint of a tall box and the
        // midpoint of its projected extent are not the same point, and using the
        // former pushed the answer towards the top or bottom edge of the view.
        let y_n = lat_to_world_y(b.north);
        let y_s = lat_to_world_y(b.south);
        let cy = (y_n + y_s) * 0.5;
        let cx = (lon_to_world_x(b.west) + lon_to_world_x(b.east)) * 0.5;
        let span_y = (y_n - y_s).abs().max(1e-6);
        let span_x = (lon_to_world_x(b.east) - lon_to_world_x(b.west)).abs().max(1e-6);
        // Fit the box into a ~900 px viewport: at zoom z the whole world is
        // 256*2^z px, so the largest span must occupy at most VIEWPORT px.
        const VIEWPORT: f64 = 860.0;
        let span = span_x.max(span_y);
        let zoom = (VIEWPORT / (256.0 * span)).log2().clamp(2.0, 17.0) as f32;
        self.cam_from = Some((self.cx, self.cy, self.zoom));
        self.cam_to = Some((cx, cy, zoom));
        self.cam_t = 0.0;
    }

    /// Jump straight to the fitted view (no animation), used on the first frame
    /// of a run so the map opens on the search area rather than flying in.
    pub fn snap_to_bounds(&mut self, b: Bounds) {
        if !b.valid() {
            return;
        }
        self.fly_to_bounds(b);
        if let Some((cx, cy, z)) = self.cam_to.take() {
            self.cx = cx;
            self.cy = cy;
            self.zoom = z;
        }
        self.cam_from = None;
        self.cam_t = 1.0;
    }

    /// Fly to a set of points (lat, lon) with padding — used to frame the
    /// candidate cluster once the solver has one, instead of leaving the answer
    /// as a speck inside the whole search box.
    pub fn fly_to_points(&mut self, pts: &[(f64, f64)]) {
        let (mut s, mut n) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut w, mut e) = (f64::INFINITY, f64::NEG_INFINITY);
        for (lat, lon) in pts {
            s = s.min(*lat);
            n = n.max(*lat);
            w = w.min(*lon);
            e = e.max(*lon);
        }
        if !s.is_finite() {
            return;
        }
        // Never zoom tighter than ~1.2 km across: a single point would otherwise
        // slam the camera to maximum zoom.
        let pad_lat = ((n - s) * 0.6).max(0.006);
        let pad_lon = ((e - w) * 0.6).max(0.009);
        self.fly_to_bounds(Bounds {
            west: w - pad_lon,
            south: s - pad_lat,
            east: e + pad_lon,
            north: n + pad_lat,
        });
    }

    /// Step the camera animation. Called once per frame from `show`.
    fn step_camera(&mut self, ctx: &egui::Context) {
        let Some((tx, ty, tz)) = self.cam_to else { return };
        let (fx, fy, fz) = self
            .cam_from
            .unwrap_or((self.cx, self.cy, self.zoom));
        let dt = ctx.input(|i| i.stable_dt).clamp(0.001, 0.1);
        self.cam_t = (self.cam_t + dt / 0.9).min(1.0);
        // Smoothstep: no jerk at either end.
        let e = self.cam_t * self.cam_t * (3.0 - 2.0 * self.cam_t);
        self.cx = fx + (tx - fx) * e as f64;
        self.cy = fy + (ty - fy) * e as f64;
        self.zoom = fz + (tz - fz) * e;
        if self.cam_t >= 1.0 {
            self.cam_from = None;
            self.cam_to = None;
        }
        ctx.request_repaint();
    }

    /// True while tile requests are still outstanding (drives the busy strip).
    pub fn tiles_pending(&self) -> bool {
        !self.pending.is_empty()
    }

    pub fn clear_heatmap(&mut self) {
        self.heat = None;
    }

    /// Take delivery of whatever the loader finished. Never blocks.
    fn drain_tiles(&mut self, ctx: &egui::Context) {
        let mut got = false;
        while let Ok((key, bytes, (w, h))) = self.tile_rx.try_recv() {
            self.pending.remove(&key);
            if bytes.is_empty() || w == 0 || h == 0 {
                self.failed.insert(key);
                continue;
            }
            let ci = egui::ColorImage::from_rgba_unmultiplied([w, h], &bytes);
            let tex = ctx.load_texture(
                format!("tile_{}_{}_{}", key.0, key.1, key.2),
                ci,
                egui::TextureOptions::LINEAR,
            );
            self.tile_cache.insert(key, tex);
            got = true;
        }
        if self.tile_cache.len() > MAX_CACHE {
            // Keep it simple: drop everything and let the visible set refetch.
            self.tile_cache.clear();
        }
        if got {
            ctx.request_repaint();
        }
    }

    pub fn show(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        self.drain_tiles(ctx);
        self.step_camera(ctx);

        // Header: mode toggle + hint.
        ui.horizontal(|ui| {
            if ui
                .selectable_label(self.mode == Mode::Pan, "Pan")
                .clicked()
            {
                self.mode = Mode::Pan;
            }
            if ui
                .selectable_label(self.mode == Mode::Select, "Select region")
                .clicked()
            {
                self.mode = Mode::Select;
            }
            ui.label(
                egui::RichText::new(format!("{} rect(s)", self.rectangles.len()))
                    .color(crate::theme::MUTED),
            );
        });
        ui.separator();

        let (rect, resp) = ui.allocate_exact_size(ui.available_size(), egui::Sense::click_and_drag());
        let painter = ui.painter_at(rect);

        // Zoom via scroll.
        // raw_scroll_delta, not smooth: the smoothed value is re-delivered over
        // several frames, so one wheel notch used to zoom three or four levels.
        // Use the raw wheel delta when it is non-zero (one notch = one step) and
        // fall back to the smoothed value otherwise; the smoothed value alone is
        // re-delivered over several frames, so a single notch zoomed 3-4 levels.
        // Accumulate the smoothed scroll and spend it in fixed 40 px chunks.
        // Applying the smoothed delta directly zoomed once per frame while it
        // decayed, so a single wheel notch jumped three or four levels.
        let scroll = ui.input(|i| i.smooth_scroll_delta.y) as f64;
        let cursor = resp.hover_pos().unwrap_or(rect.center());
        if resp.hovered() {
            // Drop leftover scroll when the direction flips, otherwise the first
            // notch the other way is swallowed by the previous notch's remainder
            // (zoom in then out appeared to do nothing).
            if (scroll > 0.0 && self.zoom_accum < 0.0) || (scroll < 0.0 && self.zoom_accum > 0.0) {
                self.zoom_accum = 0.0;
            }
            self.zoom_accum += scroll;
            while self.zoom_accum >= NOTCH_PX {
                self.zoom_accum -= NOTCH_PX;
                self.zoom_at_rect(cursor, rect, 1.0);
            }
            while self.zoom_accum <= -NOTCH_PX {
                self.zoom_accum += NOTCH_PX;
                self.zoom_at_rect(cursor, rect, -1.0);
            }
        } else {
            self.zoom_accum = 0.0;
        }

        // Pan / draw interaction.
        self.handle_pointer(resp, rect, ctx);

        // Draw tiles (bounded work: cached ones paint, the rest queue).
        self.draw_tiles(&painter, rect);

        // Draw selected rectangles + in-progress rectangle.
        self.draw_overlays(&painter, rect);
    }

    fn zoom_at_rect(&mut self, cursor: Pos2, rect: Rect, dz: f32) {
        let c = rect.center();
        let dx = (cursor.x - c.x) as f64;
        let dy = (cursor.y - c.y) as f64;
        let new_zoom = (self.zoom + dz).clamp(2.0, 18.0);
        if (new_zoom - self.zoom).abs() < f32::EPSILON {
            return;
        }
        // Keep the world point under the cursor pinned there: solve for the
        // centre that puts it back at the same screen position.
        let s_old = self.scale();
        let wx = self.cx + dx / s_old;
        let wy = self.cy + dy / s_old;
        self.zoom = new_zoom;
        let s_new = self.scale();
        self.cx = wx - dx / s_new;
        self.cy = wy - dy / s_new;
    }

    fn handle_pointer(&mut self, resp: egui::Response, rect: Rect, ctx: &egui::Context) {
        let (pos, down, pressed, released) = ctx.input(|i| {
            (
                i.pointer.latest_pos(),
                i.pointer.primary_down(),
                i.pointer.primary_pressed(),
                i.pointer.primary_released(),
            )
        });

        match self.mode {
            Mode::Pan => {
                if resp.dragged_by(egui::PointerButton::Primary) {
                    if let (Some(p), Some(last)) = (pos, self.last_pointer) {
                        let d = p - last;
                        let s = self.scale();
                        self.cx -= d.x as f64 / s;
                        self.cy -= d.y as f64 / s;
                    }
                    self.last_pointer = pos;
                } else {
                    self.last_pointer = None;
                }
            }
            Mode::Select => {
                if pressed && resp.hovered() {
                    self.draw_start = pos;
                    self.draw_end = pos;
                } else if down && self.draw_start.is_some() {
                    self.draw_end = pos;
                } else if released && self.draw_start.is_some() {
                    if let (Some(a), Some(b)) = (self.draw_start, self.draw_end) {
                        if a.distance(b) > 4.0 {
                            if let Some(bounds) = self.screen_rect_to_bounds(Rect::from_two_pos(a, b), rect) {
                                self.rectangles.push(bounds);
                            }
                        }
                    }
                    self.draw_start = None;
                    self.draw_end = None;
                }
            }
        }
    }

    fn screen_rect_to_bounds(&self, r: Rect, map_rect: Rect) -> Option<Bounds> {
        let c = map_rect.center();
        let scale = self.scale();
        let corners = [
            Pos2::new(r.min.x, r.min.y),
            Pos2::new(r.max.x, r.min.y),
            Pos2::new(r.max.x, r.max.y),
            Pos2::new(r.min.x, r.max.y),
        ];
        let mut west = f64::INFINITY;
        let mut south = f64::INFINITY;
        let mut east = f64::NEG_INFINITY;
        let mut north = f64::NEG_INFINITY;
        for p in corners {
            let wx = self.cx + (p.x - c.x) as f64 / scale;
            let wy = self.cy + (p.y - c.y) as f64 / scale;
            let lon = world_x_to_lon(wx);
            let lat = world_y_to_lat(wy);
            west = west.min(lon);
            east = east.max(lon);
            south = south.min(lat);
            north = north.max(lat);
        }
        let b = Bounds { west, south, east, north };
        if b.valid() {
            Some(b)
        } else {
            None
        }
    }

    fn draw_tiles(&mut self, painter: &egui::Painter, rect: Rect) {
        // Tile level is floored; the fractional part scales the whole layer, so
        // zooming is continuous instead of snapping a whole level at a time.
        let level = self.zoom.floor().clamp(0.0, 19.0) as u32;
        let f = 2f64.powf(self.zoom as f64 - level as f64); // 1.0 .. 2.0
        let tile_px = TILE * f;
        let max_tile = (1i64 << level) as i64;

        // World-space window of the viewport, in tile units at this level.
        let s = self.scale();
        let c = rect.center();
        let wx0 = self.cx - (rect.width() as f64 / 2.0) / s;
        let wx1 = self.cx + (rect.width() as f64 / 2.0) / s;
        let wy0 = self.cy - (rect.height() as f64 / 2.0) / s;
        let wy1 = self.cy + (rect.height() as f64 / 2.0) / s;
        let tx0 = (wx0 * max_tile as f64).floor() as i64;
        let tx1 = (wx1 * max_tile as f64).floor() as i64;
        let ty0 = (wy0 * max_tile as f64).floor() as i64;
        let ty1 = (wy1 * max_tile as f64).floor() as i64;

        // Queue what is missing, nearest to the centre first.
        let acx = self.cx * max_tile as f64;
        let acy = self.cy * max_tile as f64;
        let mut missing: Vec<(i64, i64, f64)> = Vec::new();
        for ty in ty0..=ty1 {
            for tx in tx0..=tx1 {
                if tx < 0 || tx >= max_tile || ty < 0 || ty >= max_tile {
                    continue;
                }
                let key = (level, tx as u32, ty as u32);
                if self.tile_cache.contains_key(&key)
                    || self.pending.contains(&key)
                    || self.failed.contains(&key)
                {
                    continue;
                }
                let d = (tx as f64 - acx).powi(2) + (ty as f64 - acy).powi(2);
                missing.push((tx, ty, d));
            }
        }
        missing.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
        for (tx, ty, _) in missing {
            if self.pending.len() >= MAX_IN_FLIGHT {
                break;
            }
            let key = (level, tx as u32, ty as u32);
            if self.tile_tx.send(key).is_ok() {
                self.pending.insert(key);
            }
        }

        for ty in ty0..=ty1 {
            for tx in tx0..=tx1 {
                if tx < 0 || tx >= max_tile || ty < 0 || ty >= max_tile {
                    continue;
                }
                let key = (level, tx as u32, ty as u32);
                let left = (tx as f64 / max_tile as f64 - self.cx) * s + c.x as f64;
                let top = (ty as f64 / max_tile as f64 - self.cy) * s + c.y as f64;
                let screen = Rect::from_min_size(
                    Pos2::new(left as f32, top as f32),
                    egui::vec2(tile_px as f32, tile_px as f32),
                );
                let uv_full = Rect::from_min_max(Pos2::new(0.0, 0.0), Pos2::new(1.0, 1.0));

                if let Some(tex) = self.tile_cache.get(&key) {
                    painter.image(tex.id(), screen, uv_full, Color32::WHITE);
                    continue;
                }
                // Not loaded yet: draw the best ancestor we already hold, scaled
                // into this tile's slot. Coarse-but-instant beats a black hole
                // that fills in chunk by chunk.
                if let Some((atex, uv)) = self.ancestor(level, tx as u32, ty as u32) {
                    painter.image(atex, screen, uv, Color32::WHITE);
                    continue;
                }
                if !self.failed.contains(&key) {
                    // Night-tile-ish slate, not black: an empty map should read
                    // as "loading", not as a rendering failure.
                    painter.rect_filled(screen, 0.0, Color32::from_rgb(16, 26, 38));
                } else {
                    painter.rect_filled(screen, 0.0, Color32::from_rgb(12, 18, 26));
                }
            }
        }
        if !self.pending.is_empty() {
            painter.ctx().request_repaint_after(std::time::Duration::from_millis(120));
        }
    }

    /// A cached lower-resolution ancestor for this tile, plus the sub-rectangle
    /// of it that this tile occupies.
    fn ancestor(&self, level: u32, tx: u32, ty: u32) -> Option<(egui::TextureId, Rect)> {
        for up in 1..=6u32 {
            if level < up {
                break;
            }
            let az = level - up;
            let ax = tx >> up;
            let ay = ty >> up;
            if let Some(tex) = self.tile_cache.get(&(az, ax, ay)) {
                let div = (1u32 << up) as f32;
                let fx = (tx - (ax << up)) as f32 / div;
                let fy = (ty - (ay << up)) as f32 / div;
                let uv = Rect::from_min_max(
                    Pos2::new(fx, fy),
                    Pos2::new(fx + 1.0 / div, fy + 1.0 / div),
                );
                return Some((tex.id(), uv));
            }
        }
        None
    }

    fn draw_overlays(&self, painter: &egui::Painter, rect: Rect) {
        // Heatmap first (under selection rectangles).
        self.draw_heatmap(painter, rect);
        let stroke = egui::Stroke::new(2.0, crate::theme::ACCENT_BRIGHT);
        for b in &self.rectangles {
            if let Some(r) = self.bounds_to_screen_rect(*b, rect) {
                painter.rect_filled(r, 0.0, crate::theme::ACCENT.gamma_multiply(0.25));
                painter.rect_stroke(r, 0.0, stroke, egui::StrokeKind::Middle);
            }
        }
        if let (Some(a), Some(b)) = (self.draw_start, self.draw_end) {
            let r = Rect::from_two_pos(a, b);
            painter.rect_stroke(r, 0.0, stroke, egui::StrokeKind::Middle);
        }
        // ---- live HUD: peaks found so far, then the winner ----
        // Real data only: each marker is a scored candidate cell from the run.
        let mut best: Option<([f32; 3], Pos2)> = None;
        for p in &self.peaks {
            if let Some(pos) = self.latlon_screen(p[0] as f64, p[1] as f64, rect) {
                if !rect.contains(pos) {
                    continue;
                }
                let conf = p[2].clamp(0.0, 1.0);
                let r = 4.0 + 10.0 * conf;
                let col = heat_color(conf);
                painter.circle_filled(pos, r * 0.35, col.gamma_multiply(0.55));
                painter.circle_stroke(pos, r, egui::Stroke::new(1.2, col));
                if best.map(|(b, _)| p[2] > b[2]).unwrap_or(true) {
                    best = Some((*p, pos));
                }
            }
        }
        if let Some((_, pos)) = best {
            // Crosshair reticle on the current leader.
            let s = egui::Stroke::new(1.0, Color32::from_white_alpha(180));
            painter.line_segment([Pos2::new(pos.x - 16.0, pos.y), Pos2::new(pos.x - 5.0, pos.y)], s);
            painter.line_segment([Pos2::new(pos.x + 5.0, pos.y), Pos2::new(pos.x + 16.0, pos.y)], s);
            painter.line_segment([Pos2::new(pos.x, pos.y - 16.0), Pos2::new(pos.x, pos.y - 5.0)], s);
            painter.line_segment([Pos2::new(pos.x, pos.y + 5.0), Pos2::new(pos.x, pos.y + 16.0)], s);
        }
        if let Some((lat, lon)) = self.winner {
            if let Some(pos) = self.latlon_screen(lat, lon, rect) {
                let t = painter.ctx().input(|i| i.time) as f32;
                let pulse = 0.5 + 0.5 * (t * 2.2).sin();
                // Brackets [ ] around the answer, sized to breathe slightly.
                let half = 26.0 + 3.0 * pulse;
                let arm = 12.0;
                let ink = Color32::from_rgb(255, 255, 255);
                let halo = Color32::from_black_alpha(200);
                for (thick, col) in [(4.0, halo), (1.8, ink)] {
                    // Left bracket
                    for (a, b) in [
                        (Pos2::new(-half, -half), Pos2::new(-half + arm, -half)),
                        (Pos2::new(-half, -half), Pos2::new(-half, -half + arm)),
                        (Pos2::new(-half, half), Pos2::new(-half + arm, half)),
                        (Pos2::new(-half, half), Pos2::new(-half, half - arm)),
                        // Right bracket
                        (Pos2::new(half, -half), Pos2::new(half - arm, -half)),
                        (Pos2::new(half, -half), Pos2::new(half, -half + arm)),
                        (Pos2::new(half, half), Pos2::new(half - arm, half)),
                        (Pos2::new(half, half), Pos2::new(half, half - arm)),
                    ] {
                        painter.line_segment(
                            [Pos2::new(pos.x + a.x, pos.y + a.y), Pos2::new(pos.x + b.x, pos.y + b.y)],
                            egui::Stroke::new(thick, col),
                        );
                    }
                }
                // Secondary dot so the exact cell is unambiguous inside the box.
                painter.circle_filled(pos, 2.5, ink);

                // Label: inverted chip (dark text on a bright plate) so it reads
                // against any basemap underneath.
                let label = "— LOCK";
                let galley = painter.layout_no_wrap(
                    label.to_string(),
                    egui::FontId::monospace(12.0),
                    Color32::from_rgb(6, 12, 20),
                );
                let pad = egui::vec2(7.0, 3.0);
                let at = Pos2::new(pos.x + half + 10.0, pos.y - (galley.size().y + pad.y * 2.0) * 0.5);
                painter.rect_filled(
                    egui::Rect::from_min_size(at, galley.size() + pad * 2.0),
                    2.0,
                    Color32::from_rgb(190, 235, 255),
                );
                painter.galley(at + pad, galley, Color32::from_rgb(6, 12, 20));
                painter.ctx().request_repaint();
            }
        }

        // Zoom readout, bottom-left: the header row can be crowded out.
        painter.text(
            Pos2::new(rect.max.x - 10.0, rect.min.y + 24.0),
            egui::Align2::RIGHT_TOP,
            format!("Z {:.1}", self.zoom),
            egui::FontId::monospace(11.0),
            Color32::from_white_alpha(200),
        );
        // Ground truth crosshair.
        if let Some((lat, lon)) = self.marker {
            if let Some(p) = self.latlon_screen(lat, lon, rect) {
                let s = egui::Stroke::new(2.0, Color32::from_rgb(255, 255, 255));
                painter.line_segment([Pos2::new(p.x - 8.0, p.y), Pos2::new(p.x + 8.0, p.y)], s);
                painter.line_segment([Pos2::new(p.x, p.y - 8.0), Pos2::new(p.x, p.y + 8.0)], s);
                painter.circle_stroke(p, 12.0, s);
                painter.text(
                    Pos2::new(p.x + 14.0, p.y - 8.0),
                    egui::Align2::LEFT_TOP,
                    "truth",
                    egui::FontId::proportional(12.0),
                    Color32::WHITE,
                );
            }
        }
    }

    fn latlon_screen(&self, lat: f64, lon: f64, map_rect: Rect) -> Option<Pos2> {
        let c = map_rect.center();
        let scale = self.scale();
        let x = (lon_to_world_x(lon) - self.cx) * scale + c.x as f64;
        let y = (lat_to_world_y(lat) - self.cy) * scale + c.y as f64;
        Some(Pos2::new(x as f32, y as f32))
    }

    fn draw_heatmap(&self, painter: &egui::Painter, rect: Rect) {
        let Some(heat) = &self.heat else { return };
        let Some(r) = self.bounds_to_screen_rect(heat.bounds, rect) else { return };

        // Dots, not a texture: every survivor is a fixed-size screen-space point,
        // so a zoomed-out view reads as a solid field and a zoomed-in view
        // resolves into individual samples (alpha lets the basemap through).
        let mut mesh = egui::Mesh::default();
        let mut drawn = 0usize;
        for p in &heat.points {
            let Some(pos) = self.latlon_screen(p[0] as f64, p[1] as f64, rect) else { continue };
            if !rect.contains(pos) {
                continue;
            }
            let conf = p[2].clamp(0.0, 1.0);
            // Small screen-space dots (1.2-2.8 px): dense enough to read as a
            // field when zoomed out, and they resolve into separate samples when
            // zoomed in instead of growing into blocks.
            let rad = 1.2 + 1.6 * conf;
            let col = heat_color(conf);
            let a = (120.0 + 115.0 * conf) as u8;
            let fill = Color32::from_rgba_unmultiplied(col.r(), col.g(), col.b(), a);
            // Square texel per point, batched into one mesh: reads as dithering
            // when zoomed out and as discrete samples when zoomed in.
            mesh.add_colored_rect(
                Rect::from_center_size(pos, egui::vec2(rad * 2.0, rad * 2.0)),
                fill,
            );
            drawn += 1;
        }
        if drawn > 0 {
            painter.add(egui::Shape::mesh(mesh));
        }

        // Frame the searched box so its extent is legible on the basemap.
        painter.rect_stroke(
            r,
            0.0,
            egui::Stroke::new(1.0, crate::theme::ACCENT.gamma_multiply(0.6)),
            egui::StrokeKind::Middle,
        );
        // Legend: how many cells are on screen and what the scale spans. Without
        // this an empty-looking map is indistinguishable from a failed render.
        let text = format!(
            "{} cells · best rmse {:.2} Hz · bright = best fit",
            heat.survivors, heat.best_rmse
        );
        let pos = Pos2::new(rect.min.x + 10.0, rect.min.y + 24.0);
        painter.text(
            pos + egui::vec2(1.0, 1.0),
            egui::Align2::LEFT_TOP,
            &text,
            egui::FontId::monospace(11.0),
            Color32::from_black_alpha(200),
        );
        painter.text(
            pos,
            egui::Align2::LEFT_TOP,
            text,
            egui::FontId::monospace(11.0),
            crate::theme::ACCENT_BRIGHT,
        );
    }

    fn bounds_to_screen_rect(&self, b: Bounds, map_rect: Rect) -> Option<Rect> {
        let c = map_rect.center();
        let scale = self.scale();
        let x0 = (lon_to_world_x(b.west) - self.cx) * scale + c.x as f64;
        let x1 = (lon_to_world_x(b.east) - self.cx) * scale + c.x as f64;
        let y0 = (lat_to_world_y(b.north) - self.cy) * scale + c.y as f64;
        let y1 = (lat_to_world_y(b.south) - self.cy) * scale + c.y as f64;
        Some(Rect::from_min_max(
            Pos2::new(x0 as f32, y0 as f32),
            Pos2::new(x1 as f32, y1 as f32),
        ))
    }
}

/// Runs on the loader thread: blocking HTTP + PNG decode, no egui involved.
/// Basemap style. The default is a night map: a light-on-light basemap under a
/// dark UI looks like a web page pasted into the window.
fn tile_style() -> String {
    std::env::var("AVALON_TILES").unwrap_or_else(|_| "osm-night".into())
}

fn fetch_tile(z: u32, x: u32, y: u32) -> Option<(Vec<u8>, (usize, usize))> {
    let style = tile_style();
    let url = match style.as_str() {
        // CARTO's dark basemap now watermarks tiles without an API key, so it is
        // opt-in only: set AVALON_TILES=night and a key on the URL if you have one.
        "night" | "dark" => format!("https://a.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png"),
        // Same tiles as below, tinted locally at runtime.
        "osm-night" => format!("https://tile.openstreetmap.org/{z}/{x}/{y}.png"),
        _ => format!("https://tile.openstreetmap.org/{z}/{x}/{y}.png"),
    };
    let resp = ureq::get(&url)
        .set("User-Agent", "avalon (geolocation research)")
        .call()
        .ok()?;
    let mut bytes = Vec::new();
    use std::io::Read;
    resp.into_reader().read_to_end(&mut bytes).ok()?;
    let img = image::load_from_memory(&bytes).ok()?.to_rgba8();
    let (w, h) = img.dimensions();
    let mut raw = img.into_raw();

    // Local night filter for the OSM fallback: invert luminance, then tint the
    // result a cold blue so it sits under the UI's cyan accents.
    if style == "osm-night" {
        for px in raw.chunks_mut(4) {
            let lum = (0.2126 * px[0] as f32 + 0.7152 * px[1] as f32 + 0.0722 * px[2] as f32) / 255.0;
            // Compress the range so roads stay visible without glare.
            let v = (1.0 - lum).powf(0.85) * 0.72;
            px[0] = (v * 150.0 + 8.0) as u8;
            px[1] = (v * 178.0 + 12.0) as u8;
            px[2] = (v * 205.0 + 22.0) as u8;
        }
    }
    Some((raw, (w as usize, h as usize)))
}

fn lon_to_world_x(lon: f64) -> f64 {
    (lon + 180.0) / 360.0
}

fn world_x_to_lon(x: f64) -> f64 {
    x * 360.0 - 180.0
}

fn lat_to_world_y(lat: f64) -> f64 {
    let lat = lat.clamp(-85.05, 85.05).to_radians();
    (1.0 - (lat.tan() + 1.0 / lat.cos()).ln() / std::f64::consts::PI) / 2.0
}

fn world_y_to_lat(y: f64) -> f64 {
    let n = std::f64::consts::PI * (1.0 - 2.0 * y.clamp(0.0, 1.0));
    n.sinh().atan().to_degrees()
}

/// Confidence colour: 0 = searched and rejected, 1 = the most likely receiver.
/// Green (low) → teal → blue → red → purple (maximum confidence), the ramp the
/// operator expects: purple is the answer.
fn heat_color(t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    // Cold scheme to sit on the blue-on-black UI: deep navy -> steel -> cyan ->
    // ice -> near-white lilac at maximum confidence.
    let stops: [(f32, (u8, u8, u8)); 5] = [
        (0.00, (18, 44, 82)),
        (0.35, (30, 92, 150)),
        (0.62, (64, 156, 206)),
        (0.85, (140, 214, 240)),
        (1.00, (226, 242, 255)),
    ];
    for i in 0..stops.len() - 1 {
        let (t0, c0) = stops[i];
        let (t1, c1) = stops[i + 1];
        if t <= t1 {
            let f = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
            let r = (c0.0 as f32 + (c1.0 as f32 - c0.0 as f32) * f) as u8;
            let g = (c0.1 as f32 + (c1.1 as f32 - c0.1 as f32) * f) as u8;
            let b = (c0.2 as f32 + (c1.2 as f32 - c0.2 as f32) * f) as u8;
            return Color32::from_rgb(r, g, b);
        }
    }
    Color32::from_rgb(226, 242, 255)
}

