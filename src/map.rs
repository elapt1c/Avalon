//! Lightweight OpenStreetMap slippy map for egui: pan/zoom + rectangle region
//! selection. Fetches raster tiles synchronously (cached, capped per frame).

use std::collections::HashMap;

use egui::{Color32, Pos2, Rect, TextureHandle};

use crate::model::{Bounds, HeatmapCell};

/// Max new tiles fetched per frame (avoids long network stalls).
const MAX_FETCH_PER_FRAME: usize = 6;
/// Max cached tile textures before a full clear.
const MAX_CACHE: usize = 1024;
/// Number of tiles at a given zoom (longitude axis).
const TILE: f64 = 256.0;

#[derive(Clone, Copy, PartialEq)]
pub enum Mode {
    Pan,
    Select,
}

pub struct MapState {
    // Center in world-pixel coordinates at the current zoom.
    center_px: f64,
    center_py: f64,
    zoom: f32,
    pub mode: Mode,
    tile_cache: HashMap<(u32, u32, u32), TextureHandle>,
    rectangles: Vec<Bounds>,
    /// Heatmap raster + the box it covers. Drawn as one textured quad:
    /// per-cell rectangles overlap and alias into a moiré pattern.
    heatmap: Option<(Bounds, TextureHandle)>,
    /// Optional ground-truth point drawn as a crosshair (accuracy check).
    marker: Option<(f64, f64)>,
    // interaction
    last_pointer: Option<Pos2>,
    draw_start: Option<Pos2>,
    draw_end: Option<Pos2>,
    fetched_this_frame: usize,
}

impl MapState {
    pub fn new(lat: f64, lon: f64, zoom: f32) -> Self {
        let scale = TILE * 2f64.powi(zoom as i32);
        Self {
            center_px: lon_to_world_x(lon) * scale,
            center_py: lat_to_world_y(lat) * scale,
            zoom,
            mode: Mode::Pan,
            tile_cache: HashMap::new(),
            rectangles: Vec::new(),
            heatmap: None,
            marker: None,
            last_pointer: None,
            draw_start: None,
            draw_end: None,
            fetched_this_frame: 0,
        }
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
    /// the colour scale is normalised to the 2nd–98th percentile of the cells
    /// that actually matched; otherwise one outlier flattens the whole map.
    pub fn set_heatmap(&mut self, ctx: &egui::Context, cells: &[HeatmapCell], step_deg: f64) {
        if cells.is_empty() {
            self.heatmap = None;
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
        let cols = (((e - w) / step).round() as usize + 1).clamp(1, 2048);
        let rows = (((n - s) / step).round() as usize + 1).clamp(1, 2048);

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
            self.heatmap = None;
            return;
        }
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
                let t = ((v - lo) / span).clamp(0.0, 1.0);
                let weight = 0.35 + 0.65 * t;
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

        let mut pixels = Vec::with_capacity(cols * rows);
        for i in 0..cols * rows {
            let v = field[i].min(1.0);
            if v > 0.004 {
                // Opaque core, fading skirt: the confidence ramp already encodes
                // the value, so alpha only controls how solid the patch reads.
                let c = heat_color(v);
                let a = (120.0 + 135.0 * v.powf(0.5)) as u8;
                pixels.push(Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), a));
            } else if elim[i] {
                // Eliminated: green enough to read as "searched, nothing here"
                // over the whole box, translucent enough to keep roads legible.
                pixels.push(Color32::from_rgba_unmultiplied(58, 168, 96, 104));
            } else {
                pixels.push(Color32::TRANSPARENT);
            }
        }

        let image = egui::ColorImage::new([cols, rows], pixels);
        let handle = ctx.load_texture("search_heatmap", image, egui::TextureOptions::LINEAR);
        // Store the covered box, expanded by half a cell so cells sit centred.
        self.heatmap = Some((
            Bounds {
                west: w - step / 2.0,
                south: s - step / 2.0,
                east: e + step / 2.0,
                north: n + step / 2.0,
            },
            handle,
        ));
    }

    pub fn set_marker(&mut self, p: Option<(f64, f64)>) {
        self.marker = p;
    }

    pub fn clear_heatmap(&mut self) {
        self.heatmap = None;
    }

    pub fn show(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        self.fetched_this_frame = 0;

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
        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if resp.hovered() && scroll != 0.0 {
            self.zoom_at_rect(resp.hover_pos().unwrap_or(rect.center()), rect, scroll > 0.0);
        }

        // Pan / draw interaction.
        self.handle_pointer(resp, rect, ctx);

        // Draw tiles.
        self.draw_tiles(ctx, &painter, rect);

        // Draw selected rectangles + in-progress rectangle.
        self.draw_overlays(&painter, rect);
    }

    fn zoom_at_rect(&mut self, cursor: Pos2, rect: Rect, zoom_in: bool) {
        let c = rect.center();
        let dx = (cursor.x - c.x) as f64;
        let dy = (cursor.y - c.y) as f64;
        let cwx = self.center_px + dx;
        let cwy = self.center_py + dy;
        let new_zoom = if zoom_in {
            (self.zoom + 1.0).min(18.0)
        } else {
            (self.zoom - 1.0).max(2.0)
        };
        let f = 2f64.powi(new_zoom as i32 - self.zoom as i32);
        self.center_px = cwx * f - dx;
        self.center_py = cwy * f - dy;
        self.zoom = new_zoom;
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
                        self.center_px -= d.x as f64;
                        self.center_py -= d.y as f64;
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
        let scale = TILE * 2f64.powi(self.zoom as i32);
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
            let wx = self.center_px + (p.x - c.x) as f64;
            let wy = self.center_py + (p.y - c.y) as f64;
            let lon = world_x_to_lon(wx / scale);
            let lat = world_y_to_lat(wy / scale);
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

    fn draw_tiles(&mut self, ctx: &egui::Context, painter: &egui::Painter, rect: Rect) {
        let z = self.zoom.floor() as i32;
        let z = z.clamp(0, 19) as u32;
        let max_tile = (1i64 << z) as i64;
        let c = rect.center();

        let wpx0 = self.center_px - rect.width() as f64 / 2.0;
        let wpx1 = self.center_px + rect.width() as f64 / 2.0;
        let wpy0 = self.center_py - rect.height() as f64 / 2.0;
        let wpy1 = self.center_py + rect.height() as f64 / 2.0;

        let tx0 = (wpx0 / TILE).floor() as i64;
        let tx1 = (wpx1 / TILE).floor() as i64;
        let ty0 = (wpy0 / TILE).floor() as i64;
        let ty1 = (wpy1 / TILE).floor() as i64;

        for ty in ty0..=ty1 {
            for tx in tx0..=tx1 {
                if tx < 0 || tx >= max_tile || ty < 0 || ty >= max_tile {
                    continue;
                }
                let key = (z, tx as u32, ty as u32);
                let tex_id = match self.tile_cache.get(&key) {
                    Some(t) => t.id(),
                    None => {
                        if self.fetched_this_frame >= MAX_FETCH_PER_FRAME {
                            continue;
                        }
                        self.fetched_this_frame += 1;
                        match fetch_tile(ctx, z, tx as u32, ty as u32) {
                            Some(t) => {
                                self.tile_cache.insert(key, t);
                                if self.tile_cache.len() > MAX_CACHE {
                                    self.tile_cache.clear();
                                }
                                self.tile_cache.get(&key).unwrap().id()
                            }
                            None => continue,
                        }
                    }
                };

                let left = tx as f64 * TILE - self.center_px + c.x as f64;
                let right = (tx as f64 + 1.0) * TILE - self.center_px + c.x as f64;
                let top = ty as f64 * TILE - self.center_py + c.y as f64;
                let bottom = (ty as f64 + 1.0) * TILE - self.center_py + c.y as f64;
                let screen = Rect::from_min_max(
                    Pos2::new(left as f32, top as f32),
                    Pos2::new(right as f32, bottom as f32),
                );
                let uv = Rect::from_min_max(Pos2::new(0.0, 0.0), Pos2::new(1.0, 1.0));
                painter.image(tex_id, screen, uv, Color32::WHITE);
            }
        }
    }

    fn draw_overlays(&self, painter: &egui::Painter, rect: Rect) {
        // Heatmap first (under selection rectangles).
        if self.heatmap.is_some() {
            self.draw_heatmap(painter, rect);
        }
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
        let scale = TILE * 2f64.powi(self.zoom as i32);
        let x = lon_to_world_x(lon) * scale - self.center_px + c.x as f64;
        let y = lat_to_world_y(lat) * scale - self.center_py + c.y as f64;
        Some(Pos2::new(x as f32, y as f32))
    }

    fn draw_heatmap(&self, painter: &egui::Painter, rect: Rect) {
        let Some((bounds, tex)) = &self.heatmap else { return };
        let Some(r) = self.bounds_to_screen_rect(*bounds, rect) else { return };
        let uv = Rect::from_min_max(Pos2::new(0.0, 0.0), Pos2::new(1.0, 1.0));
        // Row 0 of the texture is south, so flip V: egui textures start at top.
        let uv = Rect::from_min_max(Pos2::new(uv.min.x, 1.0), Pos2::new(uv.max.x, 0.0));
        painter.image(tex.id(), r, uv, Color32::WHITE);
        // Frame the searched box so its extent is legible on the basemap.
        painter.rect_stroke(
            r,
            0.0,
            egui::Stroke::new(1.0, crate::theme::ACCENT.gamma_multiply(0.6)),
            egui::StrokeKind::Middle,
        );
    }

    fn bounds_to_screen_rect(&self, b: Bounds, map_rect: Rect) -> Option<Rect> {
        let c = map_rect.center();
        let scale = TILE * 2f64.powi(self.zoom as i32);
        let x0 = lon_to_world_x(b.west) * scale - self.center_px + c.x as f64;
        let x1 = lon_to_world_x(b.east) * scale - self.center_px + c.x as f64;
        let y0 = lat_to_world_y(b.north) * scale - self.center_py + c.y as f64;
        let y1 = lat_to_world_y(b.south) * scale - self.center_py + c.y as f64;
        Some(Rect::from_min_max(
            Pos2::new(x0 as f32, y0 as f32),
            Pos2::new(x1 as f32, y1 as f32),
        ))
    }
}

fn fetch_tile(ctx: &egui::Context, z: u32, x: u32, y: u32) -> Option<TextureHandle> {
    // tile.openstreetmap.org usage policy: fine for light interactive use.
    let url = format!("https://tile.openstreetmap.org/{z}/{x}/{y}.png");
    let resp = ureq::get(&url).set("User-Agent", "avalon").call().ok()?;
    let mut bytes = Vec::new();
    use std::io::Read;
    resp.into_reader().read_to_end(&mut bytes).ok()?;
    let img = image::load_from_memory(&bytes).ok()?.to_rgba8();
    let (w, h) = img.dimensions();
    let ci = egui::ColorImage::from_rgba_unmultiplied([w as usize, h as usize], img.as_raw());
    Some(ctx.load_texture(
        format!("tile_{z}_{x}_{y}"),
        ci,
        egui::TextureOptions::LINEAR,
    ))
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

/// Heatmap color: t=0 (best) bright, t=1 (worst) dark.
fn heat_color(t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let stops: [(f32, (u8, u8, u8)); 5] = [
        (0.00, (70, 200, 110)),  // low confidence: green
        (0.35, (240, 170, 50)),  // orange
        (0.65, (225, 70, 50)),   // red
        (0.85, (170, 60, 160)),  // magenta
        (1.00, (110, 40, 190)),  // highest: purple
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
    Color32::from_rgb(70, 25, 90)
}
