//! Terrain elevation from Terrarium tiles (AWS open data, ~30 m at z12).
//!
//! Height is stored in a normal PNG as `(R*256 + G + B/256) - 32768` metres, so
//! no GeoTIFF/GRIB stack is needed. Tiles are cached on disk, and decoded tiles
//! are kept in memory for the duration of a search.
//!
//! Why this matters: a receiver assumed to be at sea level injects up to
//! ~2000 m of systematic error into the slant-distance geometry across CONUS,
//! which biases both the CPA range and the time the sound arrives.

use std::collections::HashMap;
use std::io::Read;
use std::path::PathBuf;
use std::sync::Arc;

/// Web-mercator zoom. 12 ≈ 38 m/px at the equator, ~27 m at 45°.
const ZOOM: u32 = 12;
/// Tile size in pixels.
const TS: usize = 256;
/// Decoded tiles kept in memory before the cache is dropped.
const MAX_TILES: usize = 600;

pub struct Terrain {
    zoom: u32,
    tiles: HashMap<(u32, u32), Option<Arc<Vec<f32>>>>,
}

impl Terrain {
    /// Zoom that keeps a region's tile count manageable. One request per tile
    /// dominates runtime, so a 9x6 degree box must not need ~9 000 of them.
    pub fn for_area(span_deg: (f64, f64)) -> Self {
        let area = (span_deg.0 * span_deg.1).max(0.01);
        let mut zoom = 12u32;
        for z in [12u32, 11, 10, 9, 8] {
            let per_deg = 2f64.powi(z as i32) / 360.0;
            let tiles = area * per_deg * per_deg;
            if tiles <= 2_000.0 {
                zoom = z;
                break;
            }
        }
        Self::new(zoom)
    }

    pub fn zoom(&self) -> u32 {
        self.zoom
    }

    /// Ensure every tile covering `points` exists on disk, downloading missing
    /// ones in parallel. Sequential fetching is ~0.45 s per tile.
    pub fn prefetch(&self, points: &[(f64, f64)]) -> (usize, usize) {
        let mut keys: Vec<(u32, u32)> = points
            .iter()
            .filter(|(lat, _)| (-85.0..=85.0).contains(lat))
            .map(|(lat, lon)| {
                let n = TS as f64 * 2f64.powi(self.zoom as i32);
                let px = (lon + 180.0) / 360.0 * n;
                let latr = lat.to_radians();
                let py = (1.0 - (latr.tan() + 1.0 / latr.cos()).ln() / std::f64::consts::PI) / 2.0 * n;
                (
                    (px as i64 / TS as i64) as u32,
                    (py as i64 / TS as i64) as u32,
                )
            })
            .collect();
        keys.sort_unstable();
        keys.dedup();
        let total = keys.len();
        let missing: Vec<(u32, u32)> = keys
            .into_iter()
            .filter(|(x, y)| !disk_path(self.zoom, *x, *y).is_file())
            .collect();
        let fetched = missing.len();
        if fetched > 0 {
            use rayon::prelude::*;
            use std::sync::atomic::{AtomicUsize, Ordering};
            let z = self.zoom;
            let done = AtomicUsize::new(0);
            let step = (fetched / 10).max(25);
            crate::log::info(
                "dem",
                &format!("downloading {fetched} tiles (parallel)"),
            );
            missing.par_iter().for_each(|(x, y)| {
                let _ = fetch_to_disk(z, *x, *y);
                let n = done.fetch_add(1, Ordering::Relaxed) + 1;
                if n % step == 0 || n == fetched {
                    crate::log::info(
                        "dem",
                        &format!("  {n}/{fetched} tiles"),
                    );
                }
            });
        }
        crate::log::info(
            "dem",
            &format!("zoom {} · {total} tiles · {fetched} downloaded", self.zoom),
        );
        (total, fetched)
    }
}

impl Default for Terrain {
    fn default() -> Self {
        Self::new(ZOOM)
    }
}

impl Terrain {
    pub fn new(zoom: u32) -> Self {
        Self {
            zoom: zoom.clamp(6, 14),
            tiles: HashMap::new(),
        }
    }

    /// Terrain height in metres (WGS84 surface), or `None` when the tile could
    /// not be fetched (offline). Callers fall back to sea level.
    pub fn elevation(&mut self, lat: f64, lon: f64) -> Option<f32> {
        if !(-85.0..=85.0).contains(&lat) {
            return None;
        }
        let (px, py) = self.world_px(lat, lon);
        // Bilinear across the four surrounding texels.
        let x0 = px.floor();
        let y0 = py.floor();
        let fx = (px - x0) as f32;
        let fy = (py - y0) as f32;
        let h00 = self.texel(x0 as i64, y0 as i64)?;
        let h10 = self.texel(x0 as i64 + 1, y0 as i64)?;
        let h01 = self.texel(x0 as i64, y0 as i64 + 1)?;
        let h11 = self.texel(x0 as i64 + 1, y0 as i64 + 1)?;
        let top = h00 + (h10 - h00) * fx;
        let bot = h01 + (h11 - h01) * fx;
        Some(top + (bot - top) * fy)
    }

    /// Elevations for many points in one pass.
    ///
    /// Sampling cell by cell thrashes the tile cache: a few-hundred-kilometre
    /// region needs thousands of tiles against a 600-tile cache, so each new
    /// point re-decodes PNGs. Grouping by tile loads each tile exactly once.
    pub fn elevations(&mut self, points: &[(f64, f64)]) -> Vec<f32> {
        let mut out = vec![0.0f32; points.len()];
        let mut by_tile: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for (i, (lat, lon)) in points.iter().enumerate() {
            if !(-85.0..=85.0).contains(lat) {
                continue;
            }
            let (px, py) = self.world_px(*lat, *lon);
            let key = (
                (px as i64 / TS as i64) as u32,
                (py as i64 / TS as i64) as u32,
            );
            by_tile.entry(key).or_default().push(i);
        }
        let tiles = by_tile.len();
        let z = self.zoom;
        for (n, ((tx, ty), idxs)) in by_tile.into_iter().enumerate() {
            if n % 200 == 0 {
                crate::log::debug("dem", &format!("tile {n}/{tiles}"));
            }
            // Read straight from disk and drop after use: keeping thousands of
            // decoded tiles resident costs hundreds of megabytes.
            let Some(tile) = read_tile(z, tx, ty, false) else { continue };
            for i in idxs {
                let (px, py) = self.world_px(points[i].0, points[i].1);
                let lx = (px as i64 % TS as i64) as f32;
                let ly = (py as i64 % TS as i64) as f32;
                let x0 = lx.floor().clamp(0.0, TS as f32 - 1.0) as usize;
                let y0 = ly.floor().clamp(0.0, TS as f32 - 1.0) as usize;
                let x1 = (x0 + 1).min(TS - 1);
                let y1 = (y0 + 1).min(TS - 1);
                let fx = lx - x0 as f32;
                let fy = ly - y0 as f32;
                let h00 = tile[y0 * TS + x0];
                let h10 = tile[y0 * TS + x1];
                let h01 = tile[y1 * TS + x0];
                let h11 = tile[y1 * TS + x1];
                let top = h00 + (h10 - h00) * fx;
                let bot = h01 + (h11 - h01) * fx;
                out[i] = top + (bot - top) * fy;
            }
        }
        crate::log::info(
            "dem",
            &format!("sampled {} points over {} tiles", points.len(), tiles),
        );
        out
    }

    fn world_px(&self, lat: f64, lon: f64) -> (f64, f64) {
        let n = TS as f64 * 2f64.powi(self.zoom as i32);
        let x = (lon + 180.0) / 360.0 * n;
        let lat = lat.clamp(-85.0, 85.0).to_radians();
        let y = (1.0 - (lat.tan() + 1.0 / lat.cos()).ln() / std::f64::consts::PI) / 2.0 * n;
        (x, y)
    }

    fn texel(&mut self, gx: i64, gy: i64) -> Option<f32> {
        // Clamp in *pixel* space: the world is 256 * 2^zoom pixels across.
        let max_px = (TS as i64) * (1i64 << self.zoom) - 1;
        let gy = gy.clamp(0, max_px);
        let gx = gx.rem_euclid(max_px + 1);
        let (tx, ty) = ((gx / TS as i64) as u32, (gy / TS as i64) as u32);
        let (lx, ly) = ((gx % TS as i64) as usize, (gy % TS as i64) as usize);
        let tile = self.tile(tx, ty)?;
        tile.get(ly * TS + lx).copied()
    }

    fn tile(&mut self, tx: u32, ty: u32) -> Option<Arc<Vec<f32>>> {
        if let Some(t) = self.tiles.get(&(tx, ty)) {
            return t.clone();
        }
        let tile = load_tile(self.zoom, tx, ty);
        if self.tiles.len() >= MAX_TILES {
            self.tiles.clear();
        }
        self.tiles.insert((tx, ty), tile.clone());
        tile
    }
}

fn read_tile(z: u32, x: u32, y: u32, allow_fetch: bool) -> Option<Arc<Vec<f32>>> {
    let path = disk_path(z, x, y);
    let bytes = match std::fs::read(&path) {
        Ok(b) => b,
        Err(_) => {
            if !allow_fetch {
                return None;
            }
            fetch_to_disk(z, x, y)?
        }
    };
    let img = image::load_from_memory(&bytes).ok()?.to_rgba8();
    if img.width() as usize != TS || img.height() as usize != TS {
        return None;
    }
    let mut heights = Vec::with_capacity(TS * TS);
    for p in img.pixels() {
        let r = p[0] as f32;
        let g = p[1] as f32;
        let b = p[2] as f32;
        heights.push(r * 256.0 + g + b / 256.0 - 32768.0);
    }
    Some(Arc::new(heights))
}

/// Download one tile to the disk cache and return its bytes.
fn fetch_to_disk(z: u32, x: u32, y: u32) -> Option<Vec<u8>> {
    let url = format!("https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png");
    let resp = ureq::get(&url)
        .set("User-Agent", "avalon")
        // Without a timeout a stalled connection blocks a worker forever.
        .timeout(std::time::Duration::from_secs(45))
        .call()
        .ok()?;
    let mut buf = Vec::new();
    resp.into_reader().read_to_end(&mut buf).ok()?;
    let path = disk_path(z, x, y);
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let _ = std::fs::write(&path, &buf);
    Some(buf)
}

fn load_tile(z: u32, x: u32, y: u32) -> Option<Arc<Vec<f32>>> {
    read_tile(z, x, y, true)
}

fn disk_path(z: u32, x: u32, y: u32) -> PathBuf {
    let base = match std::env::var("AVALON_CACHE") {
        Ok(p) => PathBuf::from(p),
        Err(_) => match std::env::var("XDG_CACHE_HOME") {
            Ok(p) => PathBuf::from(p).join("avalon"),
            Err(_) => match std::env::var_os("HOME") {
                Some(h) => PathBuf::from(h).join(".cache/avalon"),
                None => std::env::temp_dir().join("avalon"),
            },
        },
    };
    base.join("dem").join(z.to_string()).join(x.to_string()).join(format!("{y}.png"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode_offset_math() {
        // Terrarium encoding: h = (R*256 + G + B/256) - 32768
        let (r, g, b) = (128u16, 0u16, 0u16); // 32768 -> 0 m
        let h = r as f32 * 256.0 + g as f32 + b as f32 / 256.0 - 32768.0;
        assert!((h - 0.0).abs() < 0.01);
    }

    #[test]
    #[ignore = "requires network"]
    fn reads_known_elevations() {
        let mut dem = Terrain::default();
        // Denver ~1609 m, Death Valley ~-86 m, Dead Sea shore ~-430 m.
        for (lat, lon, name, target) in [
            (39.7392, -104.9903, "Denver", 1609.0),
            (36.5323, -116.9325, "Death Valley", -86.0),
            (31.5, 35.5, "Dead Sea", -430.0),
        ] {
            let h = dem.elevation(lat, lon).expect("tile");
            eprintln!("{name}: {h:.0} m");
            assert!((h - target).abs() < 250.0, "{name} got {h}");
        }
    }
}
