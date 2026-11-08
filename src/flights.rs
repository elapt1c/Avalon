//! Historical flight acquisition: download ADSB.lol daily archives and parse
//! readsb traces into `Flight` trajectories.
//!
//! Scheme (verified live): GitHub release per UTC day in
//! `adsblol/globe_history_YYYY`, tag `v{YYYY.MM.DD}-planes-readsb-{provider}`.
//! Archive assets are a single `.tar` or split `.tar.aa`/`.tar.ab`; each member
//! named `trace_full_*.json` is a gzip'd readsb trace.

use std::io::{Read, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{anyhow, Context, Result};
use flate2::read::GzDecoder;
use serde_json::Value;

use crate::dateutil;
use crate::model::{Bounds, Flight, FlightPoint};

const GITHUB_API: &str = "https://api.github.com/repos";
/// Seconds of extra trajectory kept before/after the window (approach/recession).
const PAD_S: f64 = 180.0;
/// Minimum interval between progress callbacks (keeps the UI cheap).
const PROGRESS_INTERVAL_S: f64 = 0.15;

/// A real snapshot of download/parse progress.
#[derive(Debug, Clone, Default)]
pub struct FlightProgress {
    /// Release tag currently being processed.
    pub tag: String,
    /// Asset (archive part) currently streaming.
    pub asset: String,
    /// Bytes read so far for the current day's assets.
    pub bytes_read: u64,
    /// Total bytes of the current day's assets (from release metadata).
    pub bytes_total: u64,
    /// Archive members examined.
    pub traces_seen: usize,
    /// Flights kept so far.
    pub flights_kept: usize,
    /// Days completed.
    pub days_done: usize,
    pub days_total: usize,
}

/// Download and parse flights that enter `region` within `[t0, t1]` (unix s).
///
/// `on_progress` receives throttled, real progress; `cancel` aborts the download.
pub fn fetch_flights(
    t0: i64,
    t1: i64,
    region: &[Bounds],
    provider: &str,
    on_progress: &mut dyn FnMut(FlightProgress),
    cancel: &AtomicBool,
) -> Result<Vec<Flight>> {
    let mut flights: Vec<Flight> = Vec::new();
    let day0 = t0.div_euclid(86_400);
    let day1 = t1.div_euclid(86_400);
    let days_total = (day1 - day0 + 1).max(1) as usize;
    let mut days_done = 0usize;

    for day in day0..=day1 {
        if cancel.load(Ordering::Relaxed) {
            return Err(anyhow!("cancelled"));
        }
        let (y, m, d, _, _, _) = dateutil::unix_to_civil(day * 86_400);
        let (repo, tag) = repo_tag(y, m, d, provider);

        let mut prog = FlightProgress {
            tag: tag.clone(),
            days_done,
            days_total,
            ..Default::default()
        };
        on_progress(prog.clone());
        crate::log::info(
            "flights",
            &format!("day {}/{}: resolving {tag}", days_done + 1, days_total),
        );

        // Parsed cache first: tar+gzip for a day takes minutes.
        if let Some(cached) = read_flight_cache(&tag, region, t0 as f64, t1 as f64) {
            crate::log::info(
                "flights",
                &format!("day {}/{}: {tag} parsed cache hit ({} flights)", days_done + 1, days_total, cached.len()),
            );
            let n = cached.len();
            flights.extend(cached);
            days_done += 1;
            let mut done = prog.clone();
            done.days_done = days_done;
            done.flights_kept = flights.len();
            on_progress(done);
            let _ = n;
            continue;
        }
        let assets = list_assets(&repo, &tag)?;
        if assets.is_empty() {
            return Err(anyhow!("no archive assets for {tag}"));
        }
        let bytes_total: u64 = assets.iter().map(|a| a.2).sum();
        prog.bytes_total = bytes_total;
        prog.asset = assets[0].0.clone();
        on_progress(prog.clone());

        // Download to a local cache first. Streaming a multi-gigabyte archive
        // straight into the tar parser means one dropped connection aborts the
        // whole search; a verified local file makes retries cheap and re-runs
        // instant.
        let dir = cache_dir();
        std::fs::create_dir_all(&dir)
            .with_context(|| format!("creating cache dir {}", dir.display()))?;
        let mut downloaded: u64 = 0;
        for (idx, (name, url, size)) in assets.iter().enumerate() {
            if cancel.load(Ordering::Relaxed) {
                return Err(anyhow!("cancelled"));
            }
            let path = dir.join(name);
            prog.asset = format!("{name} ({}/{})", idx + 1, assets.len());
            let cached = std::fs::metadata(&path)
                .map(|m| m.len() == *size)
                .unwrap_or(false);
            if cached {
                crate::log::info(
                    "flights",
                    &format!("{name}: cached ({:.0} MB)", *size as f64 / 1e6),
                );
            } else {
                crate::log::info(
                    "flights",
                    &format!(
                        "{name}: downloading {:.0} MB{}",
                        *size as f64 / 1e6,
                        if *size > 100_000_000 { " (this takes minutes)" } else { "" }
                    ),
                );
            }
            if !cached {
                let bytes = downloaded;
                let mut last_pct = -1i32;
                let mut report = |read: u64, prog: &mut FlightProgress,
                                  on_progress: &mut dyn FnMut(FlightProgress)| {
                    prog.bytes_read = bytes + read;
                    if *size > 0 {
                        let pct = ((read as f64 / *size as f64) * 100.0) as i32;
                        if pct >= last_pct + 5 {
                            last_pct = pct;
                            crate::log::info(
                                "flights",
                                &format!(
                                    "  {pct}% ({:.0}/{:.0} MB) · {} traces · {} flights",
                                    read as f64 / 1e6,
                                    *size as f64 / 1e6,
                                    prog.traces_seen,
                                    prog.flights_kept
                                ),
                            );
                        }
                    }
                    on_progress(prog.clone());
                };
                ensure_asset(url, &path, *size, cancel, &mut report, &mut prog, on_progress)?;
            } else {
                downloaded += *size;
                prog.bytes_read = downloaded;
                on_progress(prog.clone());
            }
            downloaded += *size;
            prog.bytes_read = downloaded;
            on_progress(prog.clone());
        }

        // Parse from the cached files.
        let mut readers: Vec<Box<dyn Read + Send>> = Vec::with_capacity(assets.len());
        for (name, _url, _size) in &assets {
            let f = std::fs::File::open(dir.join(name))
                .with_context(|| format!("opening cached {name}"))?;
            readers.push(Box::new(f));
        }
        let chained = ChainedReader {
            readers,
            idx: 0,
        };
        let mut archive = tar::Archive::new(chained);
        let entries = archive.entries().context("reading archive")?;

        prog.asset = format!("parsing {tag}");
        let day_start = flights.len();
        let mut last_emit = std::time::Instant::now();
        for entry in entries {
            if cancel.load(Ordering::Relaxed) {
                return Err(anyhow!("cancelled"));
            }
            let mut entry = entry.context("tar entry")?;
            let name = entry.path().map(|p| p.to_string_lossy().into_owned()).unwrap_or_default();
            if !name.contains("trace_full_") || !name.ends_with(".json") {
                continue;
            }
            let mut raw = Vec::new();
            entry
                .read_to_end(&mut raw)
                .context("reading trace member")?;
            prog.traces_seen += 1;
            if let Some(flight) = parse_trace(&raw, region, t0 as f64, t1 as f64) {
                flights.push(flight);
            }
            prog.flights_kept = flights.len();
            if last_emit.elapsed().as_secs_f64() >= PROGRESS_INTERVAL_S {
                last_emit = std::time::Instant::now();
                on_progress(prog.clone());
            }
        }
        // Persist everything parsed for this day (unfiltered by region) so any
        // later search over these days loads in seconds.
        let day_flights = &flights[day_start..];
        match () {
            _ if !day_flights.is_empty() => {
                let n = day_flights.len();
                match write_flight_cache(&tag, region, t0, t1, day_flights) {
                    Ok(()) => crate::log::info(
                        "flights",
                        &format!("{tag}: cached {n} parsed flights"),
                    ),
                    Err(e) => crate::log::warn("flights", &format!("cache write failed: {e}")),
                }
            }
            _ => {}
        }
        days_done += 1;
        prog.days_done = days_done;
        on_progress(prog);
    }
    Ok(flights)
}

/// Cache path for the *parsed* flights of a release tag.
fn parsed_path(tag: &str, region: &[Bounds], t0: i64, t1: i64) -> PathBuf {
    // The cache holds flights already filtered to this region/window, so the
    // key must include both.
    let mut h: u64 = 1469598103934665603;
    for b in region {
        for v in [b.west, b.south, b.east, b.north] {
            for byte in ((v * 1e4) as i64).to_le_bytes() {
                h ^= byte as u64;
                h = h.wrapping_mul(1099511628211);
            }
        }
    }
    for v in [t0, t1] {
        for byte in v.to_le_bytes() {
            h ^= byte as u64;
            h = h.wrapping_mul(1099511628211);
        }
    }
    let dir = cache_dir()
        .parent()
        .map(|p| p.join("parsed"))
        .unwrap_or_else(|| PathBuf::from("."));
    dir.join(format!("{tag}-{h:016x}.bin"))
}

const PARSED_MAGIC: &[u8; 5] = b"AOSF1";

/// Persist parsed flights so repeat searches skip tar+gzip (minutes per day).
fn write_flight_cache(
    tag: &str,
    region: &[Bounds],
    t0: i64,
    t1: i64,
    flights: &[Flight],
) -> std::io::Result<()> {
    use std::io::BufWriter;
    let path = parsed_path(tag, region, t0, t1);
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let tmp = path.with_extension("bin.tmp");
    {
        let f = std::fs::File::create(&tmp)?;
        let mut w = BufWriter::new(f);
        w.write_all(PARSED_MAGIC)?;
        w.write_all(&(flights.len() as u32).to_le_bytes())?;
        for fl in flights {
            let mut icao = [0u8; 8];
            for (i, b) in fl.icao.as_bytes().iter().take(8).enumerate() {
                icao[i] = *b;
            }
            w.write_all(&icao)?;
            w.write_all(&(fl.points.len() as u32).to_le_bytes())?;
            for p in &fl.points {
                w.write_all(&p.t_s.to_le_bytes())?;
                w.write_all(&(p.lat as f32).to_le_bytes())?;
                w.write_all(&(p.lon as f32).to_le_bytes())?;
                w.write_all(&(p.alt_m as f32).to_le_bytes())?;
            }
        }
    }
    std::fs::rename(&tmp, &path)
}

/// Load and filter cached parsed flights for a tag.
fn read_flight_cache(
    tag: &str,
    region: &[Bounds],
    t0: f64,
    t1: f64,
) -> Option<Vec<Flight>> {
    let bytes = std::fs::read(parsed_path(tag, region, t0 as i64, t1 as i64)).ok()?;
    if bytes.len() < 9 || &bytes[..5] != PARSED_MAGIC {
        return None;
    }
    let mut pos = 5usize;
    let count = u32::from_le_bytes(bytes.get(pos..pos + 4)?.try_into().ok()?) as usize;
    pos += 4;
    let mut out = Vec::new();
    for _ in 0..count {
        let icao_raw = bytes.get(pos..pos + 8)?;
        pos += 8;
        let icao = String::from_utf8_lossy(icao_raw)
            .trim_end_matches('\0')
            .to_string();
        let n = u32::from_le_bytes(bytes.get(pos..pos + 4)?.try_into().ok()?) as usize;
        pos += 4;
        let mut points = Vec::with_capacity(n);
        let mut min_t = f64::INFINITY;
        let mut max_t = f64::NEG_INFINITY;
        let mut in_region = false;
        for _ in 0..n {
            let t = f64::from_le_bytes(bytes.get(pos..pos + 8)?.try_into().ok()?);
            pos += 8;
            let lat = f32::from_le_bytes(bytes.get(pos..pos + 4)?.try_into().ok()?) as f64;
            pos += 4;
            let lon = f32::from_le_bytes(bytes.get(pos..pos + 4)?.try_into().ok()?) as f64;
            pos += 4;
            let alt = f32::from_le_bytes(bytes.get(pos..pos + 4)?.try_into().ok()?) as f64;
            pos += 4;
            if t < min_t {
                min_t = t;
            }
            if t > max_t {
                max_t = t;
            }
            if t >= t0 - PAD_S && t <= t1 + PAD_S && region.iter().any(|b| b.contains(lat, lon)) {
                in_region = true;
            }
            points.push(FlightPoint { t_s: t, lat, lon, alt_m: alt });
        }
        if in_region && max_t >= t0 - PAD_S && min_t <= t1 + PAD_S {
            out.push(Flight { icao, points });
        }
    }
    Some(out)
}

/// Where daily archives are cached.
fn cache_dir() -> PathBuf {
    if let Ok(p) = std::env::var("AVALON_CACHE") {
        return PathBuf::from(p);
    }
    if let Ok(x) = std::env::var("XDG_CACHE_HOME") {
        return PathBuf::from(x).join("avalon/globe_history");
    }
    if let Some(home) = std::env::var_os("HOME") {
        return PathBuf::from(home).join(".cache/avalon/globe_history");
    }
    std::env::temp_dir().join("avalon/globe_history")
}

/// Download one archive part, resuming a partial `.part` file when the server
/// honours a Range request. Retries transient failures.
fn ensure_asset(
    url: &str,
    path: &PathBuf,
    size: u64,
    cancel: &AtomicBool,
    on_bytes: &mut dyn FnMut(u64, &mut FlightProgress, &mut dyn FnMut(FlightProgress)),
    prog: &mut FlightProgress,
    on_progress: &mut dyn FnMut(FlightProgress),
) -> Result<()> {
    let file_name = path
        .file_name()
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_else(|| "asset".into());
    let part = path.with_file_name(format!("{file_name}.part"));

    let mut attempts = 0;
    loop {
        attempts += 1;
        let mut have = std::fs::metadata(&part).map(|m| m.len()).unwrap_or(0);
        if have > size {
            let _ = std::fs::remove_file(&part);
            have = 0;
        }
        if have == size && size > 0 {
            std::fs::rename(&part, path)
                .with_context(|| format!("finalising {file_name}"))?;
            return Ok(());
        }
        if have > 0 {
            // Surface the resumed offset so progress does not appear to restart.
            on_bytes(have, prog, on_progress);
        }

        let attempt = (|| -> Result<bool> {
            let mut req = ureq::get(url).set("User-Agent", "avalon");
            if have > 0 {
                req = req.set("Range", &format!("bytes={have}-"));
            }
            let resp = req
                .call()
                .with_context(|| format!("download {file_name}"))?;
            // Only append when the server actually resumed.
            let resumed = have > 0 && resp.status() == 206;
            if have > 0 && !resumed {
                have = 0;
            }
            let mut file = std::fs::OpenOptions::new()
                .create(true)
                .write(true)
                .append(resumed)
                .truncate(!resumed)
                .open(&part)?;
            let mut reader = resp.into_reader();
            let mut buf = vec![0u8; 1 << 20];
            let mut read = have;
            let mut last_emit = std::time::Instant::now();
            loop {
                if cancel.load(Ordering::Relaxed) {
                    return Err(anyhow!("cancelled"));
                }
                let n = match reader.read(&mut buf) {
                    Ok(0) => break,
                    Ok(n) => n,
                    // Keep what we have; the retry resumes from here.
                    Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                    Err(e) => return Err(anyhow!(e)),
                };
                file.write_all(&buf[..n])?;
                read += n as u64;
                if last_emit.elapsed().as_secs_f64() >= PROGRESS_INTERVAL_S {
                    last_emit = std::time::Instant::now();
                    on_bytes(read, prog, on_progress);
                }
            }
            drop(file);
            Ok(true)
        })();

        match attempt {
            Ok(_) => {}
            Err(e) => {
                if cancel.load(Ordering::Relaxed) || e.to_string().contains("cancelled") {
                    return Err(anyhow!("cancelled"));
                }
                if attempts >= 4 {
                    return Err(e.context(format!(
                        "download {file_name} failed after {attempts} attempts (partial data kept for resume)"
                    )));
                }
                std::thread::sleep(std::time::Duration::from_millis(800 * attempts as u64));
                continue;
            }
        }

        let got = std::fs::metadata(&part).map(|m| m.len()).unwrap_or(0);
        if got == size {
            std::fs::rename(&part, path)
                .with_context(|| format!("finalising {file_name}"))?;
            return Ok(());
        }
        if attempts >= 4 {
            return Err(anyhow!(
                "incomplete download for {file_name}: {got} of {size} bytes"
            ));
        }
        // Truncated stream: loop and resume.
    }
}

/// Compute the GitHub repo + release tag for a UTC day.
fn repo_tag(y: i64, m: u32, d: u32, provider: &str) -> (String, String) {
    let repo = format!("adsblol/globe_history_{y}");
    let tag = format!("v{y}.{m:02}.{d:02}-planes-readsb-{provider}");
    (repo, tag)
}

/// List the archive assets (name, download URL, size) for a release tag.
fn list_assets(repo: &str, tag: &str) -> Result<Vec<(String, String, u64)>> {
    let url = format!("{GITHUB_API}/{repo}/releases/tags/{tag}");
    let body = ureq::get(&url)
        .set("User-Agent", "avalon")
        .call()
        .with_context(|| format!("release metadata for {tag}"))?
        .into_string()
        .context("reading release metadata")?;
    let v: Value = serde_json::from_str(&body).context("parsing release metadata")?;
    if let Some(msg) = v.get("message").and_then(|m| m.as_str()) {
        return Err(anyhow!("{tag}: {msg}"));
    }
    let mut out = Vec::new();
    if let Some(arr) = v.get("assets").and_then(|a| a.as_array()) {
        for a in arr {
            let name = a.get("name").and_then(|n| n.as_str()).unwrap_or("");
            if !name.ends_with(".tar") && !name.contains(".tar.") {
                continue;
            }
            let url = a.get("browser_download_url").and_then(|u| u.as_str()).unwrap_or("");
            let size = a.get("size").and_then(|s| s.as_u64()).unwrap_or(0);
            if !url.is_empty() {
                out.push((name.to_string(), url.to_string(), size));
            }
        }
    }
    out.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(out)
}

/// Parse a gzip'd readsb trace into a `Flight` if it touches the region in-window.
pub fn parse_trace(bytes: &[u8], region: &[Bounds], t0: f64, t1: f64) -> Option<Flight> {
    let decoded = if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut gz = GzDecoder::new(bytes);
        let mut out = Vec::new();
        gz.read_to_end(&mut out).ok()?;
        out
    } else {
        bytes.to_vec()
    };
    let v: Value = serde_json::from_slice(&decoded).ok()?;
    let icao = v.get("icao")?.as_str()?.to_string();
    let timestamp = v.get("timestamp")?.as_f64()?;
    let trace = v.get("trace")?.as_array()?;

    let mut in_region = false;
    let mut points: Vec<FlightPoint> = Vec::new();
    for p in trace {
        let arr = p.as_array()?;
        let dt = arr.get(0)?.as_f64()?;
        let lat = arr.get(1)?.as_f64()?;
        let lon = arr.get(2)?.as_f64()?;
        let alt_ft = arr.get(3)?.as_f64()?;
        // Prefer geometric altitude (p[10]) when present.
        let geom = arr.get(10).and_then(|x| x.as_f64()).filter(|x| x.is_finite());
        let alt_m = geom.unwrap_or(alt_ft) * 0.3048;
        let t = timestamp + dt;
        if t >= t0 - PAD_S && t <= t1 + PAD_S {
            points.push(FlightPoint { t_s: t, lat, lon, alt_m });
        }
        if t >= t0 && t <= t1 && region.iter().any(|b| b.contains(lat, lon)) {
            in_region = true;
        }
    }

    if in_region && points.len() >= 2 {
        points.sort_by(|a, b| a.t_s.partial_cmp(&b.t_s).unwrap_or(std::cmp::Ordering::Equal));
        Some(Flight { icao, points })
    } else {
        None
    }
}

struct ChainedReader {
    readers: Vec<Box<dyn Read + Send>>,
    idx: usize,
}

impl Read for ChainedReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        while self.idx < self.readers.len() {
            let n = self.readers[self.idx].read(buf)?;
            if n > 0 {
                return Ok(n);
            }
            self.idx += 1;
        }
        Ok(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    #[test]
    fn parses_readsb_trace() {
        let json = serde_json::json!({
            "icao": "ABCDEF",
            "timestamp": 1000.0,
            "trace": [
                [0.0, 39.0, -98.0, 30000.0, 250.0, 0.0, 0, null, {}, null, null],
                [1.0, 39.01, -98.0, 30000.0, 250.0, 0.0, 0, null, {}, null, null],
                [2.0, 39.02, -98.0, 30000.0, 250.0, 0.0, 0, null, {}, null, null]
            ]
        });
        let raw = serde_json::to_vec(&json).unwrap();
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(&raw).unwrap();
        let gz = enc.finish().unwrap();

        let region = vec![Bounds { west: -99.0, south: 38.0, east: -97.0, north: 39.5 }];
        let flight = parse_trace(&gz, &region, 1000.0, 1005.0).unwrap();

        assert_eq!(flight.icao, "ABCDEF");
        assert_eq!(flight.points.len(), 3);
        assert!((flight.points[0].alt_m - 30000.0 * 0.3048).abs() < 1e-3);
    }
}
