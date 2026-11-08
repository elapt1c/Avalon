//! Minimal stderr logging with timestamps, plus RAM introspection.
//!
//! The GUI swallows failures into its own log pane; when a run dies from the
//! terminal there is otherwise nothing to read. `AVALON_LOG=debug|info|warn|error`
//! selects the level (default `info`).

use std::sync::OnceLock;
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum Level {
    Error = 0,
    Warn = 1,
    Info = 2,
    Debug = 3,
}

static LEVEL: OnceLock<Level> = OnceLock::new();

pub fn init() {
    let lvl = match std::env::var("AVALON_LOG").unwrap_or_default().to_lowercase().as_str() {
        "error" => Level::Error,
        "warn" => Level::Warn,
        "debug" | "trace" => Level::Debug,
        _ => Level::Info,
    };
    let _ = LEVEL.set(lvl);
    // A panic in a worker thread is otherwise silent in a GUI process.
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        log(Level::Error, "panic", &info.to_string());
        default_hook(info);
    }));
    log(
        Level::Info,
        "startup",
        &format!(
            "logging at {:?}; {:.0} MB available",
            lvl,
            available_bytes() as f64 / 1e6
        ),
    );
}

#[inline]
pub fn enabled(l: Level) -> bool {
    *LEVEL.get_or_init(|| Level::Info) >= l
}

pub fn log(level: Level, scope: &str, msg: &str) {
    if !enabled(level) {
        return;
    }
    let tag = match level {
        Level::Error => "ERROR",
        Level::Warn => "WARN ",
        Level::Info => "INFO ",
        Level::Debug => "DEBUG",
    };
    eprintln!("[{}] {} {}: {}", timestamp(), tag, scope, msg);
}

pub fn info(scope: &str, msg: &str) {
    log(Level::Info, scope, msg);
}
pub fn warn(scope: &str, msg: &str) {
    log(Level::Warn, scope, msg);
}
pub fn error(scope: &str, msg: &str) {
    log(Level::Error, scope, msg);
}
pub fn debug(scope: &str, msg: &str) {
    log(Level::Debug, scope, msg);
}

fn timestamp() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let (_, _, _, hh, mm, ss) = crate::dateutil::unix_to_civil(secs);
    format!("{hh:02}:{mm:02}:{ss:02}")
}

/// Bytes the OS says are available without swapping (`MemAvailable`).
pub fn available_bytes() -> u64 {
    if let Ok(s) = std::fs::read_to_string("/proc/meminfo") {
        for line in s.lines() {
            if let Some(rest) = line.strip_prefix("MemAvailable:") {
                let kb: u64 = rest
                    .split_whitespace()
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
                if kb > 0 {
                    return kb * 1024;
                }
            }
        }
    }
    512 * 1024 * 1024
}

/// Budget for one stage: half of what is currently available, with a floor so a
/// tiny machine still makes progress in small steps.
pub fn budget_bytes() -> u64 {
    (available_bytes() / 2).max(64 * 1024 * 1024)
}

/// Resident set size in bytes (Linux `/proc/self/statm`).
pub fn rss_bytes() -> u64 {
    if let Ok(s) = std::fs::read_to_string("/proc/self/statm") {
        if let Some(pages) = s.split_whitespace().nth(1).and_then(|v| v.parse::<u64>().ok()) {
            return pages * 4096;
        }
    }
    0
}
