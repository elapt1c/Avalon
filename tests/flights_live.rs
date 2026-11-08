//! Live network check for the flight telemetry plumbing. Ignored by default:
//! it needs internet and talks to the real ADSB.lol release API.
//!
//! Run with: cargo test --release --test flights_live -- --ignored --nocapture

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use avalon::flights;
use avalon::model::Bounds;

#[test]
#[ignore = "requires network"]
fn reports_real_download_progress_and_cancels() {
    // Isolated cache so the test neither pollutes nor depends on the user cache.
    let cache = std::env::temp_dir().join(format!("aos-test-cache-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&cache);
    // SAFETY: single-threaded test setup, no other thread reads the environment.
    unsafe { std::env::set_var("AVALON_CACHE", &cache) };
    // A single UTC day well inside the archive range.
    let t0 = 86_400 * 20_200; // 2025-04-25
    let t1 = t0 + 600;
    let region = vec![Bounds {
        west: -98.0,
        south: 38.0,
        east: -97.0,
        north: 39.0,
    }];

    let cancel = Arc::new(AtomicBool::new(false));
    let cancel_cb = cancel.clone();
    let mut seen_bytes = 0u64;
    let mut seen_total = 0u64;
    let mut seen_tag = String::new();
    let mut updates = 0usize;

    let mut on_progress = |p: flights::FlightProgress| {
        updates += 1;
        seen_tag = p.tag.clone();
        seen_bytes = seen_bytes.max(p.bytes_read);
        seen_total = seen_total.max(p.bytes_total);
        // Abort once real bytes are flowing, so the test does not pull the
        // whole multi-gigabyte archive.
        if p.bytes_read > 4_000_000 {
            cancel_cb.store(true, Ordering::Relaxed);
        }
    };

    let result = flights::fetch_flights(t0, t1, &region, "staging-0", &mut on_progress, &cancel);

    eprintln!(
        "tag={seen_tag} updates={updates} bytes={seen_bytes} total={seen_total} result={result:?}",
        result = result.as_ref().map(|f| f.len()).map_err(|e| e.to_string())
    );

    assert!(updates > 0, "expected progress callbacks");
    assert!(!seen_tag.is_empty(), "expected a resolved release tag");
    assert!(
        seen_total > 100_000_000,
        "release metadata should report the archive size, got {seen_total}"
    );
    assert!(
        seen_bytes > 4_000_000,
        "expected real streamed bytes, got {seen_bytes}"
    );
    assert!(
        result.is_err(),
        "cancellation should abort the fetch once triggered"
    );

    // The partial download must be kept so a retry resumes instead of restarting.
    let kept_partial = std::fs::read_dir(&cache)
        .map(|rd| {
            rd.flatten()
                .any(|e| e.file_name().to_string_lossy().ends_with(".part")
                    && e.metadata().map(|m| m.len() > 1_000_000).unwrap_or(false))
        })
        .unwrap_or(false);
    assert!(kept_partial, "expected a retained .part file in {cache:?}");

    // Second attempt must resume from the partial rather than starting at zero.
    let cancel2 = Arc::new(AtomicBool::new(false));
    let cancel2_cb = cancel2.clone();
    let mut first_bytes = None;
    let mut resumed_seen = 0u64;
    let mut on_progress2 = |p: flights::FlightProgress| {
        if first_bytes.is_none() && p.bytes_read > 0 {
            first_bytes = Some(p.bytes_read);
        }
        resumed_seen = resumed_seen.max(p.bytes_read);
        if p.bytes_read > 20_000_000 {
            cancel2_cb.store(true, Ordering::Relaxed);
        }
    };
    let _ = flights::fetch_flights(
        t0,
        t1,
        &region,
        "staging-0",
        &mut on_progress2,
        &cancel2,
    );
    eprintln!("resume: first_bytes={first_bytes:?} max_bytes={resumed_seen}");
    assert!(
        first_bytes.unwrap_or(0) > 1_000_000,
        "expected the retry to resume from the retained partial, got {first_bytes:?}"
    );

    let _ = std::fs::remove_dir_all(&cache);
}
