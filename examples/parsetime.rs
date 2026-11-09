//! Serial vs parallel trace parsing, on real archive members.
//!
//! Usage: parsetime [path.tar] [members]
//! Prints wall time for parsing the same members one-at-a-time and batched on
//! the rayon pool, which is exactly the change made in flights.rs.

use std::fs::File;
use std::io::Read;
use std::time::Instant;

use avalon::flights::parse_trace;
use avalon::model::Bounds;

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "/root/projects/aosint-cache/v2026.09.08-planes-readsb-staging-0.tar.aa".into());
    let want: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(400);

    let region = vec![Bounds { west: -74.0, south: 40.4, east: -73.0, north: 41.4 }];
    let (t0, t1) = (1788890400.0, 1788910200.0);

    // Pull `want` members into memory once, so the measurement is parse-only.
    let mut members: Vec<Vec<u8>> = Vec::with_capacity(want);
    {
        let mut archive = tar::Archive::new(File::open(&path)?);
        for entry in archive.entries()? {
            let mut entry = entry?;
            let name = entry.path().map(|p| p.to_string_lossy().into_owned()).unwrap_or_default();
            if !name.contains("trace_full_") || !name.ends_with(".json") {
                continue;
            }
            let mut raw = Vec::new();
            entry.read_to_end(&mut raw)?;
            members.push(raw);
            if members.len() >= want {
                break;
            }
        }
    }
    let bytes: usize = members.iter().map(|m| m.len()).sum();
    println!("{want} members, {:.1} MB compressed", bytes as f64 / 1e6);

    let t = Instant::now();
    let mut kept = 0usize;
    for m in &members {
        if parse_trace(m, &region, t0, t1).is_some() {
            kept += 1;
        }
    }
    let serial = t.elapsed().as_secs_f64();

    use rayon::prelude::*;
    let t = Instant::now();
    let n_par: usize = members
        .par_iter()
        .filter_map(|m| parse_trace(m, &region, t0, t1))
        .count();
    let parallel = t.elapsed().as_secs_f64();

    println!("serial   {serial:6.2} s  ({:.0}/s, kept {kept})", want as f64 / serial);
    println!(
        "parallel {parallel:6.2} s  ({:.0}/s, kept {n_par})  -> {:.1}x on {} threads",
        want as f64 / parallel,
        serial / parallel.max(1e-9),
        rayon::current_num_threads()
    );
    Ok(())
}
