//! Headless end-to-end benchmark: analyse the clip once, then run the locate
//! stage over one or more regions and report distance to a known truth point.
//!
//! Usage:
//!   bench [--clip PATH] [--analysis CACHE.json] [--time UNIX0:UNIX1]
//!         [--truth LAT,LON] [--jsonl OUT] LABEL:W,S,E,N [...]
//!
//! The analysis result is cached as JSON so repeated region sweeps skip the
//! (expensive) decode + ridge stage.

use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicBool;
use std::time::Instant;

use avalon::model::{
    AnalysisSettings, Bounds, DetectedEvent, Metadata, Ridge, SearchRequest, UtcWindow,
};
use avalon::{dateutil, pipeline, search};

const CLIP: &str = "/home/tpg/Downloads/Unknown location - raw audio - colsto (128p).mp4";
/// Truth for the colsto test footage (New York area).
const TRUTH: (f64, f64) = (40.87, -73.52);

fn haversine(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6_371_000.0f64;
    let (la1, la2) = (a.0.to_radians(), b.0.to_radians());
    let dla = (b.0 - a.0).to_radians();
    let dlo = (b.1 - a.1).to_radians();
    let h = (dla * 0.5).sin().powi(2) + la1.cos() * la2.cos() * (dlo * 0.5).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

fn peak_rss_mb() -> f64 {
    avalon::log::rss_bytes() as f64 / 1e6
}

#[derive(serde::Serialize, serde::Deserialize)]
struct Cached {
    events: Vec<DetectedEvent>,
    ridges: Vec<Ridge>,
}

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let mut clip = PathBuf::from(CLIP);
    let mut cache: Option<PathBuf> = None;
    let mut time: Option<(i64, i64)> = None;
    let mut truth = TRUTH;
    let mut jsonl: Option<PathBuf> = None;
    let mut provider = "staging-0".to_string();
    let mut regions: Vec<(String, Bounds)> = Vec::new();
    while let Some(a) = args.next() {
        match a.as_str() {
            "--clip" => clip = PathBuf::from(args.next().unwrap_or_default()),
            "--analysis" => cache = Some(PathBuf::from(args.next().unwrap_or_default())),
            "--jsonl" => jsonl = Some(PathBuf::from(args.next().unwrap_or_default())),
            "--provider" => provider = args.next().unwrap_or_default(),
            "--time" => {
                let s = args.next().unwrap_or_default();
                let (a, b) = s.split_once(':').unwrap_or((&s, &s));
                time = Some((a.parse()?, b.parse()?));
            }
            "--truth" => {
                let s = args.next().unwrap_or_default();
                let (a, b) = s.split_once(',').unwrap_or(("0", "0"));
                truth = (a.parse()?, b.parse()?);
            }
            _ => {
                let (label, boxs) = a.split_once(':').unwrap_or((a.as_str(), ""));
                let v: Vec<f64> = boxs.split(',').filter_map(|x| x.parse().ok()).collect();
                anyhow::ensure!(v.len() == 4, "region must be LABEL:W,S,E,N (got {a})");
                regions.push((
                    label.to_string(),
                    Bounds { west: v[0], south: v[1], east: v[2], north: v[3] },
                ));
            }
        }
    }
    let win = time.unwrap_or_else(|| {
        (
            dateutil::civil_to_unix(2026, 10, 30, 11, 0, 0),
            dateutil::civil_to_unix(2026, 10, 31, 23, 0, 0),
        )
    });

    // ---- analysis (cached) ------------------------------------------------
    let t0 = Instant::now();
    let (events, ridges) = match cache.as_ref().filter(|p| p.is_file()) {
        Some(p) => {
            let c: Cached = serde_json::from_slice(&std::fs::read(p)?)?;
            println!("analysis cache hit: {} events, {} ridges", c.events.len(), c.ridges.len());
            (c.events, c.ridges)
        }
        None => {
            let settings = AnalysisSettings::default();
            let meta = Metadata {
                region: regions.iter().map(|r| r.1).collect(),
                time_window: Some(UtcWindow { start_unix_s: win.0, end_unix_s: win.1 }),
                ..Default::default()
            };
            println!("analysing {clip:?}");
            let a = pipeline::run(Path::new(&clip), &settings, &meta, None)?;
            println!(
                "analysis: {} ridges, {} events in {:.0}s",
                a.ridges.len(),
                a.events.len(),
                t0.elapsed().as_secs_f64()
            );
            // Same selection the GUI autorun makes: plausible reps, skipping
            // the first (lowest) event.
            let mut ev = Vec::new();
            let mut rg = Vec::new();
            for e in &a.events {
                if e.id == 0 {
                    continue;
                }
                let Some(rep) = a.ridges.iter().find(|r| r.id == e.representative_ridge_id) else {
                    continue;
                };
                if !rep.plausible {
                    continue;
                }
                ev.push(e.clone());
                rg.push(rep.clone());
            }
            println!("selected {} event(s)", ev.len());
            for e in &ev {
                if let Some(f) = e.doppler {
                    println!(
                        "  event {}: {:.1}-{:.1}s f0 {:.1} Hz v {:.0} m/s slant {:.0} m cpa {:.1}s rmse {:.2}",
                        e.id, e.start_s, e.end_s, f.f0_hz, f.speed_mps,
                        f.slant_distance_m, f.t_cpa_s, f.rmse_hz
                    );
                }
            }
            if let Some(p) = cache.as_ref() {
                std::fs::write(
                    p,
                    serde_json::to_vec(&Cached { events: ev.clone(), ridges: rg.clone() })?,
                )?;
                println!("wrote analysis cache {}", p.display());
            }
            (ev, rg)
        }
    };
    anyhow::ensure!(!events.is_empty(), "no plausible events selected");

    // ---- locate per region ------------------------------------------------
    for (label, region) in &regions {
        println!("\n=== {label}: W{:.2} S{:.2} E{:.2} N{:.2}", region.west, region.south, region.east, region.north);
        let req = SearchRequest {
            region: vec![*region],
            time_window: UtcWindow { start_unix_s: win.0, end_unix_s: win.1 },
            event_ids: events.iter().map(|e| e.id).collect(),
            source_path: clip.to_string_lossy().into_owned(),
        };
        let cancel = AtomicBool::new(false);
        let mut last = String::new();
        let mut prog = |s: search::Snapshot| {
            let line = format!(
                "[{:.0}%] {} {} cand {}/{} best {:.2}@{:.3},{:.3}",
                s.fraction * 100.0, s.stage, s.detail, s.candidates_done,
                s.candidates_total, s.best_score_hz, s.best_lat, s.best_lon
            );
            if line != last {
                eprintln!("  {line}");
                last = line;
            }
        };
        let mut log = |m: String| println!("  {m}");
        let t = Instant::now();
        let res = search::run_search(&req, &events, &ridges, &provider, &mut prog, &mut log, &cancel);
        let elapsed = t.elapsed().as_secs_f64();
        match res {
            Ok(r) => {
                let survivors = r.cells.iter().filter(|c| c.score_hz > 0.0).count();
                let pairs = r.cells.iter().filter(|c| c.pair.contains('+')).count();
                let best = r.best.first().cloned();
                let (blat, blon, bscore, brmse, bpair, bdist) = match &best {
                    Some(b) => (
                        b.lat, b.lon, b.score_hz, b.rmse_hz, b.pair.clone(),
                        haversine((b.lat, b.lon), truth) / 1000.0,
                    ),
                    None => (f64::NAN, f64::NAN, 0.0, f32::NAN, String::new(), f64::NAN),
                };
                let near = r
                    .cells
                    .iter()
                    .min_by(|a, b| {
                        haversine((a.lat, a.lon), truth)
                            .partial_cmp(&haversine((b.lat, b.lon), truth))
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .map(|c| (c.lat, c.lon, c.score_hz, haversine((c.lat, c.lon), truth) / 1000.0));
                println!(
                    "RESULT {label} best={:.4},{:.4} score={:.3} rmse={:.2} pair={bpair} dist_km={:.1} \
                     survivors={survivors} pairs={pairs} candidates={} step_m={:.0} elapsed_s={elapsed:.0} rss_mb={:.0}",
                    blat, blon, bscore, brmse, bdist,
                    r.n_candidates,
                    r.step_deg * 111_320.0,
                    peak_rss_mb()
                );
                if let Some((la, lo, sc, d)) = near {
                    println!(
                        "  nearest-to-truth cell {la:.4},{lo:.4} score {sc:.3} ({d:.1} km off)"
                    );
                }
                println!("  modes (distinct near-tied solutions):");
                for m in r.modes.iter().take(6) {
                    let d = haversine((m.lat, m.lon), truth) / 1000.0;
                    println!(
                        "    {:.4},{:.4} score {:.3} rmse {:.2} Hz {:>11} {:>5} cells r~{:.0} km  -> {:.1} km from truth",
                        m.lat, m.lon, m.score, m.rmse_hz, m.pair, m.cells, m.radius_km, d
                    );
                }
                println!("  top 5:");
                for c in r.best.iter().take(5) {
                    println!(
                        "    {:.4},{:.4} score {:.3} rmse {:.2} {:>11} {:>7.1} km",
                        c.lat, c.lon, c.score_hz, c.rmse_hz, c.pair,
                        haversine((c.lat, c.lon), truth) / 1000.0
                    );
                }
                if let Some(p) = jsonl.as_ref() {
                    use std::io::Write;
                    let mut f = std::fs::OpenOptions::new().create(true).append(true).open(p)?;
                    writeln!(
                        f,
                        "{}",
                        serde_json::json!({
                            "label": label,
                            "region": [region.west, region.south, region.east, region.north],
                            "best": [blat, blon],
                            "score": bscore,
                            "rmse": brmse,
                            "pair": bpair,
                            "dist_km": bdist,
                            "nearest": near,
                            "survivors": survivors,
                            "pairs": pairs,
                            "candidates": r.n_candidates,
                            "step_m": r.step_deg * 111_320.0,
                            "elapsed_s": elapsed,
                            "rss_mb": peak_rss_mb(),
                        })
                    )?;
                }
            }
            Err(e) => println!("RESULT {label} ERROR {e}"),
        }
    }
    Ok(())
}
