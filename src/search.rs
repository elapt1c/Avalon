//! Candidate search: prune flights per candidate by closest-approach bands,
//! pair flyovers by time, and rank locations with matched-filter brightness.
//!
//! Ranking on traced-ridge RMSE smears into a ring: the metric is nearly flat
//! out to tens of kilometres. This scores the *modelled curve* against a
//! whitened spectrogram instead, sweeping f0 (±8 %) and one shared clock offset,
//! with the receiver at its real terrain elevation.

use std::sync::atomic::AtomicBool;
use std::time::Instant;

use anyhow::{anyhow, Result};
use rayon::prelude::*;

use crate::flights::FlightProgress;
use crate::model::{
    Bounds, DetectedEvent, Flight, HeatmapCell, Mode, Ridge, SearchRequest, SearchResult,
};
use crate::score::{self, Atmos, EventSpec, Track};

/// Target candidate spacing (metres) for a small region; the level-0 grid is
/// never finer than this.
const GRID_M: f64 = 200.0;
/// Cells in the level-0 grid over the whole region. Refinement levels add only
/// a few hundred cells each, so this is the only real cost knob.
const COARSE_CELLS: usize = 60_000;
/// Sub-division factor per refinement level (children per parent side).
///
/// Must be EVEN: the children are centred on the parent with `LEVEL_SPLIT/2`
/// steps either side, so an even split lands exactly on ±half a parent cell and
/// tiles the parent without gaps. An odd split leaves strips of the parent
/// cell un-sampled, and a true location in one of those strips is lost no
/// matter how many levels are added.
const LEVEL_SPLIT: usize = 6;
/// Parents refined per level (after the first).
const LEVEL_PARENTS: usize = 12;
/// Parents refined at the first level: wide on purpose, because the level-0
/// score is the cheapest and least trustworthy ranking in the whole search.
const LEVEL_PARENTS_FIRST: usize = 200;
/// Stop refining once cells are this close together.
const TARGET_STEP_M: f64 = 60.0;
/// Preferred level-0 spacing: score a candidate every 500 m.
const DENSE_STEP_M: f64 = 500.0;
/// Refinement stops before children get finer than this.
const MIN_REFINE_STEP_M: f64 = 50.0;
/// Cell budget for level 0. 500 m spacing is honoured wherever the region
/// allows it; past this budget the grid has to stretch so a run stays finite.
const MAX_LEVEL0_CELLS: usize = 400_000;
/// Hard stop on refinement depth (5^6 = 15625 x per side).
const MAX_LEVELS: usize = 6;
/// Flights per event kept as pair candidates in one cell. Not a ranking:
/// the pair test below rejects almost everything on the clock constraint.
const MAX_EVENT_TRACKS: usize = 600;
/// Flight pairs scored per cell.
const CELL_PAIRS: usize = 4;
/// Track sampling for CPA search (s).
const CPA_DT: f64 = 5.0;
/// Track sampling for scoring (s).
const SCORE_DT: f64 = 2.0;
/// Clock-offset sweep half width (s) and coarse/fine step counts.
const SHIFT_HALF_S: f64 = 10.0;

/// Half-width (s) of the shared-start tolerance when pairing two flyovers.
/// Distinct from the f0/clock sweep: this only decides whether two candidate
/// flights are *allowed* to be scored as one recording. Tunable because the
/// reference pair sits right on the edge of it.
fn shift_half_s() -> f64 {
    std::env::var("AOS_SHIFT_HALF_S")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(SHIFT_HALF_S)
}

fn pair_tol_s() -> f64 {
    std::env::var("AOS_PAIR_TOL_S")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(SHIFT_HALF_S + 2.0)
}
const SHIFT_STEPS_COARSE: usize = 3;
const SHIFT_STEPS_FINE: usize = 17;
/// f0 sweep steps.
const F0_STEPS_COARSE: usize = 5;
const F0_STEPS_FINE: usize = 25;
/// Frame subsampling for the coarse pass.
const FRAME_STEP_COARSE: usize = 16;

/// Live, real telemetry from the search.
#[derive(Debug, Clone, Default)]
pub struct Snapshot {
    pub stage: &'static str,
    pub detail: String,
    pub fraction: f32,
    pub elapsed_s: f64,
    pub eta_s: f64,
    pub bytes_read: u64,
    pub bytes_total: u64,
    pub traces_seen: usize,
    pub flights_kept: usize,
    pub days_done: usize,
    pub days_total: usize,
    pub candidates_done: usize,
    pub candidates_total: usize,
    pub cells_per_s: f64,
    /// Best (highest) brightness score found so far.
    pub best_score_hz: f32,
    pub best_lat: f64,
    pub best_lon: f64,
    pub best_flight: String,
    pub cells_matched: usize,
    /// Flight pairs scored so far.
    pub pairs_scored: usize,
}

struct Prepared {
    spec: EventSpec,
    /// Fitted closest-approach slant distance (bounds the search radius).
    slant_m: f64,
    /// Fitted source speed (m/s), used to reject flights flying the wrong speed.
    speed_mps: f64,
    /// Measured ridge curve, recording-relative: (seconds, Hz) — full detail.
    curve: Vec<(f64, f32)>,
    /// ~16-point version used by the coarse prefilter pass.
    coarse_curve: Vec<(f64, f32)>,
}

/// Run the search over the selected events.
pub fn run_search(
    req: &SearchRequest,
    events: &[DetectedEvent],
    ridges: &[Ridge],
    provider: &str,
    on_progress: &mut dyn FnMut(Snapshot),
    on_log: &mut dyn FnMut(String),
    cancel: &AtomicBool,
) -> Result<SearchResult> {
    let t_start = Instant::now();
    let mut snap = Snapshot {
        stage: "Fetching historical flights",
        ..Default::default()
    };
    on_progress(snap.clone());

    let mut download = |p: FlightProgress| {
        snap.detail = format!("{} · {}", p.tag, p.asset);
        snap.bytes_read = p.bytes_read;
        snap.bytes_total = p.bytes_total;
        snap.traces_seen = p.traces_seen;
        snap.flights_kept = p.flights_kept;
        snap.days_done = p.days_done;
        snap.days_total = p.days_total;
        snap.elapsed_s = t_start.elapsed().as_secs_f64();
        snap.fraction = if p.days_total > 0 {
            let day_frac = if p.bytes_total > 0 {
                (p.bytes_read as f64 / p.bytes_total as f64) as f32
            } else {
                0.0
            };
            ((p.days_done as f32 + day_frac) / p.days_total as f32).min(1.0)
        } else {
            0.0
        };
        snap.eta_s = if snap.fraction > 0.01 {
            snap.elapsed_s * (1.0 - snap.fraction as f64) / snap.fraction as f64
        } else {
            0.0
        };
        on_progress(snap.clone());
    };

    let flights = crate::flights::fetch_flights(
        req.time_window.start_unix_s,
        req.time_window.end_unix_s,
        &req.region,
        provider,
        &mut download,
        cancel,
    )?;
    on_log(format!(
        "{} flights intersect the region ({:.0} MB scanned)",
        flights.len(),
        snap.bytes_read as f64 / 1e6
    ));
    if flights.is_empty() {
        return Err(anyhow!("no flights found in the region/time window"));
    }

    // ---- audio -> whitened brightness per event -------------------------
    snap.stage = "Preparing audio";
    snap.fraction = 0.0;
    on_progress(snap.clone());
    let mut prepared: Vec<Prepared> = Vec::new();
    for ev in events {
        let Some(fit) = ev.doppler else {
            return Err(anyhow!("event {} has no Doppler fit", ev.id));
        };
        // Harmonics actually measured for this flyover.
        let mut harmonics: Vec<f32> = ev
            .harmonic_ratios
            .iter()
            .copied()
            .filter(|r| *r >= 0.9)
            .collect();
        if harmonics.is_empty() {
            harmonics.push(1.0);
        }
        harmonics.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        harmonics.dedup_by(|a, b| (*a - *b).abs() < 0.15);
        harmonics.truncate(4);

        let half = ((ev.end_s - ev.start_s) * 0.5).clamp(45.0, 180.0);
        let spec = EventSpec {
            f0: fit.f0_hz,
            harmonics: harmonics.clone(),
            cpa_s: fit.t_cpa_s,
            half_width_s: half,
        };
        // Measured curve for this event: the traced ridge, in seconds relative
        // to its CPA. Scoring compares this against the simulated curve.
        let Some(ridge) = ridges.iter().find(|r| r.id == ev.representative_ridge_id) else {
            return Err(anyhow!("event {}: ridge {} missing", ev.id, ev.representative_ridge_id));
        };
        // Search only on the high-contrast core of the ridge. The extension
        // tail chases faint neighbouring structure (visible as the trace
        // peeling away from the bright line) and drags f0 and the CPA with it.
        let median_w = crate::dsp::median_f32(
            &ridge.samples.iter().map(|s| s.weight).collect::<Vec<_>>(),
        );
        let w_thresh = median_w.max(1.0);
        let core: Vec<&crate::model::RidgeSample> = ridge
            .samples
            .iter()
            .filter(|s| s.weight >= w_thresh)
            .collect();
        let step = (core.len() / 120).max(1);
        let curve: Vec<(f64, f32)> = core
            .iter()
            .step_by(step)
            .map(|sm| (sm.t_s, sm.freq_hz))
            .collect();
        if curve.len() < 8 {
            return Err(anyhow!("event {}: curve too short ({} samples)", ev.id, curve.len()));
        }
        on_log(format!(
            "event {}: f0 {:.1} Hz, harmonics {:?}, {} curve samples {:.0}-{:.0} Hz",
            ev.id,
            spec.f0,
            harmonics,
            curve.len(),
            ridge.min_freq_hz,
            ridge.max_freq_hz
        ));
        let cstep = (ridge.samples.len() / 16).max(1);
        let coarse_curve: Vec<(f64, f32)> = ridge
            .samples
            .iter()
            .step_by(cstep)
            .map(|sm| (sm.t_s, sm.freq_hz))
            .collect();
        prepared.push(Prepared {
            spec,
            slant_m: fit.slant_distance_m as f64,
            speed_mps: fit.speed_mps as f64,
            curve,
            coarse_curve,
        });
    }

    // ---- flights -> tracks ----------------------------------------------
    //
    // Build tracks ONLY for flights that are airborne during an event. Doing
    // this for every flight in the day allocates a 1 s ECEF track per aircraft
    // (tens of KB each, thousands of aircraft) and is what exhausted memory
    // before scoring even started.
    snap.stage = "Building tracks";
    on_progress(snap.clone());
    // The recording happened SOMEWHERE inside the user's window, so any flight
    // airborne within it can be the one heard; the offset is fitted per
    // candidate rather than assumed. Flights entirely outside the window (plus
    // one clip length of slack) cannot have produced the audio.
    let win0 = req.time_window.start_unix_s as f64;
    let win1 = req.time_window.end_unix_s as f64;
    let clip_len = prepared
        .iter()
        .map(|p| 2.0 * p.spec.half_width_s)
        .fold(0.0f64, f64::max);
    let slack = clip_len + 120.0;
    if let Ok(list) = std::env::var("AVALON_DEBUG_ICAO") {
        for name in list.split(',').map(|x| x.trim()).filter(|x| !x.is_empty()) {
            match flights.iter().find(|f| f.icao.eq_ignore_ascii_case(name)) {
                Some(f) => crate::log::info(
                    "trace",
                    &format!(
                        "{} present in archive: {} points, t {:.0}..{:.0}, alt {:.0}..{:.0} m",
                        name,
                        f.points.len(),
                        f.points.first().map(|p| p.t_s).unwrap_or(0.0),
                        f.points.last().map(|p| p.t_s).unwrap_or(0.0),
                        f.points.iter().map(|p| p.alt_m).fold(f64::INFINITY, f64::min),
                        f.points.iter().map(|p| p.alt_m).fold(0.0, f64::max)
                    ),
                ),
                None => crate::log::info("trace", &format!("{name} NOT in parsed archive")),
            }
        }
    }

    let time_relevant: Vec<&Flight> = flights
        .iter()
        .filter(|f| {
            let (a, b) = (f.points[0].t_s, f.points[f.points.len() - 1].t_s);
            b >= win0 - slack && a <= win1 + slack
        })
        .collect();
    on_log(format!(
        "{} of {} flights airborne inside the search window (offset fitted per flight)",
        time_relevant.len(),
        flights.len()
    ));
    crate::log::info(
        "search",
        &format!(
            "{} flights, {} time-relevant, rss {:.0} MB before tracks",
            flights.len(),
            time_relevant.len(),
            crate::log::rss_bytes() as f64 / 1e6
        ),
    );

    let bounds = union_bounds(&req.region)?;
    // Anything further than this from the region centre cannot reach a
    // candidate cell within the 40 km scoring gate.
    let diag_km = {
        let rlat = (bounds.north - bounds.south) * 111.32;
        let rlon = (bounds.east - bounds.west) * 111.32 * ((bounds.south + bounds.north) * 0.5).to_radians().cos();
        (rlat * rlat + rlon * rlon).sqrt()
    };
    let radius_m = diag_km * 500.0 + 45_000.0;
    let (clat, clon) = bounds.center();
    // Track building is pure per flight (clip the trace, resample twice), so it
    // parallelises cleanly and was previously a single-core stage.
    let tracks: Vec<(String, Track, Track)> = time_relevant
        .par_iter()
        .filter_map(|f| {
            let fine = build_track(f, SCORE_DT, (clat, clon), radius_m)?;
            let coarse = build_track(f, CPA_DT, (clat, clon), radius_m)?;
            Some((f.icao.clone(), coarse, fine))
        })
        .collect();
    if let Ok(list) = std::env::var("AVALON_DEBUG_ICAO") {
        for name in list.split(',').map(|x| x.trim()).filter(|x| !x.is_empty()) {
            match tracks.iter().find(|t| t.0.eq_ignore_ascii_case(name)) {
                Some(t) => crate::log::info(
                    "trace",
                    &format!("{} has a usable track: {} coarse / {} fine samples", name, t.1.pos.len(), t.2.pos.len()),
                ),
                None => crate::log::warn(
                    "trace",
                    &format!("{name} dropped during track building (region radius {:.0} km)", radius_m / 1000.0),
                ),
            }
        }
    }
    on_log(format!("{} usable tracks", tracks.len()));

    // Everything the scoring stage needs now lives in `tracks`; the parsed
    // flights (≈1 GB of points for a continental window) can go.
    let n_flights_total = flights.len();
    drop(time_relevant);
    drop(flights);

    // Spatial hash of track passages. Scanning every track for every cell is
    // cells x tracks x samples; a 10 km bucket grid turns it into a handful of
    // lookups per cell.
    const BUCKET_KM: f64 = 10.0;
    let bucket_deg = BUCKET_KM / 111.32;
    let mut index: std::collections::HashMap<(i32, i32), Vec<u32>> =
        std::collections::HashMap::new();
    for (ti, (_, coarse, _)) in tracks.iter().enumerate() {
        let stride = ((30.0 / coarse.dt.max(1e-3)) as usize).max(1);
        let mut i = 0usize;
        while i < coarse.pos.len() {
            let (lat, lon) = ecef_to_latlon(coarse.point(i));
            index.entry(bucket_of(lat, lon, bucket_deg)).or_default().push(ti as u32);
            i += stride;
        }
    }
    // Largest plausible slant across the selected events bounds how far a
    // candidate cell can be from a contributing track.
    let gate_m = prepared
        .iter()
        .map(|p| (3.0 * p.slant_m).clamp(5_000.0, 40_000.0))
        .fold(0.0f64, f64::max);
    crate::log::info(
        "search",
        &format!(
            "spatial index: {} buckets, gate {:.0} km",
            index.len(),
            gate_m / 1000.0
        ),
    );
    let candidates_for = |lat: f64, lon: f64| -> Vec<u32> {
        let lat_deg = gate_m / 111_320.0;
        let lon_deg = gate_m / (111_320.0 * lat.to_radians().cos().abs().max(0.15));
        let (i0, j0) = bucket_of(lat - lat_deg, lon - lon_deg, bucket_deg);
        let (i1, j1) = bucket_of(lat + lat_deg, lon + lon_deg, bucket_deg);
        let mut out: Vec<u32> = Vec::new();
        for bi in i0..=i1 {
            for bj in j0..=j1 {
                if let Some(v) = index.get(&(bi, bj)) {
                    out.extend_from_slice(v);
                }
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    };
    crate::log::info(
        "search",
        &format!(
            "{} tracks built, rss {:.0} MB",
            tracks.len(),
            crate::log::rss_bytes() as f64 / 1e6
        ),
    );
    if tracks.is_empty() {
        return Err(anyhow!("no flight tracks cover the event times"));
    }

    // ---- hierarchical candidate search -----------------------------------
    //
    // One uniform grid cannot serve both a 200 km box and a continent. The
    // absolute candidate cap forced ~15 km spacing over CONUS, and the old
    // "refine" stage only re-scored those same cells at a finer *time*
    // resolution — the spatial answer could never be better than that 15 km
    // grid, which is why a continental search landed states away. Instead:
    // score a coarse grid across the whole region, then repeatedly sub-divide
    // the best cells and re-score the children at progressively higher
    // fidelity. Final spacing is `level0 / LEVEL_SPLIT^levels`, tens of metres
    // regardless of how large the region is.
    let mid_lat = (bounds.south + bounds.north) * 0.5;
    let lat_extent_m = ((bounds.north - bounds.south).abs() * 111_320.0).max(1.0);
    let lon_extent_m =
        ((bounds.east - bounds.west).abs() * 111_320.0 * mid_lat.to_rational_cos()).max(1.0);
    // Optional floor on the level-0 spacing, for experiments that measure how
    // coarse the first grid may be before the score kernel is missed entirely.
    let min_step = std::env::var("AOS_MIN_STEP_M")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(0.0);
    let coarse_cells = std::env::var("AOS_COARSE_CELLS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(COARSE_CELLS)
        .max(64);
    // Level-0 spacing: 500 m by default, stretched only when the region would
    // blow past the cell budget. The old behaviour derived spacing purely from
    // a fixed cell count, so a mid-size box silently started at 800 m and a
    // continental one at 15 km — the map showed the coarse grid, not the answer.
    let auto_step = level0_step_m(lat_extent_m, lon_extent_m, coarse_cells);
    let dense_step = (lat_extent_m * lon_extent_m / MAX_LEVEL0_CELLS as f64)
        .sqrt()
        .max(DENSE_STEP_M);
    let target = std::env::var("AOS_TARGET_STEP_M")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(TARGET_STEP_M);
    let mut step_m = auto_step
        .min(dense_step)
        .max(GRID_M)
        .max(min_step);

    let mut terrain = crate::dem::Terrain::for_area((
        bounds.east - bounds.west,
        bounds.north - bounds.south,
    ));

    snap.stage = "Scoring candidates";
    snap.fraction = 0.0;
    snap.elapsed_s = t_start.elapsed().as_secs_f64();
    on_progress(snap.clone());
    let score_start = Instant::now();

    // Every cell ever scored, kept for the heat map. Cells that were later
    // sub-divided are flagged `superseded` so their cheap-and-noisy score cannot
    // outrank the children that replaced them.
    let mut cells: Vec<Cell> = Vec::new();
    let level0_step_m = step_m;
    let mut frontier = grid_cells(&bounds, step_m);
    let mut level = 0usize;
    let mut scored_total = 0usize;
    while !frontier.is_empty() {
        let pts: Vec<(f64, f64)> = frontier.iter().map(|c| (c.lat, c.lon)).collect();
        let t_terrain = Instant::now();
        terrain.prefetch(&pts);
        let elevs = terrain.elevations(&pts);
        for (c, e) in frontier.iter_mut().zip(elevs.iter()) {
            c.obs = ecef(c.lat, c.lon, *e as f64 + 2.0);
        }
        let quality = match level {
            0 => Quality::Coarse,
            1 => Quality::Medium,
            _ => Quality::Fine,
        };
        on_log(format!(
            "level {level}: {} cells @ {:.0} m spacing, terrain {:.1}s (zoom {})",
            frontier.len(),
            step_m,
            t_terrain.elapsed().as_secs_f64(),
            terrain.zoom()
        ));
        crate::log::info(
            "search",
            &format!(
                "level {level}: {} cells @ {:.0} m, quality {quality:?}, rss {:.0} MB",
                frontier.len(),
                step_m,
                crate::log::rss_bytes() as f64 / 1e6
            ),
        );
        score_cells(&mut frontier, &candidates_for, &prepared, &tracks, quality, (win0, win1));
        scored_total += frontier.len();
        snap.candidates_done = scored_total;
        snap.candidates_total = scored_total;
        snap.cells_per_s =
            scored_total as f64 / score_start.elapsed().as_secs_f64().max(1e-6);
        snap.elapsed_s = t_start.elapsed().as_secs_f64();
        snap.cells_matched =
            frontier.iter().filter(|c| c.score.is_finite()).count();
        snap.detail = format!(
            "level {level} @ {:.0} m: {}/{} matched",
            step_m,
            snap.cells_matched,
            frontier.len()
        );
        snap.fraction = match level {
            0 => 0.5,
            l => (0.5 + 0.5 * (l as f32 / (MAX_LEVELS as f32))).min(0.95),
        };
        on_progress(snap.clone());

        if level == 0 {
            // Visibility into the coarse ranking: this is what decides whether
            // the right area survives to be refined at all.
            let mut top: Vec<(f32, f64, f64)> = frontier
                .iter()
                .filter(|c| c.score.is_finite() && c.score > 0.0)
                .map(|c| (c.score, c.lat, c.lon))
                .collect();
            top.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
            for (i, (sc, la, lo)) in top.iter().take(8).enumerate() {
                crate::log::info(
                    "search",
                    &format!("  coarse #{i} {sc:.3} @ {la:.3},{lo:.3}"),
                );
            }
            crate::log::info(
                "search",
                &format!("  coarse pass: {} of {} cells produced a hypothesis", top.len(), frontier.len()),
            );
        }
        let base = cells.len();
        cells.extend(frontier.drain(..));
        if step_m <= target || level >= MAX_LEVELS {
            break;
        }
        let n_parents = if level == 0 { LEVEL_PARENTS_FIRST } else { LEVEL_PARENTS };
        let parents = pick_clusters(&cells[base..], n_parents, step_m * 1.5);
        if parents.is_empty() {
            on_log(format!("level {level}: no cell to refine further"));
            break;
        }
        let child_step = step_m / LEVEL_SPLIT as f64;
        // The score kernel is ~1 km wide, so sub-50 m cells add nothing but
        // cost. On a dense level-0 grid this is what used to spend minutes
        // scoring an 11 m grid to move the answer by metres.
        if child_step < MIN_REFINE_STEP_M {
            on_log(format!(
                "stopping refinement: next level would be {child_step:.0} m (< {MIN_REFINE_STEP_M:.0} m adds nothing)"
            ));
            break;
        }
        let children = children_of(&cells[base..], &parents, child_step, &bounds);
        crate::log::info(
            "search",
            &format!("level {level}: {} parents -> {} children", parents.len(), children.len()),
        );
        for &pi in &parents {
            cells[base + pi].superseded = true;
        }
        if children.len() < 3 {
            break;
        }
        frontier = children;
        step_m = child_step;
        level += 1;
    }
    cells.retain(|c| !c.superseded);
    on_log(format!(
        "hierarchical search: {level} refinement level(s), final spacing {:.0} m, {} cells kept",
        step_m,
        cells.len()
    ));

    // ---- selectivity ------------------------------------------------------
    // How small a residual is achievable depends on how cleanly the audio can be
    // measured, so the accept threshold comes from the grid itself: keep the
    // cells whose residual is within a small margin of the best few percent.
    // An absolute constant either floods the map or empties it.
    let mut finite: Vec<f32> = cells.iter().map(|c| c.rmse_hz).filter(|r| r.is_finite()).collect();
    finite.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // Anchor on the best achievable residual, not on a percentile: with tens of
    // thousands of hypotheses a low percentile is dominated by flukes from
    // single-event matching and keeps a third of the grid.
    let floor = finite.first().copied().unwrap_or(f32::INFINITY);
    // Margin over the best residual achievable *anywhere on the grid*. A tight
    // multiple hides real but weaker solutions: over a continent the minimum is
    // a best-of-thousands coincidence, so the honest cut has to sit above it.
    // `AOS_CUT` overrides the multiple for experiments.
    let mult = std::env::var("AOS_CUT")
        .ok()
        .and_then(|v| v.parse::<f32>().ok())
        .unwrap_or(2.5);
    let cut = (mult * floor).max(1.0);
    crate::log::info(
        "search",
        &format!(
            "{} cells scored, residual floor {:.2} Hz -> keep <= {:.2} Hz",
            finite.len(),
            floor,
            cut
        ),
    );
    for c in cells.iter_mut() {
        if !(c.rmse_hz <= cut) {
            c.score = 0.0;
            c.rmse_hz = f32::INFINITY;
        }
    }

    let mut out: Vec<HeatmapCell> = cells
        .iter()
        .map(|c| HeatmapCell {
            lat: c.lat,
            lon: c.lon,
            score_hz: if c.score.is_finite() { c.score.max(0.0) } else { 0.0 },
            rmse_hz: c.rmse_hz,
            pair: c.pair.clone(),
        })
        .collect();
    out.sort_by(|a, b| b.score_hz.partial_cmp(&a.score_hz).unwrap_or(std::cmp::Ordering::Equal));
    let best: Vec<HeatmapCell> = out.iter().take(50).cloned().collect();

    // ---- distinct modes ---------------------------------------------------
    // Over a continent the leading hypotheses come in near-ties that are
    // hundreds of kilometres apart: two unrelated flight pairs can fit the same
    // two measured curves to within a few hundredths of a Hz. Collapsing the
    // survivors into spatially coherent clusters and reporting ALL of them is
    // the honest answer — the single best cell is over-confident by
    // construction. (The reference implementation hits the same wall: its
    // national run ranked a hypothesis 1,046 km from truth first, with the true
    // pair fourth on a near-tied score.)
    const MODE_LINK_KM: f64 = 15.0;
    let mut modes: Vec<Mode> = Vec::new();
    for c in out.iter().filter(|c| c.score_hz > 0.0) {
        let mut hit: Option<usize> = None;
        for (i, m) in modes.iter().enumerate() {
            let d = haversine_km((c.lat, c.lon), (m.lat, m.lon));
            if d <= MODE_LINK_KM + m.radius_km as f64 {
                hit = Some(i);
                break;
            }
        }
        match hit {
            Some(i) => {
                let m = &mut modes[i];
                let n = m.cells as f64;
                m.lat = (m.lat * n + c.lat) / (n + 1.0);
                m.lon = (m.lon * n + c.lon) / (n + 1.0);
                m.cells += 1;
                let d = haversine_km((c.lat, c.lon), (m.lat, m.lon)) as f32;
                m.radius_km = m.radius_km.max(d);
                if c.score_hz > m.score {
                    m.score = c.score_hz;
                    m.rmse_hz = c.rmse_hz;
                    m.pair = c.pair.clone();
                }
            }
            None => modes.push(Mode {
                lat: c.lat,
                lon: c.lon,
                score: c.score_hz,
                rmse_hz: c.rmse_hz,
                pair: c.pair.clone(),
                cells: 1,
                radius_km: 0.0,
            }),
        }
    }
    modes.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));
    modes.truncate(50);
    if let (Some(a), Some(b)) = (modes.first(), modes.get(1)) {
        let margin = (a.score - b.score) / a.score.max(1e-6);
        on_log(format!(
            "{} distinct mode(s); leader beats the runner-up by {:.0}% \
             (a margin under ~20% means both are equally supported)",
            modes.len(),
            margin * 100.0
        ));
    } else if let Some(a) = modes.first() {
        on_log(format!(
            "single mode at {:.4},{:.4} (score {:.3}, residual {:.2} Hz)",
            a.lat, a.lon, a.score, a.rmse_hz
        ));
    }

    // Simulated curves for the winning cell (visual verification aid).
    let (best_cell, model_curves) = match best.first() {
        Some(b) => {
            let obs = ecef(b.lat, b.lon, 0.0);
            let cand = candidates_for(b.lat, b.lon);
            let curves = model_curves_at(
                obs,
                &cand,
                &prepared,
                &tracks,
                (win0, win1),
                &score::shift_grid(shift_half_s(), SHIFT_STEPS_FINE),
                F0_STEPS_FINE,
            );
            // Report how closely the simulated curve tracks the measured one.
            for (i, mc) in curves.iter().enumerate() {
                if let Some(p) = prepared.get(i) {
                    let mut sum = 0.0f64;
                    let mut n = 0usize;
                    for (ti, (t, f)) in p.curve.iter().enumerate() {
                        if let Some((_, mf)) = mc.get(ti) {
                            let _ = t;
                            sum += (*mf as f64 - *f as f64).abs();
                            n += 1;
                        }
                    }
                    if n > 0 {
                        crate::log::info(
                            "search",
                            &format!(
                                "overlay event {i}: {} pts, mean |sim-measured| {:.2} Hz",
                                mc.len(),
                                sum / n as f64
                            ),
                        );
                    }
                }
            }
            (Some((b.lat, b.lon)), curves)
        }
        None => (None, Vec::new()),
    };

    // ---- honesty about the regime ----------------------------------------
    // The score kernel is a few hundred metres to a few kilometres wide, so a
    // level-0 grid coarser than that can step straight over the true location.
    // Say so; do not present a screening pass as a fix.
    let warning = if level0_step_m > 4_000.0 {
        Some(format!(
            "Region is too large for a reliable fix: the search grid is {:.1} km at its \
             coarsest level, but the score kernel is only ~1 km wide, so the true \
             location can be missed entirely between cells. Results are a screening \
             pass only — narrow the region to a few hundred kilometres (or narrow the \
             time window and re-run) before trusting a point.",
            level0_step_m / 1000.0
        ))
    } else if modes.len() > 1 {
        let (a, b) = (&modes[0], &modes[1]);
        let margin = (a.score - b.score) / a.score.max(1e-6);
        if margin < 0.20 {
            Some(format!(
                "Near-tie: the top two solutions are {:.0} km apart ({} vs {}) and differ \
                 by only {:.0}%. Two unrelated flight pairs can fit both measured curves \
                 this well; use other evidence to choose.",
                haversine_km((a.lat, a.lon), (b.lat, b.lon)),
                a.pair,
                b.pair,
                margin * 100.0
            ))
        } else {
            None
        }
    } else {
        None
    };
    if let Some(w) = &warning {
        on_log(format!("WARNING: {w}"));
        crate::log::warn("search", w);
    }

    Ok(SearchResult {
        cells: out,
        best,
        warning,
        modes,
        n_flights: n_flights_total,
        n_candidates: cells.len(),
        bounds: Some(bounds),
        step_deg: step_m / 111_320.0,
        best_cell,
        model_curves,
    })
}

struct Cell {
    lat: f64,
    lon: f64,
    obs: [f64; 3],
    score: f32,
    /// Mean curve residual in Hz for the winning hypothesis (INFINITY = none).
    rmse_hz: f32,
    pair: String,
    /// Replaced by finer children during hierarchical refinement; excluded from
    /// the final ranking so a cheap coarse score cannot outrank its own children.
    superseded: bool,
}

impl Cell {
    fn new(lat: f64, lon: f64) -> Self {
        Self {
            lat,
            lon,
            obs: ecef(lat, lon, 0.0),
            score: f32::NEG_INFINITY,
            rmse_hz: f32::INFINITY,
            pair: String::new(),
            superseded: false,
        }
    }
}

/// Scoring fidelity. Coarse is ~16x cheaper per cell than fine and is used for
/// searching, never for the final ranking.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Quality {
    Coarse,
    Medium,
    Fine,
}

impl Quality {
    fn frame_step(self) -> usize {
        match self {
            Quality::Coarse => FRAME_STEP_COARSE,
            Quality::Medium => 4,
            Quality::Fine => 1,
        }
    }
    fn f0_steps(self) -> usize {
        match self {
            Quality::Coarse => F0_STEPS_COARSE,
            Quality::Medium => 9,
            Quality::Fine => F0_STEPS_FINE,
        }
    }
    fn shift_steps(self) -> usize {
        match self {
            Quality::Coarse => SHIFT_STEPS_COARSE,
            Quality::Medium => 9,
            Quality::Fine => SHIFT_STEPS_FINE,
        }
    }
}

/// Spacing of the level-0 grid for a region of the given extent.
fn level0_step_m(lat_extent_m: f64, lon_extent_m: f64, cells: usize) -> f64 {
    (lat_extent_m.max(1.0) * lon_extent_m.max(1.0) / cells.max(64) as f64).sqrt()
}

/// Uniform grid covering `b`, spaced `step_m` metres.
fn grid_cells(b: &Bounds, step_m: f64) -> Vec<Cell> {
    let mid = (b.south + b.north) * 0.5;
    let lat_step = step_m / 111_320.0;
    let lon_step = step_m / (111_320.0 * mid.to_rational_cos());
    let n_lat = (((b.north - b.south) / lat_step).floor() as usize) + 1;
    let n_lon = (((b.east - b.west) / lon_step).floor() as usize) + 1;
    let mut out = Vec::with_capacity(n_lat.min(4096) * n_lon.min(4096));
    for i in 0..n_lat {
        for j in 0..n_lon {
            let lat = b.south + i as f64 * lat_step;
            let lon = b.west + j as f64 * lon_step;
            if lat > b.north || lon > b.east {
                continue;
            }
            out.push(Cell::new(lat, lon));
        }
    }
    out
}

/// Indices of the best `n` scored cells, at least `sep_m` apart.
fn pick_clusters(cells: &[Cell], n: usize, sep_m: f64) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..cells.len())
        .filter(|i| cells[*i].score.is_finite() && cells[*i].score > 0.0)
        .collect();
    idx.sort_by(|a, b| {
        cells[*b]
            .score
            .partial_cmp(&cells[*a].score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut out: Vec<usize> = Vec::new();
    for i in idx {
        if out.len() >= n {
            break;
        }
        let c = &cells[i];
        let coslat = c.lat.to_rational_cos();
        let far = out.iter().all(|&j| {
            let d = &cells[j];
            let dy = (c.lat - d.lat) * 111_320.0;
            let dx = (c.lon - d.lon) * 111_320.0 * coslat;
            (dx * dx + dy * dy).sqrt() >= sep_m
        });
        if far {
            out.push(i);
        }
    }
    out
}

/// A `LEVEL_SPLIT` x `LEVEL_SPLIT` sub-grid centred on each parent cell.
fn children_of(cells: &[Cell], parents: &[usize], child_step: f64, b: &Bounds) -> Vec<Cell> {
    let half = (LEVEL_SPLIT / 2) as i32;
    let mut out: Vec<Cell> = Vec::new();
    let mut seen: std::collections::HashSet<(i64, i64)> = std::collections::HashSet::new();
    for &pi in parents {
        let anchor = &cells[pi];
        let coslat = anchor.lat.to_rational_cos();
        for di in -half..=half {
            for dj in -half..=half {
                let lat = anchor.lat + di as f64 * child_step / 111_320.0;
                let lon = anchor.lon + dj as f64 * child_step / (111_320.0 * coslat);
                if lat < b.south || lat > b.north || lon < b.west || lon > b.east {
                    continue;
                }
                let key = (
                    (lat * 1e7).round() as i64,
                    ((lon + 360.0) * 1e7).round() as i64,
                );
                if !seen.insert(key) {
                    continue;
                }
                out.push(Cell::new(lat, lon));
            }
        }
    }
    out
}

/// Score a whole level in parallel.
fn score_cells(
    cells: &mut [Cell],
    candidates_for: &(dyn Fn(f64, f64) -> Vec<u32> + Sync),
    prepared: &[Prepared],
    tracks: &[(String, Track, Track)],
    quality: Quality,
    window: (f64, f64),
) {
    // Chunked so each worker amortises scheduling and keeps its scratch hot.
    const CHUNK: usize = 32;
    cells.par_chunks_mut(CHUNK).for_each(|chunk| {
        for c in chunk {
            let cand = candidates_for(c.lat, c.lon);
            let (sc, rm, pair) =
                score_cell(c.obs, (c.lat, c.lon), &cand, prepared, tracks, quality, window);
            c.score = sc;
            c.rmse_hz = rm;
            c.pair = pair;
        }
    });
}

trait CosExt {
    fn to_rational_cos(self) -> f64;
}
impl CosExt for f64 {
    fn to_rational_cos(self) -> f64 {
        self.to_radians().cos().max(0.15)
    }
}

/// Build a uniform ECEF track, keeping only the part of the trajectory that can
/// matter: an ADS-B trace covers the aircraft's whole day (up to 24 h), so
/// sampling it entirely costs ~1 MB per flight. A closest approach within the
/// search radius can only happen while the aircraft is near the region, so the
/// track is clipped to that passage.
fn build_track(
    f: &Flight,
    dt: f64,
    center: (f64, f64),
    radius_m: f64,
) -> Option<Track> {
    if f.points.len() < 2 {
        return None;
    }
    let mut first = None;
    let mut last = None;
    for (i, p) in f.points.iter().enumerate() {
        let dlat = (p.lat - center.0) * 111_320.0;
        let dlon = (p.lon - center.1) * 111_320.0 * center.0.to_radians().cos();
        if (dlat * dlat + dlon * dlon).sqrt() <= radius_m {
            first.get_or_insert(i);
            last = Some(i);
        }
    }
    let Some(a) = first else { return None };
    let Some(b) = last else { return None };
    // Grounded/parked transponders dwell in range for hours and can never be a
    // flyover. Judge that by *motion*, never by dwell time: inside a
    // continental search box a transcontinental jet legitimately stays for 6 h,
    // and the old `dwell > 4 h` cut silently discarded exactly the long flights
    // a coast-to-coast region is searching for.
    let slice = &f.points[a..=b];
    let dwell = slice.last().map(|p| p.t_s).unwrap_or(0.0)
        - slice.first().map(|p| p.t_s).unwrap_or(0.0);
    if dwell > 20.0 * 3600.0 {
        return None;
    }
    if slice.len() < 2 {
        return None;
    }
    // Airworthiness is judged **around the closest approach**, not over the whole
    // in-region slice. An aircraft that transits the box and then parks or taxis
    // inside it (its home field may lie in the region) has a slice that is mostly
    // on the ground, and the old mean-over-slice test threw it away — which is
    // exactly what happened to one half of the reference pair. The flyover we are
    // scoring is the close pass; judge that moment.
    let cpa_idx = {
        let mut best = (f64::INFINITY, 0usize);
        for (i, p) in slice.iter().enumerate() {
            let dlat = (p.lat - center.0) * 111_320.0;
            let dlon = (p.lon - center.1) * 111_320.0 * center.0.to_radians().cos();
            let d = (dlat * dlat + dlon * dlon).sqrt();
            if d < best.0 {
                best = (d, i);
            }
        }
        best.1
    };
    let t_cpa = slice[cpa_idx].t_s;
    let around: Vec<&crate::model::FlightPoint> = slice
        .iter()
        .filter(|p| (p.t_s - t_cpa).abs() <= 900.0)
        .collect();
    if around.len() < 2 {
        return None;
    }
    let max_alt = around.iter().map(|p| p.alt_m).fold(0.0f64, f64::max);
    if max_alt < 300.0 {
        return None;
    }
    // Speed at the pass: a parked transponder is stationary even if it sits in
    // the box all day, but a taxiing one still moves, so require real airspeed.
    let mut moved = 0.0f64;
    let mut dt_sum = 0.0f64;
    for w in around.windows(2) {
        let dt = w[1].t_s - w[0].t_s;
        if dt <= 0.0 {
            continue;
        }
        let dy = (w[1].lat - w[0].lat) * 111_320.0;
        let dx = (w[1].lon - w[0].lon) * 111_320.0 * w[0].lat.to_radians().cos();
        moved += (dx * dx + dy * dy).sqrt();
        dt_sum += dt;
    }
    if dt_sum <= 0.0 || moved / dt_sum < 30.0 {
        return None;
    }
    // One point of margin either side so interpolation has neighbours.
    let a = a.saturating_sub(1);
    let b = (b + 1).min(f.points.len() - 1);
    if b <= a + 1 {
        return None;
    }
    let mut pts: Vec<(f64, [f64; 3])> = Vec::with_capacity(b - a + 1);
    let mut bbox = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
    for p in &f.points[a..=b] {
        bbox[0] = bbox[0].min(p.lat);
        bbox[1] = bbox[1].min(p.lon);
        bbox[2] = bbox[2].max(p.lat);
        bbox[3] = bbox[3].max(p.lon);
        pts.push((p.t_s, ecef(p.lat, p.lon, p.alt_m)));
    }
    // Resample to a uniform grid so interpolation is O(1).
    let t0 = pts[0].0;
    let t1 = pts[pts.len() - 1].0;
    let n = ((t1 - t0) / dt).floor() as usize + 1;
    if n < 3 || n > 500_000 {
        return None;
    }
    let mut pos: Vec<[f32; 3]> = Vec::with_capacity(n);
    let origin = pts[0].1;
    let mut j = 0usize;
    for i in 0..n {
        let t = t0 + i as f64 * dt;
        while j + 2 < pts.len() && pts[j + 1].0 < t {
            j += 1;
        }
        let (ta, pa) = pts[j];
        let (tb, pb) = pts[(j + 1).min(pts.len() - 1)];
        let f = if (tb - ta).abs() > 1e-9 {
            ((t - ta) / (tb - ta)).clamp(0.0, 1.0)
        } else {
            0.0
        };
        pos.push([
            ((pa[0] + (pb[0] - pa[0]) * f) - origin[0]) as f32,
            ((pa[1] + (pb[1] - pa[1]) * f) - origin[1]) as f32,
            ((pa[2] + (pb[2] - pa[2]) * f) - origin[2]) as f32,
        ]);
    }
    Some(Track { t0, dt, origin, pos, bbox })
}

/// Score one candidate: prune flights by CPA band, pair events by time, then
/// maximise the matched-filter brightness over f0 and the shared clock offset.

/// Cell nominated by AVALON_DEBUG_CELL="lat,lon" for gate tracing.
fn debug_cell() -> Option<(f64, f64)> {
    let v = std::env::var("AVALON_DEBUG_CELL").ok()?;
    let (a, b) = v.split_once(',')?;
    Some((a.trim().parse().ok()?, b.trim().parse().ok()?))
}

mod gtracks {
    /// Is this ICAO one of the flights from the published ground truth?
    pub fn named(icao: &str) -> bool {
        matches!(
            std::env::var("AVALON_DEBUG_ICAO").ok().as_deref(),
            Some(list) if list.split(',').any(|x| x.trim().eq_ignore_ascii_case(icao))
        )
    }
}

fn score_cell(
    obs: [f64; 3],
    cell: (f64, f64),
    candidates: &[u32],
    prepared: &[Prepared],
    tracks: &[(String, Track, Track)],
    quality: Quality,
    window: (f64, f64),
) -> (f32, f32, String) {
    // Debug hook: print every gate decision for one nominated cell, so a cell
    // that "scores 0.000" can be explained instead of guessed at.
    let debug = debug_cell().map(|d| {
        let dd = ((cell.0 - d.0) * 111_320.0).powi(2)
            + ((cell.1 - d.1) * 111_320.0 * cell.0.to_radians().cos()).powi(2);
        (dd.sqrt(), d)
    });
    let is_debug_cell = debug.map(|(m, _)| m < 400.0).unwrap_or(false);
    if is_debug_cell {
        crate::log::info(
            "gate",
            &format!(
                "cell {:.4},{:.4}: {} candidate tracks, {} events",
                cell.0,
                cell.1,
                candidates.len(),
                prepared.len()
            ),
        );
        // Is a published ground-truth aircraft even in the candidate set?
        if let Ok(list) = std::env::var("AVALON_DEBUG_ICAO") {
            for name in list.split(',').map(|x| x.trim()).filter(|x| !x.is_empty()) {
                if let Some(ti) = tracks.iter().position(|t| t.0.eq_ignore_ascii_case(name)) {
                    let in_cands = candidates.contains(&(ti as u32));
                    let (d, _) = tracks[ti].1.cpa_stride(obs, 6).unwrap_or((f64::INFINITY, 0.0));
                    crate::log::info(
                        "gate",
                        &format!(
                            "  {name}: track #{ti} {} candidate list, cpa {:.1} km at this cell",
                            if in_cands { "IS in" } else { "is NOT in" },
                            d / 1000.0
                        ),
                    );
                }
            }
        }
    }

    // Per event: flights whose closest approach implies a recording start that
    // falls inside the user's window. `implied = t_cpa_flight - cpa_audio`.
    // (track index, cpa distance, ground speed, implied recording start)
    let mut per_event: Vec<Vec<(usize, f64, f64, f64)>> = Vec::with_capacity(prepared.len());
    for p in prepared {
        let mut keep = Vec::new();
        for &ti in candidates {
            let coarse_track = &tracks[ti as usize].1;
            let _ = cell;
            let stride = ((30.0 / coarse_track.dt.max(1e-3)) as usize).max(1);
            // Use the *acoustic* crossing (zero radial velocity) for the implied
            // start, and the geometric distance only for the gates. Pairing on
            // the geometric CPA mixed two different instants and cost ~10 s per
            // event on a turning aircraft.
            let geo = coarse_track.cpa_stride(obs, stride);
            let (d, t_cpa) = match (geo, coarse_track.radial_zero_near(obs, geo.map(|g| g.1))) {
                (Some((dg, tg)), Some((_dr, tr))) if (tr - tg).abs() < 120.0 => (dg, tr),
                (Some(g), _) => g,
                (None, _) => continue,
            };
            // Loose geometric gate: any close pass is worth scoring.
            if d > 40_000.0 {
                if is_debug_cell && gtracks::named(&tracks[ti as usize].0) {
                    crate::log::info("gate", &format!("  {} rejected: cpa {:.1} km > 40", tracks[ti as usize].0, d / 1000.0));
                }
                continue;
            }
            let implied = t_cpa - p.spec.cpa_s;
            if implied < window.0 - 120.0 || implied > window.1 + 120.0 {
                if is_debug_cell && gtracks::named(&tracks[ti as usize].0) {
                    crate::log::info(
                        "gate",
                        &format!(
                            "  {} rejected: implied start {:.0} outside window",
                            tracks[ti as usize].0, implied
                        ),
                    );
                }
                continue;
            }
            // The aircraft must be in the right speed class, otherwise every
            // busy-airspace flight "matches" and the map just reproduces flight
            // density. The band is wide because the fitted speed is itself
            // degenerate for a distant pass (calibrated at -28 % at 12 km); the
            // simulated-curve residual is what separates the real match.
            let gs = coarse_track
                .vel(t_cpa)
                .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
                .unwrap_or(0.0);
            if gs < 0.5 * p.speed_mps || gs > 2.0 * p.speed_mps {
                if is_debug_cell && gtracks::named(&tracks[ti as usize].0) {
                    crate::log::info(
                        "gate",
                        &format!("  {} rejected: ground speed {:.0} m/s vs fitted {:.0}", tracks[ti as usize].0, gs, p.speed_mps),
                    );
                }
                continue;
            }
            // The fitted slant is degenerate (it trades against speed), so allow
            // a wide factor. Evidence: at the true receiver the real pass was
            // 10.7 km while the audio fit said 5.9 km.
            if d < 0.1 * p.slant_m || d > 8.0 * p.slant_m {
                if is_debug_cell && gtracks::named(&tracks[ti as usize].0) {
                    crate::log::info(
                        "gate",
                        &format!("  {} rejected: cpa {:.1} km vs fitted slant {:.1} km", tracks[ti as usize].0, d / 1000.0, p.slant_m / 1000.0),
                    );
                }
                continue;
            }
            keep.push((ti as usize, d, gs, implied));
        }
        // Rank by how close the CPA is to the fitted slant, then keep only the
        // few most plausible: scoring cost is candidates x f0 x shift x frames.
        // Rank on speed agreement only. The fitted slant is unidentifiable past
        // ~7 km, so weighting it here dropped the true flight from the shortlist.
        // Do NOT rank by speed agreement and then truncate. The audio's fitted
        // speed is degenerate (documented at -28 % for a distant pass), so the
        // true jet can sort behind a turboprop that happens to fly at the fitted
        // speed and be discarded before the pair test ever sees it. With two
        // events the physical prune is "one recording start explains both",
        // applied below as a hash lookup.
        keep.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        keep.truncate(MAX_EVENT_TRACKS);
        per_event.push(keep);
    }

    let shift_grid: Vec<f64> = score::shift_grid(shift_half_s(), quality.shift_steps());
    let f0_steps = quality.f0_steps();
    let frame_step = quality.frame_step();

    // 0.0 = hard-filtered (rendered green); > 0 = survived and scored.
    let mut best = 0.0f32;
    let mut best_rmse = f32::INFINITY;
    let mut best_pair = String::new();
    let mut pair_scored = false;

    // Two or more events: pair flights across events by the time gap.
    if prepared.len() >= 2 {
        let (e0, e1) = (&prepared[0], &prepared[1]);
        // Both events were cut from ONE recording, so the same absolute
        // recording start must explain both flyovers. Bucket the second event's
        // candidates by implied start so that constraint is a hash lookup rather
        // than an all-pairs scan — this is what makes a continental region
        // tractable, and it is far stronger evidence than speed similarity.
        let tol = pair_tol_s();
        let mut by_start: std::collections::HashMap<i64, Vec<(usize, f64)>> =
            std::collections::HashMap::new();
        for &(i1, _d1, _v1, base1) in &per_event[1] {
            by_start.entry(base1.floor() as i64).or_default().push((i1, base1));
        }
        let mut cands: Vec<(usize, usize, f64)> = Vec::new();
        for &(i0, _d0, _v0, base0) in &per_event[0] {
            for key in (base0 - tol).floor() as i64..=(base0 + tol).ceil() as i64 {
                let Some(v) = by_start.get(&key) else { continue };
                for &(i1, _base1) in v {
                    if i1 != i0 {
                        cands.push((i0, i1, base0));
                    }
                }
            }
        }
        if is_debug_cell {
            if let Ok(list) = std::env::var("AVALON_DEBUG_ICAO") {
                let names: Vec<String> = list.split(',').map(|x| x.trim().to_string()).collect();
                for (ev, keep) in per_event.iter().enumerate() {
                    for (ti, d, gs, implied) in keep {
                        if names.iter().any(|n| n.eq_ignore_ascii_case(&tracks[*ti].0)) {
                            crate::log::info(
                                "gate",
                                &format!(
                                    "  event {ev} {}: cpa {:.1} km, gs {:.0} m/s, implied start {:.1}",
                                    tracks[*ti].0, d / 1000.0, gs, implied
                                ),
                            );
                        }
                    }
                }
                let truth_pair = names.len() >= 2
                    && cands.iter().any(|(i0, i1, _)| {
                        tracks[*i0].0.eq_ignore_ascii_case(&names[0])
                            && tracks[*i1].0.eq_ignore_ascii_case(&names[1])
                    });
                crate::log::info(
                    "gate",
                    &format!("  truth pair present in candidate pairs: {truth_pair} ({} pairs)", cands.len()),
                );
            }
        }

        // Rank those pairs by a cheap residual of the *subsampled* curves before
        // paying for the full f0 / clock sweep on the best few.
        let mut pairs: Vec<(f64, usize, usize, f64)> = Vec::with_capacity(cands.len());
        for (i0, i1, base0) in cands {
            let mut cost = 0.0f64;
            let mut ok = true;
            for (k, ti) in [(0usize, i0), (1usize, i1)] {
                let spec = &prepared[k].spec;
                // Sweep f0 here, not just at the audio's own estimate. The true
                // blade-pass rate can be several percent off the audio fit (the
                // reference pair sits at 113.6 and 126.3 Hz against a fit of
                // 119.6 / 136.6), and ranking on a single f0 then discards the
                // correct pair before the full sweep ever sees it.
                let mut best_k = f64::INFINITY;
                // Sweep f0 AND a small clock offset here. Ranking every pair at
                // shift 0 penalised exactly the pairs whose audio CPA is a few
                // seconds off — which is all of them, and the reference pair
                // worst of all (its two events imply starts 20 s apart).
                for f0 in score::f0_grid_rel(spec.f0 as f64, 3, 0.12) {
                    for sh in [-20.0, 0.0, 20.0] {
                        if let Some(v) = score::curve_rmse_offset(
                            &prepared[k].coarse_curve,
                            &tracks[ti].2,
                            obs,
                            f0,
                            base0,
                            sh,
                            0.0,
                            Atmos::default(),
                            40.0,
                        ) {
                            if v < best_k {
                                best_k = v;
                            }
                        }
                    }
                }
                if best_k.is_finite() {
                    cost += best_k;
                } else {
                    ok = false;
                    break;
                }
            }
            if ok {
                pairs.push((cost, i0, i1, base0));
            }
        }
        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let cell_pairs = std::env::var("AOS_CELL_PAIRS")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(CELL_PAIRS);
        pairs.truncate(cell_pairs);
        let mut scored: Vec<f32> = Vec::with_capacity(pairs.len());
        let mut top = f32::NEG_INFINITY;
        let mut top_name = String::new();
        let mut top_ids = (0usize, 0usize);
        let mut top_rmse = f32::INFINITY;
        for &(_, i0, i1, base0) in &pairs {
            let (s, rm) = score_pair(
                obs, e0, e1, &tracks[i0].2, &tracks[i1].2, base0, &shift_grid, f0_steps,
                frame_step,
            );
            scored.push(s);
            if s > top {
                top = s;
                top_rmse = rm;
                top_name = format!("{}+{}", tracks[i0].0, tracks[i1].0);
                top_ids = (i0, i1);
            }
        }
        // Silence test: the recording contains exactly these two passes, so a
        // candidate that would ALSO have heard another aircraft between them is
        // disfavoured. This is what stops busy corridors from matching.
        let mut silence = 1.0f64;
        if top.is_finite() {
            let stride = ((30.0 / tracks[top_ids.0].1.dt.max(1e-3)) as usize).max(1);
            let a = tracks[top_ids.0].1.cpa_stride(obs, stride);
            let b = tracks[top_ids.1].1.cpa_stride(obs, stride);
            if let (Some((d0, t0c)), Some((d1, t1c))) = (a, b) {
                let recv0 = t0c + d0 / crate::score::C0;
                let recv1 = t1c + d1 / crate::score::C0;
                let (lo, hi) = if recv0 < recv1 { (recv0, recv1) } else { (recv1, recv0) };
                // Count only passes that would have been as loud as the two we
                // traced. Doppler shift falls off with distance, so a pass twice
                // as far away sits ~6 dB lower and is masked by the recording's
                // own noise floor; treating any aircraft within 12 km as audible
                // vetoed every candidate in a busy terminal area.
                let audible = (d0.max(d1) * 2.0).max(4_000.0);
                let mut extra = 0usize;
                for &tj in candidates {
                    if tj as usize == top_ids.0 || tj as usize == top_ids.1 {
                        continue;
                    }
                    let tr = &tracks[tj as usize].1;
                    let st = ((30.0 / tr.dt.max(1e-3)) as usize).max(1);
                    let Some((d, t_cpa)) = tr.cpa_stride(obs, st) else { continue };
                    let recv = t_cpa + d / crate::score::C0;
                    if d <= audible && recv > lo + 5.0 && recv < hi - 5.0 {
                        extra += 1;
                    }
                }
                // Rank down rather than delete: the trace count is evidence, not
                // a hard geometric fact, and a hard cut left the map empty.
                silence = 1.0 / (1.0 + extra as f64);
            }
        }
        if top.is_finite() {
            // Contrast against the competing hypotheses in this cell: a dense
            // area where many pairs fit almost as well should not outrank a
            // quiet area with one clean match.
            let rest: Vec<f32> = scored.iter().copied().filter(|v| *v < top).collect();
            let margin = if rest.is_empty() {
                0.35 * top
            } else {
                let mean = rest.iter().sum::<f32>() / rest.len() as f32;
                (top - mean).max(0.0)
            };
            // The gap between the two modelled passes must match the audio's.
            let mut weight = 1.0f64;
            if let (Some((_, t0c)), Some((_, t1c))) = (
                tracks[top_ids.0].1.cpa(obs),
                tracks[top_ids.1].1.cpa(obs),
            ) {
                let modelled = (t1c - t0c).abs();
                let audio = (e1.spec.cpa_s - e0.spec.cpa_s).abs();
                let err = (modelled - audio).abs();
                weight = (-(err / 6.0).powi(2)).exp();
            }
            let score = (0.5 * top + 0.5 * margin) * weight as f32 * silence as f32;
            if score > best {
                best = score;
                best_rmse = top_rmse;
                best_pair = top_name;
                pair_scored = true;
            }
        }
    }

    // Single-event scoring: a fallback, and always ranked below a pair. One
    // measured curve with a free rest frequency, clock and slant can be matched
    // almost anywhere, so on its own it is density, not evidence — the user's
    // rule is that a witness heard *two* passes with a fixed offset.
    const SINGLE_WEIGHT: f32 = 0.25;
    let _ = pair_scored;
    for (ei, p) in prepared.iter().enumerate().take(if prepared.len() >= 2 { 0 } else { usize::MAX }) {
        let Some(cands) = per_event.get(ei) else { continue };
        for &(idx, _, _, base) in cands {
            let mut local = f32::NEG_INFINITY;
            let mut local_rmse = f32::INFINITY;
            for &f0 in &score::f0_grid(p.spec.f0 as f64, f0_steps) {
                for &sh in &shift_grid {
                    let c = if frame_step > 1 { &p.coarse_curve } else { &p.curve };
                    if let Some(v) = score::curve_rmse(
                        c,
                        &tracks[idx].2,
                        obs,
                        f0,
                        base,
                        sh,
                        Atmos::default(),
                        25.0,
                    ) {
                        let s = (1.0 / (1.0 + v)) as f32 * SINGLE_WEIGHT;
                        if s > local {
                            local = s;
                            local_rmse = v as f32;
                        }
                    }
                }
            }
            if local > best {
                best = local;
                best_rmse = local_rmse;
                best_pair = tracks[idx].0.clone();
            }
        }
    }

    (best, best_rmse, best_pair)
}

/// Joint score for a flight pair: one shared clock offset maximises S1 + S2.
#[allow(clippy::too_many_arguments)]
fn score_pair(
    obs: [f64; 3],
    e0: &Prepared,
    e1: &Prepared,
    track0: &Track,
    track1: &Track,
    base_abs: f64,
    shifts: &[f64],
    f0_steps: usize,
    frame_step: usize,
) -> (f32, f32) {
    // Rank on curve agreement in Hz, not brightness: the simulated Doppler must
    // follow the measured ridge. Sweep f0 per event and one shared clock shift.
    const MAX_ERR_HZ: f32 = 25.0;
    let offsets: &[f64] = if frame_step > 1 {
        &[0.0]
    } else {
        &[-3.0, -1.5, 0.0, 1.5, 3.0]
    };
    let c0 = if frame_step > 1 { &e0.coarse_curve } else { &e0.curve };
    let c1 = if frame_step > 1 { &e1.coarse_curve } else { &e1.curve };
    let mut r0 = vec![f64::INFINITY; shifts.len()];
    for &f0 in &score::f0_grid(e0.spec.f0 as f64, f0_steps) {
        for (i, &sh) in shifts.iter().enumerate() {
            for &df in offsets {
                if let Some(v) = score::curve_rmse_offset(
                    c0,
                    track0,
                    obs,
                    f0,
                    base_abs,
                    sh,
                    df,
                    Atmos::default(),
                    MAX_ERR_HZ,
                ) {
                    if v < r0[i] {
                        r0[i] = v;
                    }
                }
            }
        }
    }
    let mut r1 = vec![f64::INFINITY; shifts.len()];
    for &f0 in &score::f0_grid(e1.spec.f0 as f64, f0_steps) {
        for (i, &sh) in shifts.iter().enumerate() {
            for &df in offsets {
                if let Some(v) = score::curve_rmse_offset(
                    c1,
                    track1,
                    obs,
                    f0,
                    base_abs,
                    sh,
                    df,
                    Atmos::default(),
                    MAX_ERR_HZ,
                ) {
                    if v < r1[i] {
                        r1[i] = v;
                    }
                }
            }
        }
    }
    // Report the best joint residual rather than an absolute verdict: how small
    // a residual is achievable depends on how well the audio itself can be
    // measured, so the accept/reject threshold is applied over the whole grid.
    // The two flyovers share one recording start, but each audio CPA time is
    // estimated independently and carries a few seconds of error. Forcing a
    // SINGLE shared shift rejects the correct pair whenever one estimate is
    // biased — which is exactly what the reference case does (14 s apart).
    // Allow the two events to use different shifts, paying a small penalty per
    // second of disagreement so a genuinely shared clock still wins.
    let penalty_hz_per_s = std::env::var("AOS_OFFSET_PENALTY")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(0.02);
    // Shared shift only. A per-event offset variant was tried and measured on the
    // reference case: truth-cell score 0.001 either way, and 411 s vs 286 s. The
    // mismatch is not purely a timing offset, so the extra freedom buys nothing.
    // The two flyovers share one recording start, but each audio CPA time is
    // estimated separately; on the reference case they disagree by 14 s, which a
    // single shared shift cannot represent. Score every (event0, event1) shift
    // pair and charge a small penalty per second of disagreement, so a genuinely
    // shared clock still beats a split one. This is a combination over residuals
    // already computed above — no extra curve evaluations.
    let penalty = std::env::var("AOS_OFFSET_PENALTY")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(0.02);
    let mut best = f32::NEG_INFINITY;
    let mut best_rmse = f64::INFINITY;
    for i in 0..shifts.len() {
        if !r0[i].is_finite() {
            continue;
        }
        for j in 0..shifts.len() {
            if !r1[j].is_finite() {
                continue;
            }
            let mean = (r0[i] + r1[j]) / 2.0;
            let split = (shifts[i] - shifts[j]).abs() * penalty;
            let s = (1.0 / (1.0 + mean + split)) as f32;
            if s > best {
                best = s;
                best_rmse = mean;
            }
        }
    }
    (best, best_rmse as f32)
}

fn union_bounds(region: &[Bounds]) -> Result<Bounds> {
    if region.is_empty() {
        return Err(anyhow!("no region selected"));
    }
    let mut w = f64::INFINITY;
    let mut s = f64::INFINITY;
    let mut e = f64::NEG_INFINITY;
    let mut n = f64::NEG_INFINITY;
    for b in region {
        w = w.min(b.west);
        s = s.min(b.south);
        e = e.max(b.east);
        n = n.max(b.north);
    }
    Ok(Bounds { west: w, south: s, east: e, north: n })
}

/// Great-circle distance in kilometres.
pub fn haversine_km(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6_371.0f64;
    let (la1, la2) = (a.0.to_radians(), b.0.to_radians());
    let h = ((b.0 - a.0).to_radians() * 0.5).sin().powi(2)
        + la1.cos() * la2.cos() * ((b.1 - a.1).to_radians() * 0.5).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

/// ECEF -> geodetic (degrees), good to a few metres.
fn ecef_to_latlon(p: [f64; 3]) -> (f64, f64) {
    const A: f64 = 6_378_137.0;
    const E2: f64 = 0.00669437999014;
    let lon = p[1].atan2(p[0]);
    let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
    let mut lat = p[2].atan2(r * (1.0 - E2));
    for _ in 0..4 {
        let n = A / (1.0 - E2 * lat.sin() * lat.sin()).sqrt();
        let h = r / lat.cos() - n;
        lat = p[2].atan2(r * (1.0 - E2 * n / (n + h)));
    }
    (lat.to_degrees(), lon.to_degrees())
}

/// Bucket index for a point, at `bucket_deg` resolution.
fn bucket_of(lat: f64, lon: f64, bucket_deg: f64) -> (i32, i32) {
    (
        (lat / bucket_deg).floor() as i32,
        (lon / bucket_deg).floor() as i32,
    )
}

/// WGS84 geodetic -> ECEF (metres).
pub fn ecef(lat_deg: f64, lon_deg: f64, alt_m: f64) -> [f64; 3] {
    const A: f64 = 6_378_137.0;
    const E2: f64 = 0.00669437999014;
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();
    let n = A / (1.0 - E2 * lat.sin() * lat.sin()).sqrt();
    let x = (n + alt_m) * lat.cos() * lon.cos();
    let y = (n + alt_m) * lat.cos() * lon.sin();
    let z = (n * (1.0 - E2) + alt_m) * lat.sin();
    [x, y, z]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::FlightPoint;

    fn synth_flight(lat0: f64, lon0: f64, v: f64, alt: f64, tcpa: f64) -> Flight {
        let mut points = Vec::new();
        let mut t = tcpa - 400.0;
        while t <= tcpa + 400.0 {
            let dlat = v * (t - tcpa) / 111_320.0;
            points.push(FlightPoint { t_s: t, lat: lat0 + dlat, lon: lon0, alt_m: alt });
            t += 2.0;
        }
        Flight { icao: "TEST".into(), points }
    }

    #[test]
    fn track_cpa_matches_geometry() {
        let f = synth_flight(39.0, -98.0, 100.0, 1000.0, 1000.0);
        let tr = build_track(&f, CPA_DT, (39.0, -98.0), 50_000.0).unwrap();
        let (d, _) = tr.cpa(ecef(39.0, -98.0, 0.0)).unwrap();
        assert!((d - 1000.0).abs() < 40.0, "got {d}");
    }

    #[test]
    fn curve_residual_prefers_true_receiver() {
        let (lat0, lon0) = (39.0, -98.0);
        let (v, alt, tcpa, f0) = (100.0, 1000.0, 60.0, 120.0);
        let obs = ecef(lat0, lon0, 0.0);
        let f = synth_flight(lat0, lon0, v, alt, tcpa);
        let tr = build_track(&f, SCORE_DT, (lat0, lon0), 50_000.0).unwrap();

        // Measured curve: the exact Doppler the true receiver would hear.
        let curve: Vec<(f64, f32)> = (-50..=50)
            .map(|i| {
                let t = tcpa + i as f64;
                let ratio = score::doppler_ratio(&tr, obs, t, Atmos::default()).unwrap();
                // Recording-relative time (the recording starts at 0 here).
                (t, (f0 * ratio) as f32)
            })
            .collect();
        let spec = EventSpec {
            f0: f0 as f32,
            harmonics: vec![1.0],
            cpa_s: tcpa,
            half_width_s: 50.0,
        };
        let p = Prepared {
            spec,
            slant_m: 1000.0,
            speed_mps: 100.0,
            coarse_curve: curve.clone(),
            curve,
        };
        let tracks = vec![("TEST".to_string(), tr.clone(), tr)];

        let win = (0.0, 200.0);
        let (s_true, _, _) = score_cell(
            obs, (lat0, lon0), &[0u32], std::slice::from_ref(&p), &tracks, Quality::Fine, win,
        );
        let shifted = ecef(lat0 + 0.35, lon0, 0.0);
        let (s_wrong, _, _) = score_cell(
            shifted,
            (lat0 + 0.35, lon0),
            &[0u32],
            std::slice::from_ref(&p),
            &tracks,
            Quality::Fine,
            win,
        );
        assert!(
            s_true > s_wrong * 1.2,
            "true receiver {s_true} should beat displaced {s_wrong}"
        );
    }

    /// A jet crossing a continental box spends ~6 h inside it. The old
    /// `dwell > 4 h` filter (meant for parked transponders) threw away exactly
    /// these flights, which is why a USA-wide search matched the wrong
    /// short-hop aircraft.
    #[test]
    fn transcontinental_passage_is_kept() {
        let mut points = Vec::new();
        let mut t = 0.0;
        while t <= 6.0 * 3600.0 {
            let lat = 30.0 + 250.0 * t / 111_320.0;
            points.push(FlightPoint { t_s: t, lat, lon: -90.0, alt_m: 10_000.0 });
            t += 2.0;
        }
        let f = Flight { icao: "LONG".into(), points };
        let tr = build_track(&f, SCORE_DT, (40.0, -90.0), 3_000_000.0);
        assert!(tr.is_some(), "6 h passage inside the region was discarded");
    }

    /// A continental region cannot be gridded finely enough for the score
    /// kernel, and the search says so instead of returning a confident wrong
    /// state. This pins the arithmetic behind that warning.
    #[test]
    fn continental_region_is_flagged_as_screening_only() {
        let conus = level0_step_m(24.0 * 111_320.0, 59.0 * 111_320.0 * 0.8, COARSE_CELLS);
        assert!(conus > 4_000.0, "CONUS level-0 spacing was only {conus:.0} m");
        let city = level0_step_m(160_000.0, 190_000.0, COARSE_CELLS);
        assert!(city < 1_000.0, "a 190x160 km box gridded at {city:.0} m");
    }

    /// Parked / taxiing transponders must still be rejected - by motion, since
    /// dwell time no longer distinguishes them.
    #[test]
    fn stationary_transponder_is_rejected() {
        let points: Vec<FlightPoint> = (0..(8 * 3600))
            .step_by(30)
            .map(|s| FlightPoint { t_s: s as f64, lat: 40.0, lon: -73.0, alt_m: 9000.0 })
            .collect();
        let f = Flight { icao: "PARKED".into(), points };
        assert!(
            build_track(&f, SCORE_DT, (40.87, -73.52), 200_000.0).is_none(),
            "a stationary transponder was accepted as a flyover"
        );
    }

}

/// Recompute the winning cell's pair and emit the simulated curves.
///
/// Used only for the visual check: the simulated curve is drawn over the
/// measured ridge so a modelling error is obvious by eye.
fn model_curves_at(
    obs: [f64; 3],
    candidates: &[u32],
    prepared: &[Prepared],
    tracks: &[(String, Track, Track)],
    window: (f64, f64),
    shifts: &[f64],
    f0_steps: usize,
) -> Vec<Vec<(f64, f32)>> {
    const MAX_ERR_HZ: f32 = 25.0;
    let mut best: Option<(f64, usize, usize, f64, f64)> = None; // (rmse, i0, i1, base, shift)
    let n = prepared.len().min(2);
    if n == 2 {
        for &ti0 in candidates {
            let t0 = &tracks[ti0 as usize].1;
            let Some((_, cpa0)) = t0.cpa(obs) else { continue };
            let base0 = cpa0 - prepared[0].spec.cpa_s;
            if base0 < window.0 - 120.0 || base0 > window.1 + 120.0 {
                continue;
            }
            for &ti1 in candidates {
                if ti0 == ti1 {
                    continue;
                }
                let t1 = &tracks[ti1 as usize].1;
                let Some((_, cpa1)) = t1.cpa(obs) else { continue };
                if ((cpa1 - prepared[1].spec.cpa_s) - base0).abs() > shift_half_s() + 20.0 {
                    continue;
                }
                for &sh in shifts {
                    let mut sum = 0.0;
                    let mut cnt = 0;
                    for (k, tr) in [t0, t1].iter().enumerate() {
                        let mut best_e = f64::INFINITY;
                        for &f0 in &score::f0_grid(prepared[k].spec.f0 as f64, f0_steps) {
                            if let Some(v) = score::curve_rmse(
                                &prepared[k].curve,
                                tr,
                                obs,
                                f0,
                                base0,
                                sh,
                                Atmos::default(),
                                MAX_ERR_HZ,
                            ) {
                                if v < best_e {
                                    best_e = v;
                                }
                            }
                        }
                        if best_e.is_finite() {
                            sum += best_e;
                            cnt += 1;
                        }
                    }
                    if cnt == 2 {
                        let rmse = sum / 2.0;
                        if best.map_or(true, |(b, _, _, _, _)| rmse < b) {
                            best = Some((rmse, ti0 as usize, ti1 as usize, base0, sh));
                        }
                    }
                }
            }
        }
    }
    let Some((_rmse, i0, i1, base, shift)) = best else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for (k, ti) in [i0, i1].iter().enumerate() {
        let tr = &tracks[*ti].2;
        // Pick the f0 that best explains this event, then sample the model.
        let mut best_f0 = prepared[k].spec.f0 as f64;
        let mut best_e = f64::INFINITY;
        for &f0 in &score::f0_grid(prepared[k].spec.f0 as f64, f0_steps) {
            if let Some(v) = score::curve_rmse(
                &prepared[k].curve,
                tr,
                obs,
                f0,
                base,
                shift,
                Atmos::default(),
                MAX_ERR_HZ,
            ) {
                if v < best_e {
                    best_e = v;
                    best_f0 = f0;
                }
            }
        }
        let mut pts = Vec::new();
        let mut tau = base + prepared[k].curve.first().map(|c| c.0).unwrap_or(0.0) - shift
            - 2000.0 / score::C0;
        for &(t_rel, _) in &prepared[k].curve {
            let model_t = base + t_rel - shift;
            if let Some(r) = score::doppler_ratio_warm(tr, obs, model_t, &mut tau, Atmos::default()) {
                pts.push((t_rel, (best_f0 * r) as f32));
            }
        }
        out.push(pts);
    }
    out
}
