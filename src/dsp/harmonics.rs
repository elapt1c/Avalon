//! Harmonic grouping: consolidate ridges that are integer/rational multiples of
//! one another (e.g. 2f, 3f, or 2f and 3f with an absent fundamental) into a
//! single acoustic event. The lowest measured frequency is the representative.

use std::collections::BTreeMap;

use crate::dsp::median_f32;
use crate::model::{DetectedEvent, Ridge};

/// Max ratio numerator/denominator considered.
const MAX_RATIO: u32 = 8;
/// Relative tolerance on a candidate rational ratio.
const RATIO_TOL: f32 = 0.04;

pub fn group_harmonics(ridges: &[Ridge]) -> Vec<DetectedEvent> {
    let n = ridges.len();
    if n == 0 {
        return Vec::new();
    }

    let mut parent: Vec<usize> = (0..n).collect();

    for i in 0..n {
        for j in (i + 1)..n {
            // Only merge tracks of the same class. A verified flyover curve and
            // an unverified noise track can share a time window and a plausible
            // frequency ratio by coincidence; merging them makes a single click
            // select unrelated junk.
            if ridges[i].plausible != ridges[j].plausible {
                continue;
            }
            if time_overlap(&ridges[i], &ridges[j]) {
                let fi = median_freq(&ridges[i]);
                let fj = median_freq(&ridges[j]);
                if is_harmonic_ratio(fi, fj) {
                    union(&mut parent, i, j);
                }
            }
        }
    }

    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        groups.entry(find(&mut parent, i)).or_default().push(i);
    }

    let mut events: Vec<DetectedEvent> = Vec::new();
    for (_, mut ids) in groups {
        // Representative = the best-fitting plausible ridge; otherwise the
        // best-fitting ridge; otherwise the lowest median frequency. Ranking by
        // "lowest frequency" alone can anchor an event to a window sidelobe.
        ids.sort_by(|&a, &b| {
            let ra = &ridges[a];
            let rb = &ridges[b];
            let key = |r: &crate::model::Ridge| {
                let (plausible, rmse) = match r.fit {
                    Some(f) => (!r.plausible, f.rmse_hz),
                    None => (true, f32::MAX),
                };
                (plausible, rmse as i64, median_freq(r) as i64)
            };
            key(ra).cmp(&key(rb))
        });
        let rep_id = ids[0];
        let rep_freq = median_freq(&ridges[rep_id]);
        let ratios: Vec<f32> = ids
            .iter()
            .map(|&i| {
                let f = median_freq(&ridges[i]);
                if rep_freq > 1e-3 {
                    f / rep_freq
                } else {
                    1.0
                }
            })
            .collect();
        let start_s = ids.iter().map(|&i| ridges[i].start_s).fold(f64::INFINITY, f64::min);
        let end_s = ids.iter().map(|&i| ridges[i].end_s).fold(f64::NEG_INFINITY, f64::max);
        events.push(DetectedEvent {
            id: 0, // assigned later
            ridge_ids: ids,
            representative_ridge_id: rep_id,
            harmonic_ratios: ratios,
            doppler: None,
            start_s,
            end_s,
        });
    }

    events.sort_by(|a, b| a.start_s.partial_cmp(&b.start_s).unwrap_or(std::cmp::Ordering::Equal));
    for (i, e) in events.iter_mut().enumerate() {
        e.id = i;
    }
    events
}

fn time_overlap(a: &Ridge, b: &Ridge) -> bool {
    a.start_s < b.end_s && b.start_s < a.end_s
}

fn median_freq(r: &Ridge) -> f32 {
    let freqs: Vec<f32> = r.samples.iter().map(|s| s.freq_hz).collect();
    median_f32(&freqs)
}

/// True when `a` and `b` are within `RATIO_TOL` of a small rational ratio.
fn is_harmonic_ratio(a: f32, b: f32) -> bool {
    if a <= 1.0 || b <= 1.0 {
        return false;
    }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    let r = hi / lo;
    for q in 1..=MAX_RATIO {
        for p in q..=MAX_RATIO {
            let cand = p as f32 / q as f32;
            if (cand - r).abs() / r < RATIO_TOL {
                return true;
            }
        }
    }
    false
}

fn find(parent: &mut [usize], x: usize) -> usize {
    if parent[x] != x {
        parent[x] = find(parent, parent[x]);
    }
    parent[x]
}

fn union(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (find(parent, a), find(parent, b));
    if ra != rb {
        parent[rb] = ra;
    }
}
