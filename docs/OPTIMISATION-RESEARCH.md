# Optimisation research — Avalon

Goal: recover **where** a recording was made and **when**, when neither is known,
over a search region and a time range. Accuracy target: **≤ 10 km**.

Everything below is measured on the reference recording
(`Unknown location - raw audio - colsto (128p).mp4`, 3600 s, two flyovers at
t ≈ 2149.6 s and t ≈ 3071.4 s, Δt = 921.9 s) unless stated otherwise.

## 1. What the current pipeline actually achieves

| Setup | Result |
|---|---|
| Production scorer, CONUS grid, known pair `a67ba8,a659a2` | 118 km from truth, 35.4 Hz |
| Same, search's own #1 pick | 724 km (North Carolina) |
| Continuous nuisance fit (shift, f0, offset) at the true receiver | 8.3 – 12.3 Hz |
| Dense shift sweep, 15 km cells, truth cell | **1.0 km**, 11.1 Hz |
| Competing minimum 31 km away, shift −48 s | **3.4 Hz** — better than truth |
| Production grid at truth (5 shifts × 5 f0) | 17.3 Hz, hits the −10 s wall |

So the information is present (a 1 km solution exists), but the scorer prefers a
wrong location because **clock shift trades against position**, and the coarse
grid (5 shifts, 5 f0 values) cannot resolve the trade.

## 2. How the reference implementations solve the same problem

`mjshiggins/dopplertrack` (MIT, reimplementation of the same technique) reports on
this same video: **43 m** after refinement, 200 m blind grid, true pair ranked #1
of 328.

The differences that matter, in order of impact:

1. **Propeller aircraft only.** They drop known jets using ICAO Doc 8643 engine
   classes before matching. The tones are propeller blade-pass (60–400 Hz); a
   jet cannot produce them. Filtering removes thousands of impossible candidates
   and most of the score surface.
2. **One clock, ≤ half a second.** In their known-point check the correct pair's
   curves lie on the measured ridges "with no free parameters beyond a half-second
   time shift". A 30–60 s slide is not permitted — and our ±60 s allowance is
   exactly what creates the 31 km / 3.4 Hz false minimum.
3. **Δt is the candidate constraint, not a by-product.** Candidates must be points
   that hear an event-1-like pass and *exactly Δt later* an event-2-like pass. With
   an unknown recording start, Δt is the only absolute-time anchor available.
4. **Silence filter.** Reject candidates with another prop flyover in between,
   closer than max(slant₁, slant₂) — a hearing heuristic, configurable.
5. **Search order.** Event 1's placement absorbs measurement error; event 2 is
   locked to the flight geometry. Every further event is matched once, then
   scored with event 1's shift, so all events share one clock.
6. **Grid.** 200 m coarse, 50 m refine, over candidate (point, pair) triples that
   survive the filters — not "score every cell, then refine the winners".

`SITREC` (Metabunk guide) states the same invariants: a wrong flight pair cannot be
rescued by sliding its curves independently; a single flyover needs a known
recording-start UTC, two flyovers do not; the delay between flyovers is the useful
quantity when the start time is unknown.

`DOPPLERGEIST` methodology adds: evaluate radial velocity at **emission** time with
the arrival-time equation inverted (a constant delay is not equivalent to a receiver-
time model), split flight gaps rather than extrapolating, and treat raw scores as
ranking statistics, not calibrated probabilities.

## 3. Three findings specific to Avalon (highest priority)

1. **The search window was the wrong date.** CONFIRMED with colsto: the recording is
   **2026-09-08**, Oyster Bay NY (40.876466, −73.524951), and the aircraft are
   **bushplanes/turboprops**, not jets. Our archive window was 2026-10-30/31 — seven
   weeks late — so the true flights were never in the candidate set and the scorer
   locked onto the best-fitting geometry from the wrong day. Our measured Δt of
   921.9 s matches the published 918 s within fit error, which was the tell.
   The 2026-09-08 archive (3.87 GB) is now cached; `AVALON_WIN` selects the window.
2. **The old "known pair" was two jets.** `a67ba8` = N517CF, Cessna **Citation X**;
   `a659a2` = N508WF, Gulfstream **G450**. Turbofans do not produce a clean 100 s
   blade-pass sweep at 119.6 / 136.6 Hz; turboprops do (PC-12 ≈ 113 Hz, King Air 350
   ≈ 133 Hz at typical RPM). That explains the 8–12 Hz residual floor, f0 pinned at
   its ±3 % bound and the 30–50 s shifts the fit invented. The correct pair for
   2026-09-08 is a Pilatus PC-12 (`a1dfeb`, N22NY) and a King Air 350 (`ace391`,
   N93CQ) — both turboprops, matching colsto's "bushplanes".
3. **Jet/type information is not used anywhere in candidate generation.** Even on
   the correct day, including jets inflates the candidate set by ~an order of
   magnitude and adds false minima.

## 3b. Concrete parameters worth copying (dopplertrack 2026-11-03)

| Knob | Value | Why |
|---|---|---|
| Prop-only filter | ICAO Doc 8643 engine class (Piston/Turboprop), drop Jets | tones are blade-pass |
| Δt tolerance (`tol_s`) | 6 s between heard CPAs | event pairing, unknown epoch |
| Silence filter | 6 km / 40 s margin | no intervening closer pass |
| Max slant range | 12 km | beyond this the tone is usually inaudible |
| Candidate grid | 1 km coarse → 200 m → 50 m refine | 379k points in 19 min on a laptop |
| Clock | shared, ≤ 0.5 s | one recording clock for all events |
| Speed of sound | 340–343 m/s, no wind term | 3 % of c, same order as residuals |
| Index | 1 Hz ENU segments, time buckets + per-bucket KD-tree | O(points × flights) → O(points × nearby) |
| Scoring | **ridge mode**: per reliable frame `exp(-0.5 (df/σf)²)`, σ = 1.5 % | brightness cannot separate candidates on a faint tone; ridge position can |
| Ridge reliability | ≥ 3 dB bright AND ≥ 3 dB more than the best alternative path | rejects birds/second tones |

Ridge-mode scoring is the one change that most directly attacks our 3.4 Hz false
minimum: it scores *where the curve is*, not how bright it is, so a wrong location
that fits the amplitude cannot win.

## 4. Optimisation backlog, ranked by expected accuracy gain per unit work

| # | Change | Why | Expected |
|---|---|---|---|
| 1 | Correct archive date; re-benchmark | Wrong-day flights ⇒ any answer is meaningless | resolves the dominant error |
| 2 | Prop-only candidate filter (ICAO type/engine class) | Tones are blade-pass; jets cannot match | removes most false minima |
| 3 | Shared clock bounded to ±1 s, selected by CPA-slope alignment on the measured ridge | Kills the shift/position degeneracy that produced the 31 km / 3.4 Hz minimum | 31 km → inside 10 km |
| 4 | Hard Δt consistency filter during candidate generation | Only absolute-time anchor when the start epoch is unknown | cuts (pair, point) space ~10² |
| 5 | Silence filter (no intervening closer prop pass) | Removes physically impossible placements | fewer false positives |
| 6 | Candidate-first search: enumerate (point, pair) with CPA/slant/Δt, then score only survivors | Currently scores orders of magnitude more cells than needed | ~10× throughput |
| 7 | 200 m coarse grid → 50 m local refine on survivors | Current coarse pass (15 km, 5 shifts) starves itself — at 60 km spacing the shortlist collapsed to 1 cell | better recall + speed |
| 8 | Precompute per-flight Hermite splines; time-bucketed spatial index; rayon over candidates | The inner loop is retarded-time evaluation | 5–20× throughput |
| 9 | Harmonic-comb scoring with per-event f0 bounded ±3 % | Uses 2–3× more ridge evidence per event; stops f0 absorbing geometry error | tighter residuals |
| 10 | Estimate wind + sound speed (c = 343 m/s, zero wind today) | 10 m/s wind = 3 % of c, same order as the residuals being ranked on | fewer near-ties |
| 11 | Score calibration: report near-tied regions with margins instead of a single point | Two unrelated pairs fit to hundredths of a Hz at continental scale | honest uncertainty |

## 5. Throughput targets for "region × time range" search

The user-facing goal is a search over an **area** and a **time range**, neither exact.

- Stage A (cheap): filter the archive by aircraft class and by the *Δt* between
  events alone, ignoring position. This is a 1-D problem per flight pair and can run
  over a whole day of traffic in seconds.
- Stage B (geometry): for surviving pairs, intersect the two slant-range circles
  implied by the curve fits (CPA depths), then test only those points.
- Stage C (scoring): exact retarded-time curve match at 200 m on survivors, 50 m
  refine on the best few, all events locked to one clock.
- Stage D (report): ranked regions with margins, not a single point.

With that ordering the expensive scoring touches thousands of candidates rather
than millions, which is what makes a continental area × 24 h window tractable.

## 6. Reference numbers to beat

| Metric | Reference | Avalon today |
|---|---|---|
| Position error, blind 100×150 km box | 43 m (dopplertrack) | 118 km (wrong pair/day) |
| True pair rank | #1 of 328 | not measured cleanly |
| Scoring cost | 19 min laptop for 379,758 points × 3,068 segments | 322 s CONUS after the filter fixes |
| Clock freedom | ≤ 0.5 s | ±60 s |
