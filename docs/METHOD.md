# Method — reconstructed from DOPPLERGEIST

This documents the aircraft-Doppler acoustic geolocation method, reconstructed
from a study of the [DOPPLERGEIST](https://github.com/1337isnot1337/DOPPLERGEIST)
source and mapped to this project's Rust modules. It is a *reconstruction*, not
copied code.

## Overview

A moving aircraft produces a frequency-shifted tone whose observed frequency
traces a characteristic curve as the aircraft approaches, passes, and recedes.
Given a recording, the method:

1. **Detects** moving tonal tracks (ridges) in the spectrogram.
2. **Groups** ridges that are harmonics of one flyover.
3. **Parameterizes** each flyover: rest frequency, speed, closest-approach
   distance and time.
4. *(Later phase)* **Matches** against historical ADS-B trajectories and
   **ranks** candidate microphone locations.

This version implements steps 1–3.

## 1. Acoustic front-end (`dsp/spectrogram.rs`)

- Decode to mono PCM, resample to a target rate (default 8000 Hz).
- STFT: 2 s Hann window, 0.5 s hop, 40–1600 Hz band (configurable).
- Magnitude spectrum, then **time-whitened spectral contrast**: divide each
  frequency bin by its time-median magnitude. Stationary tones and the noise
  floor are normalized away; a moving tone passing through a bin stands out as
  a contrast spike.

## 2. Ridge detection (`dsp/ridge.rs`)

- Per frame, find local maxima of the contrast spectrogram above a seed
  threshold, then apply **non-maximum suppression at the true resolution**
  (~0.75 × `1/window`). Zero-padding makes many correlated bins per Hann
  mainlobe; without NMS one tone reads as dozens of "tonals".
- A per-frame **dynamic-range floor** (−26 dB below the frame's strongest
  component) keeps Hann sidelobes from being detected as tones.
- Tracks are found by **Viterbi dynamic programming** over the sparse peak set,
  maximizing accumulated `(contrast − seed)` minus a per-frame cost
  (`4 × seed`) and a small slope penalty. Constraints: slope ≤ 400 Hz/s,
  curvature ≤ 20 Hz/frame, monotonic descent (≤ +2 Hz/frame), gaps ≤ 2 frames.
  The per-frame cost is what stops the DP from stringing unrelated loud peaks
  into one arbitrary mega-path.
- After each track is found its peaks are **soft-masked** (gain × 0.2) and the
  DP re-solves, so distinct curves are recovered while a real curve crossing a
  junk path remains reachable.
- Rejected at build time: too short (< 15 s), too little sweep (< 4 Hz or
  < 2 % of frequency), upward-trending, or low contrast.
- **Model-guided extension.** A Viterbi path stops where the tone fades relative
  to the per-frame cost, usually mid-sweep, which also pins `(v, d)` badly. Each
  plausible track is extended frame by frame along the observed tone (monotonic
  descent, local ±3 % window; the fitted curve is only a tie-breaker),
  alternating with refits, and the extension is kept only if the result is still
  a plausible flyover. This recovers the flat recession tail, where an aircraft
  spends most of its observable time.
- **Block processing.** Recordings longer than 300 s are analysed in 300 s blocks
  with 150 s overlap — enough that any sweep under 150 s sits fully inside some
  block — then de-duplicated. A one-hour spectrogram dilutes a faint curve
  twice: whole-file baselines raise the contrast floor, and a single DP pass
  never reaches it.

Physical basis for the slope limit: a jet at 250 m/s and 500 m slant gives
`|df/dt| ≈ f0·v²/(c·d) ≈ 360 Hz/s` near 1 kHz.

## 3. Harmonic grouping (`dsp/harmonics.rs`)

- Ridges whose median frequencies are within tolerance of a small rational
  ratio (e.g. 2:1, 3:1, or 2:3 with an absent fundamental) are grouped as one
  event. Overlapping time is required.
- The lowest measured frequency is the representative; the fundamental is never
  invented. A ridge without a partner remains a single-tone event.

## 4. Doppler parameterization (`dsp/doppler.rs`)

Physical model (stationary receiver, straight constant-velocity pass, retarded
sound-arrival time neglected since `v << c`):

```
f(t) = f0 · c / (c + u_r(t))
u_r(t) = v·x / r,   x = v·(t − t_cpa),   r = √(d² + x²)
```

- `f0` rest (emitted) frequency, `v` speed, `d` closest-approach slant distance,
  `t_cpa` closest-approach time, `c` sound speed.
- Initial estimates from the asymptotic extremes:
  `f0 = 2·fmax·fmin/(fmax+fmin)`, `v = c·(fmax−fmin)/(fmax+fmin)`, `t_cpa` from
  the `f0` crossing.
- Refined with a derivative-free Nelder-Mead fit (scale-tolerant).
- Reports fit RMSE and whether both approach and recession were observed
  (complete vs. partial sweep).

These estimates are the same first-stage "bound speed and slant distance"
quantities DOPPLERGEIST derives from a ridge; they are *not* calibrated location
probabilities.

## 5. Flight matching & search (`flights.rs`, `dem.rs`, `score.rs`, `search.rs`)

- **Acquire** ADSB.lol daily archives (GitHub releases) into a local cache with
  HTTP-Range resume, then parse `trace_full_*.json` members and keep flights that
  cross the selected region inside the time window.
- **Terrain.** Receiver elevation comes from Terrarium DEM tiles (~30 m at zoom
  12, `h = R*256 + G + B/256 - 32768`), sampled bilinearly, observer at
  `elev + 2 m`. Sea level would inject up to ~2 km of slant-range bias. The zoom
  is chosen from the region size (z12 for a city, z8 for a continent) so a large
  box does not need tens of thousands of tile requests.
- **Tracks.** Positions and velocities in WGS84 ECEF, Catmull-Rom interpolated on
  a uniform grid (5 s for closest-approach pruning, 2 s for scoring). Positions
  are stored as `f32` offsets from a per-track origin: a continental window needs
  tens of thousands of tracks and the track store was the dominant allocation.
- **Airborne vs. parked.** A flight is only usable if it is airborne (mean
  altitude > 150 m) *and moving* (mean ground speed > 30 m/s over its passage).
  Dwell time is deliberately **not** a filter: inside a continental box a
  transcontinental jet legitimately stays for 6 h, and rejecting long passages
  discarded exactly the flights a coast-to-coast search is looking for.
- **Candidate search is hierarchical, not a single grid.** The level-0 grid
  covers the whole region at a spacing set by a cell budget (`COARSE_CELLS`),
  never finer than 200 m. Its best cells are then sub-divided 5x per side and
  re-scored, repeatedly, with scoring fidelity rising from coarse
  (1-in-16 frames, 5 f0 steps, 3 clock steps) to fine (every frame, 25 f0 steps,
  17 clock steps). Cells replaced by children are dropped from the ranking so a
  cheap coarse score cannot outrank its own children. Final spacing is
  `level0 / 5^levels`, capped at ~60 m; for a 200 km box that is 42 m, for CONUS
  39 m. The old code computed one spacing from the region size and "refined" by
  re-scoring the same cells at finer *time* resolution, so the answer could never
  be better than that spacing (15 km for CONUS).
- **Pruning.** A flight is only considered for a cell when its closest approach
  is within 40 km, its ground speed is within 0.5–2x the fitted speed, and its
  CPA distance within 0.1–8x the fitted slant (the fit's slant is degenerate --
  a real 10.7 km pass was fitted at 5.9 km -- so this band stays wide).
- **Pairing is by clock, never by speed.** Both events were cut from one
  recording, so the same absolute recording start must explain both flyovers.
  Candidates for the second event are bucketed by implied start and joined to the
  first event inside a +-12 s window; the surviving pairs are ranked by a cheap
  residual on the subsampled curves before the full sweep. Speed-ranked
  truncation was removed: the fitted speed is degenerate, so ranking by it put
  the true flights 122nd of 172 and discarded them before scoring.
- **Scoring (curve residual).** The measured ridge core is compared with the
  modelled curve, `f0 * (c - w.d)/(c + v.d - w.d)` with retarded time solved by
  fixed point. Sweeps: `f0` (default +-3 %, `AOS_F0_HALF` overrides), a shared
  clock offset +-10 s, and a per-event constant frequency offset of +-3 Hz.
  Residuals are capped at 25 Hz per sample and a hypothesis is rejected if more
  than 20 % of the curve cannot be placed or if the modelled frequency excursion
  disagrees with the measured one by more than 25 %.
- **Selectivity.** Only cells whose residual is within 1.6x the best residual on
  the grid are kept, so the threshold adapts to how cleanly the audio could be
  measured.
- **Modes, not a point.** Surviving cells are clustered spatially (15 km link)
  into modes, and every mode is reported with its best cell, residual, flight
  pair, cell count and radius. See "continental limit" below for why.

## Continental limit — measured, not assumed

A 60 000-cell budget over CONUS is ~15 km spacing, and a uniform score grid
there cannot converge: the score kernel is only a few hundred metres to a few
kilometres wide, so a 15-24 km grid usually misses the true location *entirely*
(measured: 149 of 24 198 level-0 cells produced any hypothesis, and the true cell
was among neither). Raising the budget does not help -- resolving a continent at
the kernel width needs ~10^6-10^7 cells.

The deeper problem is identifiability. With two events and a free rest
frequency, unrelated flight pairs fit the same two measured curves to within
hundredths of a Hz, hundreds of kilometres apart. Measured on the reference
clip: the true pair reaches ~6 Hz at its best position in the New York area,
while a coincidental pair over the Gulf of Mexico reaches 1.4 Hz. The
DOPPLERGEIST benchmark reports the same failure at national scale (top
hypothesis 1 046 km from truth, true pair 4th on a near-tied score).

So: a small or medium region (up to a few hundred km) localises, and the
hierarchy turns that into metres. A continental region is a *screening* pass
only -- the tool reports the near-tied modes and the margin between them, and
the user must narrow the region (or supply more flyovers) to decide.

## Key assumptions & limits

Constant sound speed, stationary microphone, stable source tone, adequate
sample-rate/timestamp accuracy. Wind, refraction, terrain obstruction, and
changing RPM can break matches. Raw scores rank, they are not probabilities.
