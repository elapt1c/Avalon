# Overnight log — 2026-11-09

Working notes for the accuracy push. The goal stated by the operator: one
obvious point with no credible alternative, not a 6 km plateau.

## Where the accuracy actually goes missing

Chain, with evidence from the reference case (colsto clip, 2026-09-08, receiver
40.876466, −73.524951):

1. **Both published aircraft are in the archive and both overfly the receiver.**
   `probeicao` (examples/probeicao.rs) prints, from the archive alone:
   - `a1dfeb` (N22NY) closest **2.8 km**, altitude 5.5 km, at the first flyover time
   - `ace391` (N93CQ) closest **3.5 km**, altitude 1.1 km, at the second
   Their flyovers are **942 s apart**; our two measured audio ridges are
   **921.9 s apart**. That 20 s gap is the thread to pull.

2. **One of them used to be thrown away before scoring** — *fixed this session*.
   Track building rejected a flight when its mean altitude over the whole
   in-region slice was under 150 m. A King Air that transits the box and then
   parks or taxis inside it averages well under that, so it vanished. Judged
   around the closest approach instead, both survive, and the pair is now
   present in the candidate list at the truth cell (debug output: "IS in
   candidate list").

3. **The pair is scored but fits badly.** At the truth cell:
   - `CELL_PAIRS=4` (old default): score **0.000** — the pair never reaches the
     full sweep, it is truncated by the cheap pre-rank.
   - `CELL_PAIRS=40`: score **0.001** — scored, but 100× weaker than the
     0.104 winner (a spurious pair 16 km away).

   So the remaining loss is not gating and not pairing. It is the **fit**.

4. **Why the fit is bad: the audio CPA time is wrong by 15–20 s.**
   Implied recording start from each aircraft, at the truth cell:

   | aircraft | event | implied start |
   |---|---|---|
   | a1dfeb | 0 | 1788898636.0 |
   | ace391 | 1 | 1788898650.0 |

   They must be equal — it is one recording. The 14 s disagreement means one of
   the two audio CPA estimates (2149.56 s / 3071.41 s from the ridge extractor)
   is off, or the ADS-B CPA time is (it is searched with a 30 s stride, then
   refined, so a few seconds of quantisation are expected).

## Ranked next experiments

1. **Precise audio CPA.** Fit the measured ridge to the hyperbola model and take
   its vertex, rather than using the ridge's minimum-|df/dt| frame. Then check
   whether the implied starts for `a1dfeb`/`ace391` converge. This is the
   highest-value fix: everything downstream is fighting a 15 s time error.
2. **Per-event clock offset with a prior.** The physical model allows one shared
   start, but the estimator has seconds of error. Score pairs with
   `off0, off1` independent inside ±20 s, penalised by |off|, instead of one
   shared shift. Keep the shared-shift path as the fast pre-filter.
3. **Raise CELL_PAIRS for the final refinement level only.** 4 → 40 costs ~15 %
   more time at 200 m and is what let the true pair be scored at all. Make it
   quality-dependent (coarse 4, fine 40) rather than a global constant.
4. **Silence test tuning.** The current audible radius is `max(d0,d1)*2` ≥ 4 km.
   At the true location over Long Island Sound this may still veto the right
   answer; check the count of "extra audible passes" at the truth cell.
5. **Blade-pass f0 for these airframes.** PC-12 ≈ 113 Hz and B350 ≈ 133 Hz were
   assumed; the measured ridges sit at 119.6 and 136.6 Hz. If the f0 sweep range
   does not comfortably bracket the true value while the shift is wrong, the fit
   lands in a wrong basin.

## Reproduce

```bash
cd /root/projects/AmazingOSINT/AmazingOSINT
AVALON_CACHE=/root/projects/aosint-cache \
AVALON_DEBUG_CELL="40.8765,-73.5250" AVALON_DEBUG_ICAO="a1dfeb,ace391" \
AOS_PAIR_TOL_S=40 AOS_SHIFT_HALF_S=20 AOS_CELL_PAIRS=40 \
./target/release/examples/bench --analysis /root/projects/analysis-default.json \
  --time 1788890400:1788910200 --truth 40.876466,-73.524951 \
  "tiny:-73.72,40.72,-73.42,41.02"

AVALON_DEBUG_ICAO="a1dfeb,ace391" ./target/release/examples/probeicao \
  "40.876466,-73.524951"
```

## Harness knobs

| env var | meaning |
|---|---|
| `AVALON_DEBUG_CELL` | `lat,lon` — trace every gate decision at that cell |
| `AVALON_DEBUG_ICAO` | comma list — report where those aircraft enter/leave |
| `AOS_PAIR_TOL_S` | shared-start tolerance when pairing two flyovers |
| `AOS_SHIFT_HALF_S` | half-width of the f0/clock sweep |
| `AOS_CELL_PAIRS` | flight pairs scored per cell (4 old, 40 diagnostic) |
| `AOS_BLOCK_WORKERS` | analysis block concurrency |
| `AOS_LEVEL0_STEP_M` / `AOS_TARGET_STEP_M` | candidate grid spacing |
| `AVALON_PARSE_SERIAL` | restore serial trace parsing |
