# Power-grid (ENF) timing — research note

Goal: use the **mains hum** in a recording to pin the absolute time (and coarsely the
region), so the Doppler search starts from a time prior instead of scanning a whole
window.

## What the technique is

**Electric Network Frequency (ENF) analysis.** The grid's actual frequency is never exactly
50.000 / 60.000 Hz — it wanders by tens of millihertz as load and generation move, and the
*pattern* of that wander is effectively a fingerprint of that grid over that day. Anything
mains-powered (or a mic near wiring) picks the hum up, so an audio/video recording carries
a timestamped copy of it.

Two separate claims, with very different maturity:

| Claim | Status | Accuracy |
|---|---|---|
| **Time of recording** — match the hum trace against a reference grid-frequency log | Established forensic technique (ENFSI best-practice manuals, commercial + open tooling) | **~1 second** with a clean trace and matching reference data |
| **Grid identification** — which interconnection (50 Hz vs 60 Hz; 9 grids in the SP Cup 2016 set) | Published, ~92 % classification on that benchmark | country / multi-state scale |
| **Intra-grid localisation** — position *inside* one grid | Research (Garg, Hajj-Ahmad & Wu, arXiv:2105.00668) — correlation of the high-passed trace decays with distance, trilateration against anchor recordings | "a certain geographical resolution"; needs strong ENF traces and anchor sites |

Implementation references: `openENF/openenf` (open-source CLI, adaptive Goertzel extraction,
single-second timestamping — note its **reference grid data is not in the repo**, it is
licensed separately), FNET/GridEye (UT Knoxville) for US frequency data, and the SP Cup
2016 dataset for grid classification work.

## Why it matters for Avalon

Our hardest unresolved parameter is **the recording clock**: with only two flyovers we
have to allow a shared time offset of tens of seconds, and that freedom is exactly what
produced the 31 km / 3.4 Hz false minimum. ENF would supply the absolute time directly:

- a ±1 s time anchor removes the shift parameter entirely;
- even a ±60 s anchor narrows the candidate window enough to make pairing cheap;
- grid identification (60 Hz + Eastern Interconnection) shrinks the search region to a
  continent-scale prior before any Doppler work.

## Measured on the reference recording — not usable as-is

Extracted the 40–200 Hz band from the reference clip and tracked the narrow peak near
60 Hz (10 s and 30 s windows, parabolic interpolation on the log spectrum):

| Window | Median SNR of the 60 Hz line | Frequency scatter | Windows good enough (SNR > 6) |
|---|---|---|---|
| 10 s | 2.7 | ±806 mHz | 0 of 59 |
| 30 s | 3.2 | ±555 mHz | 1 of 19 |

A real ENF trace wanders by **20–50 mHz** and stands well above the noise floor; ours
scatters by **half a hertz** and the estimator loses the line entirely in a third of the
windows (58.6 Hz, 61.3 Hz — pure noise). Occasional runs do sit within ±30 mHz of each
other, so a faint trace is probably present, but this is a 128p YouTube transcode: lossy
compression plus resampling has smeared a hum that was likely already 60+ dB below the
programme material.

**Conclusion: the concept is sound and worth building, but this file cannot support it.**
A high-bitrate original (or a recording made with the mic plugged in, near wiring) is
required.

## What it would take to ship

1. **Source audio at full quality.** Hum needs to be resolvable to a few mHz over tens of
   seconds; that is not recoverable from a 128p transcode.
2. **Reference grid data for the date and grid in question.** This is the real blocker:
   public, queryable archives are thin (FNET/GridEye is primarily a research system;
   openENF's data is separately licensed). Without a reference trace there is no match,
   only a 50/60 Hz test.
3. **A time-prior stage in the pipeline.** Extract ENF → estimate the trace → correlate
   against the reference to produce `(absolute time, uncertainty)` → weight the Doppler
   candidate windows instead of scanning them uniformly.
4. **Grid identification before geometry.** 50 Hz vs 60 Hz alone removes half the world;
   the interconnection narrows it further.

## Honest expectation

Even with perfect ENF, it gives **time** very well and **region** coarsely — not a position.
It composes with the Doppler method rather than replacing it: ENF kills the clock
ambiguity and prunes the region, Doppler then localises to metres-to-kilometres. Treat the
intra-grid trilateration result as research, not a product feature.
