# ENF reference data — sources and measured limits

Companion to `POWER-GRID-ENF.md`. Question answered here: **can we actually get
reference grid-frequency data, and how much does a degraded trace cost us?**

## Verified sources

| Source | Coverage | Resolution | History | Access | Verified |
|---|---|---|---|---|---|
| **NESO (National Grid ESO, GB)** `api.neso.energy` dataset `fnew-YYYY-M.csv` | Great Britain, one synchronous area | **1 s** | monthly files back to **Jan 2020** | free, no key, direct HTTP | **yes — downloaded Sept 2026: 69 MB, 2,592,000 rows, range 49.744–50.249 Hz, std 81.7 mHz, 1 mHz quantisation** |
| **Power-Grid Frequency Database** `power-grid-frequency.org` (OSF + arXiv:2006.01771) | 12 synchronous areas; Continental Europe 2019, Baltic 2025, plus standalone recordings | 0.2–1 s | fixed campaign periods, not continuous | free, citation required | listed, not fetched |
| **AEMO (Australia)** annual/quarterly frequency files | NEM regions | 4 s | multi-year | free CSV | listed, not fetched |
| **Fingrid (Finland)** open data API, dataset 177 | Nordic | 3 min | API, historical | free API key | listed, not fetched |
| **SP Cup 2016 (IEEE)** via IEEE DataPort | 9 grids, for grid-ID classification | high | 2014–2015 | free account | listed, not fetched |
| **FNET/GridEye (US)** UT Knoxville | North America, 300+ FDRs | 0.1 s class | research archive | **request to the lab** (`liu@utk.edu`) | gap confirmed |

### The US is the gap

There is **no free, queryable, 1-second historical US frequency archive.** PJM and ERCOT
publish load and telemetry, not frequency; EIA has no frequency product; FNET/GridEye is a
research system whose archive is obtained by request. For a North American target the
realistic options are:

1. an FNET/GridEye data request (needs a relationship with the lab, and coverage of the
   specific date must be confirmed),
2. running our own GPS-disciplined frequency recorder (works only forward in time),
3. a commercial PMU/synchrophasor data provider.

GB, Continental Europe, the Nordics and Australia are all straightforward — the technique
can be built and validated today on any of those grids, it just cannot be *applied* to a
US recording until reference data exists for that date.

## How much trace do we actually need?

Experiment on the real GB series (whole month of 1-second data, 2.59 M samples). Cut a
slice out of it, degrade it the way acoustic extraction would (quantise to 10 mm, add
Gaussian noise), then cross-correlate against the entire month. "Margin" is the score of
the winner minus the best *wrong* position, in percent.

| Slice length | Noise | Result | Margin |
|---|---|---|---|
| 10 min | 10 mHz | **wrong by 15.8 days** | 0.20 % |
| 20 min | 10 mHz | **wrong by 6.4 days** | 0.02 % |
| 30 min | 10 mHz | **wrong by 14.1 days** | 2.15 % |
| 45 min | 10 mHz | **wrong by 2.8 days** | 1.70 % |
| **60 min** | 10 mHz | **exact second** | 9.83 % |
| **60 min** | 30 mHz | **exact second** | 9.34 % |
| **60 min** | 50 mHz | **exact second** | 8.5 % |

Two things matter here:

1. **An hour of trace pins the absolute time to the second**, and keeps doing it with
   ~50 mHz of extraction noise and 10 mm quantisation — i.e. the hum only has to be
   measured to a tenth of a hertz, not to millihertz.
2. **Below ~1 hour the answer is confidently wrong.** At 30 and 45 minutes the margin is
   still a couple of percent, but the reported time is days off. The margin is *not* a
   validity check in that regime — the same failure mode as our Doppler near-tie, where a
   plausible peak was 31 km out.

Consequence for Avalon: the reference recording is ~56 minutes long, which is right at
the threshold. A single-hour ENF trace is exactly the shape of data that works — but only
if the hum is recoverable, which the 128p transcode does not permit.

## Recipe if we get a usable source audio

1. Confirm the grid (50 vs 60 Hz; then interconnection) from the hum frequency.
2. Pull the monthly reference file for that grid and month.
3. Extract ENF from the full-quality audio in ≥60 min windows, coarse quantisation is fine.
4. Sliding normalised cross-correlation against the reference → best offset = absolute
   time; report the margin, and require ≥5 % before trusting it.
5. Cross-check with a second, non-overlapping window from the same recording — two
   independent windows agreeing on the same second is a real lock; one window is a guess.
6. Feed `(time, ±uncertainty)` into the Doppler search as a prior instead of scanning.

## Bottom line

The data problem is solvable **outside the US** today with a free 1-second archive and no
special access. Inside the US it needs either a research data request or our own receiver.
And the method has a hard requirement we did not expect: **a full hour of hum**, not a few
minutes — shorter traces don't degrade gracefully, they lie.
