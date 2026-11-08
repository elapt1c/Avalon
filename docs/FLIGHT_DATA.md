# Historical flight data — sources & acquisition plan

For the **matching phase** (compare detected Doppler signatures against real
trajectories), we need historical ADS-B/aircraft track data for the recording's
time window.

## Primary: ADSB.lol (recommended)

Free, no auth, effectively no rate limit on the heavy downloads. This is the
source DOPPLERGEIST uses.

- Data lives as **GitHub releases**, one per UTC day, in repos
  `adsblol/globe_history_YYYY` (2023…2026, rolling).
- Release tag format:
  `v{YYYY.MM.DD}-planes-readsb-{provider}` where `provider ∈ {staging-0, prod-0, mlatonly-0}`.
- Metadata lookup (1 request/day): `GET https://api.github.com/repos/adsblol/globe_history_{year}/releases/tags/{tag}`.
  GitHub API anonymous limit is 60 req/hr — irrelevant here (one lookup per day).
- Archive assets: a single `.tar` (mlatonly) or split parts `.tar.aa`/`.tar.ab`.
  A full `prod-0` day ≈ **3.5 GB** (split into ~2 GB parts). Downloaded from the
  asset CDN (`browser_download_url` → release-assets.githubusercontent.com) which
  has **no rate limit**.
- Contents: gzip'd [readsb](https://github.com/wiedehopf/readsb) JSON traces
  (one JSON object per line/member, with per-sample lat/lon/alt/vel/squawk).
- Aircraft type designators: `https://raw.githubusercontent.com/wiedehopf/tar1090-db/master/db/icao_aircraft_types.js`.
- License: **ODbL** (open database) — usable, with attribution/share-alike.

Verified live 2026-11-04: `v2026.11.03-planes-readsb-prod-0` → `-prod-0.tar.aa`
(2000 MB) + `-prod-0.tar.ab` (1531.6 MB).

### Why no proxy bypass is needed

The only GitHub-API call is the per-day release-metadata lookup (60/hr anonymous
is ample). The gigabytes come over the asset CDN, which is unthrottled. Proxy
rotation is only relevant if we later hit GitHub API aggressively (e.g. backfilling
many days), where a GH token (5000 req/hr) is cleaner than proxies.

## Fallback: OpenSky Network

- Historical API: `https://opensky-network.org/api/flights/aircraft?icao24=…&begin=…&end=…`
  and `/api/flights/interval`.
- Anonymous: ~100–400 credits/day (rate-limited); authenticated account gets more.
- Bulk access via the "Impala" shell requires academic/contributor membership.
- Fine as a spot-fallback, but not suitable for bulk day-scale search without an
  account. IP-based limits make proxy rotation a ToS-grey area.

## Acquisition design (implemented)

1. Given a recording's UTC window, determine the covering UTC day(s).
2. For each day, GET the release metadata, then download the `.tar.aa`/`.tar.ab`
   parts into a local cache, resuming partial files with HTTP Range requests.
3. Parse the cached files (gzip'd readsb traces) and keep the flights that cross
   the region inside the time window.

Cache location: `$AVALON_CACHE`, else `$XDG_CACHE_HOME/avalon/globe_history`,
else `~/.cache/avalon/globe_history`. One day is ~2.9 GB; re-running a
search for an already-cached day parses locally in seconds. Incomplete parts are
kept as `<asset>.part` and resumed on retry; a part is only used once its size
matches the release metadata.

**Why cache first:** streaming a 2.9 GB archive straight into the tar parser made
a single dropped connection fatal (`tar entry: …`) and restarted retries at zero.

Terrain (later, for geometry): Copernicus GLO-90 + EGM2008 geoid.

## Matching metrics (next phase)

For each candidate receiver location × candidate flight, score agreement between
the modeled Doppler curve and the observed ridge:

- **Frequency drop** — total downward sweep magnitude (→ speed, via `v = c·(fmax−fmin)/(fmax+fmin)`).
- **Shape of the shift** — deviation from the constant-velocity S-curve reveals
  **turns/acceleration**; the exact retarded-time curve is computed from the
  trajectory's radial velocity (emission-time, not reception-time).
- **Slant distance from each point** — the ridge slope at CPA encodes `d`
  (`d = f0·v²/(c·|df/dt|_CPA)`); matches per-candidate CPA range.
- **Time between two Dopplers** — multiple flyovers lock one shared recording
  clock; the flight-derived separation between events must equal the observed gap.
- **Harmonic consistency** — all harmonics of one flyover must map to the same
  flight, one rest frequency.

Scores rank candidates; they are not probabilities. Abstain when evidence is
insufficient (DOPPLERGEIST's confidence model).
