#!/usr/bin/env bash
# Overnight accuracy matrix for Avalon.
#
# Runs the reference case under a grid of scoring settings, appending one line
# per run to /root/work/overnight/matrix.jsonl. Safe to run unattended: every
# invocation is a separate process, the parsed archive cache is warm, and a
# failure only costs its own line.
set -u
cd /root/projects/AmazingOSINT/AmazingOSINT || exit 1

OUT=/root/work/overnight
mkdir -p "$OUT"
JSONL="$OUT/matrix.jsonl"

export AVALON_CACHE=/root/projects/aosint-cache
BENCH=./target/release/examples/bench
BOX="tiny:-73.72,40.72,-73.42,41.02"
WIN="1788890400:1788910200"
TRUTH="40.876466,-73.524951"

record() {
  local label="$1"; shift
  local start end
  start=$(date +%s)
  local line
  line=$(env "$@" "$BENCH" --analysis /root/projects/analysis-default.json \
          --time "$WIN" --truth "$TRUTH" "$BOX" 2>&1 \
        | grep -E "^RESULT|nearest-to-truth" | tr '\n' ' ')
  end=$(date +%s)
  python3 - "$label" "$line" "$((end - start))" <<'PY' >> "$JSONL"
import json, re, sys
label, line, secs = sys.argv[1], sys.argv[2], int(sys.argv[3])
rec = {"label": label, "secs": secs}
m = re.search(r"best=(-?[\d.]+),(-?[\d.]+).*?rmse=([\d.]+).*?pair=(\S+).*?dist_km=([\d.]+).*?survivors=(\d+)", line)
if m:
    rec.update(lat=float(m.group(1)), lon=float(m.group(2)), rmse=float(m.group(3)),
               pair=m.group(4), dist_km=float(m.group(5)), survivors=int(m.group(6)))
m = re.search(r"nearest-to-truth cell (-?[\d.]+),(-?[\d.]+) score ([\d.]+)", line)
if m:
    rec["truth_score"] = float(m.group(3))
print(json.dumps(rec))
PY
  echo "[$(date +%H:%M)] $label -> $line"
}

# Baseline and the two knobs that changed the result tonight.
record "baseline"
record "pairs40"        AOS_CELL_PAIRS=40
record "tol40"          AOS_PAIR_TOL_S=40
record "pairs40_tol40"  AOS_CELL_PAIRS=40 AOS_PAIR_TOL_S=40
record "pairs40_shift20" AOS_CELL_PAIRS=40 AOS_SHIFT_HALF_S=20
record "pairs40_tol40_shift20" AOS_CELL_PAIRS=40 AOS_PAIR_TOL_S=40 AOS_SHIFT_HALF_S=20
record "pairs100_tol40_shift20" AOS_CELL_PAIRS=100 AOS_PAIR_TOL_S=40 AOS_SHIFT_HALF_S=20

# Wider grid: does the truth cell ever overtake the winner?
record "wide" AOS_CELL_PAIRS=60 AOS_PAIR_TOL_S=60 AOS_SHIFT_HALF_S=30

echo "done -> $JSONL"
