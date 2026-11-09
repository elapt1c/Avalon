"""How much ENF extraction error can the time match still survive?

Takes real 1-second grid frequency (GB, NESO), cuts a 10-minute slice, degrades
it the way an acoustic ENF extraction would (coarse quantisation + noise), and
cross-correlates the slice against the whole month. Reports whether the peak is
at the true position and how far it beats the best wrong position.
"""
import numpy as np
from scipy.signal import correlate

CSV = "/root/work/fnew-2026-9.csv"
f = np.loadtxt(CSV, delimiter=",", skiprows=1, usecols=1, dtype=np.float64)
print(f"series: {len(f)} s ({len(f)/86400:.2f} days), mean {f.mean():.4f} Hz, "
      f"std {f.std()*1000:.1f} mHz, range {f.min():.3f}..{f.max():.3f}")

fz = f - f.mean()
# per-offset energy via cumulative sum (O(N), exact)
cs = np.concatenate(([0.0], np.cumsum(fz * fz)))

def best_offsets(seg, top=3):
    s = seg - seg.mean()
    c = correlate(fz, s, mode="valid", method="fft")
    L = len(s)
    energy = cs[L:] - cs[:-L]          # sum of fz[i:i+L]^2 for every i
    score = c / np.sqrt(np.maximum(energy, 1e-12) * np.sum(s * s) + 1e-12)
    order = np.argsort(score)[::-1]
    picked = []
    for i in order:
        if all(abs(int(i) - p) > 600 for p in picked):
            picked.append(int(i))
        if len(picked) >= top:
            break
    return picked, score

def clock(idx):
    day, sod = divmod(idx, 86400)
    return f"day {day+1:02d} {sod//3600:02d}:{(sod%3600)//60:02d}:{sod%60:02d}"

truth = 17 * 86400 + 11 * 3600 + 40 * 60
L = 600
seg = f[truth:truth + L]
print(f"\nslice: 600 s starting {clock(truth)} (sample {truth})")

print(f"\n{'quant mm':>9} {'noise mHz':>10} {'peak at':>17} {'err s':>7} {'corr':>6} {'2nd':>6} {'margin':>7}  verdict")
rng0 = np.random.default_rng(11)
for q_mm in (1, 5, 10, 20):
    for sigma_mhz in (0, 2, 5, 10, 20, 50):
        rng = np.random.default_rng(11)
        y = np.round(seg / (q_mm / 1000.0)) * (q_mm / 1000.0) + rng.normal(0, sigma_mhz / 1000.0, L)
        picked, score = best_offsets(y)
        best = picked[0]
        err = best - truth
        margin = (score[best] - score[picked[1]]) * 100
        verdict = "LOCK" if err == 0 and margin > 5 else ("weak" if err == 0 else "LOST")
        print(f"{q_mm:>9} {sigma_mhz:>10} {clock(best):>17} {err:>7} {score[best]:>6.3f} "
              f"{score[picked[1]]:>6.3f} {margin:>6.2f}%  {verdict}")

# how does a 1-hour slice do at the same degradations? (longer window = more evidence)
print("\n--- same test with a 3600 s slice ---")
L = 3600
seg = f[truth:truth + L]
for q_mm in (5, 10):
    for sigma_mhz in (10, 20, 50):
        rng = np.random.default_rng(11)
        y = np.round(seg / (q_mm / 1000.0)) * (q_mm / 1000.0) + rng.normal(0, sigma_mhz / 1000.0, L)
        picked, score = best_offsets(y)
        best = picked[0]
        err = best - truth
        margin = (score[best] - score[picked[1]]) * 100
        verdict = "LOCK" if err == 0 and margin > 5 else ("weak" if err == 0 else "LOST")
        print(f"{q_mm:>9} {sigma_mhz:>10} {clock(best):>17} {err:>7} {score[best]:>6.3f} "
              f"{score[picked[1]]:>6.3f} {margin:>6.2f}%  {verdict}")
