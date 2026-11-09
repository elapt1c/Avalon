"""Is there a mains-hum (ENF) trace in this recording?

For every step, take a long window, band-pass the 40-200 Hz region, and find the
narrow spectral peak near 60 Hz (and its harmonics). A genuine ENF trace is a
*stable narrow line that wanders by tens of millihertz* — not a fixed 60.000 Hz
and not a broad smudge. Reports the estimate series so we can see which it is.
"""
import sys
import numpy as np
from scipy.io import wavfile
from scipy.signal import get_window

path = sys.argv[1] if len(sys.argv) > 1 else "/root/work/enf-test.wav"
sr, x = wavfile.read(path)
x = x.astype(np.float64)
if x.ndim > 1:
    x = x.mean(axis=1)
x /= (np.abs(x).max() + 1e-12)
print(f"file {path}: {len(x)/sr:.1f} s @ {sr} Hz")

def peak_near(seg, sr, f0, span=1.5, zero_pad=8):
    """Estimate the frequency of the narrow peak closest to f0."""
    n = len(seg)
    w = get_window("hann", n, fftbins=True)
    seg = (seg - seg.mean()) * w
    nfft = n * zero_pad
    spec = np.abs(np.fft.rfft(seg, nfft))
    freqs = np.fft.rfftfreq(nfft, 1.0 / sr)
    lo, hi = f0 - span, f0 + span
    m = (freqs >= lo) & (freqs <= hi)
    if not m.any():
        return None, 0.0
    idx = np.argmax(spec[m])
    f_peak = freqs[m][idx]
    band = spec[m]
    # parabolic interpolation on the log magnitude for sub-bin accuracy
    if 0 < idx < len(band) - 1:
        a, b, c = np.log(band[idx - 1] + 1e-20), np.log(band[idx] + 1e-20), np.log(band[idx + 1] + 1e-20)
        denom = (a - 2 * b + c)
        if abs(denom) > 1e-12:
            delta = 0.5 * (a - c) / denom
            f_peak += delta * (freqs[1] - freqs[0])
    snr = band[idx] / (np.median(band) + 1e-20)
    return f_peak, snr

win_s = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
step_s = float(sys.argv[3]) if len(sys.argv) > 3 else 10.0
N = int(win_s * sr)
S = int(step_s * sr)

for target, name in ((60.0, "60 Hz"), (120.0, "120 Hz"), (180.0, "180 Hz")):
    est = []
    snrs = []
    for start in range(0, max(1, len(x) - N), S):
        f, s = peak_near(x[start:start + N], sr, target)
        if f is not None and np.isfinite(f):
            est.append(f)
            snrs.append(s)
    if not est:
        print(f"{name}: no data")
        continue
    e = np.array(est)
    s = np.array(snrs)
    good = s > 6.0
    print(f"\n{name}: windows={len(e)}  median SNR={np.median(s):.1f}  windows SNR>6: {good.sum()}")
    print(f"  freq mean {e.mean():.4f} Hz  std {e.std()*1000:.1f} mHz  "
          f"range {e.min():.4f}..{e.max():.4f} ({(e.max()-e.min())*1000:.1f} mHz)")
    if good.sum() >= 5:
        eg = e[good]
        print(f"  high-SNR subset: mean {eg.mean():.4f} std {eg.std()*1000:.1f} mHz")
    print("  first 12 estimates:", " ".join(f"{v:.3f}" for v in e[:12]))
