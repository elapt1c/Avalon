//! Digital signal processing: spectrogram, ridge detection, harmonic grouping,
//! and Doppler parameterization.

pub mod doppler;
pub mod harmonics;
pub mod ridge;
pub mod spectrogram;

/// Resample mono `f32` audio to `out_sr` Hz.
///
/// Downsampling applies a windowed-sinc anti-alias lowpass before decimation so
/// out-of-band content does not alias into the analysis band, then linearly
/// resamples the (smooth, band-limited) result. Upsampling uses linear
/// interpolation (images fall above the band of interest).
pub fn resample_mono(samples: &[f32], in_sr: u32, out_sr: u32) -> Vec<f32> {
    if samples.is_empty() || out_sr == in_sr {
        return samples.to_vec();
    }
    let ratio = out_sr as f64 / in_sr as f64;
    let out_len = ((samples.len() as f64) * ratio).round() as usize;
    if out_len == 0 {
        return Vec::new();
    }

    let mut out = Vec::with_capacity(out_len);
    if ratio < 1.0 {
        // Anti-alias then decimate. The convolution is evaluated only at the
        // output positions: filtering every input sample costs `1/ratio` times
        // more work for no gain (a 1 h clip at 44.1 kHz -> 4 kHz is 15 G MACs
        // the naive way versus 1.4 G here).
        let half = 48usize;
        let h = sinc_lowpass(0.45 * ratio, half);
        for i in 0..out_len {
            let center = ((i as f64) / ratio).round() as isize;
            let mut acc = 0.0f64;
            for (k, hk) in h.iter().enumerate() {
                let idx = center + k as isize - half as isize;
                if idx >= 0 && (idx as usize) < samples.len() {
                    acc += hk * samples[idx as usize] as f64;
                }
            }
            out.push(acc as f32);
        }
    } else {
        for i in 0..out_len {
            let pos = (i as f64) / ratio;
            let i0 = pos.floor() as usize;
            let i1 = (i0 + 1).min(samples.len() - 1);
            let frac = (pos - i0 as f64) as f32;
            out.push(samples[i0] * (1.0 - frac) + samples[i1] * frac);
        }
    }
    out
}

/// Normalized Hann-windowed sinc lowpass. `cutoff` is cycles/sample (0..=0.5).
fn sinc_lowpass(cutoff: f64, half: usize) -> Vec<f64> {
    use std::f64::consts::PI;

    let n = half;
    let mut h = vec![0.0f64; 2 * n + 1];
    for (i, hk) in h.iter_mut().enumerate() {
        let k = i as isize - n as isize;
        let kk = k as f64;
        *hk = if k == 0 {
            2.0 * cutoff
        } else {
            (2.0 * cutoff * PI * kk).sin() / (PI * kk)
        };
        *hk *= 0.5 + 0.5 * (PI * kk / n as f64).cos();
    }
    let sum: f64 = h.iter().sum();
    for hi in h.iter_mut() {
        *hi /= sum;
    }
    h
}

/// Median of a slice (returns the lower-middle element for even lengths).
/// Uses selection rather than a full sort: called once per frequency bin, over
/// every frame, so O(n log n) sorting dominates long recordings.
pub fn median_f32(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut v: Vec<f32> = values.to_vec();
    let n = v.len();
    let mid = n / 2;
    let lt = |a: &f32, b: &f32| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal);
    let (_, mid_val, _) = v.select_nth_unstable_by(mid, lt);
    let hi = *mid_val;
    if n % 2 == 1 {
        hi
    } else {
        let lo = v[..mid]
            .iter()
            .copied()
            .fold(f32::NEG_INFINITY, |a, b| if b > a { b } else { a });
        (lo + hi) / 2.0
    }
}
