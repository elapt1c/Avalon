//! Short-time Fourier transform with time-whitened spectral contrast.
//!
//! Mirrors DOPPLERGEIST's acoustic front-end: a fixed-length Hann window with a
//! fixed hop, magnitude spectrum, then a per-frequency-bin contrast that removes
//! stationary tonal baselines and normalizes the noise floor. A moving tone that
//! passes through a bin therefore stands out as a contrast spike.

use rustfft::num_complex::Complex;
use rustfft::FftPlanner;

use crate::model::Spectrogram;

pub struct StftParams {
    pub window_s: f32,
    pub hop_s: f32,
    pub min_freq_hz: f32,
    pub max_freq_hz: f32,
    /// Zero-padding multiplier for the FFT. The true resolution is set by the
    /// window length (1/window_s); oversampling interpolates the spectrum so
    /// bin peaks (and therefore ridge frequency estimates) are finer-grained.
    pub fft_pad: u32,
}

pub fn compute_spectrogram(samples: &[f32], sample_rate: u32, params: &StftParams) -> Spectrogram {
    let win_len = (params.window_s * sample_rate as f32).round() as usize;
    if win_len < 4 || samples.len() < win_len {
        return Spectrogram {
            time_axis_s: Vec::new(),
            freq_axis_hz: Vec::new(),
            contrast: Vec::new(),
            resolution_hz: 0.0,
        };
    }
    let hop_len = (params.hop_s * sample_rate as f32).round().max(1.0) as usize;
    let fft_size = (win_len * params.fft_pad.max(1) as usize).next_power_of_two();
    let n_frames = (samples.len() - win_len) / hop_len + 1;

    // Hann window.
    let window: Vec<f32> = (0..win_len)
        .map(|i| {
            0.5 - 0.5 * (2.0 * std::f32::consts::PI * i as f32 / (win_len as f32 - 1.0)).cos()
        })
        .collect();

    let bin_hz = sample_rate as f32 / fft_size as f32;
    let min_bin = ((params.min_freq_hz / bin_hz).floor() as usize).min(fft_size / 2);
    let max_bin = ((params.max_freq_hz / bin_hz).ceil() as usize)
        .min(fft_size / 2 - 1)
        .max(min_bin);
    let n_freq = max_bin - min_bin + 1;

    let mut mag: Vec<Vec<f32>> = vec![vec![0.0f32; n_frames]; n_freq];

    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(fft_size);
    let mut buf = vec![Complex::new(0.0f32, 0.0f32); fft_size];

    for f in 0..n_frames {
        let start = f * hop_len;
        buf.fill(Complex::new(0.0, 0.0));
        for (i, w) in window.iter().enumerate() {
            buf[i] = Complex::new(samples[start + i] * w, 0.0);
        }
        fft.process(&mut buf);
        for k in min_bin..=max_bin {
            // Magnitude (scale cancels in the contrast ratio).
            mag[k - min_bin][f] = buf[k].norm();
        }
    }

    // Per-frame peak magnitude, used as a dynamic-range gate: a zero-padded FFT
    // of a clean tone shows Hann sidelobes across the whole band that would
    // otherwise read as dozens of separate tonals.
    let mut frame_max = vec![0.0f32; n_frames];
    for f in 0..n_frames {
        let mut m = 0.0f32;
        for k in 0..n_freq {
            m = m.max(mag[k][f]);
        }
        frame_max[f] = m;
    }

    // Per-frequency-bin median baseline, then contrast = mag / (median + eps),
    // gated to 26 dB below that frame's strongest component. Written back into
    // `mag` in place: a one-hour clip at full band is ~0.5 GB per matrix.
    for k in 0..n_freq {
        let baseline = crate::dsp::median_f32(&mag[k]).max(1e-6);
        for f in 0..n_frames {
            mag[k][f] = if mag[k][f] < 0.05 * frame_max[f] {
                0.0
            } else {
                mag[k][f] / baseline
            };
        }
    }
    let contrast = mag;

    let time_axis_s: Vec<f64> = (0..n_frames)
        .map(|f| (f * hop_len + win_len / 2) as f64 / sample_rate as f64)
        .collect();
    let freq_axis_hz: Vec<f32> = (min_bin..=max_bin)
        .map(|k| k as f32 * bin_hz)
        .collect();

    Spectrogram {
        time_axis_s,
        freq_axis_hz,
        contrast,
        resolution_hz: sample_rate as f32 / win_len as f32,
    }
}
