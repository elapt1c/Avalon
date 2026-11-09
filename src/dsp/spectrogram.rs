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
    let t_all = std::time::Instant::now();
    if std::env::var("AVALON_TIMING").is_ok() {
        crate::log::info(
            "spec",
            &format!("frames {n_frames} x bins {n_freq} (fft {fft_size}, win {win_len}, hop {hop_len})"),
        );
    }

    // Frames are independent, so the FFT pass runs on the whole pool. Output is
    // collected frame-major (cache-friendly for the workers) and scattered into
    // the freq-major matrix afterwards; one memcpy pass costs a second, the FFT
    // pass used to cost the entire analysis.
    {
        use rayon::prelude::*;
        const CHUNK_FRAMES: usize = 256;
        let mut chunk_buf = vec![0.0f32; CHUNK_FRAMES.min(n_frames.max(1)) * n_freq];
        for c0 in (0..n_frames).step_by(CHUNK_FRAMES) {
            let c1 = (c0 + CHUNK_FRAMES).min(n_frames);
            let cf = c1 - c0;
            let buf = &mut chunk_buf[..cf * n_freq];
            buf.par_chunks_mut(n_freq).enumerate().for_each_init(
                || {
                    let mut planner = FftPlanner::<f32>::new();
                    (
                        planner.plan_fft_forward(fft_size),
                        vec![Complex::new(0.0f32, 0.0f32); fft_size],
                    )
                },
                |(fft, scratch), (i, out)| {
                    let f = c0 + i;
                    let start = f * hop_len;
                    scratch.fill(Complex::new(0.0, 0.0));
                    for (j, w) in window.iter().enumerate() {
                        scratch[j] = Complex::new(samples[start + j] * w, 0.0);
                    }
                    fft.process(scratch);
                    for k in min_bin..=max_bin {
                        out[k - min_bin] = scratch[k].norm();
                    }
                },
            );
            for i in 0..cf {
                let f = c0 + i;
                let row = &buf[i * n_freq..(i + 1) * n_freq];
                for (k, v) in row.iter().enumerate() {
                    mag[k][f] = *v;
                }
            }
        }
    }

    let t_fft = t_all.elapsed();
    let t_gate = std::time::Instant::now();
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
    {
        use rayon::prelude::*;
        mag.par_iter_mut().for_each(|row| {
            let baseline = crate::dsp::median_f32(row).max(1e-6);
            for (f, v) in row.iter_mut().enumerate() {
                *v = if *v < 0.05 * frame_max[f] { 0.0 } else { *v / baseline };
            }
        });
    }
    if std::env::var("AVALON_TIMING").is_ok() {
        crate::log::info(
            "spec",
            &format!(
                "fft+scatter {:.1}s, gate+contrast {:.1}s, total {:.1}s",
                t_fft.as_secs_f64(),
                t_gate.elapsed().as_secs_f64(),
                t_all.elapsed().as_secs_f64()
            ),
        );
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
