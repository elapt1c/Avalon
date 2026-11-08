//! Render a spectrogram region to PNG for visual inspection.
//!
//! Usage: spectro <file> <out.png> [--start S] [--dur S] [--max HZ] [--rate HZ]

use avalon::dsp::spectrogram::{compute_spectrogram, StftParams};
use avalon::model::AnalysisSettings;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 2 {
        eprintln!("usage: spectro <file> <out.png> [--start S] [--dur S] [--max HZ] [--rate HZ]");
        std::process::exit(2);
    }
    let path = &args[0];
    let out = &args[1];
    let mut start_s = 0.0f64;
    let mut dur_s = 0.0f64;
    let mut settings = AnalysisSettings::default();
    let mut i = 2;
    while i + 1 < args.len() {
        let v = args[i + 1].parse::<f32>().unwrap_or(f32::NAN);
        match args[i].as_str() {
            "--start" => start_s = v as f64,
            "--dur" => dur_s = v as f64,
            "--max" => settings.max_freq_hz = v,
            "--min" => settings.min_freq_hz = v,
            "--rate" => settings.target_sample_rate = v as u32,
            "--pad" => settings.fft_pad = v as u32,
            _ => {}
        }
        i += 2;
    }

    let (audio, backend) = avalon::decode::decode(std::path::Path::new(path))?;
    let mut samples = if audio.sample_rate == settings.target_sample_rate {
        audio.samples.clone()
    } else {
        avalon::dsp::resample_mono(
            &audio.samples,
            audio.sample_rate,
            settings.target_sample_rate,
        )
    };
    let sr = settings.target_sample_rate;
    if dur_s > 0.0 {
        let a = (start_s * sr as f64) as usize;
        let b = ((start_s + dur_s) * sr as f64) as usize;
        samples = samples[a.min(samples.len())..b.min(samples.len())].to_vec();
        start_s = 0.0;
    }

    let spec = compute_spectrogram(
        &samples,
        sr,
        &StftParams {
            window_s: settings.window_s,
            hop_s: settings.hop_s,
            min_freq_hz: settings.min_freq_hz,
            max_freq_hz: settings.max_freq_hz,
            fft_pad: settings.fft_pad,
        },
    );
    let (nf, nt) = (spec.freq_axis_hz.len(), spec.time_axis_s.len());
    eprintln!(
        "backend {backend} sr {sr} -> {nf} bins x {nt} frames, res {:.3} Hz",
        spec.resolution_hz
    );
    if nf == 0 || nt == 0 {
        anyhow::bail!("empty spectrogram");
    }

    // Robust ceiling.
    let mut vals: Vec<f32> = Vec::with_capacity(nf * nt);
    for row in &spec.contrast {
        vals.extend_from_slice(row);
    }
    let idx = (((vals.len() as f64) * 0.995) as usize).min(vals.len() - 1);
    let ceil = *vals
        .select_nth_unstable_by(idx, |a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .1;
    let ceil = ceil.max(3.0);

    let mut img = image::RgbImage::new(nt as u32, nf as u32);
    for (k, row) in spec.contrast.iter().enumerate() {
        let y = (nf - 1 - k) as u32;
        for (f, &v) in row.iter().enumerate() {
            let t = ((v - 1.0) / (ceil - 1.0)).clamp(0.0, 1.0);
            img.put_pixel(f as u32, y, image::Rgb(colormap(t)));
        }
    }
    img.save(out)?;
    println!("wrote {out} ({nt}x{nf})");
    Ok(())
}

fn colormap(t: f32) -> [u8; 3] {
    let stops: [(f32, [u8; 3]); 5] = [
        (0.0, [8, 8, 12]),
        (0.25, [60, 18, 92]),
        (0.5, [150, 30, 90]),
        (0.75, [230, 110, 40]),
        (1.0, [250, 230, 140]),
    ];
    for i in 0..stops.len() - 1 {
        let (t0, c0) = stops[i];
        let (t1, c1) = stops[i + 1];
        if t <= t1 {
            let f = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
            let mut o = [0u8; 3];
            for k in 0..3 {
                o[k] = (c0[k] as f32 + (c1[k] as f32 - c0[k] as f32) * f) as u8;
            }
            return o;
        }
    }
    [250, 230, 140]
}
