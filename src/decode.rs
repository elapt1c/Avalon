//! Media decoding.
//!
//! Primary backend is `symphonia` (pure Rust, no system dependencies) which
//! demuxes WAV/MP3/FLAC/OGG/AAC/ALAC/AIFF/CAF and MP4/MKV containers. When
//! symphonia cannot open a container/codec, we fall back to a bundled or
//! PATH-resolved `ffmpeg`/`ffprobe` for full codec coverage.

use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result};

use crate::model::DecodedAudio;

/// Decode a media file to mono `f32` PCM at its native sample rate.
/// Returns the decoded audio plus the name of the backend that produced it.
pub fn decode_audio(path: &Path) -> Result<(DecodedAudio, &'static str)> {
    match symphonia_decode(path) {
        Ok(a) => Ok((a, "symphonia")),
        Err(sym_err) => match ffmpeg_decode(path) {
            Ok(a) => Ok((a, "ffmpeg")),
            Err(ff_err) => Err(anyhow::anyhow!(
                "decode failed. symphonia: {sym_err} | ffmpeg: {ff_err}"
            )),
        },
    }
}

/// Decode straight to mono at `target_sr`, streaming from ffmpeg.
///
/// Symphonia accumulates the whole file at its native rate (an hour of 44.1 kHz
/// stereo is ~635 MB of f32 before any resampling). Streaming at the analysis
/// rate holds only the final mono buffer, which is ~8x smaller.
pub fn decode_mono_at(path: &Path, target_sr: u32) -> Result<(DecodedAudio, &'static str)> {
    match ffmpeg_stream_at(path, target_sr) {
        Ok(a) => Ok((a, "ffmpeg (streaming)")),
        Err(stream_err) => match decode_audio(path) {
            Ok(mut pair) => {
                if pair.0.sample_rate != target_sr {
                    pair.0.samples = crate::dsp::resample_mono(
                        &pair.0.samples,
                        pair.0.sample_rate,
                        target_sr,
                    );
                    pair.0.duration_s = pair.0.samples.len() as f64 / target_sr as f64;
                }
                crate::log::debug(
                    "decode",
                    &format!("streaming decode unavailable ({stream_err}); used fallback"),
                );
                Ok(pair)
            }
            Err(e) => Err(anyhow::anyhow!(
                "decode failed. streaming: {stream_err} | fallback: {e}"
            )),
        },
    }
}

fn ffmpeg_stream_at(path: &Path, target_sr: u32) -> Result<DecodedAudio, String> {
    use std::io::Read as _;
    use std::process::Stdio;
    let ffmpeg = locate_binary("ffmpeg").ok_or("ffmpeg binary not found")?;
    let mut child = Command::new(&ffmpeg)
        .args(["-hide_banner", "-loglevel", "error", "-i"])
        .arg(path)
        .args([
            "-vn",
            "-ac",
            "1",
            "-ar",
            &target_sr.to_string(),
            "-f",
            "f32le",
            "pipe:1",
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| e.to_string())?;
    let mut out = child.stdout.take().ok_or("no stdout")?;
    let mut samples: Vec<f32> = Vec::new();
    let mut buf = vec![0u8; 1 << 16];
    loop {
        let n = out.read(&mut buf).map_err(|e| e.to_string())?;
        if n == 0 {
            break;
        }
        for chunk in buf[..n].chunks_exact(4) {
            samples.push(f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]));
        }
    }
    let status = child.wait().map_err(|e| e.to_string())?;
    if !status.success() || samples.is_empty() {
        return Err("ffmpeg produced no audio".into());
    }
    Ok(DecodedAudio {
        duration_s: samples.len() as f64 / target_sr as f64,
        samples,
        sample_rate: target_sr,
    })
}

// ---------------------------------------------------------------------------
// symphonia backend
// ---------------------------------------------------------------------------

fn symphonia_decode(path: &Path) -> Result<DecodedAudio, String> {
    use symphonia::core::audio::sample::Sample;
    use symphonia::core::codecs::audio::AudioDecoderOptions;
    use symphonia::core::errors::Error as SymError;
    use symphonia::core::formats::probe::Hint;
    use symphonia::core::formats::{FormatOptions, TrackType};
    use symphonia::core::io::MediaSourceStream;
    use symphonia::core::meta::MetadataOptions;

    let file = std::fs::File::open(path).map_err(|e| e.to_string())?;
    let mss = MediaSourceStream::new(Box::new(file), Default::default());

    let mut hint = Hint::new();
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        hint.with_extension(ext);
    }

    let mut format = symphonia::default::get_probe()
        .probe(
            &hint,
            mss,
            FormatOptions::default(),
            MetadataOptions::default(),
        )
        .map_err(|e| format!("unsupported format: {e}"))?;

    let track = format
        .default_track(TrackType::Audio)
        .ok_or("no audio track")?;
    let params = track
        .codec_params
        .as_ref()
        .and_then(|p| p.audio())
        .ok_or("audio codec parameters unavailable")?;
    let sample_rate = params.sample_rate.ok_or("missing sample rate")?;

    let mut decoder = symphonia::default::get_codecs()
        .make_audio_decoder(params, &AudioDecoderOptions::default())
        .map_err(|e| format!("unsupported codec: {e}"))?;
    let track_id = track.id;

    let mut mono: Vec<f32> = Vec::new();
    loop {
        let packet = match format.next_packet() {
            Ok(Some(p)) => p,
            Ok(None) => break,
            Err(SymError::ResetRequired) => break,
            Err(e) => return Err(format!("read error: {e}")),
        };
        if packet.track_id != track_id {
            continue;
        }
        match decoder.decode(&packet) {
            Ok(audio_buf) => {
                let channels = audio_buf.spec().channels().count();
                let n = audio_buf.samples_interleaved();
                let mut interleaved = vec![f32::MID; n];
                audio_buf.copy_to_slice_interleaved(&mut interleaved);
                let frames = n / channels.max(1);
                for i in 0..frames {
                    let mut acc = 0.0f32;
                    for ch in 0..channels {
                        acc += interleaved[i * channels + ch];
                    }
                    mono.push(acc / channels.max(1) as f32);
                }
            }
            Err(SymError::IoError(_)) | Err(SymError::DecodeError(_)) => continue,
            Err(_) => break,
        }
    }

    if mono.is_empty() {
        return Err("decoded no audio samples".into());
    }
    let duration_s = mono.len() as f64 / sample_rate as f64;
    Ok(DecodedAudio {
        samples: mono,
        sample_rate,
        duration_s,
    })
}

// ---------------------------------------------------------------------------
// ffmpeg fallback backend
// ---------------------------------------------------------------------------

fn ffmpeg_decode(path: &Path) -> Result<DecodedAudio, String> {
    let ffmpeg = locate_binary("ffmpeg").ok_or("ffmpeg binary not found")?;
    let ffprobe = locate_binary("ffprobe").ok_or("ffprobe binary not found")?;

    // Probe native sample rate of the first audio stream.
    let probe = Command::new(&ffprobe)
        .args([
            "-v", "error", "-select_streams", "a:0", "-show_entries",
            "stream=sample_rate", "-of", "csv=p=0",
        ])
        .arg(path)
        .output()
        .map_err(|e| e.to_string())?;
    let rate_text = String::from_utf8_lossy(&probe.stdout);
    let sample_rate: u32 = rate_text
        .trim()
        .split(',')
        .next()
        .unwrap_or("")
        .trim()
        .parse()
        .map_err(|_| format!("could not parse sample rate from ffprobe: {rate_text}"))?;

    // Decode to mono f32le on stdout.
    let out = Command::new(&ffmpeg)
        .args(["-hide_banner", "-loglevel", "error", "-i"])
        .arg(path)
        .args(["-vn", "-ac", "1", "-f", "f32le", "pipe:1"])
        .output()
        .map_err(|e| e.to_string())?;
    if !out.status.success() {
        return Err(format!(
            "ffmpeg failed: {}",
            String::from_utf8_lossy(&out.stderr).trim()
        ));
    }

    let bytes = out.stdout;
    let n = bytes.len() / 4;
    let mut mono = vec![0.0f32; n];
    for (i, chunk) in bytes.chunks_exact(4).enumerate() {
        mono[i] = f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
    }
    if mono.is_empty() {
        return Err("ffmpeg produced no audio".into());
    }
    let duration_s = mono.len() as f64 / sample_rate as f64;
    Ok(DecodedAudio {
        samples: mono,
        sample_rate,
        duration_s,
    })
}

/// Locate a bundled binary: env override, next to the executable, project
/// `resources/bin`, then `PATH`.
fn locate_binary(name: &str) -> Option<PathBuf> {
    // 1. Explicit override (useful for testing).
    let env_key = match name {
        "ffmpeg" => "AVALON_FFMPEG",
        "ffprobe" => "AVALON_FFPROBE",
        _ => return None,
    };
    if let Ok(p) = std::env::var(env_key) {
        let p = PathBuf::from(p);
        if p.is_file() {
            return Some(p);
        }
    }

    // 2. Next to the executable (and `bin/` subdir).
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let cands = [dir.join(name), dir.join("bin").join(name)];
            for c in cands {
                if c.is_file() {
                    return Some(c);
                }
            }
        }
    }

    // 3. Project resources (dev builds).
    let res = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("resources/bin").join(name);
    if res.is_file() {
        return Some(res);
    }

    // 4. PATH.
    let path_var = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path_var) {
        let c = dir.join(name);
        if c.is_file() {
            return Some(c);
        }
    }
    None
}

/// Helper kept public for diagnostics.
pub fn decoder_backend_available() -> Vec<&'static str> {
    let mut v = vec!["symphonia (built-in)"];
    if locate_binary("ffmpeg").is_some() && locate_binary("ffprobe").is_some() {
        v.push("ffmpeg (fallback)");
    }
    v
}

/// Convenience wrapper so callers get a context-rich error.
pub fn decode(path: &Path) -> Result<(DecodedAudio, &'static str)> {
    decode_audio(path).with_context(|| format!("failed to decode {}", path.display()))
}
