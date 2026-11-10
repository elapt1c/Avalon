//! Matched-filter scoring of candidate receiver positions.
//!
//! Ranking on ridge-fit RMSE measures how well a *traced line* fits a curve; it
//! is dominated by tracer error on the flat tails and cannot absorb an f0 bias.
//! Instead we score the modelled curve directly against a whitened spectrogram
//! ("brightness"): every frame and harmonic contributes the local brightness at
//! the frequency the model predicts, so a correct geometry lights up even when
//! no single traced ridge is clean.
//!
//! Reference points (digitalis-p/doppler-locate): whitening steps, the
//! `(c_e - w·d)/(c_e + v·d - w·d)` ratio, `Σ B /(nt·nh)`, f0 sweep ±8 % and one
//! shared clock offset.

use std::f64::consts::PI;

/// Mean sea-level sound speed, used when no atmosphere is supplied.
pub const C0: f64 = 343.0;

/// Per-candidate atmosphere: path-mean sound speed and wind (ECEF, m/s).
#[derive(Debug, Clone, Copy)]
pub struct Atmos {
    pub c: f64,
    pub wind: [f64; 3],
}

impl Default for Atmos {
    fn default() -> Self {
        Self { c: C0, wind: [0.0; 3] }
    }
}

/// A whitened brightness image over one event window.
pub struct Brightness {
    pub f_min: f32,
    pub df: f32,
    pub f_bins: usize,
    /// Absolute time of frame 0 and the frame spacing.
    pub t0: f64,
    pub dt: f64,
    pub frames: usize,
    /// data[f * frames + t]
    data: Vec<f32>,
}

impl Brightness {
    /// Build from mono audio, whitening exactly as the reference does:
    /// subtract the per-bin time median, per-frame ±15 Hz median, a 20 s
    /// per-bin median of what remains, then normalise by 1.4826·MAD and clip.
    pub fn build(
        audio: &[f32],
        sr: u32,
        audio_t0_abs: f64,
        band: (f32, f32),
        win_s: f32,
        hop_s: f32,
    ) -> Self {
        let win = (win_s * sr as f32).round().max(32.0) as usize;
        if audio.len() < win {
            return Self::empty(band);
        }
        let fft_size = win.next_power_of_two();
        let hop = (hop_s * sr as f32).round().max(1.0) as usize;
        let frames = (audio.len() - win) / hop + 1;
        let bin_hz = sr as f32 / fft_size as f32;
        let f_min = band.0.max(5.0);
        let f_max = band.1.max(f_min + bin_hz);
        let b0 = (f_min / bin_hz).floor() as usize;
        let b1 = ((f_max / bin_hz).ceil() as usize).min(fft_size / 2 - 1);
        let f_bins = b1.saturating_sub(b0) + 1;

        let window: Vec<f32> = (0..win)
            .map(|i| 0.5 - 0.5 * (2.0 * PI as f32 * i as f32 / (win as f32 - 1.0)).cos())
            .collect();
        let mut planner = rustfft::FftPlanner::<f32>::new();
        let fft = planner.plan_fft_forward(fft_size);
        let mut buf = vec![rustfft::num_complex::Complex::new(0.0f32, 0.0); fft_size];

        // Magnitude in dB.
        let mut img = vec![0.0f32; f_bins * frames];
        for t in 0..frames {
            buf.iter_mut().for_each(|c| *c = Default::default());
            let start = t * hop;
            for (i, w) in window.iter().enumerate() {
                buf[i] = rustfft::num_complex::Complex::new(audio[start + i] * w, 0.0);
            }
            fft.process(&mut buf);
            for f in 0..f_bins {
                let m = buf[b0 + f].norm();
                img[f * frames + t] = 20.0 * (m + 1e-9).log10();
            }
        }
        if frames == 0 || f_bins == 0 {
            return Self::empty(band);
        }

        // 1. subtract per-bin median over time.
        for f in 0..f_bins {
            let row = &img[f * frames..(f + 1) * frames];
            let med = median(row);
            for v in &mut img[f * frames..(f + 1) * frames] {
                *v -= med;
            }
        }
        // 2. subtract per-frame median over ±15 Hz.
        let half = (15.0 / bin_hz).round().max(1.0) as usize;
        let mut col = vec![0.0f32; 2 * half + 1];
        for t in 0..frames {
            let lo = 0;
            let hi = f_bins;
            let _ = (lo, hi);
            for f in 0..f_bins {
                let a = f.saturating_sub(half);
                let b = (f + half).min(f_bins - 1);
                let n = b - a + 1;
                for (i, k) in (a..=b).enumerate() {
                    col[i] = img[k * frames + t];
                }
                let med = median(&col[..n]);
                img[f * frames + t] -= med;
            }
        }
        // 3. subtract a 20 s per-bin median (steady hums / lines).
        let win_frames = ((20.0 / hop_s.max(1e-3)) as usize).max(1);
        for f in 0..f_bins {
            let row: Vec<f32> = img[f * frames..(f + 1) * frames].to_vec();
            let mut scratch: Vec<f32> = Vec::with_capacity(win_frames);
            for t in 0..frames {
                let a = t.saturating_sub(win_frames / 2);
                let b = (t + win_frames / 2 + 1).min(frames);
                scratch.clear();
                scratch.extend_from_slice(&row[a..b]);
                let med = median(&scratch).max(0.0);
                img[f * frames + t] -= med;
            }
        }
        // 4. normalise by robust sigma and clip.
        let mad = {
            let mut v = img.clone();
            let med = median(&v);
            for x in v.iter_mut() {
                *x = (*x - med).abs();
            }
            median(&v)
        };
        let sigma = (1.4826 * mad).max(1e-6);
        for v in img.iter_mut() {
            *v = (*v / sigma).clamp(0.0, 6.0);
        }
        // 5. ±1 bin max filter tolerates sub-bin model error.
        let mut out = img.clone();
        for f in 0..f_bins {
            for t in 0..frames {
                let a = f.saturating_sub(1);
                let b = (f + 1).min(f_bins - 1);
                let mut m = img[f * frames + t];
                for k in a..=b {
                    m = m.max(img[k * frames + t]);
                }
                out[f * frames + t] = m;
            }
        }

        Self {
            f_min,
            df: bin_hz,
            f_bins,
            t0: audio_t0_abs,
            dt: hop_s as f64,
            frames,
            data: out,
        }
    }

    fn empty(band: (f32, f32)) -> Self {
        Self {
            f_min: band.0,
            df: 1.0,
            f_bins: 0,
            t0: 0.0,
            dt: 1.0,
            frames: 0,
            data: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.frames == 0 || self.f_bins == 0
    }

    /// Bilinear brightness at (frequency Hz, frame index).
    #[inline]
    fn at(&self, f_hz: f32, frame: usize) -> f32 {
        let x = (f_hz - self.f_min) / self.df;
        if x < 0.0 || x > (self.f_bins - 1) as f32 {
            return 0.0;
        }
        let x0 = x.floor() as usize;
        let x1 = (x0 + 1).min(self.f_bins - 1);
        let fx = x - x0 as f32;
        let a = self.data[x0 * self.frames + frame];
        let b = self.data[x1 * self.frames + frame];
        a + (b - a) * fx
    }

    pub fn frame_time(&self, frame: usize) -> f64 {
        self.t0 + frame as f64 * self.dt
    }
}

fn median(v: &[f32]) -> f32 {
    if v.is_empty() {
        return 0.0;
    }
    let mut c = v.to_vec();
    let mid = c.len() / 2;
    let (_, m, _) = c.select_nth_unstable_by(mid, |a, b| {
        a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
    });
    *m
}

/// Cubic Hermite (Catmull-Rom) trajectory in ECEF.
///
/// Positions are stored as `f32` offsets from `origin` rather than absolute
/// ECEF triples. A continental search builds tens of thousands of tracks and
/// the track store was the single dominant allocation (24 B/point); offsets
/// halve that at ~0.25 m resolution, far below anything the acoustic model
/// cares about.
#[derive(Clone)]
pub struct Track {
    pub t0: f64,
    pub dt: f64,
    /// Absolute ECEF position the offsets are measured from.
    pub origin: [f64; 3],
    pub pos: Vec<[f32; 3]>,
    /// lat_min, lon_min, lat_max, lon_max of the retained passage.
    pub bbox: [f64; 4],
}

impl Track {
    /// Absolute ECEF position of sample `i`.
    #[inline]
    pub fn point(&self, i: usize) -> [f64; 3] {
        let p = self.pos[i];
        [
            self.origin[0] + p[0] as f64,
            self.origin[1] + p[1] as f64,
            self.origin[2] + p[2] as f64,
        ]
    }

    /// Observer position expressed in this track's offset space.
    #[inline]
    fn obs_offset(&self, obs: [f64; 3]) -> [f64; 3] {
        [obs[0] - self.origin[0], obs[1] - self.origin[1], obs[2] - self.origin[2]]
    }

    pub fn from_points(points: &[(f64, [f64; 3])]) -> Option<Self> {
        if points.len() < 2 {
            return None;
        }
        let t0 = points[0].0;
        let dt = points[1].0 - points[0].0;
        if dt <= 0.0 {
            return None;
        }
        // Uniform spacing keeps interpolation O(1); resample if needed.
        let n = ((points[points.len() - 1].0 - t0) / dt).round() as usize + 1;
        if n < 2 || n > 2_000_000 {
            return None;
        }
        let origin = points[0].1;
        let mut pos = Vec::with_capacity(n);
        let mut j = 0usize;
        for i in 0..n {
            let t = t0 + i as f64 * dt;
            while j + 1 < points.len() - 1 && points[j + 1].0 < t {
                j += 1;
            }
            let (ta, pa) = points[j];
            let (tb, pb) = points[(j + 1).min(points.len() - 1)];
            let f = if (tb - ta).abs() > 1e-9 {
                ((t - ta) / (tb - ta)).clamp(0.0, 1.0)
            } else {
                0.0
            };
            pos.push([
                ((pa[0] + (pb[0] - pa[0]) * f) - origin[0]) as f32,
                ((pa[1] + (pb[1] - pa[1]) * f) - origin[1]) as f32,
                ((pa[2] + (pb[2] - pa[2]) * f) - origin[2]) as f32,
            ]);
        }
        Some(Self { t0, dt, origin, pos, bbox: [0.0; 4] })
    }

    #[inline]
    pub fn span(&self) -> (f64, f64) {
        (self.t0, self.t0 + (self.pos.len() - 1) as f64 * self.dt)
    }

    /// Catmull-Rom position at time `t` (uniform spacing).
    pub fn at(&self, t: f64) -> Option<[f64; 3]> {
        let x = (t - self.t0) / self.dt;
        if x < 0.0 || x > (self.pos.len() - 1) as f64 {
            return None;
        }
        let i = x.floor() as usize;
        if i + 1 >= self.pos.len() {
            return Some(self.point(self.pos.len() - 1));
        }
        let f = (x - i as f64) as f64;
        let p0 = self.point(i.saturating_sub(1));
        let p1 = self.point(i);
        let p2 = self.point(i + 1);
        let p3 = self.point((i + 2).min(self.pos.len() - 1));
        let mut out = [0.0; 3];
        for k in 0..3 {
            let m1 = 0.5 * (p2[k] - p0[k]);
            let m2 = 0.5 * (p3[k] - p1[k]);
            let f2 = f * f;
            let f3 = f2 * f;
            out[k] = (2.0 * f3 - 3.0 * f2 + 1.0) * p1[k]
                + (f3 - 2.0 * f2 + f) * m1
                + (-2.0 * f3 + 3.0 * f2) * p2[k]
                + (f3 - f2) * m2;
        }
        Some(out)
    }

    /// Velocity by central difference (m/s).
    pub fn vel(&self, t: f64) -> Option<[f64; 3]> {
        let h = self.dt.min(1.0);
        let a = self.at(t - h)?;
        let b = self.at(t + h)?;
        let span = 2.0 * h;
        Some([(b[0] - a[0]) / span, (b[1] - a[1]) / span, (b[2] - a[2]) / span])
    }

    /// Time of **zero radial velocity** — the instant the aircraft crosses the
    /// receiver's perpendicular, which is what an acoustic recording actually
    /// measures (the Doppler shift passes through f0 there).
    ///
    /// This is NOT the closest-approach time for a turning aircraft: the
    /// geometric CPA can sit seconds away from the radial zero, and the search
    /// was pairing flyovers with the geometric version.
    pub fn radial_zero(&self, obs: [f64; 3]) -> Option<(f64, f64)> {
        self.radial_zero_near(obs, None)
    }

    /// As `radial_zero`, but prefers the crossing closest to `hint_t`. A track
    /// that passes the receiver twice has several zero crossings; the first one
    /// is often the wrong pass.
    pub fn radial_zero_near(&self, obs: [f64; 3], hint_t: Option<f64>) -> Option<(f64, f64)> {
        let n = self.pos.len();
        if n < 3 {
            return None;
        }
        let o = self.obs_offset(obs);
        let mut crossings: Vec<(f64, f64)> = Vec::new(); // (time, distance)
        let mut prev: Option<(f64, f64)> = None; // (vr, time)
        for i in 1..n - 1 {
            let p = self.pos[i];
            let a = self.pos[i - 1];
            let b = self.pos[i + 1];
            // Relative position and velocity in the track's offset space.
            let rx = p[0] as f64 - o[0];
            let ry = p[1] as f64 - o[1];
            let rz = p[2] as f64 - o[2];
            let range = (rx * rx + ry * ry + rz * rz).sqrt();
            if range < 1e-6 {
                continue;
            }
            let vx = (b[0] - a[0]) as f64;
            let vy = (b[1] - a[1]) as f64;
            let vz = (b[2] - a[2]) as f64;
            let vr = (vx * rx + vy * ry + vz * rz) / range;
            let t = self.t0 + i as f64 * self.dt;
            if let Some((vr0, t0)) = prev {
                if (vr0 < 0.0) != (vr < 0.0) {
                    // Linear interpolation of the sign change.
                    let w = if (vr - vr0).abs() > 1e-12 {
                        (-vr0 / (vr - vr0)).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    let tc = t0 + (t - t0) * w;
                    let idx = ((tc - self.t0) / self.dt).clamp(0.0, (n - 1) as f64);
                    let i0 = idx.floor() as usize;
                    let i1 = (i0 + 1).min(n - 1);
                    let f = idx - i0 as f64;
                    let mut d2 = 0.0;
                    for k in 0..3 {
                        let v = self.pos[i0][k] as f64
                            + (self.pos[i1][k] - self.pos[i0][k]) as f64 * f
                            - o[k];
                        d2 += v * v;
                    }
                    crossings.push((tc, d2.max(0.0).sqrt()));
                }
            }
            prev = Some((vr, t));
        }
        if crossings.is_empty() {
            return None;
        }
        let pick = match hint_t {
            Some(h) => crossings
                .into_iter()
                .min_by(|a, b| {
                    (a.0 - h).abs().partial_cmp(&(b.0 - h).abs()).unwrap_or(std::cmp::Ordering::Equal)
                })?,
            None => crossings[0],
        };
        Some((pick.1, pick.0))
    }

    /// Coarse closest approach: scan every `stride`-th point, then refine in a

    /// local window. Scoring one cell against thousands of flights must not
    /// touch every sample, or the grid loop dominates runtime.
    pub fn cpa_stride(&self, obs: [f64; 3], stride: usize) -> Option<(f64, f64)> {
        let n = self.pos.len();
        if n < 2 {
            return None;
        }
        // Work in offset space: converting the observer once is cheaper than
        // converting every track sample to absolute ECEF.
        let o = self.obs_offset(obs);
        let step = stride.max(1);
        let mut best = f64::INFINITY;
        let mut bi = 0usize;
        let mut i = 0usize;
        while i < n {
            let d = sq_dist_f32(self.pos[i], o);
            if d < best {
                best = d;
                bi = i;
            }
            i += step;
        }
        // Refine around the coarse minimum.
        let lo = bi.saturating_sub(step);
        let hi = (bi + step).min(n - 1);
        for (i, p) in self.pos.iter().enumerate().take(hi + 1).skip(lo) {
            let d = sq_dist_f32(*p, o);
            if d < best {
                best = d;
                bi = i;
            }
        }
        let d2 = |i: usize| sq_dist_f32(self.pos[i.min(n - 1)], o);
        let a = d2(bi.saturating_sub(1));
        let c = d2(bi + 1);
        let denom = a - 2.0 * best + c;
        let mut off = 0.0f64;
        let mut d_sq = best;
        if denom.abs() > 1e-9 {
            off = (0.5 * (a - c) / denom).clamp(-1.0, 1.0);
            d_sq = best - 0.25 * (a - c) * off;
        }
        Some((d_sq.max(0.0).sqrt(), self.t0 + (bi as f64 + off) * self.dt))
    }

    /// Closest approach to `obs`: (distance m, receive time). Coarse scan then
    /// parabolic refinement — cheap enough for the pruning stage.
    pub fn cpa(&self, obs: [f64; 3]) -> Option<(f64, f64)> {
        let o = self.obs_offset(obs);
        let mut best = f64::INFINITY;
        let mut bi = 0usize;
        for (i, p) in self.pos.iter().enumerate() {
            let d = sq_dist_f32(*p, o);
            if d < best {
                best = d;
                bi = i;
            }
        }
        let d2 = |i: usize| sq_dist_f32(self.pos[i.min(self.pos.len() - 1)], o);
        let a = d2(bi.saturating_sub(1));
        let c = d2(bi + 1);
        let denom = a - 2.0 * best + c;
        let mut off = 0.0f64;
        let mut d_sq = best;
        if denom.abs() > 1e-9 {
            off = (0.5 * (a - c) / denom).clamp(-1.0, 1.0);
            d_sq = best - 0.25 * (a - c) * off;
        }
        Some((d_sq.max(0.0).sqrt(), self.t0 + (bi as f64 + off) * self.dt))
    }
}

/// Squared distance between a track offset (f32) and an observer offset (f64,
/// both relative to the track origin).
#[inline]
fn sq_dist_f32(a: [f32; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] as f64 - b[0];
    let dy = a[1] as f64 - b[1];
    let dz = a[2] as f64 - b[2];
    dx * dx + dy * dy + dz * dz
}

/// One flyover as measured from the audio.
#[derive(Debug, Clone)]
pub struct EventSpec {
    /// Rest frequency from the fit (swept during scoring).
    pub f0: f32,
    /// Measured harmonic ratios (e.g. [1.0, 2.0]).
    pub harmonics: Vec<f32>,
    /// Recording-relative closest-approach time.
    pub cpa_s: f64,
    /// Half the scoring window (s).
    pub half_width_s: f64,
}

impl EventSpec {
    pub fn band(&self) -> (f32, f32) {
        let lo = 0.6 * self.f0 * self.harmonics.first().copied().unwrap_or(1.0);
        let hi = 1.6 * self.f0 * self.harmonics.last().copied().unwrap_or(1.0);
        (lo.max(15.0), hi.max(lo * 1.2 + 1.0))
    }
}

/// Model ratio `f/f0` for a receiver at `obs` and a receive time `t_recv`.
/// Solves the retarded emission time by fixed point and includes wind aloft.
#[inline]
pub fn doppler_ratio(track: &Track, obs: [f64; 3], t_recv: f64, atm: Atmos) -> Option<f64> {
    let mut tau = t_recv - 2000.0 / atm.c;
    for _ in 0..12 {
        let p = track.at(tau)?;
        let d = [p[0] - obs[0], p[1] - obs[1], p[2] - obs[2]];
        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        if r <= 0.0 {
            return None;
        }
        let dh = [d[0] / r, d[1] / r, d[2] / r];
        let w_d = atm.wind[0] * dh[0] + atm.wind[1] * dh[1] + atm.wind[2] * dh[2];
        let next = t_recv - r / (atm.c - w_d);
        if (next - tau).abs() < 1e-6 {
            tau = next;
            break;
        }
        tau = next;
    }
    let p = track.at(tau)?;
    let v = track.vel(tau)?;
    let d = [p[0] - obs[0], p[1] - obs[1], p[2] - obs[2]];
    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-6);
    let dh = [d[0] / r, d[1] / r, d[2] / r];
    let w_d = atm.wind[0] * dh[0] + atm.wind[1] * dh[1] + atm.wind[2] * dh[2];
    let v_d = v[0] * dh[0] + v[1] * dh[1] + v[2] * dh[2];
    Some((atm.c - w_d) / (atm.c + v_d - w_d))
}

/// Retarded-time ratio with a warm start.
///
/// The cold-start fixed point contracts at v/c per iteration, so for distant
/// geometry it never reaches 1e-6 and always burns the full cap. Adjacent
/// frames are 0.25 s apart, so carrying `tau` between frames converges in 1-2
/// steps; 1e-4 s of timing error is ~3 cm of path.
#[inline]
pub fn doppler_ratio_warm(
    track: &Track,
    obs: [f64; 3],
    t_recv: f64,
    tau: &mut f64,
    atm: Atmos,
) -> Option<f64> {
    for _ in 0..4 {
        let p = track.at(*tau)?;
        let d = [p[0] - obs[0], p[1] - obs[1], p[2] - obs[2]];
        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        if r <= 0.0 {
            return None;
        }
        let dh = [d[0] / r, d[1] / r, d[2] / r];
        let w_d = atm.wind[0] * dh[0] + atm.wind[1] * dh[1] + atm.wind[2] * dh[2];
        let next = t_recv - r / (atm.c - w_d);
        if (next - *tau).abs() < 1e-4 {
            *tau = next;
            break;
        }
        *tau = next;
    }
    let p = track.at(*tau)?;
    let v = track.vel(*tau)?;
    let d = [p[0] - obs[0], p[1] - obs[1], p[2] - obs[2]];
    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-6);
    let dh = [d[0] / r, d[1] / r, d[2] / r];
    let w_d = atm.wind[0] * dh[0] + atm.wind[1] * dh[1] + atm.wind[2] * dh[2];
    let v_d = v[0] * dh[0] + v[1] * dh[1] + v[2] * dh[2];
    Some((atm.c - w_d) / (atm.c + v_d - w_d))
}

/// Mean whitened brightness along the modelled curve for one event.
///
/// The denominator is `frames * harmonics`, so samples that fall out of band or
/// off the track count as zero — a curve that wanders out of the analysed band
/// is penalised rather than scoring only on its visible portion.
pub fn score_event(
    b: &Brightness,
    track: &Track,
    obs: [f64; 3],
    ev: &EventSpec,
    f0: f64,
    base_abs: f64,
    shift_s: f64,
    atm: Atmos,
    frame_step: usize,
) -> f64 {
    if b.is_empty() {
        return 0.0;
    }
    let nh = ev.harmonics.len().max(1);
    let mut acc = 0.0f32;
    let mut counted = 0usize;
    let step = frame_step.max(1);
    // Warm start carried across frames.
    let mut tau = base_abs + b.frame_time(0) - shift_s - 2000.0 / atm.c;
    for frame in (0..b.frames).step_by(step) {
        // `b` uses recording-relative time; `base_abs` is the (fitted) absolute
        // time of the recording start for this flight.
        let model_t = base_abs + b.frame_time(frame) - shift_s;
        let Some(ratio) = doppler_ratio_warm(track, obs, model_t, &mut tau, atm) else {
            tau = model_t - 2000.0 / atm.c;
            counted += nh;
            continue;
        };
        for h in &ev.harmonics {
            let f = (f0 * *h as f64 * ratio) as f32;
            acc += b.at(f, frame);
            counted += 1;
        }
    }
    if counted == 0 {
        0.0
    } else {
        acc as f64 / counted as f64
    }
}

/// RMSE (Hz) between the measured ridge curve and the modelled one.
///
/// This is the quantity the method actually ranks on: the simulated Doppler
/// curve for this candidate/flight must follow the measured ridge. Brightness
/// sums are nearly flat between competing flights in busy airspace; a residual
/// in Hz is not.
pub fn curve_rmse(
    samples: &[(f64, f32)],
    track: &Track,
    obs: [f64; 3],
    f0: f64,
    base_abs: f64,
    shift_s: f64,
    atm: Atmos,
    max_err_hz: f32,
) -> Option<f64> {
    curve_rmse_offset(samples, track, obs, f0, base_abs, shift_s, 0.0, atm, max_err_hz)
}

/// Same residual with a constant frequency offset.
///
/// A traced ridge can sit a bin or two off the true tone for its whole length
/// (visible as a vertical bias against the spectrogram), and a pure offset is
/// not absorbable by scale: at 180 Hz, 3 Hz is 1.7 %.
#[allow(clippy::too_many_arguments)]
/// Runtime switch for the scoring backend. Defaults to the DOPPLERGEIST-style
/// rules (pseudo-Huber loss, silence as a reject against the fitted slant).
pub fn dg_mode() -> bool {
    static CACHED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    // Production default: the DOPPLERGEIST-style rules. Set AOS_METHOD=legacy to
    // fall back to the old squared-error / multiplicative-silence backend.
    *CACHED.get_or_init(|| {
        std::env::var("AOS_METHOD")
            .map(|v| v != "legacy")
            .unwrap_or(true)
    })
}

pub fn curve_rmse_offset(
    samples: &[(f64, f32)],
    track: &Track,
    obs: [f64; 3],
    f0: f64,
    base_abs: f64,
    shift_s: f64,
    f_offset: f64,
    atm: Atmos,
    max_err_hz: f32,
) -> Option<f64> {
    if samples.len() < 6 {
        return None;
    }
    let mut sum = 0.0f64;
    let mut n = 0usize;
    let mut unresolved = 0usize;
    let mut tau = base_abs + samples[0].0 - shift_s - 2000.0 / atm.c;
    for &(t_rel, f_meas) in samples {
        let model_t = base_abs + t_rel - shift_s;
        let Some(ratio) = doppler_ratio_warm(track, obs, model_t, &mut tau, atm) else {
            // A sample the model cannot even place is a *failure*, not something
            // to skip: skipping let a geometry that matched a handful of lucky
            // samples report a sub-Hz residual.
            tau = model_t - 2000.0 / atm.c;
            unresolved += 1;
            sum += (max_err_hz as f64).powi(2);
            n += 1;
            continue;
        };
        let f_model = f0 * ratio + f_offset;
        let e = f_model - f_meas as f64;
        // DG method: pseudo-Huber instead of squared error, so a handful of bad
        // ridge samples cannot dominate the fit. Scale is the declared loss
        // scale (1 Hz default in their configuration).
        if dg_mode() {
            let delta = 1.0f64;
            let r = e.abs();
            sum += 2.0 * delta * delta * ((1.0 + (r / delta).powi(2)).sqrt() - 1.0);
            n += 1;
            continue;
        }
        if e.abs() > max_err_hz as f64 {
            // Grossly wrong sample: penalise but keep going so a single bad
            // ridge point cannot make the whole curve look perfect.
            sum += (max_err_hz as f64).powi(2);
        } else {
            sum += e * e;
        }
        n += 1;
    }
    if n < 6 {
        return None;
    }
    // Require the model to place most of the curve, not a lucky subset.
    if unresolved as f64 > 0.2 * n as f64 {
        return None;
    }
    // Range degeneracy guard: a distant pass has a nearly flat curve, which can
    // be matched at almost any location once f0 is free. The measured *total
    // frequency excursion* pins the geometry, so require the model's excursion
    // to agree with it.
    let mut f_min = f32::INFINITY;
    let mut f_max = f32::NEG_INFINITY;
    for &(_, f) in samples {
        f_min = f_min.min(f);
        f_max = f_max.max(f);
    }
    let measured_drop = (f_max - f_min) as f64;
    if measured_drop > 2.0 {
        let mut m_min = f64::INFINITY;
        let mut m_max = f64::NEG_INFINITY;
        let mut tau2 = base_abs + samples[0].0 - shift_s - 2000.0 / atm.c;
        for &(t_rel, _) in samples {
            let model_t = base_abs + t_rel - shift_s;
            if let Some(r) = doppler_ratio_warm(track, obs, model_t, &mut tau2, atm) {
                let f = f0 * r + f_offset;
                m_min = m_min.min(f);
                m_max = m_max.max(f);
            }
        }
        if m_max.is_finite() {
            let model_drop = m_max - m_min;
            if (model_drop - measured_drop).abs() > 0.25 * measured_drop {
                return None;
            }
        }
    }
    Some((sum / n as f64).sqrt())
}

/// f0 sweep grid (±8 %, as in the reference: straight-line fits are biased by
/// several percent when the aircraft descends or slows).
pub fn f0_grid(f0: f64, steps: usize) -> Vec<f64> {
    f0_grid_rel(f0, steps, f0_half_width())
}

/// Half-width of the rest-frequency sweep, as a fraction of the fitted `f0`.
///
/// The reference implementation sweeps ±8 %, which is enough rope for a wrong
/// flight to bend its curve onto a measured one: over a continent there are
/// thousands of candidate pairs, and the freedom to rescale the tone by 8 %
/// buys several Hz of curve-fitting error. The default is tightened to a value
/// consistent with the audio fit's own accuracy; `AOS_F0_HALF` overrides it for
/// experiments.
pub fn f0_half_width() -> f64 {
    static CACHED: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *CACHED.get_or_init(|| {
        std::env::var("AOS_F0_HALF")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .filter(|v| *v > 0.0 && *v <= 0.25)
            .unwrap_or(0.08)
    })
}

pub fn f0_grid_rel(f0: f64, steps: usize, half: f64) -> Vec<f64> {
    let n = steps.max(3);
    (0..n)
        .map(|i| {
            let t = i as f64 / (n - 1) as f64 * 2.0 - 1.0;
            f0 * (1.0 + half * t)
        })
        .collect()
}

/// Clock-offset sweep grid (seconds).
pub fn shift_grid(half_s: f64, steps: usize) -> Vec<f64> {
    let n = steps.max(3);
    (0..n)
        .map(|i| (i as f64 / (n - 1) as f64 * 2.0 - 1.0) * half_s)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn straight_track() -> Track {
        // Northbound at 100 m/s at 1 km altitude, over 200 s.
        let mut pts = Vec::new();
        for i in 0..=200 {
            let t = i as f64;
            let lat = 39.0 + 100.0 * (t - 100.0) / 111_320.0;
            pts.push((t, crate::search::ecef(lat, -98.0, 1000.0)));
        }
        Track::from_points(&pts).unwrap()
    }

    #[test]
    fn ratio_is_unity_at_cpa_and_below_after() {
        let tr = straight_track();
        let obs = crate::search::ecef(39.0, -98.0, 0.0);
        let (d, t_cpa) = tr.cpa(obs).unwrap();
        let before = doppler_ratio(&tr, obs, t_cpa - 20.0, Atmos::default()).unwrap();
        let after = doppler_ratio(&tr, obs, t_cpa + 20.0, Atmos::default()).unwrap();
        assert!(before > 1.0, "approaching should be blue-shifted: {before}");
        assert!(after < 1.0, "receding should be red-shifted: {after}");
        // The tone is unshifted when the *emission* happens at CPA, i.e. one
        // propagation delay later in receive time — not at t_cpa itself.
        let at_arrival = t_cpa + d / C0;
        let at = doppler_ratio(&tr, obs, at_arrival, Atmos::default()).unwrap();
        assert!((at - 1.0).abs() < 0.02, "ratio at arrival should be ~1: {at}");
    }

    #[test]
    fn wind_shifts_the_ratio() {
        let tr = straight_track();
        let obs = crate::search::ecef(39.0, -98.0, 0.0);
        let (_, t_cpa) = tr.cpa(obs).unwrap();
        let still = doppler_ratio(&tr, obs, t_cpa + 20.0, Atmos::default()).unwrap();
        // Wind blowing along the line of sight changes the effective speed.
        let atm = Atmos { c: 343.0, wind: [0.0, 15.0, 0.0] };
        let windy = doppler_ratio(&tr, obs, t_cpa + 20.0, atm).unwrap();
        assert!((still - windy).abs() > 1e-4, "wind must change the ratio");
    }

    #[test]
    fn cpa_distance_matches_geometry() {
        let tr = straight_track();
        let obs = crate::search::ecef(39.0, -98.0, 0.0);
        let (d, _) = tr.cpa(obs).unwrap();
        // 1 km altitude directly overhead.
        assert!((d - 1000.0).abs() < 30.0, "got {d}");
    }
}
