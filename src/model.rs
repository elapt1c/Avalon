//! Shared data types for analysis, search, and results.

use serde::{Deserialize, Serialize};

/// Sound speed used by the Doppler model (m/s). Constant air value; wind/refraction
/// are unmodeled sources of uncertainty (as in DOPPLERGEIST).
pub const SOUND_SPEED: f64 = 343.0;

/// A rectangular geographic region (degrees).
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct Bounds {
    pub west: f64,
    pub south: f64,
    pub east: f64,
    pub north: f64,
}

impl Bounds {
    pub fn valid(&self) -> bool {
        self.west < self.east && self.south < self.north
    }
    /// Centre of the box as `(lat, lon)` — the order every caller expects. It
    /// previously returned `(lon, lat)`, which put the heatmap in Antarctica.
    pub fn center(&self) -> (f64, f64) {
        ((self.south + self.north) / 2.0, (self.west + self.east) / 2.0)
    }
    pub fn contains(&self, lat: f64, lon: f64) -> bool {
        lat >= self.south && lat <= self.north && lon >= self.west && lon <= self.east
    }
    pub fn intersects(&self, o: &Bounds) -> bool {
        self.west <= o.east && o.west <= self.east && self.south <= o.north && o.south <= self.north
    }
}

/// Optional intelligence the operator can supply to narrow the search.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Metadata {
    /// Selected geographic region(s). Empty = whole world (very slow).
    pub region: Vec<Bounds>,
    /// Believed UTC window of the first audio sample.
    pub time_window: Option<UtcWindow>,
    /// Frequency band hint (Hz) constraining spectral analysis. Overrides defaults.
    pub freq_band_hz: Option<(f32, f32)>,
    /// Free-form notes.
    pub notes: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct UtcWindow {
    /// Unix seconds.
    pub start_unix_s: i64,
    pub end_unix_s: i64,
}

/// A search request: region + mandatory time window + selected flyover events.
#[derive(Debug, Clone)]
pub struct SearchRequest {
    pub region: Vec<Bounds>,
    pub time_window: UtcWindow,
    pub event_ids: Vec<usize>,
    /// Source media, re-decoded to score against the spectrogram.
    pub source_path: String,
}

/// A single flight trajectory point.
#[derive(Debug, Clone, Copy)]
pub struct FlightPoint {
    /// Unix seconds.
    pub t_s: f64,
    pub lat: f64,
    pub lon: f64,
    pub alt_m: f64,
}

/// A flight trajectory (one aircraft).
#[derive(Debug, Clone)]
pub struct Flight {
    pub icao: String,
    pub points: Vec<FlightPoint>,
}

impl Flight {
    /// Slice of the trajectory covering `[t0, t1]` (padded by `pad_s`).
    pub fn window(&self, t0: f64, t1: f64, pad_s: f64) -> Vec<FlightPoint> {
        self.points
            .iter()
            .filter(|p| p.t_s >= t0 - pad_s && p.t_s <= t1 + pad_s)
            .copied()
            .collect()
    }
    pub fn covers(&self, t0: f64, t1: f64) -> bool {
        if self.points.len() < 2 {
            return false;
        }
        let a = self.points.first().unwrap().t_s;
        let b = self.points.last().unwrap().t_s;
        t1 >= a && t0 <= b
    }
}

/// A heatmap cell: candidate receiver location + score (lower is better).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HeatmapCell {
    pub lat: f64,
    pub lon: f64,
    /// Composite ranking score (higher is better).
    pub score_hz: f32,
    /// Mean curve residual of the winning hypothesis, Hz (INFINITY = none).
    pub rmse_hz: f32,
    /// Flights that explain this cell ("icao+icao").
    pub pair: String,
}

/// One spatially coherent candidate solution.
///
/// Over a continent the top cells come in near-ties that are kilometres or
/// states apart: two flight pairs can fit the same two measured curves to
/// within a few hundredths of a Hz. Reporting the single best cell is therefore
/// over-confident by construction — the honest output is the set of modes that
/// are statistically indistinguishable, each with its own evidence.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mode {
    /// Score-weighted centre of the cluster.
    pub lat: f64,
    pub lon: f64,
    /// Best cell score in the cluster.
    pub score: f32,
    /// Residual of that best cell (Hz).
    pub rmse_hz: f32,
    /// Flights explaining the best cell ("icao+icao").
    pub pair: String,
    /// Cells in the cluster.
    pub cells: usize,
    /// Radius of the cluster (km).
    pub radius_km: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchResult {
    pub cells: Vec<HeatmapCell>,
    /// Top candidate locations, best first.
    pub best: Vec<HeatmapCell>,
    /// A human-readable caveat about this search, e.g. that the region was too
    /// large for the score kernel to be sampled. `None` when the run is sound.
    #[serde(default)]
    pub warning: Option<String>,
    /// Distinct near-tied solutions, best first. With two events and a free
    /// rest frequency this is normally longer than one; treat entries within a
    /// small score margin of the leader as equally supported.
    #[serde(default)]
    pub modes: Vec<Mode>,
    pub n_flights: usize,
    pub n_candidates: usize,
    /// Inclusive bounding box of the searched region (degrees).
    pub bounds: Option<Bounds>,
    /// Candidate grid spacing (degrees), for rendering the heatmap.
    pub step_deg: f64,
    /// Winning cell, for overlaying the simulated curve on the spectrogram.
    pub best_cell: Option<(f64, f64)>,
    /// Simulated Doppler curves at the winning cell, in recording-relative
    /// seconds: one entry per scored event.
    #[serde(default)]
    pub model_curves: Vec<Vec<(f64, f32)>>,
}

/// Decoded mono PCM at its native sample rate.
#[derive(Debug, Clone)]
pub struct DecodedAudio {
    pub samples: Vec<f32>,
    pub sample_rate: u32,
    pub duration_s: f64,
}

/// One (time, frequency, weight) measurement of a detected ridge.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct RidgeSample {
    pub t_s: f64,
    pub freq_hz: f32,
    /// Spectral contrast at this point (dimensionless).
    pub weight: f32,
}

/// A detected moving-tone ridge: one harmonic of one flyover, traced over time.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ridge {
    pub id: usize,
    pub samples: Vec<RidgeSample>,
    pub start_s: f64,
    pub end_s: f64,
    pub mean_contrast: f32,
    pub min_freq_hz: f32,
    pub max_freq_hz: f32,
    /// max - min frequency (positive for the downward Doppler sweep of a flyover).
    pub net_drop_hz: f32,
    /// Linear-regression slope of frequency vs time (Hz/s).
    pub slope_hz_s: f32,
    /// Fraction of samples moving downward (Doppler-plausibility metric, 0..=1).
    pub downward_fraction: f32,
    /// Best Doppler fit for this ridge (filled in by the pipeline).
    pub fit: Option<DopplerFit>,
    /// Whether the fit passes the physical plausibility gate.
    pub plausible: bool,
}

/// Constant-velocity Doppler parameterization of a ridge.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DopplerFit {
    /// Rest (emitted) frequency of this tone, Hz.
    pub f0_hz: f32,
    /// Estimated source speed, m/s.
    pub speed_mps: f32,
    /// Estimated closest-approach slant distance, m.
    pub slant_distance_m: f32,
    /// Time of closest approach, recording-relative seconds.
    pub t_cpa_s: f64,
    /// Root-mean-square residual of the fit, Hz.
    pub rmse_hz: f32,
    /// Fraction of ridge samples consistent with the fitted curve (0..=1).
    pub inlier_fraction: f32,
    pub n_samples: usize,
    /// True when both approach and recession extremes were observed.
    pub complete_sweep: bool,
}

/// A detected acoustic event: one flyover, possibly with several harmonics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetectedEvent {
    pub id: usize,
    /// Ridge ids belonging to this event.
    pub ridge_ids: Vec<usize>,
    /// Representative (lowest-frequency) ridge id.
    pub representative_ridge_id: usize,
    /// Harmonic ratio of each ridge relative to the representative.
    pub harmonic_ratios: Vec<f32>,
    /// Doppler fit on the representative ridge, if a fit was possible.
    pub doppler: Option<DopplerFit>,
    pub start_s: f64,
    pub end_s: f64,
}

/// Spectrogram (spectral contrast) with axes for plotting.
#[derive(Debug, Clone)]
pub struct Spectrogram {
    /// Center time of each frame, seconds.
    pub time_axis_s: Vec<f64>,
    /// Center frequency of each bin, Hz.
    pub freq_axis_hz: Vec<f32>,
    /// contrast[freq_bin][time_frame], dimensionless (>= ~0).
    pub contrast: Vec<Vec<f32>>,
    /// True frequency resolution (Hz): `sample_rate / window_length`. Bin spacing
    /// can be finer than this when the FFT is zero-padded.
    pub resolution_hz: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisSettings {
    pub min_freq_hz: f32,
    pub max_freq_hz: f32,
    /// Analysis window length (s). Sets true frequency resolution: 1/window_s.
    pub window_s: f32,
    /// Frame hop (s). Smaller = denser time sampling of the ridges.
    pub hop_s: f32,
    /// FFT zero-padding multiplier (finer spectral bin spacing).
    pub fft_pad: u32,
    /// Seed contrast threshold for a ridge candidate.
    pub seed_contrast: f32,
    pub target_sample_rate: u32,
}

impl Default for AnalysisSettings {
    fn default() -> Self {
        Self {
            min_freq_hz: 40.0,
            max_freq_hz: 1600.0,
            // 3 s window -> 0.33 Hz true resolution; 0.25 s hop -> 4 frames/s.
            window_s: 3.0,
            hop_s: 0.25,
            // FFT is 4x the window -> ~0.08 Hz bin spacing at 8 kHz.
            fft_pad: 4,
            seed_contrast: 2.5,
            target_sample_rate: 8000,
        }
    }
}

/// Full analysis result (in-memory). The spectrogram is excluded from JSON export.
#[derive(Debug, Clone)]
pub struct AnalysisResult {
    pub source_path: String,
    pub sample_rate: u32,
    pub duration_s: f64,
    /// Decoder backend that produced the audio ("symphonia" or "ffmpeg").
    pub backend: &'static str,
    pub spectrogram: Spectrogram,
    pub ridges: Vec<Ridge>,
    pub events: Vec<DetectedEvent>,
    pub settings: AnalysisSettings,
}
