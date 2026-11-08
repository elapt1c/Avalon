//! Main application UI.

use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::thread;

use egui::{Color32, RichText};
use egui_plot::{Line, Plot, PlotImage, PlotPoint, PlotPoints};

use crate::dateutil;
use crate::dsp::doppler;
use crate::map::{MapState, Mode};
use crate::model::{
    AnalysisResult, AnalysisSettings, Bounds, DetectedEvent, DopplerFit, Metadata, SearchRequest,
    SearchResult, Spectrogram, UtcWindow,
};
use crate::pipeline::{self, Progress};
use crate::search;
use crate::theme;

const PROVIDERS: [&str; 3] = ["staging-0", "prod-0", "mlatonly"];

/// Which year the test window refers to (override with AVALON_YEAR).
fn default_year() -> i64 {
    std::env::var("AVALON_YEAR")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2026)
}

/// Ground truth of the test recording (lat, lon) — used to score the result.
const TRUTH: (f64, f64) = (40.87, -73.52);

/// Great-circle distance in metres.
fn haversine_m(a: (f64, f64), b: (f64, f64)) -> f64 {
    let r = 6_371_000.0;
    let (la1, lo1) = (a.0.to_radians(), a.1.to_radians());
    let (la2, lo2) = (b.0.to_radians(), b.1.to_radians());
    let dlat = la2 - la1;
    let dlon = lo2 - lo1;
    let h = (dlat / 2.0).sin().powi(2) + la1.cos() * la2.cos() * (dlon / 2.0).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

enum AnalysisState {
    Idle,
    Running { fraction: f32, stage: String },
    Done,
    Error(String),
}

enum SearchState {
    Idle,
    Running { stage: String },
    Done,
    Error(String),
}

enum SearchMsg {
    Progress(Box<search::Snapshot>),
    Log(String),
    Done(SearchResult),
    Error(String),
}

#[derive(Clone, Copy, PartialEq)]
enum CalTarget {
    Start,
    End,
}

#[derive(Clone, Copy, PartialEq)]
enum MapPurpose {
    Region,
    Heatmap,
}

pub struct App {
    selected_path: Option<PathBuf>,

    // Intelligence.
    region: Vec<Bounds>,
    map: MapState,
    map_open: bool,
    map_purpose: MapPurpose,
    time_start_unix: Option<i64>,
    time_end_unix: Option<i64>,
    cal_open: bool,
    cal_target: CalTarget,
    cal_year: i64,
    cal_month: u32,
    cal_day: u32,
    cal_hour: u32,
    cal_min: u32,
    cal_sec: u32,
    freq_enabled: bool,
    freq_min: String,
    freq_max: String,
    notes: String,
    provider: String,

    settings: AnalysisSettings,

    // Analysis.
    analysis: AnalysisState,
    result: Option<AnalysisResult>,
    rx: Option<Receiver<Progress>>,
    selected_events: HashSet<usize>,
    spectro_texture: Option<(String, egui::TextureHandle)>,
    log: Vec<String>,
    last_export: Option<String>,

    // Search.
    search: SearchState,
    search_result: Option<SearchResult>,
    search_rx: Option<Receiver<SearchMsg>>,
    search_snapshot: Option<search::Snapshot>,
    search_log: Vec<String>,
    score_history: Vec<f32>,
    search_started: Option<std::time::Instant>,
    /// Keep the locate dashboard in the main area until the user dismisses it.
    show_processing: bool,
    /// Test hook: run the locate automatically once analysis finishes.
    auto_run: bool,
    /// Clip to open and analyse on the first frame (test hook).
    autostart: Option<PathBuf>,
    cancel: std::sync::Arc<std::sync::atomic::AtomicBool>,
}

impl App {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        theme::apply(&cc.egui_ctx);
        let auto_run = std::env::var("AVALON_AUTORUN").is_ok();
        // Test hook: open a clip and start analysis on the first frame, so a run
        // needs no clicking at all.
        let autostart = if auto_run {
            std::env::var("AVALON_CLIP")
                .ok()
                .map(PathBuf::from)
                .or_else(|| {
                    let p = PathBuf::from(
                        "/home/tpg/Downloads/Unknown location - raw audio - colsto (128p).mp4",
                    );
                    p.is_file().then_some(p)
                })
        } else {
            None
        };
        Self {
            auto_run,
            autostart,
            selected_path: None,
            // Default test region: ~50 miles around the known true location of
            // the test footage, so accuracy is measurable.
            region: vec![Bounds {
                west: TRUTH.1 - 0.9555,
                south: TRUTH.0 - 0.7226,
                east: TRUTH.1 + 0.9555,
                north: TRUTH.0 + 0.7226,
            }],
            map: MapState::new(39.0, -98.0, 4.0),
            map_open: false,
            map_purpose: MapPurpose::Region,
            // Daytime in New York. The test footage's two passes resolve to
            // 2026-10-30 15:34Z (both implied starts agree inside 7 s), so the
            // default spans the 30th and the 31st rather than the 31st alone.
            time_start_unix: Some(dateutil::civil_to_unix(default_year(), 10, 30, 11, 0, 0)),
            time_end_unix: Some(dateutil::civil_to_unix(default_year(), 10, 31, 23, 0, 0)),
            cal_open: false,
            cal_target: CalTarget::Start,
            cal_year: 2026,
            cal_month: 1,
            cal_day: 1,
            cal_hour: 0,
            cal_min: 0,
            cal_sec: 0,
            freq_enabled: false,
            freq_min: String::new(),
            freq_max: String::new(),
            notes: String::new(),
            provider: PROVIDERS[0].to_string(),
            settings: AnalysisSettings::default(),
            analysis: AnalysisState::Idle,
            result: None,
            rx: None,
            selected_events: HashSet::new(),
            spectro_texture: None,
            log: Vec::new(),
            last_export: None,
            search: SearchState::Idle,
            search_result: None,
            search_rx: None,
            search_snapshot: None,
            search_log: Vec::new(),
            score_history: Vec::new(),
            search_started: None,
            show_processing: false,
            cancel: std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false)),
        }
    }

    fn now_unix() -> i64 {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs() as i64)
            .unwrap_or(0)
    }

    fn build_metadata(&self) -> Metadata {
        let time_window = match (self.time_start_unix, self.time_end_unix) {
            (Some(a), Some(b)) => Some(UtcWindow { start_unix_s: a.min(b), end_unix_s: a.max(b) }),
            _ => None,
        };
        let freq_band_hz = if self.freq_enabled {
            match (self.freq_min.trim().parse::<f32>(), self.freq_max.trim().parse::<f32>()) {
                (Ok(a), Ok(b)) if a > 0.0 && b > a => Some((a, b)),
                _ => None,
            }
        } else {
            None
        };
        Metadata {
            region: self.region.clone(),
            time_window,
            freq_band_hz,
            notes: self.notes.clone(),
        }
    }

    // -- analysis -----------------------------------------------------------

    fn start_analysis(&mut self) {
        let Some(path) = self.selected_path.clone() else {
            self.analysis = AnalysisState::Error("Select a source file first".into());
            return;
        };
        let settings = self.settings.clone();
        let meta = self.build_metadata();
        let (tx, rx) = mpsc::channel();
        self.rx = Some(rx);
        self.analysis = AnalysisState::Running { fraction: 0.0, stage: "Starting…".into() };
        self.result = None;
        self.selected_events.clear();
        self.log.clear();
        self.last_export = None;
        self.spectro_texture = None;

        thread::spawn(move || match pipeline::run(&path, &settings, &meta, Some(tx.clone())) {
            Ok(result) => {
                let _ = tx.send(Progress::Done(result));
            }
            Err(e) => {
                let _ = tx.send(Progress::Error(e.to_string()));
            }
        });
    }

    fn poll_analysis(&mut self, ctx: &egui::Context) {
        let Some(rx) = &self.rx else { return };
        // The worker runs on another thread; without an explicit repaint the UI
        // sits idle until some input arrives, so progress only appeared to move
        // when the pointer touched the window.
        ctx.request_repaint_after(std::time::Duration::from_millis(40));
        let mut autorun_ids: Option<Vec<usize>> = None;
        loop {
            match rx.try_recv() {
                Ok(Progress::Stage { name, fraction }) => {
                    self.analysis = AnalysisState::Running { fraction, stage: name.to_string() };
                }
                Ok(Progress::Log(msg)) => self.log.push(msg),
                Ok(Progress::Done(result)) => {
                    let plausible: Vec<usize> = result
                        .events
                        .iter()
                        .filter(|e| {
                            result
                                .ridges
                                .iter()
                                .find(|r| r.id == e.representative_ridge_id)
                                .map(|r| r.plausible)
                                .unwrap_or(false)
                        })
                        .map(|e| e.id)
                        // Event 0 is left out of the autorun selection by request.
                        .filter(|id| *id != 0)
                        .collect();
                    self.result = Some(result);
                    self.analysis = AnalysisState::Done;
                    if self.auto_run && !plausible.is_empty() {
                        // Deferred: the receiver borrow above forbids &mut self.
                        autorun_ids = Some(plausible);
                    }
                }
                Ok(Progress::Error(e)) => self.analysis = AnalysisState::Error(e),
                Err(TryRecvError::Empty) => break,
                Err(TryRecvError::Disconnected) => {
                    self.rx = None;
                    break;
                }
            }
        }
        if let Some(ids) = autorun_ids {
            self.selected_events = ids.into_iter().collect();
            crate::log::info(
                "app",
                &format!(
                    "autorun: locating with {} event(s)",
                    self.selected_events.len()
                ),
            );
            self.start_search();
        }
    }

    // -- search -------------------------------------------------------------

    fn search_prereqs(&self) -> Result<(SearchRequest, Vec<DetectedEvent>, Vec<crate::model::Ridge>), String> {
        let (Some(start), Some(end)) = (self.time_start_unix, self.time_end_unix) else {
            return Err("Set a time window (mandatory)".into());
        };
        if start >= end {
            return Err("Time window start must be before end".into());
        }
        if self.region.is_empty() {
            return Err("Select a region (whole-world search is impractical)".into());
        }
        let Some(result) = &self.result else {
            return Err("Run analysis first to detect Doppler events".into());
        };
        if self.selected_events.is_empty() {
            return Err("Toggle at least one Doppler event (curve) to locate".into());
        }
        let mut ids: Vec<usize> = self.selected_events.iter().copied().collect();
        ids.sort_unstable();
        let mut events = Vec::new();
        let mut ridges = Vec::new();
        for ev_id in &ids {
            let Some(event) = result.events.get(*ev_id) else { continue };
            if event.doppler.is_none() {
                continue;
            }
            let Some(ridge) = result
                .ridges
                .iter()
                .find(|r| r.id == event.representative_ridge_id)
            else {
                continue;
            };
            events.push(event.clone());
            ridges.push(ridge.clone());
        }
        if events.is_empty() {
            return Err("Toggled events have no usable Doppler fit".into());
        }
        let req = SearchRequest {
            region: self.region.clone(),
            time_window: UtcWindow { start_unix_s: start, end_unix_s: end },
            event_ids: ids,
            source_path: result.source_path.clone(),
        };
        Ok((req, events, ridges))
    }

    fn start_search(&mut self) {
        let (req, events, ridges) = match self.search_prereqs() {
            Ok(x) => x,
            Err(e) => {
                self.search = SearchState::Error(e);
                return;
            }
        };
        let provider = self.provider.clone();
        let (tx, rx) = mpsc::channel();
        self.search_rx = Some(rx);
        self.search = SearchState::Running { stage: "Starting…".into() };
        self.search_result = None;
        self.search_snapshot = None;
        self.search_log.clear();
        self.score_history.clear();
        self.search_started = Some(std::time::Instant::now());
        self.show_processing = true;
        let cancel = self.cancel.clone();
        cancel.store(false, std::sync::atomic::Ordering::Relaxed);

        thread::spawn(move || {
            let mut progress = |s: search::Snapshot| {
                let _ = tx.send(SearchMsg::Progress(Box::new(s)));
            };
            let mut log = |msg: String| {
                let _ = tx.send(SearchMsg::Log(msg));
            };
            match search::run_search(
                &req,
                &events,
                &ridges,
                &provider,
                &mut progress,
                &mut log,
                &cancel,
            ) {
                Ok(result) => {
                    let _ = tx.send(SearchMsg::Done(result));
                }
                Err(e) => {
                    let _ = tx.send(SearchMsg::Error(e.to_string()));
                }
            }
        });
    }

    fn poll_search(&mut self, ctx: &egui::Context) {
        let Some(rx) = &self.search_rx else { return };
        ctx.request_repaint_after(std::time::Duration::from_millis(60));
        let mut got_progress = false;
        let mut open_map = false;
        loop {
            match rx.try_recv() {
                Ok(SearchMsg::Progress(snap)) => {
                    got_progress = true;
                    if matches!(self.search, SearchState::Running { .. }) {
                        self.search = SearchState::Running { stage: snap.stage.to_string() };
                    }
                    // Track the real best score over time for the live sparkline.
                    if snap.best_score_hz > 0.0 {
                        self.score_history.push(snap.best_score_hz);
                        if self.score_history.len() > 512 {
                            self.score_history.remove(0);
                        }
                    }
                    self.search_snapshot = Some(*snap);
                }
                Ok(SearchMsg::Log(msg)) => {
                    self.search_log.push(msg);
                    if self.search_log.len() > 200 {
                        self.search_log.remove(0);
                    }
                }
                Ok(SearchMsg::Done(result)) => {
                    if let Some(best) = result.best.first() {
                        let err = haversine_m((best.lat, best.lon), TRUTH) / 1000.0;
                        crate::log::info(
                            "result",
                            &format!(
                                "best {:.4},{:.4} score {:.2} · {:.1} km from truth · {} of {} cells survived",
                                best.lat,
                                best.lon,
                                best.score_hz,
                                err,
                                result.cells.iter().filter(|c| c.score_hz > 0.0).count(),
                                result.cells.len()
                            ),
                        );
                    }
                    self.search_result = Some(result);
                    self.search = SearchState::Done;
                    self.map.set_marker(Some(TRUTH));
                    open_map = true;
                }
                Ok(SearchMsg::Error(e)) => {
                    // Keep the dashboard up: the reason belongs on screen, not
                    // hidden in a side panel while the view jumps away.
                    self.search_log.push(format!("error: {e}"));
                    if e.contains("cancelled") {
                        self.search = SearchState::Idle;
                        self.search_log.push("stopped by user".into());
                    } else {
                        self.search = SearchState::Error(e);
                    }
                    ctx.request_repaint();
                }
                Err(TryRecvError::Empty) => break,
                Err(TryRecvError::Disconnected) => {
                    self.search_rx = None;
                    break;
                }
            }
        }
        if open_map {
            let c = ctx.clone();
            self.open_heatmap(&c);
        }
        if got_progress {
            ctx.request_repaint();
        }
    }

    // -- calendar -----------------------------------------------------------

    fn open_calendar(&mut self, target: CalTarget) {
        self.cal_target = target;
        let existing = match target {
            CalTarget::Start => self.time_start_unix,
            CalTarget::End => self.time_end_unix,
        }
        .unwrap_or_else(Self::now_unix);
        let (y, m, d, hh, mm, ss) = dateutil::unix_to_civil(existing);
        self.cal_year = y;
        self.cal_month = m;
        self.cal_day = d;
        self.cal_hour = hh;
        self.cal_min = mm;
        self.cal_sec = ss;
        self.cal_open = true;
    }

    fn calendar_ui(&mut self, ctx: &egui::Context) {
        if !self.cal_open {
            return;
        }
        let mut open = true;
        let mut apply = false;
        let mut close = false;
        egui::Window::new("Select date/time (UTC)")
            .collapsible(false)
            .resizable(false)
            .open(&mut open)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    if ui.button("‹").clicked() {
                        self.cal_shift_month(-1);
                    }
                    ui.label(format!("{} {}", dateutil::month_name(self.cal_month), self.cal_year));
                    if ui.button("›").clicked() {
                        self.cal_shift_month(1);
                    }
                });
                ui.separator();

                let first_wd = dateutil::weekday(dateutil::days_from_civil(
                    self.cal_year,
                    self.cal_month as i64,
                    1,
                )) as usize;
                let dim = dateutil::days_in_month(self.cal_year, self.cal_month);

                egui::Grid::new("cal_grid").num_columns(7).spacing([4.0, 4.0]).show(ui, |ui| {
                    for w in ["Su", "Mo", "Tu", "We", "Th", "Fr", "Sa"] {
                        ui.label(RichText::new(w).color(theme::MUTED).small());
                    }
                    ui.end_row();
                    let mut placed = 0u32;
                    for _ in 0..first_wd {
                        ui.label("");
                        placed += 1;
                    }
                    for d in 1..=dim {
                        let selected = d == self.cal_day;
                        let txt = if selected {
                            RichText::new(format!("{d}")).color(theme::ACCENT_BRIGHT).strong()
                        } else {
                            RichText::new(format!("{d}"))
                        };
                        if ui.button(txt).clicked() {
                            self.cal_day = d;
                        }
                        placed += 1;
                        if placed % 7 == 0 {
                            ui.end_row();
                        }
                    }
                    while placed % 7 != 0 {
                        ui.label("");
                        placed += 1;
                    }
                });

                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.label("HH");
                    ui.add(egui::DragValue::new(&mut self.cal_hour).range(0..=23));
                    ui.label("MM");
                    ui.add(egui::DragValue::new(&mut self.cal_min).range(0..=59));
                    ui.label("SS");
                    ui.add(egui::DragValue::new(&mut self.cal_sec).range(0..=59));
                });
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    if ui.button("Apply").clicked() {
                        apply = true;
                    }
                    if ui.button("Cancel").clicked() {
                        close = true;
                    }
                });
            });
        if apply {
            let unix = dateutil::civil_to_unix(
                self.cal_year,
                self.cal_month,
                self.cal_day,
                self.cal_hour,
                self.cal_min,
                self.cal_sec,
            );
            match self.cal_target {
                CalTarget::Start => self.time_start_unix = Some(unix),
                CalTarget::End => self.time_end_unix = Some(unix),
            }
            close = true;
        }
        self.cal_open = open && !close;
    }

    fn cal_shift_month(&mut self, delta: i64) {
        let m = self.cal_month as i64 + delta;
        self.cal_year += (m - 1).div_euclid(12);
        self.cal_month = ((m - 1).rem_euclid(12) + 1) as u32;
        let dim = dateutil::days_in_month(self.cal_year, self.cal_month);
        if self.cal_day > dim {
            self.cal_day = dim;
        }
    }

    // -- layout -------------------------------------------------------------

    fn layout(&mut self, ui: &mut egui::Ui) {
        // First frame: open the clip and start analysis automatically.
        if let Some(path) = self.autostart.take() {
            crate::log::info(
                "app",
                &format!("autostart: analysing {}", path.display()),
            );
            self.selected_path = Some(path);
            self.start_analysis();
        }
        self.source_panel(ui);
        self.intel_panel(ui);
        self.status_bar(ui);
        self.central_panel(ui);
    }

    fn source_panel(&mut self, ui: &mut egui::Ui) {
        egui::Panel::left("source_panel")
            .resizable(true)
            .default_size(300.0)
            .min_size(200.0)
            .max_size(420.0)
            .show(ui, |ui| {
                ui.heading(RichText::new("Avalon").color(theme::ACCENT_BRIGHT));
                ui.label(RichText::new("acoustic Doppler geolocation").color(theme::MUTED).small());
                ui.add_space(4.0);

                ui.label(RichText::new("Source").strong());
                if ui.button("Browse file…").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("Media", &["wav", "mp3", "flac", "ogg", "m4a", "aac", "mp4", "mkv", "mov", "avi", "webm", "aiff", "caf"])
                        .pick_file()
                    {
                        self.selected_path = Some(p);
                    }
                }
                if let Some(p) = &self.selected_path {
                    ui.label(
                        RichText::new(p.file_name().unwrap_or_default().to_string_lossy())
                            .color(theme::ACCENT_BRIGHT),
                    );
                }
                ui.label(
                    RichText::new("or drag & drop a file onto the window")
                        .color(theme::MUTED)
                        .small(),
                );

                let running = matches!(self.analysis, AnalysisState::Running { .. });
                let btn = ui.add_enabled(
                    !running && self.selected_path.is_some(),
                    egui::Button::new("▶ Analyze (detect ridges)"),
                );
                if btn.clicked() {
                    self.start_analysis();
                }

                ui.add_space(6.0);
                ui.collapsing("Analysis settings", |ui| {
                    egui::Grid::new("settings_grid").num_columns(2).show(ui, |ui| {
                        ui.label("Min freq (Hz)");
                        ui.add(egui::DragValue::new(&mut self.settings.min_freq_hz).range(1.0..=2000.0).speed(1.0));
                        ui.end_row();
                        ui.label("Max freq (Hz)");
                        ui.add(egui::DragValue::new(&mut self.settings.max_freq_hz).range(2.0..=12000.0).speed(10.0));
                        ui.end_row();
                        ui.label("Seed contrast");
                        ui.add(egui::DragValue::new(&mut self.settings.seed_contrast).range(1.0..=10.0).speed(0.1));
                        ui.end_row();
                        ui.label("Window (s)");
                        ui.add(egui::DragValue::new(&mut self.settings.window_s).range(0.5..=8.0).speed(0.1));
                        ui.end_row();
                        ui.label("Hop (s)");
                        ui.add(egui::DragValue::new(&mut self.settings.hop_s).range(0.1..=2.0).speed(0.05));
                        ui.end_row();
                        ui.label("FFT oversample");
                        ui.add(
                            egui::DragValue::new(&mut self.settings.fft_pad)
                                .range(1..=16)
                                .speed(1.0),
                        );
                        ui.end_row();
                    });
                    let true_res = 1.0 / self.settings.window_s.max(0.01);
                    let bin_hz = self.settings.target_sample_rate as f32
                        / ((self.settings.window_s * self.settings.target_sample_rate as f32)
                            .round() as usize
                            * self.settings.fft_pad.max(1) as usize)
                            .next_power_of_two()
                            .max(1) as f32;
                    ui.label(
                        RichText::new(format!(
                            "resolution {true_res:.3} Hz · bin spacing {bin_hz:.3} Hz · {:.0} frames/s",
                            1.0 / self.settings.hop_s.max(0.01)
                        ))
                        .color(theme::MUTED)
                        .small(),
                    );
                });

                ui.collapsing("Log", |ui| {
                    for line in &self.log {
                        ui.label(RichText::new(line).color(theme::MUTED).small());
                    }
                });
            });
    }

    fn intel_panel(&mut self, ui: &mut egui::Ui) {
        egui::Panel::right("intel_panel")
            .resizable(true)
            .default_size(340.0)
            .min_size(260.0)
            // Hard ceiling: the panel used to widen itself to fit long event
            // lines, squeezing the spectrogram and overriding the stored width.
            .max_size(420.0)
            .show(ui, |ui| {
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                    ui.heading(RichText::new("Intelligence").color(theme::ACCENT_BRIGHT));

                    // Region.
                    ui.add_space(4.0);
                    ui.label(RichText::new("Geographic region").strong());
                    ui.horizontal(|ui| {
                        if ui.button("Set on map").clicked() {
                            self.map.set_rectangles(self.region.clone());
                            self.map.clear_heatmap();
                            self.map.mode = Mode::Select;
                            self.map_purpose = MapPurpose::Region;
                            self.map_open = true;
                        }
                        if ui.button("Clear").clicked() {
                            self.region.clear();
                        }
                    });
                    if self.region.is_empty() {
                        ui.label(
                            RichText::new("No region — searches the whole world (very inefficient)")
                                .color(theme::WARN),
                        );
                    } else {
                        ui.label(format!("{} region(s)", self.region.len()));
                        for b in &self.region {
                            ui.label(
                                RichText::new(format!(
                                    "  W {:.2}  S {:.2}  E {:.2}  N {:.2}",
                                    b.west, b.south, b.east, b.north
                                ))
                                .color(theme::MUTED)
                                .small(),
                            );
                        }
                    }

                    ui.separator();

                    // Time window (mandatory).
                    ui.label(RichText::new("Time window (UTC, mandatory)").strong());
                    ui.horizontal(|ui| {
                        ui.label("Start");
                        let s = self.time_start_unix.map(fmt_civil).unwrap_or_else(|| "—".into());
                        if ui.button(s).clicked() {
                            self.open_calendar(CalTarget::Start);
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("End");
                        let s = self.time_end_unix.map(fmt_civil).unwrap_or_else(|| "—".into());
                        if ui.button(s).clicked() {
                            self.open_calendar(CalTarget::End);
                        }
                    });
                    ui.label(
                        RichText::new("Searching all of history is infeasible — a window is required.")
                            .color(theme::MUTED)
                            .small(),
                    );

                    ui.separator();

                    // Provider.
                    ui.label(RichText::new("Flight data").strong());
                    egui::ComboBox::from_id_salt("provider")
                        .selected_text(self.provider.clone())
                        .show_ui(ui, |ui| {
                            for p in PROVIDERS {
                                ui.selectable_value(&mut self.provider, p.to_string(), p);
                            }
                        });

                    // Frequency band hint.
                    ui.checkbox(&mut self.freq_enabled, "Frequency band hint (Hz)");
                    if self.freq_enabled {
                        ui.horizontal(|ui| {
                            ui.label("Min");
                            ui.text_edit_singleline(&mut self.freq_min);
                            ui.label("Max");
                            ui.text_edit_singleline(&mut self.freq_max);
                        });
                    }

                    ui.separator();
                    ui.label("Notes");
                    ui.text_edit_multiline(&mut self.notes);

                    ui.add_space(8.0);

                    // Events.
                    ui.label(RichText::new("Detected events").strong());
                    self.events_ui(ui);

                    ui.separator();

                    // Process + search results.
                    self.process_ui(ui);
                });
            });
    }

    fn events_ui(&mut self, ui: &mut egui::Ui) {
        struct Row {
            id: usize,
            start: f64,
            end: f64,
            n_ridges: usize,
            fit: Option<DopplerFit>,
            plausible: bool,
        }
        let rows: Vec<Row> = self
            .result
            .as_ref()
            .map(|r| {
                r.events
                    .iter()
                    .map(|e| {
                        let rep = r.ridges.iter().find(|x| x.id == e.representative_ridge_id);
                        Row {
                            id: e.id,
                            start: e.start_s,
                            end: e.end_s,
                            n_ridges: e.ridge_ids.len(),
                            fit: e.doppler,
                            plausible: rep.map(|x| x.plausible).unwrap_or(false),
                        }
                    })
                    .collect()
            })
            .unwrap_or_default();
        if rows.is_empty() {
            ui.label(RichText::new("no Doppler events yet").color(theme::MUTED).small());
            return;
        }
        let n_plausible = rows.iter().filter(|r| r.plausible).count();
        ui.label(
            RichText::new(format!(
                "{n_plausible} flyover curve(s) · click a line to select",
            ))
            .color(theme::MUTED)
            .small(),
        );
        for row in rows {
            let selected = self.selected_events.contains(&row.id);
            let marker = if selected { "●" } else { "○" };
            let mut txt = RichText::new(format!(
                "{marker} Event {} · {} tone(s) · {:.1}–{:.1}s",
                row.id, row.n_ridges, row.start, row.end
            ));
            txt = if selected {
                txt.color(theme::ACCENT_BRIGHT).strong()
            } else if row.plausible {
                txt.color(theme::TEXT)
            } else {
                txt.color(theme::MUTED)
            };
            ui.label(txt);

            if let Some(f) = row.fit {
                let detail = format!(
                    "    f0 {:.0} Hz · {:.0} m/s · {:.0} m · {:.2} Hz · {:.0}%",
                    f.f0_hz,
                    f.speed_mps,
                    f.slant_distance_m,
                    f.rmse_hz,
                    f.inlier_fraction * 100.0
                );
                ui.add(
                    egui::Label::new(RichText::new(detail).color(theme::MUTED).small())
                        .truncate(),
                );
            }
            if !row.plausible {
                ui.label(
                    RichText::new("    weak fit — not a flyover curve")
                        .color(theme::WARN)
                        .small(),
                );
            }
        }
    }

    fn process_ui(&mut self, ui: &mut egui::Ui) {
        let running = matches!(self.search, SearchState::Running { .. });
        let btn = ui.add_enabled(!running, egui::Button::new("⟳ Process (locate)"));
        if btn.clicked() {
            self.start_search();
        }

        match &self.search {
            SearchState::Running { stage } => {
                let mut cancel = false;
                ui.horizontal(|ui| {
                    ui.add(egui::ProgressBar::new(0.0).text(stage));
                    if ui.button("Stop").clicked() {
                        cancel = true;
                    }
                });
                if cancel {
                    self.cancel
                        .store(true, std::sync::atomic::Ordering::Relaxed);
                }
            }
            SearchState::Error(e) => {
                ui.label(RichText::new(e).color(theme::ERROR));
            }
            _ => {}
        }

        if let Some(res) = &self.search_result {
            ui.add_space(6.0);
            ui.label(RichText::new("Results").strong());
            ui.label(format!(
                "{} flights · {} candidates · step {:.4}°",
                res.n_flights, res.n_candidates, res.step_deg
            ));
            for m in res.modes.iter().take(6) {
                ui.label(
                    RichText::new(format!(
                        "  {:.3} @ ({:.4}, {:.4}) · {} cells · {:.1} Hz · {}",
                        m.score, m.lat, m.lon, m.cells, m.rmse_hz, m.pair
                    ))
                    .color(theme::TEXT)
                    .small(),
                );
            }
            // Near-ties are the normal continental failure mode: two unrelated
            // flight pairs can fit both measured curves to within a few
            // hundredths of a Hz. Say so instead of presenting one point as fact.
            if let Some(w) = &res.warning {
                ui.add_space(4.0);
                ui.label(RichText::new(w).color(theme::WARN).small());
            }
            if let (Some(a), Some(b)) = (res.modes.first(), res.modes.get(1)) {
                let margin = (a.score - b.score) / a.score.max(1e-6);
                if margin < 0.20 && res.warning.is_none() {
                    ui.label(
                        RichText::new(format!(
                            "⚠ near-tie: {:.0}% margin between the top {} modes — pick using other evidence",
                            margin * 100.0,
                            res.modes.len()
                        ))
                        .color(theme::WARN),
                    );
                }
            }
            if ui.button("Show heatmap on map").clicked() {
                let ctx = ui.ctx().clone();
                self.open_heatmap(&ctx);
            }
        }
    }

    /// Open the map with the current search heatmap.
    fn open_heatmap(&mut self, ctx: &egui::Context) {
        let Some(res) = self.search_result.clone() else { return };
        // Zoom so the searched box fills the map: pick the zoom whose span in
        // degrees is closest to the box size.
        if let Some(b) = res.bounds {
            let (lat, lon) = b.center();
            let span = (b.north - b.south).max(b.east - b.west).max(0.02);
            let zoom = (360.0 / span / 3.0).log2().clamp(3.0, 11.0) as f32;
            self.map = MapState::new(lat, lon, zoom);
        }
        self.map.set_heatmap(ctx, &res.cells, res.step_deg);
        self.map.set_marker(Some(TRUTH));
        self.map.mode = Mode::Pan;
        self.map_purpose = MapPurpose::Heatmap;
        self.map_open = true;

        if let Some(best) = res.best.first() {
            let err = haversine_m((best.lat, best.lon), TRUTH) / 1000.0;
            crate::log::info(
                "result",
                &format!(
                    "best cell {:.4},{:.4} (score {:.2}) is {:.1} km from truth; {} of {} cells survived",
                    best.lat,
                    best.lon,
                    best.score_hz,
                    err,
                    res.cells.iter().filter(|c| c.score_hz > 0.0).count(),
                    res.cells.len()
                ),
            );
        }
    }

    /// Live processing dashboard: only measured values, throttled updates.
    fn processing_view(&mut self, ui: &mut egui::Ui) {
        // Outcome banner: running / complete / failed. A failed locate used to
        // drop straight back to the spectrogram, which looked like nothing
        // happened.
        let banner: Option<(Color32, String)> = match &self.search {
            SearchState::Running { .. } => None,
            SearchState::Done => Some((
                theme::ACCENT_BRIGHT,
                format!(
                    "Locate complete — {} flights, {} candidates",
                    self.search_result.as_ref().map(|r| r.n_flights).unwrap_or(0),
                    self.search_result.as_ref().map(|r| r.n_candidates).unwrap_or(0),
                ),
            )),
            SearchState::Error(e) => Some((theme::ERROR, format!("Locate failed: {e}"))),
            SearchState::Idle => Some((theme::WARN, "Stopped".into())),
        };

        if let Some((color, text)) = banner {
            let mut want_heatmap = false;
            ui.horizontal(|ui| {
                ui.label(RichText::new(text).color(color).strong());
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("Back to spectrogram").clicked() {
                        self.show_processing = false;
                    }
                    if ui.button("Retry").clicked() {
                        self.start_search();
                    }
                    if matches!(self.search, SearchState::Done) && ui.button("Show heatmap").clicked() {
                        want_heatmap = true;
                    }
                });
            });
            if want_heatmap {
                let ctx = ui.ctx().clone();
                self.open_heatmap(&ctx);
            }
            ui.separator();
        }

        let Some(snap) = self.search_snapshot.clone() else {
            ui.add_space(ui.available_height() * 0.3);
            ui.vertical_centered(|ui| {
                ui.heading(RichText::new("Starting locate…").color(theme::TEXT));
            });
            return;
        };

        ui.horizontal(|ui| {
            ui.heading(RichText::new("Locating").color(theme::ACCENT_BRIGHT));
            ui.label(
                RichText::new(format!("· {}", snap.stage)).color(theme::TEXT).strong(),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let mut cancel = false;
                if ui.button("■ Stop").clicked() {
                    cancel = true;
                }
                if cancel {
                    self.cancel
                        .store(true, std::sync::atomic::Ordering::Relaxed);
                }
                ui.label(
                    RichText::new(format!(
                        "elapsed {}   eta {}",
                        fmt_duration(snap.elapsed_s),
                        if snap.eta_s > 0.5 { fmt_duration(snap.eta_s) } else { "—".into() }
                    ))
                    .color(theme::MUTED)
                    .monospace(),
                );
            });
        });

        ui.add_space(4.0);
        ui.add(
            egui::ProgressBar::new(snap.fraction)
                .desired_height(18.0)
                .text(
                    egui::RichText::new(format!("{:.1}%", snap.fraction * 100.0))
                        .monospace()
                        .size(11.0),
                ),
        );
        ui.label(RichText::new(&snap.detail).color(theme::MUTED).monospace().small());

        ui.add_space(8.0);

        // Real counters, two columns of tiles.
        let tiles: [(&str, String); 8] = [
            (
                "archive read",
                if snap.bytes_total > 0 {
                    format!(
                        "{:.0} / {:.0} MB",
                        snap.bytes_read as f64 / 1e6,
                        snap.bytes_total as f64 / 1e6
                    )
                } else {
                    format!("{:.0} MB", snap.bytes_read as f64 / 1e6)
                },
            ),
            ("traces parsed", format!("{}", snap.traces_seen)),
            ("flights in region", format!("{}", snap.flights_kept)),
            (
                "days",
                format!("{}/{}", snap.days_done.max(1), snap.days_total.max(1)),
            ),
            (
                "candidates",
                format!("{} / {}", snap.candidates_done, snap.candidates_total),
            ),
            ("throughput", format!("{:.0} cells/s", snap.cells_per_s)),
            ("cells matched", format!("{}", snap.cells_matched)),
            (
                "best residual",
                if snap.best_score_hz > 0.0 {
                    format!(
                        "{:.2} @ {:.3},{:.3}",
                        snap.best_score_hz, snap.best_lat, snap.best_lon
                    )
                } else {
                    "—".into()
                },
            ),
        ];
        egui::Grid::new("proc_tiles")
            .num_columns(4)
            .spacing([18.0, 6.0])
            .show(ui, |ui| {
                for (i, (label, value)) in tiles.iter().enumerate() {
                    tile(ui, label, value);
                    if i % 4 == 3 {
                        ui.end_row();
                    }
                }
            });

        if !snap.best_flight.is_empty() {
            ui.label(
                RichText::new(format!("best match ICAO {}", snap.best_flight))
                    .color(theme::ACCENT_BRIGHT)
                    .monospace()
                    .small(),
            );
        }

        ui.add_space(10.0);
        ui.label(RichText::new("Best residual over grid progress").color(theme::MUTED).small());
        // Sparkline of the real running best; cheap (<=512 points, no axes).
        let pts: Vec<[f64; 2]> = self
            .score_history
            .iter()
            .enumerate()
            .map(|(i, v)| [i as f64, *v as f64])
            .collect();
        let height = (ui.available_height() * 0.32).clamp(70.0, 180.0);
        if pts.len() > 1 {
            egui_plot::Plot::new("score_spark")
                .height(height)
                .allow_drag(false)
                .allow_zoom(false)
                .allow_scroll(false)
                .allow_boxed_zoom(false)
                .show_axes([false, true])
                .show_grid(false)
                .y_axis_formatter(|m, _| format!("{:.1}", m.value))
                .show(ui, |plot_ui| {
                    plot_ui.line(
                        egui_plot::Line::new("best", egui_plot::PlotPoints::new(pts))
                            .color(theme::ACCENT)
                            .width(2.0),
                    );
                });
        } else {
            ui.allocate_space(egui::vec2(ui.available_width(), height));
        }

        ui.add_space(6.0);
        ui.label(RichText::new("Log").color(theme::MUTED).small());
        egui::ScrollArea::vertical()
            .max_height(ui.available_height().max(40.0))
            .stick_to_bottom(true)
            .show(ui, |ui| {
                for line in &self.search_log {
                    ui.label(RichText::new(line).color(theme::MUTED).monospace().small());
                }
            });
    }

    fn central_panel(&mut self, ui: &mut egui::Ui) {
        egui::CentralPanel::default().show(ui, |ui| {
            // The locate stage is where the real work happens; keep it in the
            // main area for the whole run *and* afterwards, so a failure shows
            // its reason here instead of silently reverting to the spectrogram.
            if self.show_processing {
                self.processing_view(ui);
                return;
            }
            match &self.analysis {
                AnalysisState::Running { fraction, stage } => {
                    ui.add_space(ui.available_height() * 0.4);
                    ui.vertical_centered(|ui| {
                        ui.heading(RichText::new(stage).color(theme::TEXT));
                        ui.add(egui::ProgressBar::new(*fraction).desired_width(360.0).text(format!("{:.0}%", fraction * 100.0)));
                    });
                }
                AnalysisState::Error(e) => {
                    ui.add_space(ui.available_height() * 0.4);
                    ui.vertical_centered(|ui| {
                        ui.label(RichText::new("Analysis failed").color(theme::ERROR).strong());
                        ui.label(RichText::new(e).color(theme::ERROR));
                    });
                }
                _ if self.result.is_none() => {
                    ui.add_space(ui.available_height() * 0.35);
                    ui.vertical_centered(|ui| {
                        ui.label(
                            RichText::new("Drop a video/audio file here, or browse on the left")
                                .color(theme::MUTED),
                        );
                        ui.label(
                            RichText::new("Analyze detects flyover Doppler signatures; then set region + time and Process to locate")
                                .color(theme::MUTED)
                                .small(),
                        );
                    });
                }
                _ => {
                    self.render_spectrogram(ui);
                }
            }
        });
    }

    fn status_bar(&mut self, ui: &mut egui::Ui) {
        egui::Panel::bottom("status_bar").show(ui, |ui| {
            ui.horizontal(|ui| {
                match &self.result {
                    Some(r) => {
                        ui.label(format!(
                            "{} · {} Hz · {:.1}s · {} events · {}",
                            r.source_path, r.sample_rate, r.duration_s, r.events.len(), r.backend
                        ));
                    }
                    None => {
                        ui.label("idle");
                    }
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let backends = crate::decode::decoder_backend_available();
                    ui.label(RichText::new(backends.join(" · ")).color(theme::MUTED).small());
                });
            });
        });
    }

    // -- map modal ----------------------------------------------------------

    fn map_modal(&mut self, ctx: &egui::Context) {
        if !self.map_open {
            return;
        }
        let mut open = true;
        let mut apply = false;
        let mut clear = false;
        let mut close = false;
        let title = match self.map_purpose {
            MapPurpose::Region => "Select search region",
            MapPurpose::Heatmap => "Candidate heatmap",
        };
        egui::Window::new(title)
            .default_size([840.0, 620.0])
            .min_size([400.0, 300.0])
            .open(&mut open)
            .show(ctx, |ui| {
                // Action bar first, so it is always visible.
                ui.horizontal(|ui| {
                    match self.map_purpose {
                        MapPurpose::Region => {
                            if ui.button("✔ Apply region").clicked() {
                                apply = true;
                            }
                            if ui.button("Clear").clicked() {
                                clear = true;
                            }
                        }
                        MapPurpose::Heatmap => {
                            let survivors = self
                                .search_result
                                .as_ref()
                                .map(|r| r.cells.iter().filter(|c| c.score_hz > 0.0).count())
                                .unwrap_or(0);
                            let total = self
                                .search_result
                                .as_ref()
                                .map(|r| r.cells.len())
                                .unwrap_or(0);
                            ui.label(
                                RichText::new(format!(
                                    "{survivors} of {total} cells passed the hard filters · green = eliminated → orange → red → purple = highest"
                                ))
                                .color(theme::MUTED),
                            );
                        }
                    }
                    if ui.button("Close").clicked() {
                        close = true;
                    }
                });
                ui.separator();
                self.map.show(ui, ctx);
            });
        if clear {
            self.map.clear_rectangles();
        }
        if apply {
            self.region = self.map.rectangles().to_vec();
            close = true;
        }
        self.map_open = open && !close;
    }

    fn handle_drop(&mut self, ctx: &egui::Context) {
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        if let Some(f) = dropped.first() {
            self.selected_path = Some(f.path().to_path_buf());
        }
    }

    // -- spectrogram --------------------------------------------------------

    fn render_spectrogram(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let dims = self
            .result
            .as_ref()
            .map(|r| (r.spectrogram.time_axis_s.len(), r.spectrogram.freq_axis_hz.len()));
        if let Some((w, h)) = dims {
            let key = format!("{w}x{h}");
            let rebuild = self.spectro_texture.as_ref().map(|(k, _)| k != &key).unwrap_or(true);
            if rebuild {
                let image = spectrogram_to_image(&self.result.as_ref().unwrap().spectrogram);
                let handle = ctx.load_texture("spectrogram", image, egui::TextureOptions::LINEAR);
                self.spectro_texture = Some((key, handle));
            }
        }

        let Some(result) = self.result.as_ref() else { return };
        let sim_curves: Vec<Vec<(f64, f32)>> = self
            .search_result
            .as_ref()
            .map(|r| r.model_curves.clone())
            .unwrap_or_default();
        let Some((_, tex)) = self.spectro_texture.as_ref() else { return };
        let tex_id = tex.id();

        let spec = &result.spectrogram;
        let n_frames = spec.time_axis_s.len();
        let n_freq = spec.freq_axis_hz.len();
        if n_frames == 0 || n_freq == 0 {
            return;
        }
        let min_f = spec.freq_axis_hz.first().copied().unwrap_or(0.0);
        let max_f = spec.freq_axis_hz.last().copied().unwrap_or(0.0);
        let duration = spec.time_axis_s.last().copied().unwrap_or(0.0).max(1e-6);

        let selected_events: Vec<&DetectedEvent> = result
            .events
            .iter()
            .filter(|e| self.selected_events.contains(&e.id))
            .collect();
        // Doppler fit curve for each toggled event.
        let selected_fits: Vec<(DopplerFit, (f64, f64))> = selected_events
            .iter()
            .filter_map(|ev| ev.doppler.map(|f| (f, (ev.start_s, ev.end_s))))
            .collect();

        // The spectrogram is a time series: give it a landscape aspect that
        // grows with clip length instead of filling a square panel.
        let avail = ui.available_size();
        let aspect = (result.duration_s as f32 / 15.0).clamp(1.4, 3.4);
        let plot_height = (avail.x / aspect).min(avail.y).max(120.0);
        let spare = (avail.y - plot_height).max(0.0);
        if spare > 0.0 {
            ui.add_space(spare * 0.18);
        }

        let plot_response = Plot::new("spectrogram")
            .height(plot_height)
            .show_x(true)
            .show_y(true)
            .allow_drag(true)
            .allow_zoom(true)
            // Wheel is handled as zoom below; letting the plot scroll it too
            // would pan the view instead (egui only zooms on ctrl+wheel).
            .allow_scroll(false)
            .allow_boxed_zoom(true)
            .allow_double_click_reset(true)
            .x_axis_formatter(|m, _| format!("{:.0}s", m.value))
            .y_axis_formatter(|m, _| format!("{:.0}Hz", m.value))
            .show(ui, |plot_ui| {
                // Plain wheel zooms around the cursor.
                let scroll = plot_ui.ctx().input(|i| i.smooth_scroll_delta.y);
                if scroll != 0.0 && plot_ui.response().hovered() {
                    let factor = (scroll * 0.0025).exp();
                    plot_ui.zoom_bounds_around_hovered(egui::vec2(factor, factor));
                }
                plot_ui.image(PlotImage::new(
                    "spectrogram-image",
                    tex_id,
                    PlotPoint::new(duration / 2.0, (min_f as f64 + max_f as f64) / 2.0),
                    egui::Vec2::new(duration as f32, (max_f - min_f) as f32),
                ));

                // Draw every ridge: plausible flyover curves bold, weak fits dim
                // but still clickable, selected ones green and thick.
                for ridge in &result.ridges {
                    let pts: Vec<[f64; 2]> =
                        ridge.samples.iter().map(|s| [s.t_s, s.freq_hz as f64]).collect();
                    let owner = result
                        .events
                        .iter()
                        .find(|e| e.ridge_ids.contains(&ridge.id))
                        .map(|e| e.id);
                    let selected = owner
                        .map(|id| self.selected_events.contains(&id))
                        .unwrap_or(false);
                    let (color, width) = if selected {
                        (theme::ACCENT_BRIGHT, 4.0)
                    } else if ridge.plausible {
                        (Color32::from_rgb(120, 190, 200), 2.5)
                    } else {
                        (Color32::from_rgb(78, 96, 104), 1.5)
                    };
                    plot_ui.line(
                        Line::new(format!("ridge{}", ridge.id), PlotPoints::new(pts))
                            .color(color)
                            .width(width),
                    );
                }

                // Simulated curve at the winning cell, over the measured ridges:
                // a visual check that the model follows the real Doppler.
                for (i, curve) in sim_curves.iter().enumerate() {
                    if curve.len() < 2 {
                        continue;
                    }
                    let pts: Vec<[f64; 2]> =
                        curve.iter().map(|(t, f)| [*t, *f as f64]).collect();
                    plot_ui.line(
                        Line::new(format!("sim-{i}"), PlotPoints::new(pts))
                            .color(Color32::from_rgb(90, 220, 255))
                            .width(2.5),
                    );
                }

                for (i, (fit, (t0, t1))) in selected_fits.iter().enumerate() {
                    let n = 240;
                    let pts: Vec<[f64; 2]> = (0..=n)
                        .map(|j| {
                            let t = t0 + (t1 - t0) * j as f64 / n as f64;
                            [t, doppler::evaluate_fit(fit, t)]
                        })
                        .collect();
                    plot_ui.line(
                        Line::new(format!("fit-{i}"), PlotPoints::new(pts))
                            .color(theme::ACCENT)
                            .width(3.0),
                    );
                }
            });

        // Pixel-space hit test: pick the curve nearest the click within a
        // generous radius (thin-line hit boxes are unreliable).
        let mut toggle_event: Option<usize> = None;
        if plot_response.response.clicked() {
            if let Some(pointer) = plot_response.response.interact_pointer_pos() {
                let mut best: Option<(f32, usize)> = None;
                for ev in &result.events {
                    for rid in &ev.ridge_ids {
                        let Some(ridge) = result.ridges.iter().find(|r| r.id == *rid) else {
                            continue;
                        };
                        for s in &ridge.samples {
                            let p = plot_response
                                .transform
                                .position_from_point(&PlotPoint::new(s.t_s, s.freq_hz as f64));
                            let d = p.distance(pointer);
                            if d <= 12.0 && best.map_or(true, |(bd, _)| d < bd) {
                                best = Some((d, ev.id));
                            }
                        }
                    }
                }
                toggle_event = best.map(|(_, id)| id);
            }
        }
        if let Some(ev_id) = toggle_event {
            if !self.selected_events.remove(&ev_id) {
                self.selected_events.insert(ev_id);
            }
        }
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.poll_analysis(&ctx);
        self.poll_search(&ctx);
        self.handle_drop(&ctx);
        self.layout(ui);
        self.map_modal(&ctx);
        self.calendar_ui(&ctx);
    }
}

// ---------------------------------------------------------------------------
// Spectrogram -> texture
// ---------------------------------------------------------------------------

fn spectrogram_to_image(spec: &Spectrogram) -> egui::ColorImage {
    let n_frames = spec.time_axis_s.len();
    let n_freq = spec.freq_axis_hz.len();
    if n_frames == 0 || n_freq == 0 {
        return egui::ColorImage::filled([1, 1], theme::BG);
    }
    let (max_w, max_h) = (2048usize, 1024usize);
    let sx = (n_frames + max_w - 1) / max_w;
    let sy = (n_freq + max_h - 1) / max_h;
    let w = ((n_frames + sx - 1) / sx).max(1);
    let h = ((n_freq + sy - 1) / sy).max(1);

    // Defensive: rows can be shorter than the time axis if a producer ever
    // emits a ragged matrix. Never index past a row.
    let row_len = spec.contrast.iter().map(|r| r.len()).min().unwrap_or(0);
    if row_len == 0 {
        return egui::ColorImage::filled([1, 1], theme::BG);
    }
    let mut vals: Vec<f32> = Vec::with_capacity(w * h);
    for k in (0..n_freq).step_by(sy.max(1)) {
        for c in (0..n_frames.min(row_len)).step_by(sx.max(1)) {
            vals.push(spec.contrast[k][c]);
        }
    }
    if vals.is_empty() {
        return egui::ColorImage::filled([1, 1], theme::BG);
    }
    let idx = (((vals.len() as f64) * 0.99) as usize).min(vals.len() - 1);
    let p99 = *vals
        .select_nth_unstable_by(idx, |a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .1;
    let p99 = p99.max(2.0);

    let mut pixels = Vec::with_capacity(w * h);
    for r in 0..h {
        let k = n_freq - 1 - (r * sy).min(n_freq - 1);
        let row = &spec.contrast[k];
        for c in 0..w {
            let frame = (c * sx).min(row_len - 1);
            let val = row[frame];
            let t = ((val - 1.0) / (p99 - 1.0)).clamp(0.0, 1.0);
            pixels.push(colormap(t));
        }
    }
    egui::ColorImage::new([w, h], pixels)
}

fn colormap(t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let stops: [(f32, (u8, u8, u8)); 5] = [
        (0.00, (8, 8, 12)),
        (0.25, (60, 18, 92)),
        (0.50, (150, 30, 90)),
        (0.75, (230, 110, 40)),
        (1.00, (250, 230, 140)),
    ];
    for i in 0..stops.len() - 1 {
        let (t0, c0) = stops[i];
        let (t1, c1) = stops[i + 1];
        if t <= t1 {
            let f = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
            let r = (c0.0 as f32 + (c1.0 as f32 - c0.0 as f32) * f) as u8;
            let g = (c0.1 as f32 + (c1.1 as f32 - c0.1 as f32) * f) as u8;
            let b = (c0.2 as f32 + (c1.2 as f32 - c0.2 as f32) * f) as u8;
            return Color32::from_rgb(r, g, b);
        }
    }
    Color32::from_rgb(250, 230, 140)
}

fn fmt_civil(unix: i64) -> String {
    let (y, m, d, hh, mm, ss) = dateutil::unix_to_civil(unix);
    format!("{y:04}-{m:02}-{d:02} {hh:02}:{mm:02}:{ss:02}")
}

/// Compact duration: `4.2s`, `1m 07s`, `1h 12m`.
fn fmt_duration(secs: f64) -> String {
    let s = secs.max(0.0);
    if s < 60.0 {
        format!("{s:.1}s")
    } else if s < 3600.0 {
        format!("{}m {:02}s", (s / 60.0) as u64, (s % 60.0) as u64)
    } else {
        format!("{}h {:02}m", (s / 3600.0) as u64, ((s % 3600.0) / 60.0) as u64)
    }
}

/// One labelled statistic tile.
fn tile(ui: &mut egui::Ui, label: &str, value: &str) {
    ui.vertical(|ui| {
        ui.label(RichText::new(label).color(theme::MUTED).small());
        ui.label(
            RichText::new(value)
                .color(theme::ACCENT_BRIGHT)
                .monospace()
                .size(15.0),
        );
    });
}
