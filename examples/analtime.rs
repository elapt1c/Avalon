//! Time the analysis stage alone (decode → spectrogram → ridges → events).
//!
//! Usage: analtime [clip]

use std::path::Path;
use std::time::Instant;

use avalon::model::{AnalysisSettings, Metadata};
use avalon::pipeline;

fn main() -> anyhow::Result<()> {
    let clip = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "/root/projects/Unknown location - raw audio - colsto (128p).mp4".into());
    let settings = AnalysisSettings { target_sample_rate: 4000, ..Default::default() };
    let meta = Metadata::default();
    let t = Instant::now();
    let a = pipeline::run(Path::new(&clip), &settings, &meta, None)?;
    println!(
        "ANALTIME {:.2} s · {} ridges · {} events · {:.0}s audio",
        t.elapsed().as_secs_f64(),
        a.ridges.len(),
        a.events.len(),
        a.duration_s
    );
    Ok(())
}
