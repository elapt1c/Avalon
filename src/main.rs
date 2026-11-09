//! Avalon — acoustic Doppler geolocation from video/audio sources.

use avalon::app::App;

fn main() -> eframe::Result {
    avalon::log::init();
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1600.0, 1000.0])
            .with_min_inner_size([1000.0, 640.0])
            .with_title("AVALON — Acoustic Signal Geolocation"),
        ..Default::default()
    };
    eframe::run_native(
        "Avalon",
        options,
        Box::new(|cc| Ok(Box::new(App::new(cc)))),
    )
}
