//! Dark, QGIS-inspired application theme.

use egui::{Color32, CornerRadius, Stroke};

// Palette.
pub const BG: Color32 = Color32::from_rgb(27, 29, 31); // canvas / deepest background
pub const PANEL: Color32 = Color32::from_rgb(38, 41, 44); // side panels
pub const PANEL_2: Color32 = Color32::from_rgb(47, 51, 55); // raised widgets
pub const BORDER: Color32 = Color32::from_rgb(58, 63, 68);
pub const TEXT: Color32 = Color32::from_rgb(214, 214, 214);
pub const MUTED: Color32 = Color32::from_rgb(140, 145, 150);
pub const ACCENT: Color32 = Color32::from_rgb(86, 176, 128); // teal-green
pub const ACCENT_BRIGHT: Color32 = Color32::from_rgb(112, 210, 162);
pub const WARN: Color32 = Color32::from_rgb(214, 158, 62);
pub const ERROR: Color32 = Color32::from_rgb(214, 92, 92);

pub fn apply(ctx: &egui::Context) {
    let mut v = egui::Visuals::dark();
    v.dark_mode = true;
    v.panel_fill = PANEL;
    v.window_fill = PANEL;
    v.extreme_bg_color = BG;
    v.faint_bg_color = PANEL_2;
    v.code_bg_color = PANEL_2;
    v.text_edit_bg_color = Some(BG);
    v.override_text_color = Some(TEXT);
    v.weak_text_color = Some(MUTED);
    v.hyperlink_color = ACCENT_BRIGHT;
    v.warn_fg_color = WARN;
    v.error_fg_color = ERROR;
    v.selection.bg_fill = ACCENT.gamma_multiply(0.35);
    v.selection.stroke = Stroke::new(1.0, ACCENT_BRIGHT);
    v.window_stroke = Stroke::new(1.0, BORDER);

    // Widget fills.
    let r = CornerRadius::same(4);
    v.widgets.noninteractive.bg_fill = PANEL_2;
    v.widgets.noninteractive.corner_radius = r;
    v.widgets.inactive.bg_fill = PANEL_2;
    v.widgets.inactive.corner_radius = r;
    v.widgets.inactive.fg_stroke = Stroke::new(1.0, TEXT);
    v.widgets.hovered.bg_fill = Color32::from_rgb(58, 64, 69);
    v.widgets.hovered.corner_radius = r;
    v.widgets.hovered.fg_stroke = Stroke::new(1.0, ACCENT_BRIGHT);
    v.widgets.active.bg_fill = ACCENT;
    v.widgets.active.corner_radius = r;
    v.widgets.active.fg_stroke = Stroke::new(1.0, Color32::from_rgb(16, 20, 18));
    v.widgets.open.bg_fill = PANEL_2;
    v.widgets.open.corner_radius = r;

    ctx.set_visuals(v);
}
