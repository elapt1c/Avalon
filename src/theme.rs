//! Tactical HUD theme.
//!
//! Design language: near-black substrate, hairline borders, monospaced data,
//! one cyan signal accent with amber/red reserved for state. Everything is
//! uppercase and tightly spaced so the UI reads like an instrument panel rather
//! than a form.

use egui::{Color32, CornerRadius, Stroke};

// -- palette ---------------------------------------------------------------
pub const VOID: Color32 = Color32::from_rgb(5, 7, 9); // behind everything
pub const BG: Color32 = Color32::from_rgb(8, 11, 13); // canvas
pub const PANEL: Color32 = Color32::from_rgb(12, 17, 20); // side panels
pub const PANEL_2: Color32 = Color32::from_rgb(18, 25, 29); // raised widgets
pub const PANEL_3: Color32 = Color32::from_rgb(24, 33, 38); // hover
pub const BORDER: Color32 = Color32::from_rgb(36, 50, 57); // hairlines
pub const BORDER_HI: Color32 = Color32::from_rgb(48, 68, 76);

pub const TEXT: Color32 = Color32::from_rgb(198, 213, 219);
pub const MUTED: Color32 = Color32::from_rgb(126, 147, 155);
pub const DIM: Color32 = Color32::from_rgb(84, 103, 111);

pub const ACCENT: Color32 = Color32::from_rgb(46, 214, 196); // signal cyan
pub const ACCENT_BRIGHT: Color32 = Color32::from_rgb(120, 245, 226);
pub const ACCENT_DIM: Color32 = Color32::from_rgb(26, 122, 112);
pub const WARN: Color32 = Color32::from_rgb(255, 176, 32);
pub const ERROR: Color32 = Color32::from_rgb(255, 77, 79);
pub const OK: Color32 = Color32::from_rgb(94, 219, 138);

// -- metrics ---------------------------------------------------------------
pub const HAIR: f32 = 1.0;
pub const PAD: f32 = 12.0;
pub const RADIUS: u8 = 2;

/// Uppercase monospaced label — the typographic base of the whole UI.
pub fn mono(text: impl Into<String>) -> egui::RichText {
    egui::RichText::new(text.into().to_uppercase())
        .monospace()
        .color(TEXT)
}

pub fn label(text: impl Into<String>) -> egui::RichText {
    egui::RichText::new(text.into().to_uppercase())
        .monospace()
        .color(MUTED)
        .size(11.0)
}

pub fn value(text: impl Into<String>) -> egui::RichText {
    egui::RichText::new(text.into())
        .monospace()
        .color(TEXT)
        .size(12.0)
}

/// Section header: `▌ LABEL ──────────────────`
pub fn section(ui: &mut egui::Ui, title: &str) {
    ui.add_space(10.0);
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(ui.available_width(), 16.0),
        egui::Sense::hover(),
    );
    let p = ui.painter();
    p.rect_filled(
        egui::Rect::from_min_size(rect.min, egui::vec2(3.0, 10.0)),
        CornerRadius::same(1),
        ACCENT,
    );
    p.text(
        rect.min + egui::vec2(10.0, 1.0),
        egui::Align2::LEFT_TOP,
        title.to_uppercase(),
        egui::FontId::monospace(10.5),
        MUTED,
    );
    let text_w = title.chars().count() as f32 * 6.6 + 18.0;
    p.line_segment(
        [
            egui::pos2(rect.min.x + text_w, rect.center().y),
            egui::pos2(rect.max.x, rect.center().y),
        ],
        Stroke::new(HAIR, BORDER),
    );
    ui.add_space(4.0);
}

/// `KEY .............. VALUE` row.
pub fn kv(ui: &mut egui::Ui, key: &str, value_text: &str, color: Color32) {
    ui.add_space(1.0);
    ui.horizontal(|ui| {
        ui.label(label(key));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(
                egui::RichText::new(value_text)
                    .monospace()
                    .color(color)
                    .size(11.5),
            );
        });
    });
}

/// Full-width hairline.
pub fn rule(ui: &mut egui::Ui) {
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(ui.available_width(), HAIR),
        egui::Sense::hover(),
    );
    ui.painter().hline(rect.x_range(), rect.center().y, Stroke::new(HAIR, BORDER));
}

/// Status pill: a coloured dot plus uppercase text.
pub fn pill(ui: &mut egui::Ui, text: &str, color: Color32) {
    let galley = ui.painter().layout_no_wrap(
        text.to_uppercase(),
        egui::FontId::monospace(10.5),
        color,
    );
    let size = egui::vec2(galley.size().x + 22.0, 18.0);
    let (rect, _) = ui.allocate_exact_size(size, egui::Sense::hover());
    let p = ui.painter();
    p.rect_filled(rect, CornerRadius::same(2), PANEL_2);
    p.rect_stroke(
        rect,
        CornerRadius::same(2),
        Stroke::new(HAIR, BORDER),
        egui::StrokeKind::Inside,
    );
    p.circle_filled(egui::pos2(rect.min.x + 9.0, rect.center().y), 3.0, color);
    p.galley(
        egui::pos2(rect.min.x + 16.0, rect.center().y - galley.size().y / 2.0),
        galley,
        color,
    );
}

/// Corner brackets around a rect — the HUD framing used on the canvas.
pub fn brackets(painter: &egui::Painter, rect: egui::Rect, len: f32, color: Color32) {
    let s = Stroke::new(HAIR, color);
    let c = [
        (rect.left_top(), egui::vec2(1.0, 0.0), egui::vec2(0.0, 1.0)),
        (rect.right_top(), egui::vec2(-1.0, 0.0), egui::vec2(0.0, 1.0)),
        (rect.left_bottom(), egui::vec2(1.0, 0.0), egui::vec2(0.0, -1.0)),
        (rect.right_bottom(), egui::vec2(-1.0, 0.0), egui::vec2(0.0, -1.0)),
    ];
    for (corner, dx, dy) in c {
        painter.line_segment([corner, corner + dx * len], s);
        painter.line_segment([corner, corner + dy * len], s);
    }
}

/// Faint crosshair grid, drawn on the canvas substrate.
pub fn grid(painter: &egui::Painter, rect: egui::Rect, step: f32, color: Color32) {
    let s = Stroke::new(0.5, color);
    let mut x = rect.left();
    while x < rect.right() {
        painter.vline(x, rect.y_range(), s);
        x += step;
    }
    let mut y = rect.top();
    while y < rect.bottom() {
        painter.hline(rect.x_range(), y, s);
        y += step;
    }
}

/// Panel frame: 1px border, dark fill, no rounding.
pub fn panel_frame(fill: Color32) -> egui::Frame {
    egui::Frame::new()
        .fill(fill)
        .stroke(Stroke::new(HAIR, BORDER))
        .inner_margin(egui::Margin::symmetric(PAD as i8, 8))
}

pub fn apply(ctx: &egui::Context) {
    let mut v = egui::Visuals::dark();
    v.dark_mode = true;
    v.panel_fill = PANEL;
    v.window_fill = PANEL;
    v.extreme_bg_color = VOID;
    v.faint_bg_color = PANEL_2;
    v.code_bg_color = PANEL_2;
    v.text_edit_bg_color = Some(BG);
    v.override_text_color = Some(TEXT);
    v.weak_text_color = Some(MUTED);
    v.hyperlink_color = ACCENT_BRIGHT;
    v.warn_fg_color = WARN;
    v.error_fg_color = ERROR;
    v.selection.bg_fill = ACCENT_DIM;
    v.selection.stroke = Stroke::new(HAIR, ACCENT_BRIGHT);
    v.window_stroke = Stroke::new(HAIR, BORDER_HI);
    v.window_corner_radius = CornerRadius::same(RADIUS);
    v.menu_corner_radius = CornerRadius::same(RADIUS);

    let r = CornerRadius::same(RADIUS);
    v.widgets.noninteractive.bg_fill = PANEL_2;
    v.widgets.noninteractive.bg_stroke = Stroke::new(HAIR, BORDER);
    v.widgets.noninteractive.fg_stroke = Stroke::new(HAIR, TEXT);
    v.widgets.noninteractive.corner_radius = r;

    v.widgets.inactive.bg_fill = PANEL_2;
    v.widgets.inactive.weak_bg_fill = PANEL_2;
    v.widgets.inactive.bg_stroke = Stroke::new(HAIR, BORDER_HI);
    v.widgets.inactive.fg_stroke = Stroke::new(HAIR, TEXT);
    v.widgets.inactive.corner_radius = r;

    v.widgets.hovered.bg_fill = PANEL_3;
    v.widgets.hovered.weak_bg_fill = PANEL_3;
    v.widgets.hovered.bg_stroke = Stroke::new(HAIR, ACCENT);
    v.widgets.hovered.fg_stroke = Stroke::new(HAIR, ACCENT_BRIGHT);
    v.widgets.hovered.corner_radius = r;

    v.widgets.active.bg_fill = ACCENT_DIM;
    v.widgets.active.weak_bg_fill = ACCENT_DIM;
    v.widgets.active.bg_stroke = Stroke::new(HAIR, ACCENT);
    v.widgets.active.fg_stroke = Stroke::new(HAIR, ACCENT_BRIGHT);
    v.widgets.active.corner_radius = r;

    v.widgets.open.bg_fill = PANEL_2;
    v.widgets.open.bg_stroke = Stroke::new(HAIR, BORDER_HI);
    v.widgets.open.corner_radius = r;

    ctx.set_visuals(v);

    // Monospace everywhere: the surface is data.
    ctx.all_styles_mut(|style| {
        for (_, font) in style.text_styles.iter_mut() {
            font.family = egui::FontFamily::Monospace;
        }
        style.text_styles = [
            (egui::TextStyle::Heading, egui::FontId::monospace(17.0)),
            (egui::TextStyle::Body, egui::FontId::monospace(13.0)),
            (egui::TextStyle::Button, egui::FontId::monospace(12.5)),
            (egui::TextStyle::Small, egui::FontId::monospace(11.0)),
            (egui::TextStyle::Monospace, egui::FontId::monospace(12.5)),
        ].into();
        style.spacing.item_spacing = egui::vec2(7.0, 6.0);
        style.spacing.button_padding = egui::vec2(10.0, 5.0);
        style.spacing.interact_size.y = 22.0;
        // Instrument-style: text truncates inside its column instead of drawing
        // across the neighbouring panel.
        style.wrap_mode = Some(egui::TextWrapMode::Truncate);
    });
}
