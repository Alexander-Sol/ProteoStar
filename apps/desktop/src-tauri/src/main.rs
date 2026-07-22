// Prevents an extra console window on Windows in release builds.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod arrow_out;
mod commands;
mod state;

use state::AppState;

fn main() {
    let builder = tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .manage(AppState::default())
        .invoke_handler(tauri::generate_handler![
            commands::open_dataset,
            commands::close_dataset,
            commands::get_metadata,
            commands::get_scan_summaries,
            commands::get_nearest_scan,
            commands::get_tic_trace,
            commands::get_range_xic,
            commands::get_spectrum,
            commands::get_spectrum_at_rt,
            commands::get_ms2_spectrum,
            commands::get_ms2_for_precursor,
            commands::load_features,
            commands::load_psms,
            commands::run_feature_detection,
            commands::score_seed_ladder,
        ]);

    // Playwright E2E: mount the control server and grant its `pw_result` IPC the
    // permission it needs. Both are gated behind the `e2e-testing` feature so a
    // normal build never links the plugin or references the `playwright:default`
    // permission (which only exists when the plugin is compiled in). The
    // capability is added at runtime — the file lives outside capabilities/ so
    // the always-on ACL build never tries to resolve it. See e2e/README.md.
    #[cfg(feature = "e2e-testing")]
    let builder = builder
        .plugin(tauri_plugin_playwright::init_with_config(
            // Explicit port so Windows (TCP) and the Playwright fixture agree; on
            // unix the plugin still prefers its default socket.
            tauri_plugin_playwright::PluginConfig::new().tcp_port(6274),
        ))
        .setup(|app| {
            use tauri::Manager;
            app.add_capability(include_str!("../e2e-capabilities/playwright.json"))?;
            Ok(())
        });

    builder
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
