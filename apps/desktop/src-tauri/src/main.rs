// Prevents an extra console window on Windows in release builds.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod arrow_out;
mod commands;
mod state;

use state::AppState;

fn main() {
    tauri::Builder::default()
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
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
