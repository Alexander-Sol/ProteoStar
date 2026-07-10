// Prevents an extra console window on Windows in release builds.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod arrow_out;
mod commands;
mod state;

use state::AppState;

fn main() {
    tauri::Builder::default()
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
            commands::get_ms2_for_precursor,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
