//! Tauri command surface (IPC contract §1/§6). M0 returns REAL Arrow with a
//! FAKE (synthetic) source — every shape matches the contract; the data is
//! stubbed. M1 replaces the command bodies with flashlfq-core calls (see §7).

use std::collections::HashMap;
use std::sync::atomic::Ordering;

use serde::Serialize;
use tauri::ipc::{Channel, Response};
use tauri::State;

use crate::arrow_out;
use crate::state::{AppState, StubDataset};

// Stub dataset dimensions (contract §-M0 task spec).
const SCAN_COUNT: usize = 500;
const RT_STEP: f32 = 0.1;

// ------------------------------------------------------------------ DTOs (JSON)

#[derive(Clone, Serialize)]
pub struct NumericRange {
    pub min: f64,
    pub max: f64,
}

#[derive(Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct DatasetMetadata {
    pub file_name: String,
    pub format: String,
    pub scan_count: u32,
    pub ms1_scan_count: u32,
    pub ms_levels_present: Vec<u32>,
    pub retention_time_range: Option<NumericRange>,
    pub mz_range: Option<NumericRange>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Precursor {
    pub mz: f64,
    pub scan_index: i32,
    pub isolation_low: f64,
    pub isolation_high: f64,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ScanSummary {
    pub scan_index: u32,
    pub one_based_scan_number: u32,
    pub retention_time: f64,
    pub tic: f64,
    pub ms_level: u32,
    pub precursor: Option<Precursor>,
}

#[derive(Serialize)]
pub struct OpenResult {
    pub handle: u64,
    pub metadata: DatasetMetadata,
}

#[derive(Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ProgressEvent {
    pub phase: String,
    pub scans_done: u32,
    pub scans_total: u32,
}

/// Rejected-command payload (contract §2 ViewerError).
#[derive(Serialize)]
pub struct ViewerError {
    pub code: String,
    pub message: String,
}

impl ViewerError {
    fn handle_not_found(handle: u64) -> Self {
        ViewerError {
            code: "HANDLE_NOT_FOUND".into(),
            message: format!("no open dataset for handle {handle}"),
        }
    }
    fn internal(message: impl Into<String>) -> Self {
        ViewerError {
            code: "INTERNAL".into(),
            message: message.into(),
        }
    }
    fn scan_out_of_range(scan_index: u32) -> Self {
        ViewerError {
            code: "SCAN_OUT_OF_RANGE".into(),
            message: format!("scan index {scan_index} exceeds {SCAN_COUNT}"),
        }
    }
}

// ------------------------------------------------------------- synthetic source

/// A visibly peak-shaped chromatogram: low baseline + a few Gaussian bumps.
fn synthetic_intensity(scan_index: usize) -> f32 {
    let rt = scan_index as f32 * RT_STEP;
    let bump = |center: f32, height: f32, width: f32| {
        let d = rt - center;
        height * (-(d * d) / (2.0 * width * width)).exp()
    };
    500.0 + bump(15.0, 8000.0, 1.5) + bump(30.0, 12000.0, 2.0) + bump(40.0, 5000.0, 1.0)
}

fn basename(path: &str) -> String {
    path.rsplit(['/', '\\'])
        .next()
        .filter(|s| !s.is_empty())
        .unwrap_or(path)
        .to_string()
}

fn stub_metadata(path: &str) -> DatasetMetadata {
    DatasetMetadata {
        file_name: basename(path),
        format: "mzml".into(),
        scan_count: SCAN_COUNT as u32,
        ms1_scan_count: SCAN_COUNT as u32,
        ms_levels_present: vec![1],
        retention_time_range: Some(NumericRange { min: 0.0, max: 50.0 }),
        mz_range: Some(NumericRange {
            min: 200.0,
            max: 2000.0,
        }),
    }
}

// ------------------------------------------------------------------- commands

#[tauri::command]
pub async fn open_dataset(
    path: String,
    on_progress: Channel<ProgressEvent>,
    state: State<'_, AppState>,
) -> Result<OpenResult, ViewerError> {
    // Emit a couple of progress events so the UI wiring is exercised end-to-end.
    let _ = on_progress.send(ProgressEvent {
        phase: "reading".into(),
        scans_done: 0,
        scans_total: SCAN_COUNT as u32,
    });
    let _ = on_progress.send(ProgressEvent {
        phase: "indexing".into(),
        scans_done: SCAN_COUNT as u32,
        scans_total: SCAN_COUNT as u32,
    });

    let metadata = stub_metadata(&path);
    let handle = state.next.fetch_add(1, Ordering::SeqCst) + 1;
    {
        let mut datasets = state
            .datasets
            .lock()
            .map_err(|e| ViewerError::internal(e.to_string()))?;
        datasets.insert(
            handle,
            StubDataset {
                metadata: metadata.clone(),
            },
        );
    }

    let _ = on_progress.send(ProgressEvent {
        phase: "done".into(),
        scans_done: SCAN_COUNT as u32,
        scans_total: SCAN_COUNT as u32,
    });

    Ok(OpenResult { handle, metadata })
}

#[tauri::command]
pub async fn close_dataset(handle: u64, state: State<'_, AppState>) -> Result<(), ViewerError> {
    let mut datasets = state
        .datasets
        .lock()
        .map_err(|e| ViewerError::internal(e.to_string()))?;
    datasets.remove(&handle);
    Ok(())
}

#[tauri::command]
pub async fn get_metadata(
    handle: u64,
    state: State<'_, AppState>,
) -> Result<DatasetMetadata, ViewerError> {
    let datasets = state
        .datasets
        .lock()
        .map_err(|e| ViewerError::internal(e.to_string()))?;
    datasets
        .get(&handle)
        .map(|d| d.metadata.clone())
        .ok_or_else(|| ViewerError::handle_not_found(handle))
}

#[tauri::command]
pub async fn get_scan_summaries(
    handle: u64,
    state: State<'_, AppState>,
) -> Result<Response, ViewerError> {
    ensure_handle(&state, handle)?;

    let n = SCAN_COUNT;
    let one_based: Vec<u32> = (1..=n as u32).collect();
    let rt: Vec<f32> = (0..n).map(|i| i as f32 * RT_STEP).collect();
    let tic: Vec<f32> = (0..n).map(synthetic_intensity).collect();
    let ms_level: Vec<u8> = vec![1u8; n];
    let null_f64: Vec<Option<f64>> = vec![None; n];
    let null_i32: Vec<Option<i32>> = vec![None; n];

    let bytes = arrow_out::build_scan_summaries(
        one_based,
        rt,
        tic,
        ms_level,
        null_f64.clone(),
        null_i32,
        null_f64.clone(),
        null_f64,
    )
    .map_err(|e| ViewerError::internal(e.to_string()))?;
    Ok(Response::new(bytes))
}

#[tauri::command]
pub async fn get_nearest_scan(
    handle: u64,
    retention_time: f64,
    ms_level: Option<u32>,
    state: State<'_, AppState>,
) -> Result<Option<ScanSummary>, ViewerError> {
    ensure_handle(&state, handle)?;
    let _ = ms_level; // stub is MS1-only

    let idx = ((retention_time / RT_STEP as f64).round() as i64)
        .clamp(0, (SCAN_COUNT - 1) as i64) as u32;
    Ok(Some(ScanSummary {
        scan_index: idx,
        one_based_scan_number: idx + 1,
        retention_time: idx as f64 * RT_STEP as f64,
        tic: synthetic_intensity(idx as usize) as f64,
        ms_level: 1,
        precursor: None,
    }))
}

#[tauri::command]
pub async fn get_tic_trace(
    handle: u64,
    rt_min: Option<f64>,
    rt_max: Option<f64>,
    ms_level: u32,
    max_points: Option<u32>,
    state: State<'_, AppState>,
) -> Result<Response, ViewerError> {
    ensure_handle(&state, handle)?;
    let _ = (rt_min, rt_max, ms_level, max_points); // stub ignores reduction knobs

    let n = SCAN_COUNT;
    let rt: Vec<f32> = (0..n).map(|i| i as f32 * RT_STEP).collect();
    let intensity: Vec<f32> = (0..n).map(synthetic_intensity).collect();
    let scan_index: Vec<u32> = (0..n as u32).collect();

    let bytes = arrow_out::build_trace(rt, intensity, scan_index)
        .map_err(|e| ViewerError::internal(e.to_string()))?;
    Ok(Response::new(bytes))
}

#[tauri::command]
pub async fn get_range_xic(
    handle: u64,
    mz_low: f64,
    mz_high: f64,
    rt_min: Option<f64>,
    rt_max: Option<f64>,
    max_points: Option<u32>,
    state: State<'_, AppState>,
) -> Result<Response, ViewerError> {
    ensure_handle(&state, handle)?;
    let _ = (rt_min, rt_max, max_points);

    // Synthetic window-XIC: scale the baseline chromatogram by the window width
    // so the stub reacts to its m/z arguments without a real index.
    let width = (mz_high - mz_low).max(0.0) as f32;
    let scale = (width / 100.0).clamp(0.05, 1.0);
    let n = SCAN_COUNT;
    let rt: Vec<f32> = (0..n).map(|i| i as f32 * RT_STEP).collect();
    let intensity: Vec<f32> = (0..n).map(|i| synthetic_intensity(i) * scale).collect();
    let scan_index: Vec<u32> = (0..n as u32).collect();

    let bytes = arrow_out::build_trace(rt, intensity, scan_index)
        .map_err(|e| ViewerError::internal(e.to_string()))?;
    Ok(Response::new(bytes))
}

#[tauri::command]
pub async fn get_spectrum(
    handle: u64,
    scan_index: u32,
    mz_min: Option<f64>,
    mz_max: Option<f64>,
    max_peaks: Option<u32>,
    state: State<'_, AppState>,
) -> Result<Response, ViewerError> {
    ensure_handle(&state, handle)?;
    let _ = (mz_min, mz_max, max_peaks);
    if scan_index as usize >= SCAN_COUNT {
        return Err(ViewerError::scan_out_of_range(scan_index));
    }

    // ~20 canned peaks across the m/z range, envelope keyed to the scan's RT.
    let rt = scan_index as f32 * RT_STEP;
    let env = synthetic_intensity(scan_index as usize);
    let mut mz = Vec::with_capacity(20);
    let mut intensity = Vec::with_capacity(20);
    for k in 0..20u32 {
        let m = 250.0 + k as f64 * 85.0; // 250 .. 1865
        let phase = (k as f32 * 0.7).sin() * 0.5 + 0.5;
        mz.push(m);
        intensity.push(env * (0.15 + 0.85 * phase));
    }

    let mut metadata = HashMap::new();
    metadata.insert("scanIndex".into(), scan_index.to_string());
    metadata.insert("oneBasedScanNumber".into(), (scan_index + 1).to_string());
    metadata.insert("retentionTime".into(), rt.to_string());
    metadata.insert("msLevel".into(), "1".to_string());

    let bytes = arrow_out::build_spectrum(mz, intensity, metadata)
        .map_err(|e| ViewerError::internal(e.to_string()))?;
    Ok(Response::new(bytes))
}

#[tauri::command]
pub async fn get_ms2_for_precursor(
    handle: u64,
    mz: f64,
    ms1_scan_index: u32,
    state: State<'_, AppState>,
) -> Result<Vec<ScanSummary>, ViewerError> {
    ensure_handle(&state, handle)?;
    let _ = (mz, ms1_scan_index); // MS1-only stub: no MS2 scans
    Ok(Vec::new())
}

// --------------------------------------------------------------------- helpers

fn ensure_handle(state: &State<'_, AppState>, handle: u64) -> Result<(), ViewerError> {
    let datasets = state
        .datasets
        .lock()
        .map_err(|e| ViewerError::internal(e.to_string()))?;
    if datasets.contains_key(&handle) {
        Ok(())
    } else {
        Err(ViewerError::handle_not_found(handle))
    }
}
