//! Tauri command surface (IPC contract §1/§6).
//!
//! M1+: real data. `open_dataset` reads the MS1 scans via flashlfq-core
//! (`read_ms1_scans` → `PeakIndexingEngine`) and every query is served from the
//! in-memory `Vec<Scan>`: the TIC is per-scan summed centroided intensity,
//! spectra are the stored MS1 peak lists, and the range-XIC sums peaks in an m/z
//! window per scan. MS2 is not yet retained (the vendored reader is MS1-only), so
//! `get_ms2_for_precursor` returns empty and every scan reports `msLevel = 1`.

use std::collections::HashMap;
use std::sync::atomic::Ordering;

use std::sync::Arc;

use serde::Serialize;
use tauri::ipc::{Channel, Response};
use tauri::{Emitter, State};

use flashlfq_core::feature_refinement::{
    self, refine_feature_multi, resolve_consensus_by_apex, RefinedFeature, ResolvedFeature,
};
use flashlfq_core::peak_indexing::{
    read_ms1_scans, read_ms1_tic_metadata, PeakIndexingEngine, PeakSource, RandomAccessMs1Reader,
    Scan,
};
use flashlfq_core::spectral_averaging::SpectralAveragingParameters;
use flashlfq_core::trace_kernel::{
    detect_features, median_ms1_scan_spacing_minutes, seed_ladder_diagnostics,
    SeedLadderDiagnostics, TraceKernelParameters, FWHM_TO_SIGMA,
};

use crate::arrow_out;
use crate::state::{AppState, IndexedData, OpenDataset};

/// Proton mass (Da) for neutral-mass ↔ m/z conversions.
const PROTON_MASS: f64 = 1.007_276_466_8;

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
    /// `Some(message)` if the background index build failed. The frontend polls this so a failed
    /// index surfaces an error instead of hanging on "indexing…" forever. `None` while indexing
    /// or once indexed successfully (a non-zero `ms1_scan_count` is the success signal).
    pub index_error: Option<String>,
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
#[derive(Debug, Serialize)]
pub struct ViewerError {
    pub code: String,
    pub message: String,
}

impl ViewerError {
    fn new(code: &str, message: impl Into<String>) -> Self {
        ViewerError { code: code.into(), message: message.into() }
    }
    fn handle_not_found(handle: u64) -> Self {
        ViewerError::new("HANDLE_NOT_FOUND", format!("no open dataset for handle {handle}"))
    }
    fn internal(message: impl Into<String>) -> Self {
        ViewerError::new("INTERNAL", message)
    }
    fn scan_out_of_range(scan_index: u32, count: usize) -> Self {
        ViewerError::new(
            "SCAN_OUT_OF_RANGE",
            format!("scan index {scan_index} exceeds {count}"),
        )
    }
    /// The dataset is open and its TIC is available, but the background peak index
    /// isn't built yet. Callers that need per-scan peaks get this until it lands.
    fn indexing() -> Self {
        ViewerError::new("INDEXING", "dataset is still being indexed")
    }
}

// ---------------------------------------------------------------- dataset build

fn basename(path: &str) -> String {
    path.rsplit(['/', '\\'])
        .next()
        .filter(|s| !s.is_empty())
        .unwrap_or(path)
        .to_string()
}

/// Fast pass: open the file for on-demand reads and build a **smooth MS1-only TIC** without the
/// full peak index. The displayed TIC is read from per-scan `total ion current` metadata
/// ([`read_ms1_tic_metadata`] — no peak decode, no index), so it appears quickly and is MS1-only
/// (not the jagged native TIC, which interleaves MS2). Falls back to the native instrument TIC
/// when the file exposes no per-scan TIC metadata; an empty pair then means the served TIC fills
/// in from the background index instead. Runs on a blocking thread (`.raw` pulls a .NET runtime).
fn open_fast(path: &str) -> Result<(RandomAccessMs1Reader, Vec<f64>, Vec<f32>), ViewerError> {
    let mut reader =
        RandomAccessMs1Reader::open(path).map_err(|e| ViewerError::new("READ_ERROR", e.to_string()))?;
    let (rt, intensity) = match read_ms1_tic_metadata(path) {
        Ok(tic) if !tic.retention_times.is_empty() => (tic.retention_times, tic.intensities),
        // No per-scan TIC metadata: fall back to the native instrument TIC (may be jagged), or an
        // empty pair (served TIC then comes from the background MS1 index once it lands).
        _ => match reader.tic() {
            Some(tic) => (tic.retention_times, tic.intensities),
            None => (Vec::new(), Vec::new()),
        },
    };
    Ok((reader, rt, intensity))
}

/// Background pass: read every MS1 scan, build the peak index, and compute the per-scan
/// summed-centroid TIC. This is the heavy work the fast path defers. Runs on a blocking
/// thread (both the read and the index build are CPU-bound; `.raw` pulls a .NET runtime).
fn build_indexed(path: &str) -> Result<IndexedData, ViewerError> {
    let scans = read_ms1_scans(path).map_err(|e| ViewerError::new("READ_ERROR", e.to_string()))?;
    if scans.is_empty() {
        return Err(ViewerError::new("EMPTY_INDEX", "no MS1 scans found in file"));
    }
    let engine = PeakIndexingEngine::index_peaks(&scans)
        .ok_or_else(|| ViewerError::new("EMPTY_INDEX", "no indexable MS1 peaks"))?;

    // Per-scan TIC = summed centroided MS1 intensity. Used by scan summaries + nearest-scan;
    // the displayed TIC chromatogram is the native one read on the fast path.
    let tic: Vec<f32> = scans
        .iter()
        .map(|s| s.intensity.iter().sum::<f64>() as f32)
        .collect();

    Ok(IndexedData {
        scans: Arc::new(scans),
        tic: Arc::new(tic),
        engine: Arc::new(engine),
    })
}

/// Provisional metadata available the moment the fast TIC is read: file name, format, and
/// the RT range spanned by the native TIC. The MS1 scan count and m/z range are unknown
/// until the index is built, so they start empty and [`refine_metadata`] fills them in.
fn provisional_metadata(path: &str, tic_rt: &[f64]) -> DatasetMetadata {
    let (rt_min, rt_max) = tic_rt.iter().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(lo, hi), &t| (lo.min(t), hi.max(t)),
    );
    let format = if path.to_lowercase().ends_with(".raw") {
        "thermo_raw"
    } else {
        "mzml"
    };
    DatasetMetadata {
        file_name: basename(path),
        format: format.into(),
        scan_count: tic_rt.len() as u32, // all scans (provisional); MS1 count filled on index
        ms1_scan_count: 0,               // 0 marks "not indexed yet" to the frontend
        ms_levels_present: vec![1],
        retention_time_range: rt_min
            .is_finite()
            .then_some(NumericRange { min: rt_min, max: rt_max }),
        mz_range: None,
        index_error: None,
    }
}

/// Refine provisional metadata once the MS1 scans are read: exact RT range, m/z range, and
/// the real MS1 scan count (a non-zero `ms1_scan_count` is the frontend's "indexed" signal).
fn refine_metadata(meta: &mut DatasetMetadata, scans: &[Scan]) {
    let (mut rt_min, mut rt_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut mz_min, mut mz_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for s in scans {
        rt_min = rt_min.min(s.retention_time);
        rt_max = rt_max.max(s.retention_time);
        if let Some(&m) = s.mz.first() {
            mz_min = mz_min.min(m); // mz is ascending
        }
        if let Some(&m) = s.mz.last() {
            mz_max = mz_max.max(m);
        }
    }
    meta.ms1_scan_count = scans.len() as u32;
    if rt_min.is_finite() {
        meta.retention_time_range = Some(NumericRange { min: rt_min, max: rt_max });
    }
    if mz_min.is_finite() {
        meta.mz_range = Some(NumericRange { min: mz_min, max: mz_max });
    }
}

/// Clone the indexed read handles for a dataset, or return an `INDEXING` error if the
/// background index build hasn't finished. All three are `Arc`s, so this is cheap and the
/// caller can drop the state lock (which this releases on return) before doing real work.
fn indexed_handles(
    state: &State<'_, AppState>,
    handle: u64,
) -> Result<(Arc<Vec<Scan>>, Arc<Vec<f32>>, Arc<PeakIndexingEngine>), ViewerError> {
    let datasets = state.datasets.lock().map_err(|e| ViewerError::internal(e.to_string()))?;
    let d = datasets.get(&handle).ok_or_else(|| ViewerError::handle_not_found(handle))?;
    let guard = d.indexed.lock().map_err(|e| ViewerError::internal(e.to_string()))?;
    let idx = guard.as_ref().ok_or_else(ViewerError::indexing)?;
    Ok((idx.scans.clone(), idx.tic.clone(), idx.engine.clone()))
}

// --------------------------------------------------------------------- helpers

fn ensure_handle(state: &State<'_, AppState>, handle: u64) -> Result<(), ViewerError> {
    let datasets = state.datasets.lock().map_err(|e| ViewerError::internal(e.to_string()))?;
    if datasets.contains_key(&handle) {
        Ok(())
    } else {
        Err(ViewerError::handle_not_found(handle))
    }
}

/// Stride-decimate parallel arrays to at most `max_points`, always keeping the
/// last point so the chromatogram spans its full RT range.
fn decimate(
    rt: Vec<f32>,
    intensity: Vec<f32>,
    scan_index: Vec<u32>,
    max_points: Option<u32>,
) -> (Vec<f32>, Vec<f32>, Vec<u32>) {
    let n = rt.len();
    let cap = max_points.map(|m| m as usize).unwrap_or(usize::MAX);
    if cap == 0 || n <= cap {
        return (rt, intensity, scan_index);
    }
    let stride = n.div_ceil(cap);
    let mut or = Vec::with_capacity(cap + 1);
    let mut oi = Vec::with_capacity(cap + 1);
    let mut os = Vec::with_capacity(cap + 1);
    let mut i = 0;
    while i < n {
        or.push(rt[i]);
        oi.push(intensity[i]);
        os.push(scan_index[i]);
        i += stride;
    }
    if *os.last().unwrap() as usize != n - 1 {
        or.push(rt[n - 1]);
        oi.push(intensity[n - 1]);
        os.push(scan_index[n - 1]);
    }
    (or, oi, os)
}

/// True when a scan's RT is inside the optional `[rt_min, rt_max]` window.
fn in_rt_window(rt: f64, rt_min: Option<f64>, rt_max: Option<f64>) -> bool {
    rt_min.map_or(true, |lo| rt >= lo) && rt_max.map_or(true, |hi| rt <= hi)
}

// ------------------------------------------------------------------- commands

#[tauri::command]
pub async fn open_dataset(
    path: String,
    on_progress: Channel<ProgressEvent>,
    app: tauri::AppHandle,
    state: State<'_, AppState>,
) -> Result<OpenResult, ViewerError> {
    if !std::path::Path::new(&path).exists() {
        return Err(ViewerError::new("FILE_NOT_FOUND", format!("no such file: {path}")));
    }
    let _ = on_progress.send(ProgressEvent {
        phase: "reading".into(),
        scans_done: 0,
        scans_total: 0,
    });

    // Phase 1 — fast: open the reader + native TIC only, so the viewer can paint the TIC and
    // serve spectra by RT immediately (no peak decode / index build yet).
    let path_for_open = path.clone();
    let (reader, tic_rt, tic_intensity) =
        tauri::async_runtime::spawn_blocking(move || open_fast(&path_for_open))
            .await
            .map_err(|e| ViewerError::internal(format!("open task failed: {e}")))??;

    let metadata = provisional_metadata(&path, &tic_rt);
    let handle = state.next.fetch_add(1, Ordering::SeqCst) + 1;
    let dataset = OpenDataset {
        metadata: Arc::new(std::sync::Mutex::new(metadata.clone())),
        tic_rt,
        tic_intensity,
        reader: Arc::new(std::sync::Mutex::new(reader)),
        indexed: Arc::new(std::sync::Mutex::new(None)),
    };
    let indexed_slot = dataset.indexed.clone();
    let metadata_slot = dataset.metadata.clone();
    {
        let mut datasets = state.datasets.lock().map_err(|e| ViewerError::internal(e.to_string()))?;
        datasets.insert(handle, dataset);
    }

    // Phase 2 — background: read scans + build the peak index off the open path, then fill
    // the indexed slot, refine the metadata, and signal readiness. Commands that need
    // per-scan peaks return INDEXING until this completes.
    let path_for_index = path.clone();
    tauri::async_runtime::spawn(async move {
        let built =
            tauri::async_runtime::spawn_blocking(move || build_indexed(&path_for_index)).await;
        match built {
            Ok(Ok(indexed)) => {
                if let Ok(mut m) = metadata_slot.lock() {
                    refine_metadata(&mut m, &indexed.scans);
                }
                if let Ok(mut slot) = indexed_slot.lock() {
                    *slot = Some(indexed);
                }
                let _ = app.emit("dataset-indexed", handle);
            }
            Ok(Err(e)) => {
                // Record the failure in metadata so the frontend's poll surfaces it (rather
                // than hanging on "indexing…"), and also emit for any event-based listener.
                if let Ok(mut m) = metadata_slot.lock() {
                    m.index_error = Some(e.message.clone());
                }
                let _ = app.emit("dataset-index-error", (handle, e.message));
            }
            Err(e) => {
                let msg = format!("index task failed: {e}");
                if let Ok(mut m) = metadata_slot.lock() {
                    m.index_error = Some(msg.clone());
                }
                let _ = app.emit("dataset-index-error", (handle, msg));
            }
        }
    });

    let _ = on_progress.send(ProgressEvent {
        phase: "ready".into(),
        scans_done: 0,
        scans_total: 0,
    });

    Ok(OpenResult { handle, metadata })
}

#[tauri::command]
pub async fn close_dataset(handle: u64, state: State<'_, AppState>) -> Result<(), ViewerError> {
    let mut datasets = state.datasets.lock().map_err(|e| ViewerError::internal(e.to_string()))?;
    datasets.remove(&handle);
    Ok(())
}

#[tauri::command]
pub async fn get_metadata(
    handle: u64,
    state: State<'_, AppState>,
) -> Result<DatasetMetadata, ViewerError> {
    let datasets = state.datasets.lock().map_err(|e| ViewerError::internal(e.to_string()))?;
    let d = datasets.get(&handle).ok_or_else(|| ViewerError::handle_not_found(handle))?;
    let meta = d.metadata.lock().map_err(|e| ViewerError::internal(e.to_string()))?;
    Ok(meta.clone())
}

#[tauri::command]
pub async fn get_scan_summaries(
    handle: u64,
    state: State<'_, AppState>,
) -> Result<Response, ViewerError> {
    let (scans, tic_arc, _) = indexed_handles(&state, handle)?;

    let n = scans.len();
    let one_based: Vec<u32> = scans.iter().map(|s| s.one_based_scan_number.max(0) as u32).collect();
    let rt: Vec<f32> = scans.iter().map(|s| s.retention_time as f32).collect();
    let tic = (*tic_arc).clone();
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
    let _ = ms_level; // MS1-only
    let (scans, tic, _) = indexed_handles(&state, handle)?;

    let mut best: Option<(usize, f64)> = None;
    for (i, s) in scans.iter().enumerate() {
        let dist = (s.retention_time - retention_time).abs();
        if best.map_or(true, |(_, b)| dist < b) {
            best = Some((i, dist));
        }
    }
    Ok(best.map(|(i, _)| {
        let s = &scans[i];
        ScanSummary {
            scan_index: i as u32,
            one_based_scan_number: s.one_based_scan_number.max(0) as u32,
            retention_time: s.retention_time,
            tic: tic[i] as f64,
            ms_level: 1,
            precursor: None,
        }
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
    let _ = ms_level;
    let datasets = state.datasets.lock().map_err(|e| ViewerError::internal(e.to_string()))?;
    let d = datasets.get(&handle).ok_or_else(|| ViewerError::handle_not_found(handle))?;

    let (mut rt, mut intensity, mut scan_index) = (Vec::new(), Vec::new(), Vec::new());
    if !d.tic_rt.is_empty() {
        // Preferred: the smooth **MS1-only** TIC read from per-scan `total ion current` metadata on
        // the fast path (`read_ms1_tic_metadata`) — available immediately, no peak index needed,
        // and MS1-only (so it doesn't need to wait on the full index and isn't the jagged native
        // TIC). `scanIndex` here is the TIC point's ordinal, not an MS1 scan index; TIC clicks
        // navigate by RT, which `get_nearest_scan` resolves to a real scan once indexing is done.
        for (i, (&t, &inten)) in d.tic_rt.iter().zip(d.tic_intensity.iter()).enumerate() {
            if in_rt_window(t, rt_min, rt_max) {
                rt.push(t as f32);
                intensity.push(inten);
                scan_index.push(i as u32);
            }
        }
    } else if let Some((scans, tic)) = d
        .indexed
        .lock()
        .ok()
        .and_then(|g| g.as_ref().map(|i| (i.scans.clone(), i.tic.clone())))
    {
        // Fallback (file exposed no per-scan TIC metadata): the per-scan summed-centroid MS1 TIC,
        // available only once the index is built. Here `scanIndex` *is* the MS1 scan index.
        for (i, s) in scans.iter().enumerate() {
            if in_rt_window(s.retention_time, rt_min, rt_max) {
                rt.push(s.retention_time as f32);
                intensity.push(tic[i]);
                scan_index.push(i as u32);
            }
        }
    }
    let (rt, intensity, scan_index) = decimate(rt, intensity, scan_index, max_points);
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
    let (scans, _, _) = indexed_handles(&state, handle)?;

    let (lo, hi) = if mz_low <= mz_high { (mz_low, mz_high) } else { (mz_high, mz_low) };
    let (mut rt, mut intensity, mut scan_index) = (Vec::new(), Vec::new(), Vec::new());
    for (i, s) in scans.iter().enumerate() {
        if !in_rt_window(s.retention_time, rt_min, rt_max) {
            continue;
        }
        // mz ascending: sum intensities of peaks within [lo, hi].
        let a = s.mz.partition_point(|&m| m < lo);
        let b = s.mz.partition_point(|&m| m <= hi);
        let sum: f64 = s.intensity[a..b].iter().sum();
        rt.push(s.retention_time as f32);
        intensity.push(sum as f32);
        scan_index.push(i as u32);
    }
    let (rt, intensity, scan_index) = decimate(rt, intensity, scan_index, max_points);
    let bytes = arrow_out::build_trace(rt, intensity, scan_index)
        .map_err(|e| ViewerError::internal(e.to_string()))?;
    Ok(Response::new(bytes))
}

/// Build the Arrow spectrum `Response` for one scan: apply the m/z window, keep the
/// `max_peaks` most intense (restoring m/z order), and attach the metadata columns. `scan_index`
/// is what to report as `scanIndex` (the MS1 ordinal for indexed reads, or the file spectrum
/// index for on-demand RT reads, which have no MS1 ordinal yet).
fn build_spectrum_response(
    s: &Scan,
    scan_index: i64,
    mz_min: Option<f64>,
    mz_max: Option<f64>,
    max_peaks: Option<u32>,
) -> Result<Response, ViewerError> {
    // m/z window (mz ascending).
    let a = mz_min.map_or(0, |lo| s.mz.partition_point(|&m| m < lo));
    let b = mz_max.map_or(s.mz.len(), |hi| s.mz.partition_point(|&m| m <= hi));
    let mut mz: Vec<f64> = s.mz[a..b].to_vec();
    let mut intensity: Vec<f32> = s.intensity[a..b].iter().map(|&v| v as f32).collect();

    // Level-of-detail: keep the `max_peaks` most intense, then restore m/z order.
    if let Some(cap) = max_peaks {
        let cap = cap as usize;
        if cap > 0 && mz.len() > cap {
            let mut idx: Vec<usize> = (0..mz.len()).collect();
            idx.sort_unstable_by(|&x, &y| intensity[y].total_cmp(&intensity[x]));
            idx.truncate(cap);
            idx.sort_unstable();
            mz = idx.iter().map(|&i| mz[i]).collect();
            intensity = idx.iter().map(|&i| intensity[i]).collect();
        }
    }

    let mut metadata = HashMap::new();
    metadata.insert("scanIndex".into(), scan_index.to_string());
    metadata.insert("oneBasedScanNumber".into(), s.one_based_scan_number.max(0).to_string());
    metadata.insert("retentionTime".into(), s.retention_time.to_string());
    // Real MS level (1 for the MS1 reads; 2+ for an MS2 scan pulled by scan number).
    metadata.insert("msLevel".into(), s.msn_order.max(1).to_string());

    let bytes = arrow_out::build_spectrum(mz, intensity, metadata)
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
    let (scans, _, _) = indexed_handles(&state, handle)?;

    let si = scan_index as usize;
    if si >= scans.len() {
        return Err(ViewerError::scan_out_of_range(scan_index, scans.len()));
    }
    build_spectrum_response(&scans[si], scan_index as i64, mz_min, mz_max, max_peaks)
}

/// Fetch the MS1 spectrum nearest a retention time via an on-demand single-scan read from the
/// kept-open reader — no peak index required, so it works during (and after) indexing. This is
/// the RT-first path the viewer uses for TIC clicks and feature selection.
#[tauri::command]
pub async fn get_spectrum_at_rt(
    handle: u64,
    retention_time: f64,
    mz_min: Option<f64>,
    mz_max: Option<f64>,
    max_peaks: Option<u32>,
    state: State<'_, AppState>,
) -> Result<Response, ViewerError> {
    let reader = {
        let datasets = state.datasets.lock().map_err(|e| ViewerError::internal(e.to_string()))?;
        let d = datasets.get(&handle).ok_or_else(|| ViewerError::handle_not_found(handle))?;
        d.reader.clone()
    };
    // Read off the async runtime: a single-scan `.raw` read crosses into the .NET runtime.
    let scan = tauri::async_runtime::spawn_blocking(move || {
        let mut r = reader.lock().map_err(|e| e.to_string())?;
        Ok::<Option<Scan>, String>(r.ms1_scan_at_rt(retention_time))
    })
    .await
    .map_err(|e| ViewerError::internal(format!("spectrum read task failed: {e}")))?
    .map_err(ViewerError::internal)?;

    let Some(s) = scan else {
        return Err(ViewerError::new("NO_SCAN", "no MS1 scan near that retention time"));
    };
    // No MS1 ordinal exists pre-index; report the file spectrum index (one_based − 1).
    let file_index = (s.one_based_scan_number as i64 - 1).max(-1);
    build_spectrum_response(&s, file_index, mz_min, mz_max, max_peaks)
}

/// Fetch a specific MS2 (or MSn) spectrum by its **one-based scan number** — the `Scan Number`
/// a MetaMorpheus PSM carries. Read on demand from the kept-open reader by file index
/// (`scan_number - 1`), so it needs no peak index (which is MS1-only anyway). This is how a
/// selected PSM shows its identified fragment spectrum.
#[tauri::command]
pub async fn get_ms2_spectrum(
    handle: u64,
    scan_number: i32,
    mz_min: Option<f64>,
    mz_max: Option<f64>,
    max_peaks: Option<u32>,
    state: State<'_, AppState>,
) -> Result<Response, ViewerError> {
    if scan_number < 1 {
        return Err(ViewerError::new("SCAN_OUT_OF_RANGE", "scan number must be >= 1"));
    }
    let reader = {
        let datasets = state.datasets.lock().map_err(|e| ViewerError::internal(e.to_string()))?;
        let d = datasets.get(&handle).ok_or_else(|| ViewerError::handle_not_found(handle))?;
        d.reader.clone()
    };
    // Read off the async runtime: a single-scan `.raw` read crosses into the .NET runtime.
    let scan = tauri::async_runtime::spawn_blocking(move || {
        let mut r = reader.lock().map_err(|e| e.to_string())?;
        Ok::<Option<Scan>, String>(r.scan_by_one_based_number(scan_number))
    })
    .await
    .map_err(|e| ViewerError::internal(format!("MS2 read task failed: {e}")))?
    .map_err(ViewerError::internal)?;

    let Some(s) = scan else {
        return Err(ViewerError::new("NO_SCAN", format!("no scan with number {scan_number}")));
    };
    build_spectrum_response(&s, scan_number as i64, mz_min, mz_max, max_peaks)
}

#[tauri::command]
pub async fn get_ms2_for_precursor(
    handle: u64,
    mz: f64,
    ms1_scan_index: u32,
    state: State<'_, AppState>,
) -> Result<Vec<ScanSummary>, ViewerError> {
    ensure_handle(&state, handle)?;
    let _ = (mz, ms1_scan_index); // MS1-only reader: no MS2 scans retained
    Ok(Vec::new())
}

// ------------------------------------------------------------- feature overlay

/// One observed isotope-peak m/z for a specific charge state of a resolved feature.
#[derive(Debug, Serialize, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PerChargeMz {
    pub charge: i32,
    pub mz: f64,
}

/// A resolved (peptide/proteoform-level) feature parsed from the runner's output
/// TSV — the unit the overlay draws and the feature list browses.
#[derive(Debug, Serialize, Clone)]
#[serde(rename_all = "camelCase")]
pub struct Feature {
    pub detected_mz: f64,
    pub rt_start: f64,
    pub rt_apex: f64,
    pub rt_end: f64,
    pub charge_states: Vec<i32>,
    pub per_charge_mz: Vec<PerChargeMz>,
    pub primary_charge: i32,
    pub monoisotopic_mass: f64,
    pub mono_mz: f64,
    pub summed_intensity: f64,
    pub cross_charge_support: u32,
    pub num_members: u32,
}

/// Parse the resolved-feature TSV written by `detect_features_tsv` (columns keyed
/// by header name, tolerant of reordering). Returns a `FEATURE_PARSE`-coded error
/// if the file is unreadable or lacks the expected header.
#[tauri::command]
pub async fn load_features(path: String) -> Result<Vec<Feature>, ViewerError> {
    let text = std::fs::read_to_string(&path)
        .map_err(|e| ViewerError::new("FEATURE_READ", format!("{path}: {e}")))?;
    parse_resolved_features(&text)
}

/// Parse the resolved-feature TSV body (header + rows). Pure so it can be unit-tested.
fn parse_resolved_features(text: &str) -> Result<Vec<Feature>, ViewerError> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| ViewerError::new("FEATURE_PARSE", "empty feature file"))?;
    let cols: HashMap<&str, usize> =
        header.split('\t').enumerate().map(|(i, c)| (c.trim(), i)).collect();

    // Required columns (resolved format).
    let need = |name: &str| -> Result<usize, ViewerError> {
        cols.get(name)
            .copied()
            .ok_or_else(|| ViewerError::new("FEATURE_PARSE", format!("missing column '{name}'")))
    };
    let c_detected_mz = need("Detected m/z (primary)")?;
    let c_rt_start = need("RT Start")?;
    let c_rt_apex = need("RT Apex")?;
    let c_rt_end = need("RT End")?;
    let c_charges = need("Charge States")?;
    let c_per_charge = need("Per-Charge Detected m/z")?;
    let c_primary = need("Primary Charge")?;
    let c_mass = need("Monoisotopic Mass")?;
    let c_mono_mz = need("Mono m/z (primary)")?;
    let c_intensity = need("Summed Intensity")?;
    let c_support = need("Cross-Charge Support")?;
    let c_members = need("Num Members")?;

    let fnum = |f: &[&str], i: usize| -> f64 { f.get(i).and_then(|s| s.trim().parse().ok()).unwrap_or(0.0) };
    let inum = |f: &[&str], i: usize| -> i32 { f.get(i).and_then(|s| s.trim().parse().ok()).unwrap_or(0) };
    let unum = |f: &[&str], i: usize| -> u32 { f.get(i).and_then(|s| s.trim().parse().ok()).unwrap_or(0) };

    let mut out = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let charge_states: Vec<i32> = f
            .get(c_charges)
            .map(|s| s.split(';').filter_map(|p| p.trim().parse().ok()).collect())
            .unwrap_or_default();
        // "z2:497.2584;z3:331.8416"
        let per_charge_mz: Vec<PerChargeMz> = f
            .get(c_per_charge)
            .map(|s| {
                s.split(';')
                    .filter_map(|p| {
                        let p = p.trim().trim_start_matches('z');
                        let (z, mz) = p.split_once(':')?;
                        Some(PerChargeMz {
                            charge: z.trim().parse().ok()?,
                            mz: mz.trim().parse().ok()?,
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();
        out.push(Feature {
            detected_mz: fnum(&f, c_detected_mz),
            rt_start: fnum(&f, c_rt_start),
            rt_apex: fnum(&f, c_rt_apex),
            rt_end: fnum(&f, c_rt_end),
            charge_states,
            per_charge_mz,
            primary_charge: inum(&f, c_primary),
            monoisotopic_mass: fnum(&f, c_mass),
            mono_mz: fnum(&f, c_mono_mz),
            summed_intensity: fnum(&f, c_intensity),
            cross_charge_support: unum(&f, c_support),
            num_members: unum(&f, c_members),
        });
    }
    Ok(out)
}

/// Neutral monoisotopic mass → m/z at a given charge.
pub fn mass_to_mz(mass: f64, charge: i32) -> f64 {
    (mass + charge as f64 * PROTON_MASS) / charge.max(1) as f64
}

// ------------------------------------------------------------------------ PSMs

/// One MetaMorpheus PSM (`.psmtsv` row), reduced to the fields the selector, the
/// feature link, the XIC extraction, and the MS2 spectrum pull actually use. Ten-column
/// read cap: the reader (`psm_tsv::Identification`) now resolves ten columns; we surface
/// nine of them — dropping the redundant base sequence (the full sequence is what we
/// display and filter on) — plus the derived theoretical precursor `mz`. Linking to a
/// detected feature is done client-side (mass + RT + charge join) where the feature list lives.
#[derive(Debug, Serialize, Clone)]
#[serde(rename_all = "camelCase")]
pub struct Psm {
    /// Full (modified) sequence — the display + filter key.
    pub full_sequence: String,
    /// Peptide monoisotopic (neutral) mass in daltons.
    pub monoisotopic_mass: f64,
    /// Precursor charge state.
    pub precursor_charge: i32,
    /// Theoretical precursor m/z derived from mass + charge — the XIC fallback when the
    /// PSM links to no detected feature.
    pub precursor_mz: f64,
    /// MS2 retention time (minutes); `-1.0` when the psmtsv omits it.
    pub ms2_retention_time: f64,
    /// One-based MS2 scan number; `-1` when the psmtsv omits it. Used to pull the exact
    /// identified MS2 spectrum for display when the PSM is selected.
    pub ms2_scan_number: i32,
    /// PSM q-value.
    pub q_value: f64,
    /// PSM score.
    pub score: f64,
    /// Source spectra-file name (extension stripped), for filtering to the open run.
    pub file_name: String,
    /// Decoy flag (the `Decoy/Contaminant/Target` cell contained a `D`).
    pub is_decoy: bool,
}

/// Parse a MetaMorpheus `.psmtsv` into [`Psm`] rows. Reads the file via the parity-tested
/// `psm_tsv` reader and computes each row's theoretical precursor m/z. Independent of any
/// open dataset — PSMs are overlaid on the raw data like features.
#[tauri::command]
pub async fn load_psms(path: String) -> Result<Vec<Psm>, ViewerError> {
    let ids = flashlfq_core::psm_tsv::read_identifications(&path)
        .map_err(|e| ViewerError::new("PSM_READ", format!("{path}: {e}")))?;
    Ok(ids
        .into_iter()
        .map(|id| {
            let precursor_mz = if id.precursor_charge_state > 0 {
                mass_to_mz(id.monoisotopic_mass, id.precursor_charge_state)
            } else {
                0.0
            };
            Psm {
                full_sequence: id.modified_sequence,
                monoisotopic_mass: id.monoisotopic_mass,
                precursor_charge: id.precursor_charge_state,
                precursor_mz,
                ms2_retention_time: id.ms2_retention_time_in_minutes,
                ms2_scan_number: id.ms2_scan_number,
                q_value: id.q_value,
                score: id.score,
                file_name: id.file_name,
                is_decoy: id.is_decoy,
            }
        })
        .collect())
}

// ---------------------------------------------------------- in-app detection

/// Options for `run_feature_detection`. `max_charge` is capped at 25 by default —
/// the guardrail that avoids the full-range top-down detect crash; `min_seed_intensity`
/// trades feature count/recall against runtime.
#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DetectOptions {
    pub max_charge: Option<i32>,
    pub min_seed_intensity: Option<f64>,
}

/// Map a resolved core feature to the wire `Feature` (same fields the runner's
/// resolved TSV carries; see `detect_features_tsv::write_tsv`).
fn resolved_to_feature(r: &ResolvedFeature) -> Feature {
    let member_detected_mz = |m: &RefinedFeature| -> Option<f64> {
        m.detected
            .peaks
            .iter()
            .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
            .map(|p| p.m() as f64)
    };
    let primary_member = r
        .members
        .iter()
        .max_by(|a, b| a.detected.summed_intensity.total_cmp(&b.detected.summed_intensity));
    let primary_charge = primary_member.map(|m| m.refined_charge).unwrap_or(0);
    let mono_mz = if primary_charge != 0 {
        mass_to_mz(r.monoisotopic_mass, primary_charge)
    } else {
        0.0
    };
    let detected_mz = primary_member.and_then(member_detected_mz).unwrap_or(0.0);
    let per_charge_mz: Vec<PerChargeMz> = r
        .charge_states
        .iter()
        .map(|&z| {
            let mz = r
                .members
                .iter()
                .filter(|m| m.refined_charge == z)
                .max_by(|a, b| a.detected.summed_intensity.total_cmp(&b.detected.summed_intensity))
                .and_then(member_detected_mz)
                .unwrap_or(0.0);
            PerChargeMz { charge: z, mz }
        })
        .collect();
    Feature {
        detected_mz,
        rt_start: r.start_rt,
        rt_apex: r.apex_rt,
        rt_end: r.end_rt,
        charge_states: r.charge_states.clone(),
        per_charge_mz,
        primary_charge,
        monoisotopic_mass: r.monoisotopic_mass,
        mono_mz,
        summed_intensity: r.summed_intensity,
        cross_charge_support: r.cross_charge_support as u32,
        num_members: r.members.len() as u32,
    }
}

/// The top-down detect → multi-envelope refine → apex-consensus pipeline, run
/// in-process (mirrors the `TOPDOWN` preset of `detect_features_tsv`, minus IsoDec).
/// Emits coarse phase progress on `on_progress`. CPU-bound — call from a blocking task.
fn run_topdown_pipeline(
    scans: &[Scan],
    engine: &PeakIndexingEngine,
    max_charge: i32,
    min_seed_intensity: f64,
    progress: impl Fn(&str, u32, u32),
) -> Vec<Feature> {
    // TOPDOWN globals: recharge over the full charge range, mono-offset fit on, no
    // cross-charge off-by-one bridge (conflates near-1-Da species).
    feature_refinement::set_recharge_max_charge(max_charge);
    feature_refinement::set_td_mono_fit_kmax(3);
    feature_refinement::set_offbyone_units(0);

    let base = TraceKernelParameters {
        min_charge: 1,
        max_charge,
        max_isotopes: 60,
        min_isotopes_observed: 3,
        min_seed_intensity,
        coverage_target: 1.0,
        trace_max_half_width_minutes: 1.0,
        min_feature_scans: 2,
        ..TraceKernelParameters::default()
    };
    let params = base.with_rt_from_index(engine, 36.0);

    progress("detecting", 0, 0);
    let detected = detect_features(engine, &params);

    // Composite averaging window from the measured FWHM (as the runner derives it).
    let mut avg = SpectralAveragingParameters::default();
    let fwhm_seconds = params.rt_sigma_minutes * FWHM_TO_SIGMA * 60.0;
    let spacing_seconds = median_ms1_scan_spacing_minutes(engine.scan_info()) * 60.0;
    avg.avg_scans = feature_refinement::derived_avg_scans(fwhm_seconds, spacing_seconds);

    let total = detected.len() as u32;
    let mut refined: Vec<RefinedFeature> = Vec::with_capacity(detected.len());
    for (i, f) in detected.iter().enumerate() {
        if let Some(r) = refine_feature_multi(f, scans, &avg, 20.0, 4) {
            refined.push(r);
        }
        if (i + 1) % 2000 == 0 || i + 1 == detected.len() {
            progress("refining", (i + 1) as u32, total);
        }
    }

    progress("resolving", total, total);
    let resolved = resolve_consensus_by_apex(&refined, 15.0, 0.3, 0);

    let feats: Vec<Feature> = resolved.iter().map(resolved_to_feature).collect();
    progress("done", total, total);
    feats
}

/// Run feature finding on an already-open dataset, in-process, returning the
/// resolved features directly (no TSV round-trip). Progress streams on the channel.
#[tauri::command]
pub async fn run_feature_detection(
    handle: u64,
    options: Option<DetectOptions>,
    on_progress: Channel<ProgressEvent>,
    state: State<'_, AppState>,
) -> Result<Vec<Feature>, ViewerError> {
    // Clone cheap Arc handles and drop the lock before the multi-minute run. Detection needs
    // the peak index, so this returns INDEXING until the background build lands.
    let (scans, _, engine) = indexed_handles(&state, handle)?;
    let opts = options.unwrap_or(DetectOptions { max_charge: None, min_seed_intensity: None });
    let max_charge = opts.max_charge.unwrap_or(25).clamp(1, 60);
    let min_seed = opts.min_seed_intensity.unwrap_or(10_000.0).max(0.0);

    tauri::async_runtime::spawn_blocking(move || {
        run_topdown_pipeline(&scans, &engine, max_charge, min_seed, |phase, done, total| {
            let _ = on_progress.send(ProgressEvent {
                phase: phase.into(),
                scans_done: done,
                scans_total: total,
            });
        })
    })
    .await
    .map_err(|e| ViewerError::internal(format!("detection task failed: {e}")))
}

// ------------------------------------------------ feature-finding walkthrough

/// One isotope tooth of a charge's comb in the ladder walkthrough (JSON DTO).
#[derive(Debug, Serialize, Clone)]
#[serde(rename_all = "camelCase")]
pub struct LadderToothDto {
    pub isotope_index: usize,
    pub expected_mz: f64,
    pub weight: f64,
    pub observed_mz: Option<f64>,
    pub observed_intensity: Option<f64>,
    pub credited: bool,
}

/// One charge state's comb over the seed's RT window for the candidate mass (JSON DTO).
#[derive(Debug, Serialize, Clone)]
#[serde(rename_all = "camelCase")]
pub struct LadderChargeDto {
    pub charge: i32,
    pub mono_mz: f64,
    pub spacing: f64,
    pub response: f64,
    pub num_isotopes_observed: usize,
    pub retained: bool,
    pub teeth: Vec<LadderToothDto>,
}

/// The full top-down charge-state-ladder fit for one `(seed, z_seed)` (JSON DTO), returned by
/// `score_seed_ladder` and rendered by the walkthrough panel.
#[derive(Debug, Serialize, Clone)]
#[serde(rename_all = "camelCase")]
pub struct SeedLadderDto {
    pub seed_mz: f64,
    pub seed_scan_index: i32,
    pub seed_rt: f64,
    pub z_seed: i32,
    pub screen_teeth: usize,
    pub screen_passed: bool,
    pub i_star: usize,
    pub mono_mz: f64,
    pub mono_mass: f64,
    pub mass_in_range: bool,
    pub refined_mono_mass: f64,
    pub total_response: f64,
    pub num_charge_states: usize,
    pub accepted: bool,
    pub window_scan_count: usize,
    /// Highest charge in the ladder (`min_charge` is fixed at 1) — the range the UI's z_seed picker spans.
    pub max_charge: i32,
    pub charges: Vec<LadderChargeDto>,
}

fn map_ladder(d: SeedLadderDiagnostics, max_charge: i32) -> SeedLadderDto {
    SeedLadderDto {
        seed_mz: d.seed_mz,
        seed_scan_index: d.seed_scan_index,
        seed_rt: d.seed_rt,
        z_seed: d.z_seed,
        screen_teeth: d.screen_teeth,
        screen_passed: d.screen_passed,
        i_star: d.i_star,
        mono_mz: d.mono_mz,
        mono_mass: d.mono_mass,
        mass_in_range: d.mass_in_range,
        refined_mono_mass: d.refined_mono_mass,
        total_response: d.total_response,
        num_charge_states: d.num_charge_states,
        accepted: d.accepted,
        window_scan_count: d.window_scan_count,
        max_charge,
        charges: d
            .charges
            .into_iter()
            .map(|c| LadderChargeDto {
                charge: c.charge,
                mono_mz: c.mono_mz,
                spacing: c.spacing,
                response: c.response,
                num_isotopes_observed: c.num_isotopes_observed,
                retained: c.retained,
                teeth: c
                    .teeth
                    .into_iter()
                    .map(|t| LadderToothDto {
                        isotope_index: t.isotope_index,
                        expected_mz: t.expected_mz,
                        weight: t.weight,
                        observed_mz: t.observed_mz,
                        observed_intensity: t.observed_intensity,
                        credited: t.credited,
                    })
                    .collect(),
            })
            .collect(),
    }
}

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct LadderOptions {
    /// Top of the charge ladder (clamped to 1..=60). Defaults to the top-down 30.
    pub max_charge: Option<i32>,
}

/// Compute the top-down charge-state-ladder fit for a manually chosen seed and anchoring charge,
/// for the interactive walkthrough. Takes a retention time (not a scan index) and resolves the
/// nearest MS1 scan in the engine's own index space, so the seed lookup can't drift from the
/// displayed spectrum. Requires the peak index (returns INDEXING otherwise).
#[tauri::command]
pub async fn score_seed_ladder(
    handle: u64,
    retention_time: f64,
    seed_mz: f64,
    z_seed: i32,
    options: Option<LadderOptions>,
    state: State<'_, AppState>,
) -> Result<SeedLadderDto, ViewerError> {
    let (_, _, engine) = indexed_handles(&state, handle)?;
    let max_charge = options.and_then(|o| o.max_charge).unwrap_or(30).clamp(1, 60);

    let diag = tauri::async_runtime::spawn_blocking(move || {
        // Nearest MS1 scan to the requested RT, in the engine's zero-based index space.
        let scan_info = PeakSource::scan_info(&*engine);
        let scan_index = scan_info
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (a.retention_time - retention_time)
                    .abs()
                    .total_cmp(&(b.retention_time - retention_time).abs())
            })
            .map(|(i, _)| i as i32)
            .unwrap_or(0);

        // Top-down detection parameters (mirrors `run_topdown_pipeline`'s base), so the ladder
        // reflects what real detection would score.
        let base = TraceKernelParameters {
            min_charge: 1,
            max_charge,
            max_isotopes: 60,
            min_isotopes_observed: 3,
            coverage_target: 1.0,
            trace_max_half_width_minutes: 1.0,
            min_feature_scans: 2,
            ..TraceKernelParameters::default()
        };
        let params = base.with_rt_from_index(&engine, 36.0);
        seed_ladder_diagnostics(&engine, seed_mz, scan_index, z_seed, &params)
    })
    .await
    .map_err(|e| ViewerError::internal(format!("ladder task failed: {e}")))?;

    Ok(map_ladder(diag, max_charge))
}

#[cfg(test)]
mod tests {
    use super::*;

    // Exactly the header + a row in the resolved TSV format written by
    // `detect_features_tsv::write_tsv`.
    const RESOLVED_TSV: &str = "Detected m/z (primary)\tRT Start\tRT Apex\tRT End\tCharge States\tPer-Charge Detected m/z\tNum Charge States\tPrimary Charge\tMonoisotopic Mass\tMono m/z (primary)\tSummed Intensity\tCross-Charge Support\tNum Members\n\
497.25840\t30.1200\t30.4500\t30.9000\t2;3\tz2:497.2584;z3:331.8416\t2\t2\t992.50220\t497.25840\t1.2345e6\t1\t4\n\
331.84160\t12.0000\t12.5000\t13.0000\t3\tz3:331.8416\t1\t3\t992.50220\t331.84836\t9.8765e5\t0\t2\n";

    #[test]
    fn parses_resolved_features() {
        let feats = parse_resolved_features(RESOLVED_TSV).expect("parse ok");
        assert_eq!(feats.len(), 2);

        let f = &feats[0];
        assert!((f.detected_mz - 497.2584).abs() < 1e-4);
        assert!((f.rt_start - 30.12).abs() < 1e-4);
        assert!((f.rt_apex - 30.45).abs() < 1e-4);
        assert!((f.rt_end - 30.90).abs() < 1e-4);
        assert_eq!(f.charge_states, vec![2, 3]);
        assert_eq!(f.per_charge_mz.len(), 2);
        assert_eq!(f.per_charge_mz[0].charge, 2);
        assert!((f.per_charge_mz[0].mz - 497.2584).abs() < 1e-4);
        assert_eq!(f.per_charge_mz[1].charge, 3);
        assert_eq!(f.primary_charge, 2);
        assert!((f.monoisotopic_mass - 992.5022).abs() < 1e-3);
        assert_eq!(f.cross_charge_support, 1);
        assert_eq!(f.num_members, 4);

        assert_eq!(feats[1].charge_states, vec![3]);
        assert_eq!(feats[1].per_charge_mz.len(), 1);
    }

    #[test]
    fn missing_column_errors() {
        let bad = "Detected m/z (primary)\tRT Start\n1.0\t2.0\n";
        let err = parse_resolved_features(bad).unwrap_err();
        assert_eq!(err.code, "FEATURE_PARSE");
    }

    #[test]
    fn empty_file_errors() {
        assert_eq!(parse_resolved_features("").unwrap_err().code, "FEATURE_PARSE");
    }

    // Real detector output on the Jurkat top-down fraction (machine-specific path).
    // Run explicitly: `cargo test --bin msviewer-tauri -- --ignored parses_real`.
    #[test]
    #[ignore]
    fn parses_real_jurkat_tsv() {
        let text = std::fs::read_to_string(r"D:\JurkatTopdown\jurkat_features.tsv").unwrap();
        let feats = parse_resolved_features(&text).unwrap();
        assert_eq!(feats.len(), 466_666);
        // First data row: a ~13.77 kDa proteoform seen at z 8,9,15,17,18,19.
        let f = &feats[0];
        assert_eq!(f.charge_states, vec![8, 9, 15, 17, 18, 19]);
        assert_eq!(f.per_charge_mz.len(), f.charge_states.len());
        assert_eq!(f.primary_charge, 18);
        assert!((f.monoisotopic_mass - 13766.51078).abs() < 1e-2);
    }

    // Full in-app pipeline on the real Jurkat raw (machine-specific; needs .NET 8).
    // Run: `cargo test --bin msviewer-tauri -- --ignored in_app_detection`.
    // High seed floor keeps it quick (few strong proteoforms).
    #[test]
    #[ignore]
    fn in_app_detection_on_real_raw() {
        let ds = build_dataset(r"D:\JurkatTopdown\02-18-20_jurkat_td_rep1_fract6.raw")
            .expect("read raw");
        assert_eq!(ds.scans.len(), 2851);
        let feats = run_topdown_pipeline(&ds.scans, &ds.engine, 25, 300_000.0, |_, _, _| {});
        assert!(!feats.is_empty(), "expected some features");
        let f = &feats[0];
        assert!(f.monoisotopic_mass > 0.0);
        assert!(!f.charge_states.is_empty());
        assert_eq!(f.per_charge_mz.len(), f.charge_states.len());
        assert!(f.primary_charge >= 1 && f.primary_charge <= 25);
        eprintln!(
            "in-app detection: {} features; strongest {:.2} Da z{:?}",
            feats.len(),
            f.monoisotopic_mass,
            f.charge_states
        );
    }

    #[test]
    fn mass_to_mz_matches_manual() {
        // Neutral 992.5022 Da at z=2 → (992.5022 + 2·1.00727664)/2.
        let mz = mass_to_mz(992.5022, 2);
        assert!((mz - 497.258476).abs() < 1e-4, "got {mz}");
    }
}
