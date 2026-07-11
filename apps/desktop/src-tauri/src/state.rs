//! App state: open datasets keyed by a u64 handle.
//!
//! M1+: each dataset is a real `OpenDataset` built from flashlfq-core — the MS1
//! `Vec<Scan>` read from the raw/mzML file, a per-scan TIC (summed centroided
//! intensity), a `PeakIndexingEngine` over those peaks (kept for future indexed
//! XIC + in-process feature detection), and the derived `DatasetMetadata`.
//!
//! Scan indexing: the viewer's `scanIndex` is the zero-based position in `scans`
//! (MS1-only), which is also the row position in the scan-summary Arrow batch and
//! the argument `get_spectrum`/`get_nearest_scan` take.

use std::collections::HashMap;
use std::sync::atomic::AtomicU64;
use std::sync::{Arc, Mutex};

use flashlfq_core::peak_indexing::{PeakIndexingEngine, Scan};

use crate::commands::DatasetMetadata;

/// One open dataset: the MS1 scans plus everything derived from them.
///
/// `scans` and `engine` are `Arc`-wrapped so in-app feature detection can clone
/// cheap read handles, drop the state lock, and run the (multi-minute) pipeline on
/// a blocking thread without blocking other commands.
pub struct OpenDataset {
    pub metadata: DatasetMetadata,
    /// MS1 scans in file order; `scanIndex` indexes this vector.
    pub scans: Arc<Vec<Scan>>,
    /// Per-scan summed centroided MS1 intensity (the served TIC), parallel to `scans`.
    pub tic: Vec<f32>,
    /// Peak index over `scans` — feeds in-process feature detection and future
    /// indexed XIC. Not required by the current range-XIC path (which scans directly).
    pub engine: Arc<PeakIndexingEngine>,
}

#[derive(Default)]
pub struct AppState {
    pub datasets: Mutex<HashMap<u64, OpenDataset>>,
    pub next: AtomicU64,
}
