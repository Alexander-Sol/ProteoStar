//! App state: open datasets keyed by a u64 handle.
//!
//! Two-phase load. When a dataset opens, only the **fast** fields are populated
//! synchronously: the native instrument TIC chromatogram (`tic_rt`/`tic_intensity`,
//! read straight from the file with no peak decode) plus a provisional `metadata`.
//! That lets the viewer paint the TIC in ~1-2 s. The heavy `indexed` state — the
//! full MS1 `Vec<Scan>`, its summed-centroid TIC, and the `PeakIndexingEngine` — is
//! built on a background task; commands that need per-scan peaks (spectra, range-XIC,
//! scan summaries, in-app detection) return an `INDEXING` error until it lands, and a
//! `dataset-indexed` event fires when it does.
//!
//! Scan indexing: the viewer's `scanIndex` is the zero-based position in `scans`
//! (MS1-only), which is also the row position in the scan-summary Arrow batch and
//! the argument `get_spectrum`/`get_nearest_scan` take. The native TIC's own point
//! ordinal is *not* an MS1 scan index (it spans all MS levels); TIC navigation is by
//! retention time, resolved to a scan via `get_nearest_scan`.

use std::collections::HashMap;
use std::sync::atomic::AtomicU64;
use std::sync::{Arc, Mutex};

use flashlfq_core::peak_indexing::{PeakIndexingEngine, RandomAccessMs1Reader, Scan};

use crate::commands::DatasetMetadata;

/// One open dataset.
///
/// The fast fields are set at open; `indexed` is filled later by the background
/// indexing task. `metadata` and `indexed` sit behind locks (and `Arc`s) so that
/// task can update them without holding the datasets map, and so a long-running
/// detection can clone cheap handles and drop the lock before the multi-minute run.
pub struct OpenDataset {
    /// Provisional at open (file name/format + RT range from the TIC); refined to
    /// exact RT/m-z ranges and the real MS1 scan count once indexing completes.
    pub metadata: Arc<Mutex<DatasetMetadata>>,
    /// Native instrument-TIC retention times (minutes), available immediately.
    pub tic_rt: Vec<f64>,
    /// Native instrument-TIC intensities, parallel to `tic_rt`.
    pub tic_intensity: Vec<f32>,
    /// The file kept open for on-demand single-scan reads, so spectra can be shown by RT
    /// before the full index exists. Behind a `Mutex` (the reader needs `&mut self` per read)
    /// and an `Arc` so a command can clone it, drop the state lock, and read off-thread.
    pub reader: Arc<Mutex<RandomAccessMs1Reader>>,
    /// Full indexed state; `None` until the background indexing task fills it.
    pub indexed: Arc<Mutex<Option<IndexedData>>>,
}

/// The heavy per-scan state, built off the open path by the background indexing task.
///
/// All three are `Arc`-wrapped so a command can clone read handles and drop the state
/// lock before doing real work (indexed queries, and the multi-minute detection run).
pub struct IndexedData {
    /// MS1 scans in file order; `scanIndex` indexes this vector.
    pub scans: Arc<Vec<Scan>>,
    /// Per-scan summed centroided MS1 intensity, parallel to `scans`.
    pub tic: Arc<Vec<f32>>,
    /// Peak index over `scans` — feeds indexed XIC + in-process feature detection.
    pub engine: Arc<PeakIndexingEngine>,
}

#[derive(Default)]
pub struct AppState {
    pub datasets: Mutex<HashMap<u64, OpenDataset>>,
    pub next: AtomicU64,
}
