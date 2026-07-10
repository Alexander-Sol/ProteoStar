//! App state: open datasets keyed by a u64 handle. In M0 each dataset is a
//! `StubDataset` (synthetic metadata only); M1 replaces this with an
//! `OpenDataset { reader, index, scan_summaries }` from flashlfq-core.

use std::collections::HashMap;
use std::sync::atomic::AtomicU64;
use std::sync::Mutex;

use crate::commands::DatasetMetadata;

/// One open (stub) dataset. M1 will grow a reader + PeakIndexingEngine here.
pub struct StubDataset {
    pub metadata: DatasetMetadata,
}

#[derive(Default)]
pub struct AppState {
    pub datasets: Mutex<HashMap<u64, StubDataset>>,
    pub next: AtomicU64,
}
