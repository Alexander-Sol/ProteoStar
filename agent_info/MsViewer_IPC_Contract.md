# MsViewer — IPC Contract (v1)

The wire contract between the React front end and the Rust core over Tauri IPC. Companion to
`MsViewer_Architecture.md`. This pins down every command, its arguments, and the **exact byte
layout** of each response so both sides can be implemented independently.

---

## 0. Conventions

- **Numeric arrays are Apache Arrow, not hand-rolled bytes.** Rust builds an Arrow
  `RecordBatch` per query and serializes it to the **Arrow IPC stream** format; the command
  returns those bytes via `tauri::ipc::Response::new(bytes)`. JS `invoke` resolves to an
  `ArrayBuffer`, decoded with `apache-arrow` (`tableFromIPC`). This reuses the Rust `arrow`
  crate the core already depends on (`flashlfq-core` builds RecordBatches today) and drops all
  bespoke offset/alignment bookkeeping — column access on the JS side is `col.toArray()`,
  which yields a `TypedArray` for primitive columns.
- **Two response shapes:**
  - **JSON** — small, structured, variable-shape data (metadata, single scan, errors, progress).
  - **Arrow IPC** — the four numeric-array queries (scan summaries, TIC, XIC, spectrum).
- **Column dtypes:** `m/z → float64` (full precision, zero decode math), `intensity → float32`,
  `retentionTime → float32` (minutes), `scanIndex → uint32`. Nullable columns (precursor
  fields) use **Arrow nulls**, not sentinels. *(Fixed-point `uint32` m/z stays available as a
  payload optimization if profiling shows the m/z column dominates — but reduction caps counts,
  so float64 is the simpler default.)*
- **Scan index:** zero-based; it is the **row position** in the scan-summary batch.
- **Per-query scalar metadata** (e.g. a spectrum's scan number / precursor) rides in the Arrow
  **schema `custom_metadata`** key-value map, so each query is still a single response.
- **Version alignment:** the Rust `arrow` crate (54.x in `flashlfq-core`) and the JS
  `apache-arrow` package are versioned independently, but the **IPC stream format is stable and
  cross-compatible**; no lockstep required.
- **Handles:** `open_dataset` returns a `u64` handle. The core holds a map
  `handle → OpenDataset { reader, index, scan_summaries }`. The front end associates a handle
  with a UI slot (0 or 1); the core is slot-agnostic.

---

## 1. Command surface

| Command | Args | Response | Notes |
|---|---|---|---|
| `open_dataset` | `path: string`, `onProgress: Channel` | JSON `OpenResult` | async; builds index; streams progress |
| `close_dataset` | `handle: u64` | JSON `{}` | frees reader + index |
| `get_metadata` | `handle` | JSON `DatasetMetadata` | also embedded in `OpenResult` |
| `get_scan_summaries` | `handle` | **Arrow** `scan_summaries` | all scans, all MS levels |
| `get_nearest_scan` | `handle`, `retentionTime: f64`, `msLevel?: u32` | JSON `ScanSummary \| null` | for TIC click → select |
| `get_tic_trace` | `handle`, `rtMin?`, `rtMax?`, `msLevel: u32=1`, `maxPoints?` | **Arrow** `trace` | decimated |
| `get_range_xic` | `handle`, `mzLow: f64`, `mzHigh: f64`, `rtMin?`, `rtMax?`, `maxPoints?` | **Arrow** `trace` | window-XIC, MS1 index |
| `get_spectrum` | `handle`, `scanIndex: u32`, `mzMin?`, `mzMax?`, `maxPeaks?` | **Arrow** `spectrum` | any MS level (mzdata random access) |
| `get_ms2_for_precursor` | `handle`, `mz: f64`, `ms1ScanIndex: u32` | JSON `ScanSummary[]` | MS1→MS2 navigation |

`maxPoints` / `maxPeaks` are the **level-of-detail** knobs — the core decimates/aggregates
server-side so the payload is bounded regardless of file size (see `MsViewer_Architecture.md` §5, §6).

---

## 2. JSON schemas

### `OpenResult`
```jsonc
{
  "handle": 1,                       // u64
  "metadata": { /* DatasetMetadata */ }
}
```

### `DatasetMetadata`
```jsonc
{
  "fileName": "jurkat_rep2.raw",
  "format": "mzml" | "thermo_raw",
  "scanCount": 41235,                // all levels
  "ms1ScanCount": 8210,
  "msLevelsPresent": [1, 2],
  "retentionTimeRange": { "min": 0.01, "max": 120.0 } , // minutes, or null
  "mzRange": { "min": 350.0, "max": 1800.0 }            // MS1 index m/z span, or null
}
```

### `ScanSummary` (JSON form — used by `get_nearest_scan`, `get_ms2_for_precursor`)
```jsonc
{
  "scanIndex": 1204,                 // zero-based row position
  "oneBasedScanNumber": 1205,
  "retentionTime": 31.42,            // minutes
  "tic": 4.2e7,
  "msLevel": 2,
  "precursor": {                     // null for MS1
    "mz": 655.8123,
    "scanIndex": 1198,               // the MS1 scan this fragmented from (-1 if unknown)
    "isolationLow": 654.8,
    "isolationHigh": 656.8
  }
}
```

### `ViewerError` (rejected commands)
```jsonc
{ "code": "SCAN_OUT_OF_RANGE", "message": "scan index 99999 exceeds 41235" }
```
Codes: `FILE_NOT_FOUND · UNSUPPORTED_FORMAT · READ_ERROR · INDEX_BUILD_FAILED ·
HANDLE_NOT_FOUND · SCAN_OUT_OF_RANGE · EMPTY_INDEX · THERMO_RUNTIME_MISSING · INTERNAL`.

### `ProgressEvent` (streamed on the `open_dataset` channel)
```jsonc
{ "phase": "reading" | "indexing" | "done", "scansDone": 5000, "scansTotal": 41235 }
```

---

## 3. Arrow schemas

Each numeric-array command returns **one Arrow `RecordBatch`** serialized as an Arrow IPC
stream. Columns are named; the JS side reads them by name (`table.getChild("mz").toArray()`).
Nullable columns use Arrow's null bitmap — no sentinel values. Per-query scalars live in the
schema's `custom_metadata`.

### `scan_summaries` (one row per scan, all MS levels; `scanIndex` = row position)
| Column | Arrow type | Nullable | Notes |
|---|---|---|---|
| `oneBasedScanNumber` | `uint32` | no | |
| `retentionTime` | `float32` | no | minutes |
| `tic` | `float32` | no | per-format (§ arch 4b) |
| `msLevel` | `uint8` | no | 1, 2, … |
| `precursorMz` | `float64` | **yes** | null for MS1 |
| `precursorScanIndex` | `int32` | **yes** | null if unknown/MS1 |
| `isolationLow` | `float64` | **yes** | null for MS1 |
| `isolationHigh` | `float64` | **yes** | null for MS1 |

### `trace` (shared by `get_tic_trace` and `get_range_xic`; one row per point)
| Column | Arrow type | Nullable | Notes |
|---|---|---|---|
| `retentionTime` | `float32` | no | minutes |
| `intensity` | `float32` | no | TIC, or summed intensity over `[mzLow, mzHigh]` |
| `scanIndex` | `uint32` | no | representative scan (click → select) |

### `spectrum` (one row per peak) — columns
| Column | Arrow type | Nullable | Notes |
|---|---|---|---|
| `mz` | `float64` | no | |
| `intensity` | `float32` | no | |

**`spectrum` schema `custom_metadata`** (string values; parse on JS side):
`scanIndex`, `oneBasedScanNumber`, `retentionTime`, `msLevel`, and — MS2 only —
`precursorMz`, `precursorScanIndex`, `isolationLow`, `isolationHigh` (keys omitted for MS1).

---

## 4. TypeScript `DatasetProvider` interface

The whole UI depends only on this. The Tauri adapter (§5) implements it; unit tests supply a fake.

```ts
export interface NumericRange { min: number; max: number; }

export interface DatasetMetadata {
  fileName: string;
  format: "mzml" | "thermo_raw";
  scanCount: number;
  ms1ScanCount: number;
  msLevelsPresent: number[];
  retentionTimeRange: NumericRange | null;
  mzRange: NumericRange | null;
}

export interface Precursor {
  mz: number;
  scanIndex: number;        // -1 if unknown
  isolationLow: number;
  isolationHigh: number;
}

export interface ScanSummary {
  scanIndex: number;
  oneBasedScanNumber: number;
  retentionTime: number;
  tic: number;
  msLevel: number;
  precursor: Precursor | null;
}

export interface TicPoint { scanIndex: number; retentionTime: number; intensity: number; }
export interface XicPoint { scanIndex: number; retentionTime: number; intensity: number; }

export interface SpectrumPeak { mz: number; intensity: number; }
export interface Spectrum {
  scanIndex: number;
  oneBasedScanNumber: number;
  retentionTime: number;
  msLevel: number;
  precursor: Precursor | null;
  peaks: readonly SpectrumPeak[];
}

export interface DatasetProvider {
  getMetadata(): Promise<DatasetMetadata>;
  getScanSummaries(): Promise<readonly ScanSummary[]>;
  getNearestScan(retentionTime: number, msLevel?: number): Promise<ScanSummary | null>;
  getTicTrace(opts?: { rtRange?: NumericRange; msLevel?: number; maxPoints?: number }): Promise<readonly TicPoint[]>;
  getRangeXic(mzLow: number, mzHigh: number, opts?: { rtRange?: NumericRange; maxPoints?: number }): Promise<readonly XicPoint[]>;
  getSpectrum(scanIndex: number, opts?: { mzRange?: NumericRange; maxPeaks?: number }): Promise<Spectrum>;
  getMs2ForPrecursor(mz: number, ms1ScanIndex: number): Promise<readonly ScanSummary[]>;
}
```

`getScanSummaries` / `getTicTrace` / `getRangeXic` / `getSpectrum` decode Arrow (§5); the rest
are plain JSON. This matches the existing worker-provider shape — the RPC is just retargeted
from `postMessage` to `invoke`.

---

## 5. Transport adapter (`tauri-dataset-provider.ts`) — decode sketch

Decoding is `apache-arrow`'s `tableFromIPC`; primitive columns come back as `TypedArray` via
`.toArray()`. No manual offsets, no alignment bookkeeping.

```ts
import { invoke, Channel } from "@tauri-apps/api/core";
import { tableFromIPC, type Table } from "apache-arrow";

const toTable = (buf: ArrayBuffer): Table => tableFromIPC(new Uint8Array(buf));

export async function openDataset(
  path: string,
  onProgress?: (p: ProgressEvent) => void
): Promise<{ handle: number; provider: DatasetProvider }> {
  const channel = new Channel<ProgressEvent>();
  if (onProgress) channel.onmessage = onProgress;
  const { handle, metadata } = await invoke<OpenResult>("open_dataset", { path, onProgress: channel });

  const arrow = (cmd: string, args: object) => invoke<ArrayBuffer>(cmd, { handle, ...args }).then(toTable);

  const provider: DatasetProvider = {
    getMetadata: async () => metadata,

    getScanSummaries: async () => decodeScanSummaries(await arrow("get_scan_summaries", {})),

    getNearestScan: (retentionTime, msLevel) =>
      invoke("get_nearest_scan", { handle, retentionTime, msLevel }),

    getTicTrace: async (o = {}) =>
      decodeTrace(await arrow("get_tic_trace", {
        rtMin: o.rtRange?.min, rtMax: o.rtRange?.max, msLevel: o.msLevel ?? 1, maxPoints: o.maxPoints,
      })),

    getRangeXic: async (mzLow, mzHigh, o = {}) =>
      decodeTrace(await arrow("get_range_xic", {
        mzLow, mzHigh, rtMin: o.rtRange?.min, rtMax: o.rtRange?.max, maxPoints: o.maxPoints,
      })),

    getSpectrum: async (scanIndex, o = {}) =>
      decodeSpectrum(await arrow("get_spectrum", {
        scanIndex, mzMin: o.mzRange?.min, mzMax: o.mzRange?.max, maxPeaks: o.maxPeaks,
      })),

    getMs2ForPrecursor: (mz, ms1ScanIndex) =>
      invoke("get_ms2_for_precursor", { handle, mz, ms1ScanIndex }),
  };
  return { handle, provider };
}

// `trace` schema — shared by get_tic_trace and get_range_xic
function decodeTrace(t: Table): TicPoint[] {
  const rt = t.getChild("retentionTime")!.toArray() as Float32Array;
  const int = t.getChild("intensity")!.toArray() as Float32Array;
  const si = t.getChild("scanIndex")!.toArray() as Uint32Array;
  const out: TicPoint[] = new Array(rt.length);
  for (let i = 0; i < rt.length; i++) out[i] = { scanIndex: si[i], retentionTime: rt[i], intensity: int[i] };
  return out;
}

function decodeSpectrum(t: Table): Spectrum {
  const mz = t.getChild("mz")!.toArray() as Float64Array;
  const int = t.getChild("intensity")!.toArray() as Float32Array;
  const peaks: SpectrumPeak[] = new Array(mz.length);
  for (let i = 0; i < mz.length; i++) peaks[i] = { mz: mz[i], intensity: int[i] };

  const m = t.schema.metadata;                    // Map<string,string>
  const precMz = m.get("precursorMz");
  return {
    scanIndex: Number(m.get("scanIndex")),
    oneBasedScanNumber: Number(m.get("oneBasedScanNumber")),
    retentionTime: Number(m.get("retentionTime")),
    msLevel: Number(m.get("msLevel")),
    precursor: precMz === undefined ? null : {
      mz: Number(precMz),
      scanIndex: Number(m.get("precursorScanIndex") ?? -1),
      isolationLow: Number(m.get("isolationLow") ?? NaN),
      isolationHigh: Number(m.get("isolationHigh") ?? NaN),
    },
    peaks,
  };
}
// decodeScanSummaries: read the 8 named columns per §3; precursor fields are nullable —
// use the column's null check (row-wise) to emit `precursor: null` for MS1 rows.
```

**Future fast path (noted, not v1):** `plot-adapter` currently takes object arrays; the loops
above materialize objects from the Arrow `TypedArray`s. At very high peak counts that dominates.
Later, `plot-adapter` can consume the Arrow `TypedArray` columns directly (Plotly accepts
`Float32Array`/`Float64Array` for `x`/`y`), skipping per-point allocation entirely. Reduction
(`maxPeaks`/`maxPoints`) keeps v1 counts low enough that objects are fine to start.

---

## 6. Rust command signatures (`msviewer-tauri/commands.rs`)

```rust
use tauri::{State, ipc::{Channel, Response}};

#[derive(serde::Serialize)]
struct OpenResult { handle: u64, metadata: DatasetMetadata }

#[tauri::command]
async fn open_dataset(
    path: String,
    on_progress: Channel<ProgressEvent>,
    state: State<'_, AppState>,
) -> Result<OpenResult, ViewerError> { /* read → index (emit progress) → insert handle */ }

#[tauri::command]
async fn close_dataset(handle: u64, state: State<'_, AppState>) -> Result<(), ViewerError>;

#[tauri::command]
async fn get_metadata(handle: u64, state: State<'_, AppState>) -> Result<DatasetMetadata, ViewerError>;

#[tauri::command]                                   // Arrow: scan_summaries
async fn get_scan_summaries(handle: u64, state: State<'_, AppState>) -> Result<Response, ViewerError>;

#[tauri::command]
async fn get_nearest_scan(
    handle: u64, retention_time: f64, ms_level: Option<u32>, state: State<'_, AppState>,
) -> Result<Option<ScanSummary>, ViewerError>;

#[tauri::command]                                   // Arrow: trace
async fn get_tic_trace(
    handle: u64, rt_min: Option<f64>, rt_max: Option<f64>, ms_level: u32, max_points: Option<u32>,
    state: State<'_, AppState>,
) -> Result<Response, ViewerError>;

#[tauri::command]                                   // Arrow: trace
async fn get_range_xic(
    handle: u64, mz_low: f64, mz_high: f64, rt_min: Option<f64>, rt_max: Option<f64>,
    max_points: Option<u32>, state: State<'_, AppState>,
) -> Result<Response, ViewerError>;

#[tauri::command]                                   // Arrow: spectrum
async fn get_spectrum(
    handle: u64, scan_index: u32, mz_min: Option<f64>, mz_max: Option<f64>, max_peaks: Option<u32>,
    state: State<'_, AppState>,
) -> Result<Response, ViewerError>;

#[tauri::command]
async fn get_ms2_for_precursor(
    handle: u64, mz: f64, ms1_scan_index: u32, state: State<'_, AppState>,
) -> Result<Vec<ScanSummary>, ViewerError>;
```

`Response` is built with `Response::new(bytes)` where `bytes: Vec<u8>` is an Arrow IPC stream:
build a `RecordBatch` for the schema in §3, write it with `arrow::ipc::writer::StreamWriter`
into a `Vec<u8>`, attach any scalars via the schema's `custom_metadata`. A small `arrow_out.rs`
helper (one builder per schema) keeps each command's tail short. `AppState` holds
`Mutex<HashMap<u64, OpenDataset>>`; each `OpenDataset` owns its `mzdata` reader (for
`get_spectrum` random access) plus the built `PeakIndexingEngine` and cached scan summaries.

---

## 7. Mapping to reused `flashlfq-core` symbols

| Command | Backed by |
|---|---|
| `open_dataset` | `read_ms1_scans` (+ new all-levels reader) → `PeakIndexingEngine::index_peaks` |
| `get_scan_summaries` | `scan_info()` + per-format TIC + MS2 precursor from mzdata |
| `get_tic_trace` | scan summaries' TIC column, decimated to `maxPoints` |
| `get_range_xic` | **new** `get_range_xic` over `get_bins_in_range` |
| `get_spectrum` | **new** mzdata random-access read by `scan_index` + LOD |
| `get_nearest_scan` | linear/binary search over `scan_info()` RT |
| `get_ms2_for_precursor` | filter all-levels summaries by isolation window |

New Rust work is small and localized: the all-levels reader, `get_range_xic`, random-access
`get_spectrum`, per-format TIC, and the `arrow_out.rs` RecordBatch builders. Everything else is
reuse — and `flashlfq-core` already produces Arrow, so the builders have a head start.
