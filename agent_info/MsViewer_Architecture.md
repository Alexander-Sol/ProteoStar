# MsViewer — Architecture Design

**Status:** Design draft · **Supersedes:** the `.imsp`-based roadmap for the new tool.

A desktop mass-spec data viewer: a **Rust core** (reads raw data via `mzdata`, indexes
peaks, extracts XICs) behind a **Tauri** shell, driving the **existing React front end**
from the MsBrowser repo. No `.imsp` files — raw data is read directly.

---

## 1. Decisions locked

| Area | Decision |
|---|---|
| Deployment | **Tauri desktop app** — Rust backend + React UI in the webview, IPC via `invoke` |
| Input | **Raw data read in Rust via `mzdata`** (mzML/mzMLb/MGF; Thermo `.raw` / Bruker `.d` behind feature flags). No `.imsp`. |
| Peak index / XIC | **Reuse existing Rust indexing engine** (already written; indexes peaks, fast XIC extraction) |
| Query API | Extend the contract with **server-side reduction** — decimated TIC, aggregated XIC, level-of-detail spectra |
| MS levels | **MS1 + MS2** (§9). Index stays MS1-only; summaries & spectra span all levels; precursor→fragment navigation. |
| Core reuse | Lift reader + `PeakIndexingEngine` from `flashlfq-rust/rust/flashlfq-core` (§4b) — already `mzdata`-based and parity-tested. |
| Wire format | **Apache Arrow, not JSON** (§6a) — numeric arrays cross the IPC boundary as Arrow IPC streams; reuses the `arrow` crate the core already has. Decided (see `MsViewer_IPC_Contract.md`). |
| TIC source | **Per-format, computed in Rust** (§4b) — Thermo `.raw` uses the native `GetTic`; mzML reads `total_ion_current` / computes. Decided. |
| Window-XIC | **New `get_range_xic` core method** (§4b) — the default chromatogram XIC; not the FlashLFQ traced `get_xic`. Decided. |
| Disk cache | **Open — benchmark required** (§12) — measure disk-cached index reload vs. rebuild-on-open before committing. |

## 2. Guiding principle — keep the seam

The current front end already talks to data through **one interface**, `DatasetProvider`
(6 async methods). Every UI layer — reducer, plots, components — depends only on that
interface, never on the binary format. The whole migration is therefore:

> **Keep `DatasetProvider` as the contract. Delete the TS `.imsp` parser. Write one new
> transport adapter that fulfills the contract by calling Rust over Tauri IPC.**

Nothing in `viewer-state`, `plot-adapter`, or `ui` changes. That is the payoff of the
existing layering, and it's the invariant we protect through every milestone.

## 3. Reuse map (what happens to each existing piece)

| Existing | Fate |
|---|---|
| `packages/ui` (layout, primitives, inputs) | **Keep verbatim** |
| `packages/plot-adapter` (Plotly TIC + spectrum, DTOs) | **Keep verbatim** |
| `packages/viewer-state` (reducer, store, 2-slot compare, pin/zoom) | **Keep verbatim** |
| `DatasetProvider` interface + worker-protocol *types* | **Keep as the contract** (move to a shared `contract` package; extend with reduction ops) |
| `parseImsp`, `createImspDatasetProvider`, `IMSP_*` constants | **Delete** — Rust owns reading now |
| `apps/web/app/imsp-worker.ts`, `worker-dataset-provider.ts` | **Replace** with a Tauri adapter (the RPC shape is the template) |
| `viewer-controller.ts` `loadViewerDataset` | **Rework**: native file dialog → `invoke('open_dataset', {path})` → handle; no more `ArrayBuffer` in JS |
| `apps/web` (Next.js) | **Reshell** into a Tauri front end (Vite recommended over Next for a desktop SPA; see §8) |

## 4. Rust workspace layout

Split the pure engine from the Tauri shell so the engine stays testable and portable
(your existing code drops into `core`):

```
crates/
  msviewer-core/        # lib, NO tauri dependency
    ├─ source/          # mzdata wrappers: open, iterate, random-access spectra
    ├─ index/           # << your existing peak-index + XIC engine lives here >>
    ├─ query/           # query engine: metadata, scan summaries, TIC, XIC, spectrum LOD
    ├─ reduce/          # decimation / aggregation (server-side reduction)
    └─ model.rs         # DTOs mirroring the TS contract (serde)
  msviewer-tauri/       # bin, the shell
    ├─ commands.rs      # #[tauri::command] wrappers over core
    ├─ state.rs         # open datasets keyed by handle id (slot 0 / slot 1)
    └─ main.rs
```

`msviewer-core` is a plain Rust lib with no UI dependency — unit-testable against small
mzML fixtures, and reusable outside Tauri (CLI, tests, a future server) if ever needed.

## 4b. Reuse from `flashlfq-rust` (`rust/flashlfq-core/src/peak_indexing.rs`)

Most of the core already exists and is parity-tested. It reads raw data via `mzdata 0.65`
(mzML pure-Rust; Thermo `.raw` via `ThermoRawReader`, which bundles a .NET 8 runtime) with
format sniffing already done. The viewer lifts these directly:

| Viewer need | Existing symbol | Notes |
|---|---|---|
| Open raw file | `read_ms1_scans(path)` → `collect_ms1_scans(reader)` over `MZReader::open_path` / `ThermoRawReader` | Format sniff + centroiding + zero-intensity strip all done |
| Build m/z index | `PeakIndexingEngine::index_peaks(&scans)` / `from_spectra_file(path)` | Same `BINS_PER_DALTON = 100` binning the `.imsp` used |
| Scan summaries | `scan_info() → &[ScanInfo]` = `{one_based_scan_number, zero_based_scan_index, retention_time, msn_order}` | `msn_order` already present — MS-level ready. **Gap: no TIC field** — add it, sourced **per format in Rust** (see TIC note below) |
| XIC | `get_xic(m, rt, ppm, …)`, `get_xic_by_scan_index(…)`, `get_all_xics(…)`, `ExtractedIonChromatogram { apex_rt, apex_scan_index }` | See semantics note below |
| Point query | `get_indexed_peak(m, scanIndex, ppm)` | Nearest-peak-in-tolerance |
| All peaks | `all_peaks()` | bin-then-scan order |
| Peak model | `IndexedMassSpectralPeak { mz, intensity, zero_based_scan_index, retention_time }` (f32) | — |

**XIC — new core method (decided).** The existing `get_xic` is FlashLFQ's *traced* peak:
it walks one m/z across RT with a ppm tolerance, missed-scan allowance, and a half-width
limit — it stops at gaps and follows a single trace. The viewer's default chromatogram
wants a **simple m/z-window XIC**: sum every peak in `[mzLo, mzHi]` per scan across the whole
run, no gap logic. This is a **new `get_range_xic` core method**, built on the existing bins
(`get_bins_in_range` + per-scan accumulation). Both stay available — window-XIC for the
chromatogram view, traced `get_xic` (M6) when the user wants FlashLFQ-style peak tracing.

**TIC — per-format, computed in Rust (decided).** TIC is sourced by data type behind one
core function so the front end never sees the difference:
- **Thermo `.raw`** → the native **`GetTic`** the reader exposes (authoritative instrument TIC).
- **mzML** → the per-spectrum `total_ion_current` cvParam when present; otherwise sum
  filtered intensities.
This fills the `ScanInfo` TIC gap above and keeps `getTicTrace` format-independent on the JS side.

## 4c. Two gaps the viewer must add to the core

The index is intentionally **MS1-only** (XIC is an MS1 concept — keep it that way). The
viewer needs two things the index doesn't provide:

1. **Full-spectrum readout for any scan (MS1 *or* MS2)** — `getSpectrum(scanIndex)` must
   return the actual peak list of that scan, including MS2 fragment spectra, which are never
   indexed. Serve this from **mzdata random access** (`MZReader` as a
   `RandomAccessSpectrumIterator` — read the spectrum by index on demand), not from the
   index. This means the dataset **handle keeps the mzdata reader open** alongside the index
   (see §6b), not just the built index.
2. **All-levels scan metadata + precursor info** — `read_ms1_scans` filters to MS1. Add a
   sibling reader that retains MS2 scans' summaries and their **precursor** (isolation m/z /
   window / precursor-scan link) from mzdata, to drive MS1→MS2 navigation (§9).

## 5. The query contract (evolved)

TypeScript interface the UI keeps depending on; Rust implements it. Reduction params are
the new part — Rust does the aggregation so payloads over IPC stay small.

```ts
interface DatasetProvider {
  getMetadata(): Promise<DatasetMetadata>;            // run-level: rt range, m/z range, scan count, ms levels
  getScanSummaries(): Promise<ScanSummary[]>;         // + msLevel per scan (forward-compat for MS2)
  getNearestScan(rt: number): Promise<ScanSummary | null>;

  getTicTrace(opts?: { rtRange?; maxPoints?; }): Promise<TicPoint[]>;      // decimated to maxPoints
  getXic(mzLo, mzHi, opts?: { rtRange?; }): Promise<XicPoint[]>;          // aggregated per scan  (NEW)
  getSpectrum(scanIndex, opts?: { mzRange?; maxPeaks?; }): Promise<Spectrum>; // level-of-detail

  // getPeaksInMzRange stays available for raw/debug access, but getXic is the UI path.
}
```

Rust command surface (Tauri `invoke`) mirrors it 1:1:
`open_dataset · close_dataset · get_metadata · get_scan_summaries · get_nearest_scan ·
get_tic_trace · get_xic · get_spectrum`.

## 6. Two load-bearing implementation concerns

### 6a. Wire format for numeric arrays — Apache Arrow (decided)
Tauri v2 IPC serializes to JSON by default — fine for metadata, wasteful for the arrays
that dominate MS data (spectrum peaks, XIC/TIC points). **We pass numeric arrays as Apache
Arrow IPC streams, not JSON**, rather than hand-rolling a byte format. Rationale:
- The Rust core **already depends on `arrow`** (`flashlfq-core` builds RecordBatches today),
  so the encode side has a head start and no new bespoke format to maintain.
- JS decodes with `apache-arrow`'s `tableFromIPC`; primitive columns come back as
  `TypedArray` via `.toArray()` — no manual byte offsets, no alignment bookkeeping.
- **Server-side reduction** still caps array sizes before they cross the boundary, so the
  Arrow payloads stay small; Arrow's structure/compression buys little at those sizes but its
  ergonomics and reuse win.
- Metadata / scan summaries that aren't array-shaped stay JSON (small, structured).

Column dtypes (see `MsViewer_IPC_Contract.md` §0, §3): **intensity → f32**, **RT → f32**,
**scanIndex → u32**, **m/z → f64** (full precision, zero decode math). Nullable columns
(precursor fields) use **Arrow nulls**, not sentinels. *(The earlier fixed-point `uint32` m/z
plan is superseded — with Arrow, f64 is simpler and more precise; fixed-point stays available
as a payload optimization only if profiling shows the m/z column dominates.)*

Implementation: commands return `tauri::ipc::Response::new(bytes)` where `bytes` is an Arrow IPC
stream (`arrow::ipc::writer::StreamWriter`); per-query scalars ride in the schema
`custom_metadata`. Full schemas and the decode adapter live in `MsViewer_IPC_Contract.md`.

### 6b. Dataset handles, async load, progress
- Opening + indexing a large run (your Jurkat file is 64 MB, millions of peaks) takes
  seconds. `open_dataset` is **async**, runs on a Tokio/worker task so the UI thread never
  blocks, and streams progress via a **Tauri Channel/event**. The existing
  `DatasetLoadState` already models `loading → ready/error` — wire progress into it.
- Core holds **multiple open datasets keyed by handle id**, matching the front end's
  `datasetSlots: [slot0, slot1]` (the 2-run comparison view — the "match between runs"
  workflow). Each slot = one open reader + its index.
- Native **file dialog + fs** come from Tauri: the JS never reads an `ArrayBuffer`; it
  passes a path. This removes the browser-memory ceiling entirely — the big win of going desktop.

## 7. Data flow (single query)

```
[React UI] --dispatch--> [viewer-state reducer]        (unchanged)
   |
   v  provider.getXic(mzLo,mzHi)
[tauri-dataset-provider.ts]  --invoke('get_xic', {handle,mzLo,mzHi})-->
[commands.rs] -> [msviewer-core: index -> reduce] -> packed bytes
   |
   v  rebuild Float32Array
[plot-adapter]  renders XIC          (unchanged)
```

## 8. Front-end shell: Next.js → Vite

The current app is Next.js (SSR-oriented). A Tauri desktop app wants a **static SPA**;
**Vite + React** is the conventional Tauri front end and drops the SSR machinery you don't
need. `ui`, `plot-adapter`, and `viewer-state` are framework-agnostic and move over
untouched; only the `apps/web` shell (`page.tsx`, `layout.tsx`, worker files) is rebuilt as
a Vite entry that mounts `ViewerPage`.

## 9. MS1 + MS2 (decided)

The tool handles **both** MS levels. The split of responsibility:
- **Index = MS1 only.** XIC, m/z-range queries, and the TIC chromatogram operate on the
  MS1 peak index (reuses `PeakIndexingEngine` unchanged).
- **Scan summaries span all levels.** `ScanSummary` carries `msLevel` (already `msn_order`
  in `ScanInfo`) plus, for MS2, a `precursor { mz, isolationWindow, precursorScanIndex }`.
- **Spectrum readout spans all levels** via mzdata random access (§4c), so MS2 fragment
  spectra render in the same spectrum panel.

**MS2 navigation (the UI surface):**
- MS-level filter on the scan list / TIC (show MS1, MS2, or both).
- Click an MS1 peak (or scan) → list the MS2 scans whose isolation window contains that m/z
  → select one → its fragment spectrum loads in the spectrum panel.
- Breadcrumb back to the precursor MS1 scan via `precursorScanIndex`.
- The existing 2-slot compare + pin/zoom machinery is untouched; MS2 is additive UI.

**Core additions for MS2:** the all-levels reader (§4c #2) and the random-access spectrum
reader (§4c #1). No change to the indexing engine.

## 10. Milestones

Much of M1/M3 is *wiring existing, parity-tested code*, not new algorithms — the reuse in
§4b compresses the risky middle of the plan.

| # | Goal | Proves |
|---|---|---|
| M0 ✅ | Scaffold Tauri + Vite (`apps/desktop`), reuse `plot-adapter`/`ui`, stub commands return **real Arrow IPC** | **Done** — webview ↔ Rust round-trip verified; `cargo build` + `vite build` green; window launches |
| M1 | Lift `flashlfq-core` reader + `PeakIndexingEngine` into `msviewer-core`; `open_dataset` (async + progress), `get_metadata`, `get_scan_summaries` (+TIC), `get_tic_trace` | real mzML/`.raw` → TIC chromatogram renders |
| M2 | `get_spectrum` via mzdata random access (LOD), `get_nearest_scan`; click TIC → spectrum | single-dataset spectrum readout, all MS levels |
| M3 | `get_range_xic` (new, on existing bins) + reduction params + **binary wire format** | fast window-XIC on the 64 MB file, small payloads |
| M4 | Two handles (slot 0/1), robust async load + progress, error taxonomy | dual-run compare, large-file UX |
| M5 | MS2: all-levels reader + precursor info, MS-level filter, precursor→fragment navigation | DDA/DIA workflows |
| M6 | *(optional)* traced XIC exposed (`get_xic`) for FlashLFQ-style peak tracing | quant-adjacent power feature |

## 11. Testing strategy

- **Rust core:** unit tests against tiny hand-checkable mzML fixtures with known metadata,
  scan counts, a known spectrum, and known XIC results — the same discipline as the old
  `tiny-known.imsp`, now in mzML. Index/XIC engine tested in isolation (no Tauri).
- **Contract fixtures:** a shared JSON of expected query results, asserted by Rust tests —
  the single source of truth for what each query must return.
- **Front end:** `viewer-page.test.tsx` and friends keep running against a **fake
  `DatasetProvider`** — no Tauri, no Rust needed in the UI test loop. This is exactly why we
  preserve the interface: UI tests stay pure and fast.

## 12. Risks / things to pin down next

1. **`mzdata` random-access read-by-index** (§4c #1) — confirm `MZReader` /
   `ThermoRawReader` support efficient `get_spectrum_by_index` (RandomAccessSpectrumIterator)
   for the formats you target. This is the one genuinely new dependency the viewer leans on
   that FlashLFQ's streaming pass didn't need. **Verify first** — it gates `getSpectrum`.
2. **Reader lifetime + threading** — the dataset handle keeps a `mzdata` reader open for
   random access; its read is `&mut self`, so wrap each dataset (reader + index) behind a
   `Mutex` on its own task. Thermo `.raw` pulls in a .NET 8 runtime — confirm that's
   acceptable on target machines (mzML stays pure-Rust).
3. **Wire format** — ✅ **decided: Apache Arrow IPC** (§6a; contract doc). Column dtypes
   settled: m/z **f64**, intensity/RT **f32**, scanIndex **u32**, precursor fields nullable.
4. **TIC source** — ✅ **decided: per-format in Rust** (§4b). Thermo `.raw` → native `GetTic`;
   mzML → `total_ion_current` cvParam / computed. One core function, format-independent to JS.
5. **XIC semantics** — ✅ **decided: new `get_range_xic`** window-XIC is the default
   chromatogram (§4b); traced `get_xic` is an M6 power feature.
6. **Disk cache vs. rebuild-on-open** — ⏳ **benchmark before deciding.** Measure, on the
   medium/large/stress files: (a) rebuild the index from the raw file each open, vs. (b)
   write a compact index cache to disk on first open and mmap/reload it after. Compare
   cold-open latency, warm-reopen latency, disk footprint, and invalidation cost. Only adopt
   the cache if the reopen win is real on your typical files. (Note: the old `.imsp` format
   is effectively this cache's on-disk shape if we go this route.) Folds together with the
   index memory-footprint question on the 250 MB+ stress file.
7. **`flashlfq-core` extraction boundary** — ✅ **decided: vendored into this repo** at
   `crates/flashlfq-core/` (library source + `data/`, minus the 13 MB parity goldens;
   `Cargo.toml` de-workspaced, `arrow`/`parquet` pinned to 54.x, `thermo` feature kept).
   Builds standalone here (`cargo build --lib`, 198 crates, ~48 s cold). Not yet wired into
   `apps/desktop/src-tauri` — that happens in M1 when the stub command bodies are replaced.
   May later split just the indexing module into a leaner `msviewer-core` if the quant
   surface proves to be dead weight.
