# Spectra during indexing — attempt, reverted (deferred)

## Goal
Let the MsViewer desktop app display MS1 spectra (by retention time, via a TIC click or feature
selection) **while the background peak index is still building**, not only after it finishes. The
fast native TIC already paints first; the ask was for on-demand spectra during the multi-minute
index build.

## Why it doesn't work today
For Thermo `.raw`, both the full-file index read (`read_ms1_scans`) and the on-demand single-scan
reads (`RandomAccessMs1Reader`, used by `get_spectrum_at_rt`) cross into the same in-process **.NET
runtime** hosted by the `thermorawfilereader` crate. While the index build is iterating every scan,
on-demand reads are starved on that shared runtime, so a spectrum request effectively blocks until
indexing completes. (mzML is pure-Rust and would not have this contention.)

## What was attempted (and is being reverted)
A **streaming partial-scan buffer**: the background build streamed each MS1 scan into an in-memory
`Vec<Scan>` as it was decoded, and `get_spectrum_at_rt` served the nearest already-read scan from
that buffer during indexing, bypassing the contended reader.

**It did not solve the problem** — spectra still did not display during indexing — and it
**regressed** both time-to-first-TIC and total indexing time. Decision: revert it, accept that
spectra are unavailable until indexing finishes, and move on. The fast-path TIC must still display
first (it does; that path was never touched).

## Exactly what was reverted
All edits below were removed, restoring the code to its pre-attempt state.

### `crates/flashlfq-core/src/peak_indexing.rs`
- Removed `pub fn stream_ms1_scans(...)` and the private `stream_ms1_scans_from_iter(...)` helper.
- Restored `collect_ms1_scans` to its original standalone loop (no longer delegates to the streaming
  helper).

### `apps/desktop/src-tauri/src/state.rs`
- Removed the `partial_scans: Arc<Mutex<Vec<Scan>>>` field from `OpenDataset`.

### `apps/desktop/src-tauri/src/commands.rs`
- Import: `stream_ms1_scans` → back to `read_ms1_scans`.
- `build_index_from_scans(scans: Vec<Scan>)` → back to `build_indexed(path: &str)` (reads the scans
  itself via `read_ms1_scans`).
- `open_dataset`: removed the `partial_scans` field initialiser and the `partial_slot` clone; the
  background task is back to `spawn_blocking(move || build_indexed(&path_for_index))`.
- `get_spectrum_at_rt`: removed the partial-buffer fast path and the `nearest_scan_by_rt` helper;
  back to reader-only.

## What was NOT reverted (kept)
These are unrelated to the indexing attempt and stay:
- Walkthrough charge-ladder panel and the `score_seed_ladder` command.
- Spectrum **hover/click** fix (invisible marker hit-target over the thin bars).
- Comb overlay **detected** indicator (solid vs dotted per observed peak).
- Spectrum **zoom on charge selection** behavior.
- **Arrow-key MS1 scan navigation** (← / → step scans, zoom held constant).
- Drawer `rightInset` layout fix (plots reflow when a drawer opens).

## If revisited later
The partial-buffer idea is a dead end as implemented. Better directions to explore:
- Drive **both** the index build and on-demand reads from a **single shared reader** so there is
  one .NET consumer, and yield to pending spectrum requests between scan reads.
- Reduce index-build wall-clock so the unavailable window is short enough not to matter.
- Investigate whether the `thermorawfilereader` .NET host can be made concurrent / re-entrant.

## Verification of the revert
`cargo check` passes for `flashlfq-core` and `msviewer-tauri`; the full `flashlfq-core` suite passes
(243 tests).
