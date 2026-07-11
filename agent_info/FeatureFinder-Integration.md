# FlashLFQ feature-finder ↔ MsViewer integration — status & handoff

Integrates FlashLFQ-Rust MS1 feature detection into the MsViewer desktop app: real raw/mzML
reading, a feature-result overlay, and in-app detection. Companion to `MsViewer_Architecture.md`,
`MsViewer_IPC_Contract.md`, `MsViewer_FeatureOverlay.md`. Started 2026-07-10.

## What's done (all verified: cargo build/test, tsc, vite build; the app launches + runs)

**Real data reading (was M0 stub).** Vendored `flashlfq-core` wired into `apps/desktop/src-tauri`.
`open_dataset` → `read_ms1_scans` + `PeakIndexingEngine`; TIC/spectra/range-XIC/scan-summaries/
nearest-scan all served from the in-memory `Vec<Scan>` (no random access). MS1-only (the vendored
reader discards MS2), so `get_ms2_for_precursor` is empty and every scan is `msLevel 1`.
`OpenDataset.scans`/`engine` are `Arc` (so detection can clone handles and drop the state lock).

**Feature overlay.** `load_features(path)` parses the runner's **resolved** TSV
(`parse_resolved_features`, header-keyed, unit-tested + real-466k-file test). Frontend: `features.ts`
(loadFeatures, isotopeGrid), plot-adapter extended with `featureRug`+`regions` (TIC) and `envelope`
dotted isotope lines (spectrum). `App.tsx`: file-open + load-features dialogs (`@tauri-apps/plugin-dialog`),
a feature-list drawer, click-to-select → zoom TIC + load apex spectrum + draw predicted isotope grid.
Features sorted by intensity on load; drawer capped 800, rug capped 8000 (map is 100k+ on noisy runs).

**In-app detection (`run_feature_detection`).** Re-vendored the NEWEST `flashlfq-core` (added
`isodec.rs`, `feature_export.rs`, and the `models/` dir it `include_bytes!`s). `run_topdown_pipeline`
runs detect → `refine_feature_multi` → `resolve_consensus_by_apex` on a spawn_blocking task, streams
detect/refine/resolve progress, returns `Vec<Feature>` into the same overlay. `maxCharge`=25 default
(crash guardrail), `minSeedIntensity` default 10000. No IsoDec in-app (detector charge only). Frontend:
"Run feature finding" button. Integration-tested on the real Jurkat raw (98s release, passes).

**UI fixes (2026-07-10).** TIC click now routed through Plotly's own `onClick` (the pan drag-layer
swallows raw DOM clicks) using the line-point's RT, plus an invisible clickable marker layer on the TIC
line. Pin buttons made mutually exclusive.

## Test data
- Raw: `D:\JurkatTopdown\02-18-20_jurkat_td_rep1_fract6.raw` (2851 MS1 scans, 22.6M peaks).
- Demo resolved TSV: `D:\JurkatTopdown\jurkat_features.tsv` (466,666 features; from a capped run
  `TOPDOWN=1 ISODEC_CHARGE=0 MAX_CHARGE=25 MAX_ISOTOPES=40 MIN_SEED_INTENSITY=5000`).

## Run
`cd apps/desktop && bun run tauri dev` (MUST be from `apps/desktop` — the repo-root package.json has
no `tauri` script). Debug build → **in-app detection is ~15× slower**; use `bun run tauri dev --release`
for interactive detection. TSV-load path is fast in debug.

## Build gotchas (all fixed — re-check if they reappear)
1. `thermorawfilereader`/`dotnetrawfilereader-sys` must be pinned to **=0.7.0** in the src-tauri lock;
   0.7.1 breaks `mzdata 0.65.2` (`raw_view`/`bytes`). `cargo update -p <crate> --precise 0.7.0`.
2. Windows `tauri-build` needs `icons/icon.ico` (`bun x @tauri-apps/cli icon icons/icon.png`, then keep
   only icon.ico + the original PNGs; delete the android/ios/Square* cruft it also emits).
3. When re-vendoring the core, copy `models/` too (isodec `include_bytes!`).
4. `bunx` isn't installed here — tauri.conf.json uses `bun x` for before{Dev,Build}Command.
5. `@tauri-apps/plugin-dialog@^2` must be in apps/desktop/package.json (a failed `bun add` once left it
   in node_modules only).

## NEXT / PENDING: fast-TIC-first load (proposed, NOT built)
**Problem:** the TIC takes ~20s to appear because `open_dataset` reads every scan's full centroided
peaks + builds the index before summing the TIC. We do NOT use the native Thermo `GetTic`.

**Plan (two-phase load):**
1. Fast pass: iterate at `DetailLevel::MetadataOnly` pulling per-MS1-scan `(rt, tic, msLevel)` — no peak
   decode → TIC (the *real* instrument TIC) appears in ~1–2s. Requires adding `mzdata` as a direct
   src-tauri dep (same features as core: `mzml, miniz_oxide, thermo`); verify the 0.65 MetadataOnly +
   `total_ion_current()` API (vendored core only ever uses `DetailLevel::Full`).
2. Background index: full `read_ms1_scans` + `PeakIndexingEngine` on a task; `get_spectrum`/
   `get_range_xic`/`run_feature_detection` return an `INDEXING` state until ready, then light up. Needs a
   readiness signal to the UI (Tauri event, or poll).
3. Optional further win: `get_spectrum` does an on-demand random-access single-scan read via mzdata
   (`get_spectrum_by_index`) so spectra are available immediately, not gated on the full index — this is
   the random-access path the original architecture doc anticipated (keep the reader open behind a Mutex).

`OpenDataset` would split into always-present `fast_rt`/`fast_tic`/`fast_one_based` + optional
`Arc<Vec<Scan>>`/`Arc<PeakIndexingEngine>` + an `indexed` flag. User asked for this; awaiting go-ahead.

## Known limitations
- MS1 only (no MS2 nav; served TIC currently = summed centroids, fixed by the plan above).
- In-app detect capped at charge 25 (full 1–60 crashes the flashlfq-rust detector, exit -1 mid-detect).
- No IsoDec in-app (the `+10pp` neural charge step; `isodec` module is vendored, port is a follow-up).
- Vendored core snapshot must be re-synced from flashlfq-rust when the upstream pipeline changes.
