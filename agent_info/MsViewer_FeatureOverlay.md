# MsViewer — Feature-finding overlay

Visualize FlashLFQ-Rust feature-detection output on top of the raw data, to see where
feature finding succeeds and fails. Built on the M1 real-data path (raw/mzML read in Rust
via the vendored `flashlfq-core`).

## What it does

1. **Open a raw file** (`Open file…`) — reads the MS1 scans (`read_ms1_scans` → `PeakIndexingEngine`)
   and renders the summed-MS1-intensity chromatogram. Click anywhere to load that scan's spectrum.
2. **Load features** (`Load features…`) — parses the runner's **resolved** feature TSV
   (`detect_features_tsv`'s final `out.tsv`) and overlays them:
   - **TIC:** a marker at each feature's apex RT (the "feature rug"), coloured by charge. Click a
     marker — or a row in the **Features** drawer — to select it.
   - On select: the TIC zooms to the feature's traced elution window (shaded band), the apex-scan
     spectrum loads, and the spectrum panel draws **dotted lines at the predicted isotope m/z**
     (`monoMz + k·1.00335/z`) for every detected charge state. Misaligned dotted lines vs. real
     peaks = a mis-placed monoisotope / wrong charge — the diagnostic you want.

## Two ways to get features

1. **Load a TSV** produced by the CLI runner (`Load features…`).
2. **Run feature finding in-app** (`Run feature finding`) — the newest vendored `flashlfq-core`
   runs the top-down pipeline (detect → multi-envelope refine → apex consensus) in-process on the
   open dataset, streaming detect/refine/resolve progress, and feeds the results straight into the
   same overlay/drawer path. No TSV round-trip. Charge is capped at 25 (`maxCharge`) — the guardrail
   that avoids the full-range detect crash — and `minSeedIntensity` trades feature count vs runtime.

## Commands / contract added

- Rust: `load_features(path) -> Vec<Feature>` (parser `parse_resolved_features`, unit-tested);
  `run_feature_detection(handle, options) -> Vec<Feature>` (`run_topdown_pipeline`, integration-tested
  on the real Jurkat raw).
- TS: `Feature` / `PerChargeMz` in `contract.ts`; `loadFeatures` + `runFeatureDetection` + isotope-grid
  helpers in `features.ts`.
- plot-adapter: optional `featureRug` + `regions` on `TicPlot`, `envelope` on `SpectrumPlot`.

## Run

```
cd apps/desktop
bun tauri dev
```

Generate a resolved feature TSV from the sibling `flashlfq-rust` repo:

```
cd rust/flashlfq-core
TOPDOWN=1 cargo run --release --example detect_features_tsv -- <file.raw> <out.tsv>
```

## Limitations (current)

- **MS1 only** — the vendored reader discards MS2, so there is no precursor→fragment navigation yet
  and the served TIC is summed centroided MS1 intensity (not the vendor instrument TIC).
- **Charge cap 25** — the in-app detect runs at `maxCharge=25`; the full 1–60 range crashes the
  detector (a flashlfq-rust issue). Raise it once that's fixed.
- **No IsoDec in-app** — `run_feature_detection` uses the detector's own charge calls; the IsoDec
  neural charge re-assignment (a `+10pp` strict-recall step in the CLI runner) is not yet ported into
  the in-app path. The CLI TSV route still gets it.
