# Untargeted feature-detection — TODO

Running task list for the untargeted MS1 feature-detection work. Check items off as they land;
add a one-line result/commit note when closing one. Priority tags: **[high] / [med] / [low]**.

See also: `Feature-Detection-Design.md` (source of truth), `Detector-Improvement-Plan.md`,
`Miss-Characterization-and-Discriminator-Plan.md`.

## Test-data reference (medium + long gradients)

Existing case — **short ~10-min** (CA/Lumos): raw `D:\SP_Tutorial\Lumos\04-17-23_CA_Tryp_HCD_10min.raw`,
ref `...\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv` (~620 peaks). ~1.9 s FWHM.

**Medium ~65-min gradient** — `D:\D_Morgen_Glyco`.
- Refs: `Byonic and MSFragger ID\ConvertedBionicResults\*_Ngly_Generic.tsv`. Each ref's filename prefix
  maps to a `.raw` in `RawFiles\`. Many pairs available; start with one.
- Example pair:
  - raw: `D:\D_Morgen_Glyco\RawFiles\HFX_MB_14751_5_02062022.raw`
  - ref: `D:\D_Morgen_Glyco\Byonic and MSFragger ID\ConvertedBionicResults\HFX_MB_14751_5_02062022.raw_20230206_WC_Ngly_Generic.tsv`
- **~7–9 MS1 scans inside each peak's FWHM** (Alex's visual estimate = the ground-truth ballpark for
  scans-to-average here).
- **Gotcha:** ref reports the **scan RT at which a PSM was collected, NOT the peak apex** → use a
  **lenient RT acceptance delta** when scoring recall.
- Ref columns (differ from AllQuantifiedPeaks — needs a format adapter): `File Name`, `Base Sequence`,
  `Full Sequence`, `Peptide Monoisotopic Mass`, `Precursor Charge`, `Protein Accession`,
  `Scan Retention Time`. No apex/intensity columns. 13,179 rows = per-PSM (dedupe to unique
  peptide+charge for a peak-level ground truth).

**Long ~120-min gradient** — `D:\PXD003881_IonStar_SpikeIn` (IonStar spike-in).
- Raws at the top level of the folder. Refs: `MM_Master_noMBR_Norm\Task1-SearchTask\Individual File
  Results\*_QuantifiedPeaks.tsv`.
- Example pair:
  - raw: `D:\PXD003881_IonStar_SpikeIn\B03_19_150304_human_ecoli_B_3ul_3um_column_95_HCD_OT_2hrs_30B_9B.raw`
  - ref: `...\Individual File Results\B03_19_150304_...9B_QuantifiedPeaks.tsv`
- Ref is a standard FlashLFQ QuantifiedPeaks table (has `Peak RT Apex`, `Peak Charge`, `Peak intensity`,
  `Peak MZ`, `Peak Detection Type`, `Decoy Peptide`, `File Name`). **504,784 rows = many files pooled** →
  filter to this raw's `File Name`, drop `Decoy Peptide`, and consider excluding MBR rows via
  `Peak Detection Type`. The runner's existing compare logic mostly works after that filtering.
- **~20 MS1 scans inside each peak's FWHM** (Alex's estimate = the ground-truth ballpark for
  scans-to-average here).

Scans-to-average ground truth by gradient: short ~10-min ≈ 3 (apex±1) · medium ~65-min ≈ 7–9 ·
long ~120-min ≈ 20. The data-dependent formula should reproduce this progression.

## High priority

### Data-dependent averaging window (umbrella goal)
Replace hard-coded `MAX_SCANS_TO_AVERAGE = 3` with a scans-to-average count derived from the measured
chromatographic FWHM / scans-per-peak, so it adapts to gradient length instead of being tuned to
CA/Lumos 10-min. Broken into the tasks below.

**DONE (umbrella):** `MAX_SCANS_TO_AVERAGE = 3` is now only a floor; the operative count comes from
`derived_avg_scans(fwhm_seconds, scan_spacing_seconds)` → `averaging_params.avg_scans`
(`feature_refinement.rs:156, 683`). Validated in `FWHM-Validation.md` / `FWHM-Flow-Audit.md`.

- [x] **[high] Prep the medium + long test cases.** — DONE: medium/long pairs wired into `Benchmarking-Guide.md`; see `Longer-Gradient-Recall.md`.
  Pick one medium pair (D_Morgen_Glyco) and one long pair (IonStar) to start. Build a format adapter
  for each ref (medium `Ngly_Generic` PSM table; long `QuantifiedPeaks` filtered to one file). Dedupe
  the medium ref to a peak-level ground truth. Record Alex's ground-truth scans-to-average per case
  (medium ≈ 7–9).
- [x] **[high] Run the detector on the longer files; record run times.** — DONE: per-stage timings in `Longer-Gradient-Recall.md` / `Detector-Perf-and-Parallelization.md`.
  Full-pipeline runs on the medium (~65-min) and long (~120-min) raws. Capture per-stage timings
  (read/detect/refine/consensus) — these files are much larger than the 10-min case and will stress
  the detector (see perf/parallelization task).
- [x] **[high] Test recall on medium + long.** — DONE: `Long-Recall-vs-Coverage.md` + `Longer-Gradient-Recall.md`.
  Score recall vs each ref. Medium: **lenient RT delta** (ref RT = PSM scan time, not apex). Long:
  filter ref to the single file + non-decoy. Report recall + charge-match like the 10-min case.
- [x] **[high] Test FWHM + n-scans-to-average calculations.** — DONE: `FWHM-Validation.md`.
  Verify the measured FWHM and the derived scans-to-average on each gradient. Sanity-check against
  ground truth (short ~3, medium ~7–9, long ~20 scans/FWHM). This is where the data-dependent formula
  is validated.
- [x] **[high] Ensure FWHM info flows correctly to ALL downstream steps.** — DONE: audited in `FWHM-Flow-Audit.md`.
  Confirm the measured FWHM / scans-to-average is fed consistently into the averaging window, the
  detector RT σ / matched-filter window, and anything else that assumes a peak width — no step left
  reading a stale hard-coded value.

### Detector performance
- [x] **[high] Detector performance + parallelization.** — DONE: `Detector-Perf-and-Parallelization.md`; 2-D tiling is now the shipped default detect path (index packing + scoring-lookup reductions landed).
  Profile the detector (~76–79% of total wall-clock; the real cost, not refinement). Find
  optimizations and design a parallelization strategy — evaluate per-file-on-its-own-thread (likely
  simplest/effective) vs. intra-file parallelism. The longer files above make this urgent.

### Data-dependent seed-intensity floor
- [ ] **[high] Set `MIN_SEED_INTENSITY` data-dependently instead of the hard-coded 1000.** ← **STILL OPEN** (the main remaining high-pri item). Confirmed still `unwrap_or(1000.0)` at `detect_features_tsv.rs:401`; no per-file noise-floor / percentile derivation yet.
  Today the detector walks all peaks tallest-first and stops at a fixed global floor
  (`min_seed_intensity`, default **1000**) — that floor, not `COVERAGE_TARGET`, is what actually
  bounds an "uncapped" run. **1000 is too low for the medium/long files:** it admits a large tail of
  near-noise seeds that are almost all rejected (wasted compute), and it does not adapt per file. Set
  it from the data instead — e.g. off the estimated noise floor (`estimate_noise_floor`) or an
  intensity percentile — mirroring the FWHM-derived averaging window. Evidence from the reject-cap A/B
  (2026-07-08 reject-cap runs): at floor 1000 an uncapped run explained only
  **90.0% / 90.2% ΣTIC** on the 10-min / 65-min (2-hr 96.5%) — the missing ~10% is sub-1000 peaks
  never seeded. The per-tile reject-rate cap is a crude *adaptive* proxy for the same goal (skip the
  noise tail), but a principled per-file seed floor is the real fix. Coordinate with the reject-cap
  default-`frac` decision — they trade off against each other.

## Medium priority

- [x] **[med] Turn averaging off completely via params.** — DONE: `build_feature_slices` now takes an `average_spectra: bool`; when false the composite is not computed at all (`feature_refinement.rs:269, 1046`).
  Today `build_feature_slices` unconditionally calls `average_spectra` (`feature_refinement.rs:738`) —
  even `shift_apex` builds a composite it never reads. Add a real parameter to disable averaging so the
  composite is not computed at all when unused.
- [x] **[med] Change the default to NOT average (apex-only).** — DONE: bottom-up default `refine_method = "shift_apex"` → `use_shift_apex = true` → `average_spectra = false` (`detect_features_tsv.rs:877-897`). No composite built on the default path.
  A/B on CA/Lumos 10-min: apex-only beats the averaged composite — 97.4% vs 96.0% recall, 94.4% vs
  92.3% charge. Make apex-only the default in the library/param defaults. Coordinate with the
  turn-averaging-off item. NOTE: revisit once data-dependent averaging lands — averaging may help on
  the longer gradients (more scans/peak) even though it hurts on 10-min.
- [~] **[med] Supported output formats.** — SCOPED: format decision written up in `Feature-Output-Formats.md`; confirm the exporter is actually wired before checking off.
  Determine feature-output formats to emit (e.g. ms-align files, feature files). Must be supported by
  mzLib; ideally widely used. Pick target format(s) and wire an exporter.
- [x] **[med] Quantification — kick-off.** — DONE (kickoff/scope): `Quantification-Kickoff.md`.
  Scope quantification: choose test cases and benchmark against FlashLFQ classic. (The IonStar long
  case ships full FlashLFQ QuantifiedPeaks — a natural quant benchmark.)

## Low priority

- [x] **[low] Clean up unused / dead code.** — DONE: `Dead-Code-Sweep.md`.
  Sweep `#[allow(dead_code)]` and abandoned experiment paths.
- [x] **[low] Algorithm explainer.** — DONE: `Algorithm-Explainer.md`.
  Write up the basic workflow / how the algorithm works (index → detect → refine → resolve).

## MsViewer desktop app (`apps/desktop`)

GUI/viewer workstream, distinct from the detection-engine tasks above. See
`MsViewer_Architecture.md`. Priority tags below are a first-pass guess — reorder as needed.

- [ ] **[high] Async indexing — view spectra before indexing finishes.**
  Today spectra can't be displayed until `.raw`/mzML indexing completes — indexing blocks the thread.
  Make indexing async so a free thread stays responsive to plot clicks and can render a clicked
  spectrum on demand while indexing is still running.
- [x] **[med] Spectrum plot must fit without scrolling.** — DONE (verified 2026-07-22): `index.html` body-margin
  reset + `#root` 100vh/overflow-hidden; `layout.tsx` grid tracks `1fr`→`minmax(0,1fr)`. Tauri e2e: no page
  scroll (scrollHeight==innerHeight) and spectrum plot bottom within viewport. `msviewer-tasks.spec.ts` Task 1.
- [x] **[med] Persist spectrum zoom across scan changes.** — DONE (verified 2026-07-22): `stepScan` freezes the
  fitted y-range on the first arrow step while x-zoomed (`resolveSpectrumYRange`/`hasPersistedY` in
  `plot-adapter/viewport.ts`); "Reset zoom" refits. Tauri e2e: y-range held across a scan step, x refits on
  reset. `msviewer-tasks.spec.ts` Task 4. (NB: intentionally reverses the earlier "reset on new scan" behavior.)
- [x] **[med] Label m/z + charge for prominent MS1 features.** — DONE (verified 2026-07-22): `annotate.ts`
  (`computePeakLabels`/`estimateCharge`, isotope-spacing Δ≈1.00235/z) → `SpectrumPlot annotations`. Tauri e2e on
  the yeast apex scan renders labels like `356.19 · z1`, `433.74 · z2`. `msviewer-tasks.spec.ts` Task 2.
- [ ] **[high] Annotated MS2 spectra for PSM data (MetaDraw-style) — LARGE.**
  Port MetaMorpheus/MetaDraw annotated-MS2 display: given PSM data, render the MS2 spectrum with
  matched fragment-ion labels. Big feature — substantial work; likely warrants its own design doc and
  a sub-task breakdown before starting.
- [x] **[med] Enable resizing of the PSM panel.** — DONE (verified 2026-07-22): `ResizableDrawer` in `App.tsx`
  adds a draggable left-edge handle to the right-side drawer (`drawerWidth` state, clamp 280–760); the future
  PSM panel reuses it. Tauri e2e drags the handle → controlled width + clamps at both ends.
  `msviewer-tasks.spec.ts` Task 3. **CAVEAT:** no PSM panel exists yet — implemented on the existing drawer as
  the reusable mechanism; confirm this interpretation with the user.
- [ ] **[med] Drag-and-drop file loading.**
  Let the user drop files onto the window to have them recognized and loaded (route by type through
  the same open flow as the file picker — spectra files → dataset, PSM/results files → PSM panel).
- [ ] **[high] Fully implement pinning + associated cross-panel functionality.**
  Pin state exists (independent TIC/spectrum pins) but the cross-panel actions it should gate are not
  wired. Decomposes into at least two tasks:
  - [ ] **Spectrum pinned → RT-region select on TIC averages MS1 scans.** When the spectrum is pinned,
    selecting an RT region on the TIC plot should average all MS1 scans in that region and display the
    resulting averaged spectrum.
  - [ ] **TIC pinned → m/z-region select on spectrum adds an XIC.** When the TIC is pinned, selecting an
    m/z region of the spectrum should add a new XIC (extracted-ion chromatogram) trace to the TIC for
    the selected m/z range.
- [ ] **[high] Unified PSM (peptide-spectral-match) reader interface.** Design captured in
  `PSM-Reader-Interface-Plan.md` (not yet implemented). Define a `PeptideSpectralMatch` abstraction in
  `flashlfq-core` as two layered traits (scan-level; and scan-level + MS2 fragment ions), make the
  existing `.psmtsv` reader produce PSM members, and add readers for other search software
  (MSFragger — `psm.tsv` first, then `.pepXML`). Foundation for the annotated-MS2 / PSM-panel work.
