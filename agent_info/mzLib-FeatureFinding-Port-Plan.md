# Porting ProteoStar feature finding into mzLib — plan

Plan for bringing ProteoStar's untargeted MS1 feature detector into
`C:\Users\Alex\Source\Repos\mzLib\mzLib` as C#. Written 2026-07-23, after a full two-sided code
map (ProteoStar `crates/flashlfq-core/src` ↔ mzLib `FeatureFinding` branch).

## TL;DR

**This is a reconciliation, not a greenfield port.** mzLib is already on a `FeatureFinding` branch
that contains a substantially-complete C# feature-finder, *and* it already has every primitive
ProteoStar originally ported *from* it. The genuinely novel ProteoStar work is four things:

1. **The matched-filter detector** (`trace_kernel::detect_features`) — a sparse 2-D matched filter
   (isotope comb × RT Gaussian) with cross-charge non-max suppression and an intensity-descending,
   claim-as-you-go, coverage-driven greedy loop. It scores charge hypotheses directly against the
   *raw m/z* peak index and only deconvolutes *accepted* features. This is a different, more scalable
   seeding strategy than mzLib's current "deconvolute every scan, then trace" approach.
2. **Detect-then-refine** — refine only accepted features on a spectral-averaged composite.
3. **Cross-charge mono resolution** (`resolve_charge_state_consensus`) — the true monoisotopic mass
   appears in *every* charge state's candidate set, so cross-charge voting fixes the ±1 Da off-by-one.
4. **A decoy-envelope FDR framework** — dual-axis decoys (weight / position / RT) threaded end-to-end,
   with 2-D/3-D cosine scoring, giving a calibrated target–decoy FDR.

**Recommendation: hybrid.** *Reuse* mzLib's mature downstream (`MassTrace` → `TraceCorrector` →
`MassFeature`/`MassFeatureBuilder` → `DeconvolutionScorer` → `Ms1Feature`/`_ms1.feature` writer) and
*port* the ProteoStar matched-filter detector as an alternative **front-end** plus the decoy-FDR layer.
Add the one thing both sides lack: a **public orchestrator + parameters object**. Bottom-up first;
top-down multi-charge and IsoDec are later, optional phases.

## What exists where (ground truth from the code map)

### Already in mzLib (`FeatureFinding` branch) — reuse, do not re-port

Primitives ProteoStar took *from* mzLib (so they're native C# already):

| ProteoStar module | mzLib home |
|---|---|
| `deconvolution.rs` (ClassicDecon + Averagine) | `MassSpectrometry/Deconvolution/Algorithms/ClassicDeconvolutionAlgorithm.cs`, `.../AverageResidue/Averagine.cs`, `Deconvoluter.cs`, `IsotopicEnvelope` |
| `spectral_averaging.rs` | `SpectralAveraging/Algorithms/SpectraAveraging.cs` (`AverageSpectra`/`MzBinning`) |
| `peak_indexing.rs` (raw m/z index) | `MassSpectrometry/PeakIndexing/` — `IndexingEngine<T>`, `IndexedMassSpectralPeak`, `PeakIndexingEngine` (FlashLFQ), `MassIndexingEngine` (deconvoluted), `ExtractedIonChromatogram` (+ `CutPeak`) |
| `tolerance.rs` | `MzLibUtil/Tolerance.cs` (`PpmTolerance`, `AbsoluteTolerance`) |
| `isotopic_distribution.rs`, chemistry | `Chemistry/IsotopicDistribution.cs`, `ChemicalFormula`, `PeriodicTable`, `Constants.C13MinusC12` |

An in-progress C# feature pipeline that parallels ProteoStar's *downstream*:

- `MassSpectrometry/Deconvolution/Consensus/MassTrace.cs`, `MassTraceBuilder.cs` — per-scan
  ClassicDecon → **charge-locked, scan-order greedy** mass traces. (≈ ProteoStar per-charge traces,
  but a *different seeding order* — scan order, not intensity-descending matched filter.)
- `.../Consensus/TraceCorrector.cs` — weighted-median consensus mass, ¹³C off-by-one rescue,
  resolvable-species splitting. (≈ ProteoStar `refine` mono-correction, per trace.)
- `.../Consensus/MassFeature.cs`, `MassFeatureBuilder.cs` — cross-charge union-find stitch into a
  `MassFeature`. (≈ ProteoStar `resolve_charge_state_consensus` → `ResolvedFeature`.)
- `.../Scoring/EnvelopeScoreFeatures.cs`, `DeconvolutionScorer.cs` — logistic envelope scorer.
  **Coefficients are provisional / "RECALIBRATE AGAINST LABELLED DATA".**
- `Readers/.../Ms1Feature.cs`, `Ms1FeatureFile.cs`, `Ms1FeatureFromMassFeatureExtensions.cs` —
  `MassFeature` → TopFD-style `_ms1.feature`, round-trip tested.
- A **second front-end already exists** as a test: `JurkatMs1FeatureGeneratorViaXic.cs` uses
  `MassIndexingEngine.GetAllXics(...)` + `XicToMassTrace` → identical downstream. Proof the trace
  layer is front-end-agnostic.
- `DeconvoluteWithDecoys(...)` + `DecoyAveragine` exist, but not the full decoy-comb FDR.

**Gaps in mzLib's current pipeline** (from the map): no public orchestrator/params object (the whole
thing only exists wired inside `[Explicit]` test drivers); no threading/progress/cancellation;
`MassTrace` **drops per-envelope peak lists** (forces scoring inline at generation); scorer not
calibrated; the primary path seeds in scan order (deconvolute-everything), not intensity-descending.

### Novel in ProteoStar — the actual port work

| ProteoStar module | LOC | What it adds over mzLib |
|---|---:|---|
| `trace_kernel.rs` | 4181 | Matched-filter detector, cross-z NMS, coverage/knee/reject stop, Poisson+averagine comb weights, greedy claim-as-you-go. **The core deliverable.** |
| `feature_refinement.rs` | 2263 | Detect-then-refine (composite decon), `resolve_charge_state_consensus` cross-charge mono voting, `NeighborIndex` masking. |
| `isotope_shift_decon.rs` | 1303 | FlashLFQ-style −1/0/+1 shift placement / recharge (a refine variant). |
| `joint_fit.rs` | 380 | NNLS multi-envelope joint fit (a refine variant). |
| `deconvolution.rs` decoy funcs | — | `EnvelopeModel` + exotic decoy envelopes (Cl/Fe/B, shuffled, rotated, hybrid) + `bin_envelope_to_comb_weights`. |
| `feature_export.rs` | 205 | `.msalign` writer (mzLib can already *read* this format via `Ms1Feature`). |
| top-down multi-charge | (in trace_kernel/refine) | `detect_features_multicharge`, joint charge-ladder scoring, IsoDec charge NN (`isodec.rs`). |

## Recommended architecture — the seam

The clean insight: ProteoStar's per-charge `RefinedFeature` ≈ mzLib's `CorrectedTrace`, and
`ResolvedFeature` ≈ mzLib's `MassFeature`. So the ProteoStar detector can be dropped in **above**
mzLib's existing consensus/output layer:

```
                 ProteoStar contribution                 |     reuse mzLib as-is
  raw MS1 scans                                           |
    → PeakIndexingEngine (m/z index; already in mzLib)    |
    → [NEW] MatchedFilterDetector.DetectFeatures()        |   ← port trace_kernel (serial path)
         emits per-charge DetectedTrace (mass,z,rt,PEAKS) |
    → [NEW] FeatureRefiner.Refine()                       |   ← port detect-then-refine;
         composite via SpectraAveraging + Deconvoluter    |     uses mzLib ClassicDecon + averaging
         → CorrectedTrace-compatible                      |
                                                          |   → MassFeatureBuilder.BuildFeatures()
                                                          |   → DeconvolutionScorer (calibrate w/ FDR)
                                                          |   → Ms1FeatureFile.FromMassFeatures().Write()
```

Two front-ends then coexist behind one orchestrator + params object: (a) mzLib's existing
scan-order/`GetAllXics` trace builder, (b) the ProteoStar matched filter. Benchmark them against the
same `_ms1.feature` output and the same ground truth.

**Critical: retain peaks.** ProteoStar's `DetectedFeature.peaks: Vec<IndexedMassSpectralPeak>` carries
the claimed peaks forward; mzLib's `MassTrace` drops them. Port the detector to keep peaks on the
trace so refinement and scoring don't have to re-extract — this also removes mzLib's "score must run
inline at generation" limitation.

## Phased plan

Each phase is independently landable and testable. Phases 1–4 are the bottom-up MVP.

### Phase 0 — Reconcile & spike (small)
- Read the four mzLib Consensus files + the two Jurkat generators end-to-end; run
  `TestConsensusToMs1FeatureEndToEnd` locally to confirm the existing pipeline works and to capture a
  baseline `_ms1.feature` + timing on a shared test raw.
- Decide the seam types: extend `MassTrace` to carry peaks, or add a parallel `DetectedTrace`. (Recommend
  extend `MassTrace` with an optional `List<IIndexedPeak> Peaks` so both front-ends and the writer are shared.)

### Phase 1 — Public orchestrator + parameters (small–medium) — *value even without the detector*
- New `MassSpectrometry/FeatureFinding/` namespace. Add `Ms1FeatureFinder` (public orchestrator) and
  `Ms1FeatureFindingParameters` (deconv ppm/ratio, charge range, trace tol/gap, cross-charge ppm,
  coverage target, seed floor, front-end selector, thread count).
- Wire the *existing* pipeline behind it (the thing today's `[Explicit]` tests hand-wire): scans →
  per-scan decon → traces → correct → features → write. Add `Parallel.ForEach` over independent
  features for refine/score. This alone gives mzLib a usable public API.

### Phase 2 — Port the matched-filter detector (large) — *the core*
- Port `PpmTolerance` semantics are already present; add the two `PeakIndexingEngine` methods the
  detector needs and mzLib lacks: a **multi-tooth scan-range lookup** (`get_indexed_peaks_in_scan_range`)
  and **`tallest_unclaimed_in_mz_at_scan`**, plus a **claimed-peak `HashSet`** mechanism. `IIndexedPeak`
  already exposes `float M` — matches ProteoStar's f32 keys, so bit-compatible.
- Port `trace_kernel::detect_features_serial` (the canonical, deterministic path) as
  `MatchedFilterDetector`: comb offsets (`tooth_offsets`), Poisson + averagine comb weights
  (`Averagine` already precomputes the table — reuse it; add `AveragineCombWeights`), `score_hypothesis`,
  cross-z NMS, greedy claim-as-you-go, coverage/knee/reject/seed-floor stops. Emit `DetectedTrace`
  (mass, charge, mono m/z, apex, RT bounds, summed intensity, score, **peaks**).
- **Defer the parallel tile2d/red-black paths** — they are *not* bit-identical to serial (boundary
  claim-order). Ship serial first; optionally add `Parallel.ForEach` over the 4 colour phases later.
- Unit tests mirror ProteoStar's: synthetic z=2 envelope → correct mass/charge; z=2 beats z=1/z=4;
  single-detection via claiming; Poisson mode shift at heavy mass; scan-spacing/window derivation.

### Phase 3 — Detect-then-refine into MassFeature (medium)
- Port `refine_feature`: for each `DetectedTrace`, slice the local m/z neighbourhood over apex±N scans,
  `SpectraAveraging.AverageSpectra` → `Deconvoluter.Deconvolute(ClassicDeconvolutionParameters)` → pick
  the same-charge envelope nearest the seed; surface `monoisotopic_mass_predictions` (ProteoStar added
  this to its `DeconEnvelope`; mzLib's `IsotopicEnvelope` needs the equivalent field, populated in
  `SetMedianMonoisotopicMass`). Cap the averaging window (`MAX_SCANS_TO_AVERAGE`, ~3–7).
- Emit `CorrectedTrace`-compatible output so `MassFeatureBuilder.BuildFeatures` and the `_ms1.feature`
  writer are reused unchanged.

### Phase 4 — Cross-charge mono resolution (medium)
- Either (a) reuse `MassFeatureBuilder`'s union-find stitch and rely on `TraceCorrector` for
  off-by-one, or (b) port ProteoStar's `resolve_charge_state_consensus` /
  `resolve_mass_by_cross_charge` (pool candidate masses across grouped charges, single-linkage cluster,
  rank by distinct-charge support). Port (b) as an alternative resolver and benchmark against (a) — the
  cross-charge voting is ProteoStar's answer to the ±1 Da problem and may beat per-trace correction.
- Validation gate: PSM re-discovery recall on a shared raw (mass+RT+charge), matching ProteoStar's
  Benchmarking-Guide method.

### Phase 5 — Decoy-envelope FDR + scorer calibration (medium–large)
- Port `EnvelopeModel` + the decoy weight builders (`bin_envelope_to_comb_weights`, chlorinated,
  shuffled, rotated, Cl/Fe custom, hybrid) into mzLib's decon area (alongside `DecoyAveragine`). Thread
  the model through detection **and** refinement (the "laundering fix": a decoy-detected feature must be
  refined/scored against the *same* decoy envelope).
- Port the 2-D (apex) and 3-D (m/z×RT) cosine scoring from `decoy_score_export.rs`. Use the resulting
  target–decoy score distributions to **calibrate `DeconvolutionScorer`'s provisional coefficients**
  against real labelled data — directly closing mzLib's biggest open scoring gap.

### Phase 6 — Top-down multi-charge + IsoDec (large, optional)
- Port `detect_features_multicharge` (score a neutral mass across the whole charge ladder),
  `joint_fit` (NNLS), and the shift-decon refine variants. mzLib already has an IsoDec algorithm
  (`IsoDecAlgorithm.cs`, native `isodeclib.dll`) — prefer wiring that over porting `isodec.rs`. This is
  the path for the Jurkat top-down datasets; skip unless top-down is a target.

## Validation strategy (keep the project's parity discipline)

- **Cell-parity where a golden exists.** ClassicDecon, SpectraAveraging, Averagine, PeakIndexing are
  already mzLib code — no re-porting, no re-golden. The *ported* pieces (comb weights, matched-filter
  response, consensus) get **unit tests mirroring ProteoStar's Rust `#[test]`s** (synthetic envelopes,
  hand-computed expectations).
- **Cross-implementation diff.** Run ProteoStar's detector and the C# detector on the same raw at the
  same parameters; diff detected `(mass, z, apex_rt, intensity)`. Not expected bit-identical (f32
  ordering, tie-breaks) but should agree within tolerance on the strong features.
- **Real-data gates** (from `agent_info/Benchmarking-Guide.md`): PSM re-discovery recall, ΣTIC coverage,
  and (Phase 5) target–decoy AUC. Record results in the guide's §5 table.

## Where the code lives

**`MassSpectrometry` project** (agent-confirmed): the whole existing pipeline lives there and it depends
only on `Chemistry` + `MzLibUtil`, so `Readers`/`FlashLFQ`/MetaMorpheus can consume it. New code under
`MassSpectrometry/FeatureFinding/` (detector, refiner, orchestrator, params) and extend
`Deconvolution/Consensus/` + `Deconvolution/Scoring/`. **Not** FlashLFQ (wrong dependency direction —
entangles with ID-anchored quant). Tests under `Test/MassSpectrometryTests/Deconvolution/` (NUnit 4,
`net8.0-windows`), following the existing `[Explicit]` driver pattern for data-dependent runs.

## Risks & gotchas

- **f32 masses.** ProteoStar casts to f32 to mirror C# `(float)`; mzLib's `IIndexedPeak.M` is already
  `float`. Keep peak keys on `float` for compatible dedup.
- **Determinism.** Only the *serial* detector is deterministic; the rayon tile/red-black paths are not
  bit-identical. Port serial first; treat parallelism as an opt-in perf phase.
- **Averagine already precomputed.** mzLib's `Averagine` builds its 1500-envelope table in a static
  ctor — reuse it directly instead of porting ProteoStar's `LazyLock` table. Add comb-weight accessors
  on top (`GetAllTheoreticalIntensities` already gives the envelope).
- **Peak retention.** Must extend `MassTrace`/the new `DetectedTrace` to carry peaks, or the scorer
  stays crippled (mzLib's current limitation).
- **`MassIndexingEngine` vs the matched filter.** mzLib's `MassIndexingEngine.IndexPeaks` deconvolutes
  *every* scan up front — the exact cost the ProteoStar detector avoids. Do not route the matched filter
  through `MassIndexingEngine`; it consumes the raw-m/z `PeakIndexingEngine`.
- **Rayon → `Parallel.ForEach`.** Refine/score fan out trivially (features independent). The detector's
  greedy claiming is inherently sequential in the serial path.

## Open decisions (recommend confirming before Phase 2)

1. **Extend vs. parallel front-end.** Recommend: reuse mzLib downstream, add the matched filter as a
   selectable front-end (not a rewrite). Confirm this is the intent vs. a full swap.
2. **Scope.** Recommend bottom-up MVP (Phases 1–4) first; decoy-FDR (5) next since it fixes mzLib's
   uncalibrated scorer; top-down (6) only if Jurkat/top-down is a target.
3. **Consensus resolver.** Reuse `MassFeatureBuilder` + `TraceCorrector`, or also port
   `resolve_charge_state_consensus`? Recommend porting it as an alternative and benchmarking — it is
   ProteoStar's distinctive answer to the ±1 Da off-by-one.
4. **Upstream intent.** Is this destined for a mzLib PR (public API, CI, docs) or a local fork? Affects
   how much Phase-1 API polish / test-coverage is required.
