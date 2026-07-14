# Decoy feature detection — envelope decoys for the untargeted MS1 detector

Ported from the `flashlfq-rust` `worktree-decoy-comb-fdr` branch. Adds the ability to run the whole
detect→refine pipeline with a **decoy isotope-envelope model** — a theoretical envelope no real
peptide produces — so that "features" are found on noise/coincidence. Scoring the target run against
the real averagine and each decoy run against *its own* decoy envelope gives the score distributions a
target–decoy FDR (or a classifier) needs.

## The three decoy strategies (as requested)

| strategy | what changes | how it's selected | scored against |
|----------|--------------|-------------------|----------------|
| **shifted spacing (0.94 Da)** | tooth *positions*: the isotope comb is laid at 0.94-Da spacing instead of the ¹³C 1.00336 Da — an off-lattice comb | `DECOY_LATTICE=scaled DECOY_SPACING_SCALE=0.9368` (0.94 / 1.00336) | averagine weights on the 0.9368× lattice |
| **weird averagine (Fe+Cl)** | tooth *weights*: the averagine backbone is **replaced** by the exotic composition `P2 C1 N1 Cl1 Fe1` (⁵⁷/⁵⁸Fe + ³⁷Cl A+2 satellites) | `COMB_MODEL=custom` (`CUSTOM_AVERAGINE` overrides the composition) | the custom Fe+Cl envelope |
| **shuffled envelope** | tooth *weights*: the real averagine weights are randomly permuted (deterministic, `SHUFFLE_SEED`), teeth stay on the ¹³C lattice | `COMB_MODEL=shuffled` (`SHUFFLE_SEED` varies it) | the shuffled envelope |

All three are applied **end-to-end**: the trace-kernel comb (`CombWeightModel` / `LatticeMode`), the
refinement deconvoluter (`EnvelopeModel` + spacing scale), and the reported `Decon Score` all use the
decoy model — so a decoy is not re-anchored onto real peaks during refinement (the "laundering" fix).

Extra models also ported (not part of the three-way request, available for experiments):
`COMB_MODEL=decoy` (chlorinated averagine), `rotated` (rotate-half), `cbp` (chloro-boro-phosphate),
`hybrid` (Cl=Fe below `HYBRID_MASS`, shuffled above); `DECOY_LATTICE=mixed`; `RT_PROFILE=uniform|inverted`.

## Where the code lives (ProteoStar `crates/flashlfq-core`)

- `src/deconvolution.rs` — `EnvelopeModel` enum + all decoy envelope generators (`decoy_comb_weights`,
  `rotated_*`, `shuffled_*`, `cbp_*`, `custom_*`, `hybrid_*`) and the shared `bin_envelope_to_comb_weights`.
- `src/isotope_shift_decon.rs` — `Deconvoluter` service bound to one `EnvelopeModel` + spacing scale;
  every placement/score/recharge/walk-back method has a model-aware `_model` variant behind it.
- `src/trace_kernel.rs` — decoy `CombWeightModel` variants, `LatticeMode` (Uniform/Scaled/MixedCharge),
  `RtProfile`, `tooth_offsets()`, `rt_weight()`; threaded through `score_hypothesis`,
  `gather_extent_peaks`, `seed_rt_window` (the bottom-up production path). *Caveat:* the top-down
  `detect_features_multicharge` path was left on the physical lattice — the bottom-up benchmark does not
  use it, and the decoy lattice defaults to Uniform.
- `src/feature_refinement.rs` — `refine_feature_shift_with` / `refine_feature_shift_neighbor_with` take
  an explicit `Deconvoluter` so detection and refinement share one model.
- `examples/detect_features_tsv.rs` — env-var knobs above wired into the runner.
- `examples/decoy_score_export.rs` — **the scorer.** For each refined feature computes a 2-D (apex-scan
  m/z) and 3-D (m/z × RT window) union-grid envelope-fit cosine against a given `<model> <spacing_scale>`.

## Scoring — 2-D vs 3-D cosine

Both are the uncentered cosine `Σ T·O / (‖T‖·‖O‖)` on the **union grid** (predicted teeth + unexplained
in-window peaks, so a missing tooth hurts completeness and an unexplained peak hurts fraction-explained).

- **2-D**: single apex scan, template `T[k] = w_k` (the detector's `Decon Score` shape).
- **3-D**: over the RT window (±2σ, σ = 0.15 min), template `T[k,s] = w_k · g(s)` with `g` the elution
  Gaussian; observed intensities at every (tooth, scan). The same weighting is laid for target and decoy.

The 3-D score is lower in absolute value (it integrates off-apex noise into ‖O‖) — the question is
whether the RT dimension *separates* target from decoy better than the apex-only 2-D score.

## Reproduce

```powershell
# detection (per file, per model) — APEX_PREGATE=0 for the benchmark baseline
$env:APEX_PREGATE="0"
$det = "F:\ProteoStar\target\release\examples\detect_features_tsv.exe"
& $det <raw> out_target.tsv                                                   # target
$env:DECOY_LATTICE="scaled"; $env:DECOY_SPACING_SCALE="0.9368"; & $det <raw> out_shifted.tsv
$env:COMB_MODEL="custom";   & $det <raw> out_weird.tsv                        # Fe+Cl
$env:COMB_MODEL="shuffled"; & $det <raw> out_shuffled.tsv

# scoring (2-D + 3-D vs the run's own envelope)
$scr = "F:\ProteoStar\target\release\examples\decoy_score_export.exe"
& $scr <raw> out_target.refined.tsv   s_target.tsv   averagine 1.0
& $scr <raw> out_shifted.refined.tsv  s_shifted.tsv  averagine 0.9368
& $scr <raw> out_weird.refined.tsv    s_weird.tsv    custom    1.0
& $scr <raw> out_shuffled.refined.tsv s_shuffled.tsv shuffled  1.0

# analysis + plots (ID-matches targets to the recall GT, overlays distributions, reports AUC)
python analysis/decoy_features/plot_decoys.py
```

Full orchestration: `analysis/decoy_features/plot_decoys.py` + the run script used to produce the
figures under `analysis/decoy_features/plots/`.

## Top-down extension (2026-07-14) — decoys threaded into the multi-charge detector

The three decoys were extended to the **top-down** joint charge-ladder detector, which the original port
had left on the physical averagine lattice. `score_mass_hypothesis`, `refine_mono_offset`,
`gather_charge_extent`, `emit_mass_hypothesis`, and the seed spacing screen in `detect_features_multicharge`
now use model-aware weights (`multicharge_weights`, a mono-keyed dispatch analogous to `comb_weights`) and
`tooth_offsets` / `tooth_step` under `params.lattice_mode`. Averagine/`Uniform` (the target path) is bitwise
unchanged. `decoy_score_export` gained `SCORE_MAX_ISOTOPES` / `SCORE_RT_SIGMA` (heavy proteoforms span
>24 teeth and elute broader than 0.15 min). Full write-up + tables + plots:
`analysis/decoy_features_topdown/` (`report.md`, `analyze_topdown_decoys.py`, `run_topdown_decoys.ps1`).

Datasets: **Golden** (`golden.raw`, GT `gt_golden_all.tsv`, 930) and **Jurkat**
(`02-18-20_jurkat_td_rep1_fract6.raw`, GT `gt_jurkat_intersect.tsv`, 1541 hi-conf). Config isolates the
envelope model: `TOPDOWN=1 APEX_PREGATE=0 REFINE_METHOD=shift_apex TD_MONO_FIT=0 ISODEC_CHARGE=0` — the
multi-envelope refine, averagine mono-fit, and IsoDec are OFF because they re-anchor a decoy onto the real
averagine (laundering); `shift_apex` refine and the multi-charge `MC_MONO_KMAX` refine are decoy-aware.

Results (ROC-AUC, ID-matched target vs decoy):

| dataset | shifted 2-D/3-D | weird 2-D/3-D | shuffled 2-D/3-D |
|---------|:---------------:|:-------------:|:----------------:|
| Golden | 0.715 / 0.726 | **0.938** / 0.885 | 0.896 / 0.865 |
| Jurkat | 0.694 / 0.669 | **0.922** / 0.785 | 0.899 / 0.795 |

Findings: (1) same decoy-strength ranking as bottom-up (**weird > shuffled > shifted**); at a 5%-decoy-pass
threshold, weird retains ~76–79% of real IDs on both files. (2) **Shifted spacing is a weak top-down decoy**
(AUC 0.69–0.72, worse than bottom-up) — a 0.94-Da comb at high charge barely leaves the ¹³C lattice and the
dense high-charge spectra still hit real peaks by coincidence. (3) **2-D > 3-D on both files** (wide-window
3-D folds in off-apex noise). (4) Real-feature cosines are lower on Jurkat (heavy 13–21 kDa envelopes), so
thresholds are lower and file-specific — calibrate per file.

**Threshold terminology (correction):** the two reported operating points are **95% recall** (5th pct of
ID-matched target scores) and **95% decoy-rejection** (95th pct of decoy scores → only 5% of decoys pass).
The latter is a 5% false-positive / 95%-specificity cut (a decoy-based critical value at α=0.05), **not**
classifier precision `TP/(TP+FP)`; the two diverge sharply here because decoy features outnumber ID-matched
targets ~10:1.

## Results (2026-07-14, `decoy-features` branch)

Full 3-file × 4-model sweep, `APEX_PREGATE=0`, `MIN_SEED_INTENSITY=1000`, `COVERAGE_TARGET=1.0`. Target
feature counts reproduce the recall baseline exactly (65-min 974,162; 2-hr 1,116,158), so the target
path is unchanged. Targets split into ID-matched (mass ≤ 20 ppm + RT within the file's tolerance to a
reference PSM/quantified peak) vs unmatched. ROC-AUC = ID-matched-target (positive) vs decoy (negative);
0.5 = no separation, 1.0 = perfect. Figures: `analysis/decoy_features/plots/decoy_scores_{10,65,2h}.png`.

| file | ID-matched | decoy | AUC 2-D | AUC 3-D |
|------|-----------:|-------|:-------:|:-------:|
| 10-min | 1,617 (1.1%) | shifted spacing (0.94 Da) | **0.882** | 0.720 |
| 10-min | | weird averagine (Fe+Cl) | **0.921** | 0.795 |
| 10-min | | shuffled envelope | **0.901** | 0.772 |
| 65-min | 23,379 (2.4%) | shifted spacing (0.94 Da) | **0.810** | 0.797 |
| 65-min | | weird averagine (Fe+Cl) | **0.898** | 0.846 |
| 65-min | | shuffled envelope | **0.848** | 0.830 |
| 2-hr | 63,032 (5.6%) | shifted spacing (0.94 Da) | 0.792 | **0.809** |
| 2-hr | | weird averagine (Fe+Cl) | **0.871** | 0.848 |
| 2-hr | | shuffled envelope | 0.838 | **0.841** |

**Takeaways.**
- All three decoys separate real (ID-matched) features from decoys (AUC 0.72–0.92 everywhere): the
  ID-matched targets pile up near cosine 1.0, the decoys near 0.
- **Ranking of decoy strength is consistent:** weird averagine (Fe+Cl) separates best, then shuffled,
  then shifted spacing. The Fe+Cl envelope's A+2-dominated shape is the hardest for real peptides to mimic.
- **2-D vs 3-D:** the apex-only 2-D cosine separates clearly better on the fast **10-min** gradient
  (0.88–0.92 vs 0.72–0.80). On the longer **65-min / 2-hr** gradients the gap collapses and the 3-D score
  ties or slightly *wins* in a few cases (2-hr shifted 0.809 vs 0.792; 2-hr shuffled 0.841 vs 0.838) —
  denser, better-RT-resolved elution gives the RT axis real information, whereas on a compressed gradient
  the wide-window 3-D integration mostly folds in off-apex noise and dilutes the discriminant. Net: the
  RT dimension is not free — it helps only where chromatography is well resolved.
