# Benchmarking Guide — untargeted MS1 feature-detection recall

How the recall numbers (**98% / 92% / 94%** on the 10-min / 65-min / 2-hr files)
are computed, end to end, and how to reproduce them. Recall = the fraction of
PSM-identified peaks that the untargeted detector independently rediscovered
(mass + retention-time match).

Last run: 2026-07-08. All artifacts live under `analysis/dinosaur_bench/`
(scorer, ground-truth builders) and `analysis/longer_gradients/` (65-min/2-hr
ground-truth builder). This is analysis-branch material.

## Top-down benchmark (separate)

This guide covers the **bottom-up** recall benchmark. The **top-down** proteomics adaptation has its own
ground truth, tolerances, and `TOPDOWN=1` preset — see `analysis/topdown_bench/RESULTS.md` (builder
`build_topdown_gt.py`, scorer `score_topdown.py`). Headline: Jurkat 84.2% strict / 99.9% isotope-tolerant
against the high-confidence Classic∩IsoDec intersection GT. Design: `agent_info/TopDown-Deconvolution-Plan.md`.

## 0. TL;DR — one-line reproduction

```powershell
$env:APEX_PREGATE="0"   # REQUIRED for the baseline — see the caveat in §2
$bin = "rust/target/release/examples/detect_features_tsv.exe"
& $bin <raw> <out.tsv>                               # detect + refine + resolve
python analysis/dinosaur_bench/score_any.py <out.tsv> <ground_truth.tsv> <rt_delta...>
```

## 1. Test files, ground-truth sources, and RT tolerance

Full paths for every input. All raws and PSM references live on the `D:` data
drive; the normalized ground-truth tables they are converted into live in the repo
(see §3).

### 10-min (CA/Lumos) — 555 ref peaks, RT = real apex → ±0.5 min
- **Raw:** `D:\SP_Tutorial\Lumos\04-17-23_CA_Tryp_HCD_10min.raw`
- **Reference (FlashLFQ AllQuantifiedPeaks):**
  `D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv`

### 65-min (D_Morgen glyco) — 10,021 ref peaks, RT = PSM scan RT (not apex) → lenient ±0.5 / ±1.0 min
- **Raw:** `D:\D_Morgen_Glyco\RawFiles\HFX_MB_14751_5_02062022.raw`
- **Reference (Byonic N-glyco PSM table):**
  `D:\D_Morgen_Glyco\Byonic and MSFragger ID\ConvertedBionicResults\HFX_MB_14751_5_02062022.raw_20230206_WC_Ngly_Generic.tsv`

### 2-hr (IonStar) — 24,642 ref peaks, RT = real apex → ±0.5 min
- **Raw:** `D:\PXD003881_IonStar_SpikeIn\B03_19_150304_human_ecoli_B_3ul_3um_column_95_HCD_OT_2hrs_30B_9B.raw`
- **Reference (FlashLFQ QuantifiedPeaks, per-file):**
  `D:\PXD003881_IonStar_SpikeIn\MM_Master_noMBR_Norm\Task1-SearchTask\Individual File Results\B03_19_150304_human_ecoli_B_3ul_3um_column_95_HCD_OT_2hrs_30B_9B_QuantifiedPeaks.tsv`

The RT delta differs by source: the 65-min Byonic reference records the *scan RT at
which the PSM was collected*, not the chromatographic apex, so recall there is
scored with a lenient window. The 10-min and 2-hr references carry a real `Peak RT
Apex`, scored tighter. (The 2-hr `QuantifiedPeaks.tsv` pools many files — the
builder in §3 filters it to this raw's `File Name` and drops decoys.)

## 2. Stage 1 — run the detector (produces the feature table)

Binary: `rust/flashlfq-core/examples/detect_features_tsv.rs`, built with
`cargo build --release --example detect_features_tsv`. It runs the full pipeline —
`read MS1 (.raw or mzML)` → `detect_features` (trace kernel) → `refine_feature`
(shift-apex + envelope-fit recharge) → `resolve_charge_state_consensus` — and
writes three TSVs; recall scoring uses the **resolved** one (`<out>.tsv`).

Command: `detect_features_tsv.exe <raw> <out.tsv>`

Relevant defaults (all baked into the example / library defaults):
- `MIN_SEED_INTENSITY = 1000` (the "threshold=1000" baseline; env-overridable)
- `COVERAGE_TARGET = 1.0` (uncapped → 2-D tiling detector engages)
- 10 ppm comb tolerance, charges 1..=6, data-driven σ_RT from measured FWHM

> **CAVEAT — `APEX_PREGATE` must be 0 for these numbers.** The apex-charge pre-gate
> (prune charges lacking an apex isotope neighbour) is an in-progress feature whose
> default is `coverage_target >= 1.0` (example line ~411) — i.e. it turns **ON** in
> the standard uncapped config. With it ON, recall drops (65-min 92.2→90.4%, 2-hr
> 93.8→93.0%) and feature counts fall ~30%. The 98/92/94 baseline is the pre-gate
> **OFF** path, so set `APEX_PREGATE=0` explicitly until the default is fixed to be
> opt-in (default off).

The resolved TSV columns the scorer reads: `Monoisotopic Mass`, `RT Apex`,
`Charge States` (semicolon-separated list).

## 3. Stage 2 — build the ground truth (normalized peak table)

All three references are adapted to one schema — `mono_mass  mz  charge  rt
intensity  detection` — one row per expected peak:

- **10-min**: `analysis/dinosaur_bench/build_gt_10min.py` → `gt_10min.tsv`.
  Drops decoys, dedupes to unique (`Full Sequence`, `Peak Charge`); rt = `Peak RT Apex`.
- **65-min & 2-hr**: `analysis/longer_gradients/build_ground_truth.py` →
  `medium_ground_truth.tsv`, `long_ground_truth_msms.tsv`.
  - Medium: dedupe Byonic PSMs to unique (`Full Sequence`, `Precursor Charge`);
    mz from `Peptide Monoisotopic Mass` + charge; rt = `Scan Retention Time`.
  - Long: filter `QuantifiedPeaks` to this raw's `File Name`, drop
    `Decoy Peptide == True`, keep MSMS-detected (primary truth); rt = `Peak RT Apex`.

## 4. Stage 3 — score recall (`analysis/dinosaur_bench/score_any.py`)

For every ground-truth peak `(ref_mass, ref_rt, ref_charge)`, the scorer looks for a
detected feature that matches on BOTH mass and RT:

- **mass match**: `|feat_mass - ref_mass| / ref_mass * 1e6 <= 20 ppm`
  (features pre-sorted by mass; the ppm window is found by binary search / `bisect`)
- **RT match**: `|feat_rt_apex - ref_rt| <= rt_delta`
- first feature satisfying both counts the ref as **matched** (first-hit, not best-hit)
- **charge match** (secondary): the ref's charge is among the matched feature's
  `Charge States` list

Then:

```
recall        = (# refs matched)        / (# ref peaks)
charge-recall = (# refs matched w/ charge) / (# ref peaks)
```

`score_any.py` auto-detects the table format from the header, so the identical scorer
also scores Dinosaur output (`mass`/`rtApex`/`charge`) for the cross-tool comparison
(see `RESULTS.md`). Note our resolved features aggregate charge states into one row,
so the charge-match metric is not apples-to-apples vs one-row-per-charge tools.

## 5. Results — baseline (`MIN_SEED_INTENSITY=1000`, `APEX_PREGATE=0`), 2026-07-08

Full 24-thread runs (`rust/target/release`), scored as above.

| File | detected | resolved | %ΣTIC explained | recall (tight) | recall (lenient) | runtime |
|------|---------:|---------:|:---------------:|:--------------:|:----------------:|:-------:|
| **10-min** | 143,132 | 125,919 | 90.0% | **98.2%** (±0.5) · 96.0% (±0.2) | — | ~15 s |
| **65-min** | 974,162 | 840,344 | 90.2% | 92.2% (±0.5) | **95.6%** (±1.0) | ~63 s |
| **2-hr** | 1,116,158 | 930,499 | 96.5% | **93.8%** (±0.5) | 96.3% (±1.0) | ~128 s |

Charge-match (ref charge within the feature's charge list): 10-min 86.3% (±0.5),
65-min 77.4% (±0.5), 2-hr 75.6% (±0.5).

### Exact reproduction commands

```powershell
$env:APEX_PREGATE="0"
$bin = "F:\flashlfq-rust\rust\target\release\examples\detect_features_tsv.exe"
$sa  = "F:\flashlfq-rust\analysis\dinosaur_bench\score_any.py"

& $bin "D:\SP_Tutorial\Lumos\04-17-23_CA_Tryp_HCD_10min.raw" out10.tsv
python $sa out10.tsv "F:\flashlfq-rust\analysis\dinosaur_bench\gt_10min.tsv" 0.2 0.5

& $bin "D:\D_Morgen_Glyco\RawFiles\HFX_MB_14751_5_02062022.raw" out65.tsv
python $sa out65.tsv "F:\flashlfq-rust\analysis\longer_gradients\medium_ground_truth.tsv" 0.5 1.0

& $bin "D:\PXD003881_IonStar_SpikeIn\B03_19_150304_human_ecoli_B_3ul_3um_column_95_HCD_OT_2hrs_30B_9B.raw" out2h.tsv
python $sa out2h.tsv "F:\flashlfq-rust\analysis\longer_gradients\long_ground_truth_msms.tsv" 0.5 1.0
```
