# Feature-scoring investigation — findings (2026-07-22)

Goal: cut the untargeted detector's 10–20× feature over-count (vs Dinosaur et al.)
without losing PSM recall, by finding a per-feature score that ranks real features
above the junk tail. Branch: `feature-scoring`. Companion branch: `noise-floor-reduction`.

Method: `score_tradeoff.py` ranks resolved features by a score column, keeps the top
fraction, and scores recall against the PSM ground truth (same ±20 ppm / RT-window logic
as `score_any.py`). A good score holds recall high as the retained fraction shrinks.
Baselines run with `APEX_PREGATE=0` per the Benchmarking-Guide caveat.

## Scores evaluated

Existing resolved-TSV columns (all intensity-correlated): Summed Intensity, Num Charge
States, Cross-Charge Support, Num Members. Newly added (this branch), intended to be
intensity-orthogonal: **Decon Score** (primary member envelope-fit cosine), **Min Decon
Score** (weakest charge), **Max Num Isotopes** (envelope completeness), **PPM Spread**
(stddev of per-isotope-tooth mass errors across the apex envelope; 999 sentinel when the
apex has <2 isotopes).

## Headline result: no single score beats Summed Intensity

Recall (rt ±1.0 min) vs retained top-fraction, **65-min glyco** (baseline 95.6%, 840,344 features):

| retained | Intensity | NumCharge | Decon | PPM Spread | MaxIsotopes |
|---------:|----------:|----------:|------:|-----------:|------------:|
| 80% | 95.4 | 95.4 | 94.9 | 95.3 | 94.9 |
| 60% | 94.5 | 94.5 | 93.7 | 94.5 | 90.3 |
| 50% | 93.7 | 93.8 | 93.0 | 93.7 | 80.2 |
| 40% | 92.6 | 92.6 | 91.5 | 91.8 | 68.2 |
| 30% | 91.5 | 91.5 | 89.9 | 87.8 | 56.7 |
| 20% | 89.7 | 89.7 | 87.7 | 78.0 | 45.8 |

10-min (rich, flat-plateau file) tells the same story — intensity dominates the aggressive
regime; PPM Spread ties it only down to ~70% retained then falls off.

- **Summed Intensity and Num Charge States are the best and are tied.** The 10–20× junk is
  dominated by the low-intensity tail, which intensity already ranks best.
- **Decon fit and PPM Spread are close but never better** than intensity.
- **Max Num Isotopes is worst** (integer quantization + real low-abundance peptides have few
  isotopes → cutting by it removes real signal).
- **56% of 65-min features are single-isotope at their apex** (PPM-Spread sentinel). Prime
  junk suspects, but isotope-count gating is too blunt to remove them safely.

## Interpretation

1. Intensity (i.e. the noise floor) is already near-optimal as a **1-D** ranker. A better
   noise floor — specifically the **per-scan local floor** the seed-floor sweep recommended —
   is the highest-value lever, because it makes the winning signal adaptive (cut hard in
   dense regions, gently in sparse ones). → `noise-floor-reduction` branch.
2. PSM recall **cannot** directly measure junk removal (dropped features may be real-but-
   unidentified). Two scores that hold the same recall may differ greatly in how much true
   junk they removed. The **target-decoy** framework (decoy-features branch) is the correct
   evaluation tool for specificity.
3. **Open question the tradeoff can't answer:** do these scores *add* to intensity (a 2-D
   gate: drop only features that are BOTH low-intensity AND poor-fitting)? Single-score curves
   don't test additivity. This is the decoy-calibrated multivariate-scoring question and is
   the productive direction for a scoring system — single new columns are not.

## Semi-supervised SVM rescore (2026-07-22) — decoys as negatives DON'T beat intensity

Tried feature-level Percolator (`svm_rescore.py`): NEGATIVES = decoy-model detections (spacing =
`DECOY_LATTICE=scaled`, weird-averagine = `COMB_MODEL=custom`/`shuffled`), POSITIVES = confident
targets seeded by apex fit, re-selected each iteration by out-of-fold SVM score at a decoy-FDR
threshold, 3-fold CV.

**Result: the SVM score ranks targets WORSE than plain Summed Intensity at every retained
fraction** (10-min, e.g. 40% kept: SVM 95.9% vs intensity 97.5%; 20% kept: 82.9% vs 93.5%).

**Why (diagnosed, not guessed):**
- The semi-supervised loop collapses — every target scores above the decoys (q≈0), so the positive
  set becomes *all* targets. The SVM learns "target-like vs decoy-like", not "real vs junk".
- The learned model's standardized weights lean on `ppm_spread` (+0.30) and `decon_score` (+0.24)
  and assign `log_intensity` **+0.008 — essentially zero**. It *discards* intensity.
- Root cause: **the decoys have the same intensity distribution as targets** (log-intensity median
  target 14.12 vs shuffled 14.18 / custom 13.90 / spacing 14.50). A decoy detection claims real
  peaks, so it is a *normal-intensity* feature that fits a *wrong* template. Target junk is the
  opposite: a *low-intensity* feature that fits the *right* template. So intensity — the true
  junk axis — carries no target-vs-decoy signal and the objective zeroes it out, forcing the model
  onto shape features we already know underperform intensity.

**Conclusion:** these decoys model the wrong failure mode ("wrong envelope shape", normal
intensity), not the actual junk ("low-abundance noise coincidence", real envelope shape). Any
classifier trained against them is structurally barred from using intensity and so cannot beat it.

**The fix:** a *noise-faithful* decoy — the **real** averagine comb placed at m/z positions
**offset off the true isotope grid** so it can only claim noise. Its detections would be
low-intensity, low-persistence — the same profile as target junk — so intensity would separate
target-real from decoy-noise and the SVM could keep intensity AND add shape refinements. The
`decoy-features` branch reportedly has a "shifted" decoy of this kind; porting/enabling it is the
concrete next step. (A pure m/z offset, not `DECOY_LATTICE=scaled`, which lands teeth near real
peaks and inherits normal intensity.)

## Artifacts
- `svm_rescore.py` — semi-supervised SVM rescorer (Percolator-style, 3-fold CV, decoy-FDR).
- `score_tradeoff.py` — the ranker/curve tool.
- New resolved-TSV columns: Decon Score, Min Decon Score, Max Num Isotopes, PPM Spread
  (see `apex_ppm_spread` in `examples/detect_features_tsv.rs`).

## Peak-width cap (asked separately) — not a lever
`trace_max_half_width_minutes` default 0.5 min hard-caps claimable elution width; broad
elutions fragment (unclaimed shoulders re-seed same-m/z features). Measured on 65-min:
widening 0.5→2.0 min cut only **1.7%** of features AND cost **0.7 pp** recall (over-merges
co-eluting distinct peptides). Net-negative; the cap is already reasonably tuned.
