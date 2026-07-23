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

## intensity x fit 2-D gate (2026-07-22) — BEATS intensity alone

`combo_gate.py`: supervised logistic regression, label = PSM-match, evaluated out-of-fold (fair
generalization test), nested feature sets. Unlike the decoy-SVM, the model KEEPS intensity
(standardized weight ~+1.0) and ADDS envelope fit on top (ppm-spread ~-1.0, decon ~+0.4). Result —
recall at retained top-fraction, `int` = intensity baseline:

10-min (rt ±0.5), baseline 98.2%:

| retained | intensity | int+decon | full |
|---------:|----------:|----------:|-----:|
| 30% | 96.2% | 96.8% | **97.1%** |
| 25% | 94.6% | 96.2% | **96.8%** |
| 20% | 93.5% | 95.9% | **96.0%** |
| 15% | 90.8% | 94.6% | **95.1%** |
| 10% | 84.7% | 91.7% | **92.1%** |

65-min (rt ±1.0), baseline 95.6%: smaller but consistent, e.g. 40% retained 92.6% → 93.2%,
30% 91.5% → 92.0%, 10% 86.0% → 86.9%.

**This is the combination the single-score curves and the decoy-SVM both missed.** `int+decon`
captures most of the gain (a shippable 2-feature gate); the full model adds a little more on the
10-min. Biggest wins are in the aggressive-cut regime (10–25% retained) — exactly where we want to
push feature count down toward comparable-tool levels.

**Caveats:** (1) supervised on the abundance-biased PSM label, so this measures "retains identified
peptides better" and is a *floor* on value — low-abundance rescue is still unmeasured. (2) A shipped
detector can't use PSM labels at runtime — needs either fixed trained coefficients (coefficients are
stable across folds, so plausible) or the noise-faithful decoy to supply negatives label-free. That
decoy is now doubly motivated: it would let this same model train without PSM labels AND probe the
low-abundance regime the PSM metric can't see.

## Option A: noise-faithful decoy + within-band semi-supervised (2026-07-22) — MATCHES the supervised ceiling, label-free

Built a **noise-faithful decoy** (`NOISE_DECOY_SHIFT=<frac>`, detector): rigidly offsets the WHOLE
real-averagine comb off the seed by `frac × (¹³C/z)` so every tooth (anchor included) samples noise —
unlike the old decoys whose anchor stayed on the intense seed. Verified faithful: within every
intensity band the decoy has worse fit than targets (Q5 decon 0.388 vs 0.547, ppm 1.95 vs 1.33).

`within_band_svm.py` — label-free rescoring: negatives = noise decoy; positives = per-band top-q
targets by apex cosine, re-selected each iteration by out-of-fold model score (semi-supervised);
train a SHAPE-only model (per-band labels keep intensity out of it), then rank by
`composite = z(log_intensity) + λ·z(shape)`. PSM GT used for EVAL ONLY.

10-min, recall at retained fraction (λ=0.5, seed-q=0.5), vs intensity and the PSM-supervised ceiling:

| retained | intensity | within-band (label-free) | combo_gate (PSM-supervised) |
|---------:|----------:|-------------------------:|----------------------------:|
| 25% | 94.6% | 96.8% | 96.8% |
| 20% | 93.5% | 95.7% | 96.0% |
| 15% | 90.8% | 94.4% | 95.1% |
| 10% | 84.7% | **91.9%** | 92.1% |

**The label-free method recovers essentially the entire supervised gain** (10% retained: 91.9% vs
supervised 92.1% vs intensity 84.7%). λ≈0.5 (intensity primary, shape secondary) is the sweet spot;
higher λ over-weights shape and regresses. Shape alone still loses to intensity — the win is the blend.

**Seed fraction — counterintuitive:** top-10% positives (the "confident minority" instinct) does
WORSE than top-50% (10% retained: 88.3% vs 91.9%). Because the negatives are clean decoys (not
unlabeled targets), positive *purity* matters less than positive *quantity* — more labeled positives =
more signal, and the decoy supplies the clean contrast. So Percolator's confident-minority seeding
doesn't transfer here; a permissive per-band seed is better.

### All-three breakdown — the 10-min gain does NOT generalize (2026-07-22, revised)

within-band λ=0.5 minus intensity, at each retained fraction:

| retained | 10-min | 65-min | 2-hr |
|---------:|-------:|-------:|-----:|
| 40% | +0.2 | +0.7 | −0.3 |
| 25% | +2.2 | −0.2 | −1.3 |
| 20% | +2.2 | −0.9 | −2.0 |
| 15% | +3.6 | −2.4 | −3.5 |
| 10% | **+7.2** | **−5.4** | **−6.6** |

**The clear win is ONLY the 10-min.** On BOTH large realistic files — the 65-min glyco AND the 2-hr
IonStar (tryptic!) — the label-free gate HURTS at aggressive cuts. So the earlier "works on tryptic,
fails on glyco (glyco≠averagine)" story was WRONG: the 2-hr is tryptic and also has no usable shape
signal.

**Why — it's a signal problem, confirmed by the supervised ceiling.** The PSM-supervised combo_gate
gain is +7 pp on the 10-min but only **+1.3 pp (2-hr)** and **+0.9 pp (65-min)** at a 10% cut, and the
shape-vs-PSM within-band AUC is ~0.82 on the 10-min but **~0.5 on both big files**. So realistic,
complex runs carry little real-vs-junk shape signal at all; the 10-min (a small 555-ref, abundant
"tutorial" sample) is the outlier where shape separates strongly.

**Why label-free HURTS on the big files while the supervised ceiling merely stalls:** the supervised
model keeps intensity dominant (weight ~+1.1) and adds shape only as much as it helps, so it never
drops below intensity. The label-free score blends a near-noise (AUC≈0.5) shape score at a FIXED
λ=0.5 — injecting noise into the ranking. A data-driven λ (→0 when the decoy shows shape isn't
separating) would remove the harm, but on the big files it would then just reduce to intensity — i.e.
no gain to be had.

**Bottom line for option A (honest):** the noise-decoy + within-band pipeline is sound and recovers
whatever shape signal exists label-free. But on realistic complex samples that signal is small
(~1 pp supervised ceiling) and a fixed-λ blend actively hurts, so the impressive 10-min result does
not generalize. To ship at all it needs a data-driven λ (off when the decoy shows no within-band
separation); even then the expected gain on production-scale runs is ~0–1 pp, not the 10-min's +7.
The real headroom for feature reduction remains the intensity/noise floor (per-scan local floor),
not a shape score.

Caveat unchanged: eval is abundance-biased PSM recall, so low-abundance rescue is still unmeasured —
but the decoy now provides label-free negatives at every intensity, which is the machinery a
decoy-FDR (not recall) evaluation of the low bands would need.

## CORRECTION (2026-07-22) — "low intensity = junk" is NOT established

Earlier phrasing in this doc ("the junk is the low-intensity tail", "target junk is a low-intensity
feature") overstated the case. The `within_band.py` test shows:
- PSM-matched (real) features are 69–76% in the top intensity quintile; the lowest quintile holds
  almost none. But the PSM reference is **abundance-biased** (an ID needs MS2-level signal), so it is
  blind to genuine low-abundance features. Intensity's dominance on PSM recall is therefore partly a
  *reference artifact*, not proof that low-intensity features are noise.
- Envelope fit (decon) separates real from unmatched with **AUC 0.82 in the HIGH-intensity band** —
  so there is substantial junk even among abundant features, and fit catches it. => intensity **and**
  fit together should beat intensity alone (a 2-D gate never actually tested). Persistence separates
  real in the lowest band (Q1 AUC 0.62); other shape scores are ~0.5 (inconclusive — tiny labeled
  counts + real-unidentified contamination, NOT evidence shape fails there).

Implication for the SVM below: a useful model must **keep intensity and add fit/persistence**, not
discard intensity as the decoy-only training did.

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
