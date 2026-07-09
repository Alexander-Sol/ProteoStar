# Top-down deconvolution — design plan (proteoform-aware, multi-species)

Consolidates the design direction for the top-down monoisotope/overlap problem, from the 2026-07-09
benchmark findings (`analysis/topdown_bench/RESULTS.md`) and the session's design discussion. The
detector already *finds* ~99% of proteoforms; the work is all in **deconvolution / mass assignment**.

## Findings that shape the design
- **Off-by-one consensus bridging must be OFF for top-down** (`OFFBYONE_UNITS=0`, +16pp). It conflates
  genuinely-close co-eluting species. Notably **deamidation (+0.984 Da)** is only 0.019 Da from a
  +1.003 Da isotope step — indistinguishable by mass tolerance at high mass — so any off-by-one
  "correction" risks destroying real deamidated proteoforms.
- **Averaging during refine is not on by default** (shipped default = apex scan only). The FWHM-window
  composite gives the low isotope peaks (which discriminate the monoisotope) the SNR they need.
- Single-species averagine-cosine mono placement plateaus (~+3pp) — the metric is dominated by the big
  apex peaks and is near-blind to a ±1 shift for large masses.

## Design: proteoform-family joint deconvolution
Refine a co-eluting cluster **against a combination of envelopes that are all the same protein bearing
common/expected mass shifts**, solved jointly. Machinery already exists: `joint_fit.rs` has NNLS
(`nnls`), `joint_envelope_fit` (fit N averagine `Component`s to a window), and `joint_fit_target_shift`
(search the target mono while neighbours are fixed).

**Component generation from a PTM-shift family.** Anchor on a confidently-placed proteoform in the
co-elution window; generate sibling `Component`s at the anchor mass + each expected Δmass, at each
plausible charge. Fit all jointly by NNLS; the apportioned amplitudes separate shared/overlapping
peaks, and the family prior stops a modified form's peak from being read as an off-by-one of the base.

**Expected Δmass table (monoisotopic, Da)** — common PTMs / proteoform relationships:
- Deamidation **+0.98402** (N→D, Q→E) — the off-by-one look-alike
- Oxidation **+15.99491** (and 2× = +31.99)
- Methylation **+14.01565**, di- **+28.0313**, tri-/acetyl-adjacent
- Acetylation **+42.01057**, trimethyl **+42.04695** (near-isobaric with acetyl — 0.036 Da apart)
- Phosphorylation **+79.96633**
- Dehydration/water loss **−18.01056**, ammonia loss **−17.02655**
- N-terminal Met excision **−131.04049**; +Met **+131.04049**
- Disulfide (−2H per bond) **−2.01565**
- Consider integer multiples / combinations within a bounded total.

**Mono regularization.** When a feature's mono placement is ambiguous, prefer the placement that puts it
at a common-PTM Δmass from a confidently-placed co-eluting neighbour (a soft prior added to the
envelope-fit score), rather than the placement the cosine barely prefers.

**Two-dimensional consensus: share across the family × charge grid.** Today cross-charge consensus
shares information *within one species* (its charge ladder). But every modified form also appears at
many charges, so a co-eluting proteoform family is a 2-D grid — {base, base+Δ₁, base+Δ₂, …} ×
{z=6..25} — in which (a) family members are locked to each other by *known* Δmass and (b) each member's
mass is over-determined by its own charge ladder. Resolve the whole grid jointly:
- Cluster co-eluting features; find a base and the siblings whose neutral-mass deltas match the PTM
  table (within tight neutral-mass tolerance — tight enough to separate +0.984 deamidation from a
  +1.003 isotope step, which needs the cross-charge agreement, not per-charge ppm).
- The base's clean high-charge envelopes fix the family's mass frame; each sibling's mono is then the
  base mono + its (known, integer-PTM) Δmass, cross-checked against that sibling's own charge ladder.
- This **distinguishes a real +0.984 deamidated species from a −1 monoisotope miscount** — the very
  ambiguity that forced `OFFBYONE_UNITS=0` — because the deamidated form has its *own* consistent charge
  ladder at base+0.984, whereas an isotope miscount would not recur coherently across charges at a
  constant neutral-mass offset.

## Consensus = joint cross-charge proteoform inference (NOT merging)
The correct framing (clarified in design discussion): the cross-charge step is **not** a grouping/merging
of features. It is a **joint inference** over a co-eluting window that simultaneously determines
**(1) how many distinct proteoforms are present** and **(2) the mass of each**, gaining accuracy by
**sharing information across charge states** — every proteoform appears at many charges, so its mass is
massively over-determined.

Mechanism: pool the per-charge mass evidence across the whole co-eluting window; a candidate proteoform
mass that recurs *consistently across multiple charges* is real (its mass = the cross-charge consensus,
far more precise than any single charge), while a mass seen at only one charge is noise or a per-charge
mono slip. This is **multi-charge support**, and it is the opposite of merging: a real +0.984 deamidation
and its base are BOTH kept as well-supported proteoforms (each with its own charge ladder), instead of
being forced into one. This is why flexible-apex *grouping* was the wrong tool (it merges genuinely-close
species) — the separating power comes from requiring cross-charge support, not from a mass-tolerance link.

Concrete form: enumerate candidate masses (each charge's multi-envelope monos), cluster tightly in
NEUTRAL mass, score each cluster by distinct-charge support, keep/rank by support; the surviving clusters
ARE the proteoform set and their weighted masses ARE the answers. Distinct from the old mono-keyed
consensus in that it (a) never bridges off-by-one, (b) treats count as an output, (c) shares across
charges for mass precision. Next implementation step; measure vs the intersection GT.

## (Superseded framing) Cross-charge consensus must be reworked for multi-envelope
Multi-envelope refine gives *each charge* its own set of species monos (one per co-eluting envelope),
each from an independent shift search — so the same proteoform can land at mono X (z=12) and X±1 (z=15),
and the current **mass-keyed** consensus (group features by monoisotopic mass) fails to combine them,
especially with the off-by-one bridge disabled. Fix by **inverting what we key on**:
- **Group cross-charge by the apex (most-abundant) peak's neutral mass** — charge-independent
  (|z|·apex_mz − |z|·proton is the same physical value at every charge) and cleanly observed (the
  tallest isotope), so it reliably collects a proteoform's charge ladder even when the monos disagree.
- **Resolve the monoisotope once per group**, pooling every charge's composite (and the PTM-family
  siblings). This is where family × charge evidence combines; it distinguishes a real +0.984 deamidation
  (a consistent apex ladder at base+0.984 across charges) from a −1 monoisotope miscount (no consistent
  ladder). Invert the logic: **group on the robust quantity (apex), resolve the fragile one (mono)
  collectively.** Re-measure after — it changes how features merge.

## Work order
1. **Measure first** (per the "measure the overlap cost before building" decision): quantify what plain
   averaging (`REFINE_METHOD=shift_composite`/`classic`) and the existing `JOINT_FIT` NNLS pass buy on
   Jurkat, on top of `OFFBYONE_UNITS=0`. (In progress.)
2. If joint fit helps, **generate joint components from the PTM-family** (above) instead of / in
   addition to raw co-eluting neighbours; add the Δmass prior to component ranking.
3. Add the deamidation-aware guard so mono correction never turns a +0.984 species into its −1 neighbour.
4. Longer term: **IsoDec** (see `IsoDec-Port-Investigation.md`) — a learned charge/mono assigner that
   sidesteps the averagine-cosine weakness entirely; strongest single lever if ported.

## Benchmarks to gate against
`analysis/topdown_bench/` — Jurkat (2885 proteoforms) and Golden (930, per L/M/H tier), scored strict
(20 ppm) and isotope-tolerant by `score_topdown.py`. Iterate refine cheaply via the `LOAD_DETECTED`
cache; validate finalists with a full detect run.
