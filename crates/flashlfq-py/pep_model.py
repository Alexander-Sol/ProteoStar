"""FlashLFQ MBR posterior-error-probability (PEP) model — PLAN.md P3.3.

Replaces the C# ``Microsoft.ML`` FastTree model (``FlashLFQ.PEP.PepAnalysisEngine``)
with a scikit-learn gradient-boosted-tree classifier. It trains on the MBR feature
table emitted by ``flashlfq_py.quant(..., match_between_runs=True)`` — real transferred
peaks (targets) versus ``random_rt`` peaks (decoys) — and assigns each peak
``mbr_pep = 1 - P(target)``.

Usage (the round trip the binding performs)::

    import flashlfq_py
    from pep_model import compute_pep
    table = flashlfq_py.quant(psmtsv, mzmls, match_between_runs=True,
                              pep_model=compute_pep)
    # `table` now has a trailing `mbr_pep` column, and the Rust peaks carry the score.

``compute_pep(record_batch_or_table)`` takes the feature table (a ``pyarrow``
``RecordBatch`` or ``Table``) and returns ``list[float]`` — one PEP per row, **in
the table's row order**, which is exactly what ``apply_mbr_pep`` on the Rust side
expects.

Fidelity. The training *scheme* mirrors the C# engine: peaks are grouped by donor
peptide, ordered, split into three cross-validation partitions (round-robin plus the
target/decoy equalization swap), and the model is retrained ten times — iteration 0
selects positive examples by MBR score, iterations 1-9 by the PEP from the previous
round. The FastTree hyper-parameters (100 trees, 20 leaves, min 10 examples/leaf,
learning rate 0.2, unbalanced-set weighting, fixed seed 42) map onto
``HistGradientBoostingClassifier``. Exact numerical parity with FastTree is not
attainable across ML frameworks — and there is no C# golden for the MBR orchestration
that feeds this — so this is a faithful *structural* port, not a cell-parity one.

Divergences from C# (documented):
  * Donor grouping keys on ``donor_modified_sequence`` (the feature table's available
    proxy) rather than the donor ``Identification`` object reference.
  * The ``swappedDonors`` set in ``EqualizeDonorGroupIndices`` is never populated in C#
    (a latent bug — ``GroupSwap`` does not add to it), so its ``Where(!Contains)``
    filters never filter; replicated here (the set stays empty).
"""

import math

import numpy as np
import pyarrow as pa
from sklearn.ensemble import HistGradientBoostingClassifier

_RANDOM_SEED = 42
_NUM_GROUPS = 3
_NUM_ITERATIONS = 10  # C#: 1 first pass + 9 iterative passes
_TRAINING_FRACTION = 0.25
# Floor for intensities before log2 so a zero/absent intensity does not produce -inf
# (C# computes Math.Log2(Intensity) directly; real transferred peaks are positive).
_TINY = 1e-300


class _Peak:
    """One feature-table row: the 10 model features plus the labels the scheme reads."""

    __slots__ = ("row", "donor", "random_rt", "mbr_score", "features", "pep")

    def __init__(self, row, donor, random_rt, mbr_score, features):
        self.row = row
        self.donor = donor
        self.random_rt = random_rt
        self.mbr_score = mbr_score
        self.features = features
        self.pep = None  # mirrors C# nullable MbrPep (null until scored)


class _DonorGroup:
    """Acceptor peaks sharing one donor peptide (C# FlashLFQ.PEP.DonorGroup)."""

    __slots__ = ("donor", "targets", "decoys")

    def __init__(self, donor, targets, decoys):
        self.donor = donor
        self.targets = targets
        self.decoys = decoys

    @property
    def all(self):
        # C# enumerator: TargetAcceptors.Concat(DecoyAcceptors)
        return self.targets + self.decoys

    @property
    def best_target_score(self):
        return 0.0 if not self.targets else max(p.mbr_score for p in self.targets)


def _to_table(table):
    if isinstance(table, pa.RecordBatch):
        return pa.Table.from_batches([table])
    return table


def _extract_peaks(table):
    cols = {name: table.column(name).to_numpy(zero_copy_only=False)
            for name in table.schema.names}
    intensity = cols["intensity"].astype(np.float64)
    features = np.column_stack([
        cols["ppm_score"],                                   # PpmErrorScore
        cols["intensity_score"],                             # IntensityScore
        cols["rt_score"],                                    # RtScore
        cols["scan_count_score"],                            # ScanCountScore
        cols["isotopic_distribution_score"],                 # IsotopicDistributionScore
        np.abs(cols["mass_error"]),                          # PpmErrorRaw = |MassError|
        np.log2(np.maximum(intensity, _TINY)),              # IntensityRaw = log2(Intensity)
        np.abs(cols["rt_prediction_error"]),                 # RtPredictionErrorRaw
        cols["scan_count"].astype(np.float64),               # ScanCountRaw
        cols["isotopic_pearson_correlation"],                # IsotopicPearsonCorrelation
    ]).astype(np.float64)

    donor = cols["donor_modified_sequence"]
    random_rt = cols["random_rt"]
    mbr_score = cols["mbr_score"].astype(np.float64)
    peaks = []
    for i in range(table.num_rows):
        peaks.append(_Peak(
            row=i,
            donor=str(donor[i]),
            random_rt=bool(random_rt[i]),
            mbr_score=float(mbr_score[i]),
            features=features[i],
        ))
    return peaks


def _build_donor_groups(peaks):
    # C#: OrderByDescending(MbrScore).GroupBy(Identifications.First()).
    # Python's sort is stable, so a descending sort + first-appearance grouping matches.
    ordered = sorted(peaks, key=lambda p: p.mbr_score, reverse=True)
    groups = {}
    order = []
    for p in ordered:
        if p.donor not in groups:
            groups[p.donor] = []
            order.append(p.donor)
        groups[p.donor].append(p)
    donors = []
    for donor in order:
        plist = groups[donor]
        targets = [p for p in plist if not p.random_rt]
        decoys = [p for p in plist if p.random_rt]
        donors.append(_DonorGroup(donor, targets, decoys))
    return donors


def _order_donor_groups(donors):
    # C#: OrderByDescending(TargetAcceptors.Count).ThenByDescending(DecoyAcceptors.Count)
    #     .ThenByDescending(BestTargetMbrScore). Python stable sort on a descending key tuple.
    return sorted(
        donors,
        key=lambda d: (len(d.targets), len(d.decoys), d.best_target_score),
        reverse=True,
    )


def _get_donor_group_indices(donors, num_groups, score_cutoff):
    groups = [[] for _ in range(num_groups)]
    my_index = 0
    while my_index < len(donors):
        sub = 0
        while sub < num_groups and my_index < len(donors):
            groups[sub].append(my_index)
            sub += 1
            my_index += 1
    _equalize_donor_group_indices(donors, groups, score_cutoff, num_groups)
    return groups


def _targets_ge(donor, cutoff):
    return sum(1 for p in donor.targets if p.mbr_score >= cutoff)


def _targets_gt(donor, cutoff):
    return sum(1 for p in donor.targets if p.mbr_score > cutoff)


def _group_swap(donors, groups, donor_a, donor_b, gi_a, gi_b, cutoff,
                target_surplus, decoy_surplus):
    # Surplus is the A-minus-B difference, so a swap shifts it by twice the per-donor delta.
    target_surplus += 2 * (_targets_ge(donors[donor_b], cutoff) - _targets_ge(donors[donor_a], cutoff))
    decoy_surplus += 2 * (len(donors[donor_b].decoys) - len(donors[donor_a].decoys))
    groups[gi_a].append(donor_b)
    groups[gi_a].remove(donor_a)
    groups[gi_b].append(donor_a)
    groups[gi_b].remove(donor_b)
    return target_surplus, decoy_surplus


def _equalize_donor_group_indices(donors, groups, cutoff, num_groups=3):
    swapped = set()  # never populated (faithful to the C# GroupSwap quirk)
    for i in range(num_groups * 3 - 1):
        group_a = i % num_groups
        group_b = (i + 1) % num_groups
        targets_a = sum(_targets_ge(donors[idx], cutoff) for idx in groups[group_a])
        targets_b = sum(_targets_ge(donors[idx], cutoff) for idx in groups[group_b])
        decoys_a = sum(len(donors[idx].decoys) for idx in groups[group_a])
        decoys_b = sum(len(donors[idx].decoys) for idx in groups[group_b])

        stuck = False
        outer = 0
        if not groups[group_a]:
            continue
        min_index = min(groups[group_a])
        target_surplus = targets_a - targets_b
        decoy_surplus = decoys_a - decoys_b

        while (abs(target_surplus) > 1 or abs(decoy_surplus) > 1) and not stuck and outer < 3:
            outer += 1

            inner = 0
            while abs(target_surplus) > 1 and not stuck and inner < 3:
                inner += 1
                swapped_now = False
                for donor_a in sorted((idx for idx in groups[group_a] if idx not in swapped),
                                      reverse=True):
                    a_count = _targets_gt(donors[donor_a], cutoff)
                    if target_surplus > 0:
                        if a_count < 1:
                            continue  # C#: `continue` skips the min_index check below
                        for donor_b in sorted((idx for idx in groups[group_b] if idx not in swapped),
                                              reverse=True):
                            if _targets_gt(donors[donor_b], cutoff) < a_count:
                                target_surplus, decoy_surplus = _group_swap(
                                    donors, groups, donor_a, donor_b, group_a, group_b,
                                    cutoff, target_surplus, decoy_surplus)
                                swapped_now = True
                                break
                    else:
                        for donor_b in sorted((idx for idx in groups[group_b] if idx not in swapped),
                                              reverse=True):
                            if _targets_gt(donors[donor_b], cutoff) > a_count:
                                target_surplus, decoy_surplus = _group_swap(
                                    donors, groups, donor_a, donor_b, group_a, group_b,
                                    cutoff, target_surplus, decoy_surplus)
                                swapped_now = True
                                break
                    if donor_a == min_index:
                        stuck = True
                        break
                    if swapped_now:
                        break

            inner = 0
            while abs(decoy_surplus) > 1 and not stuck and inner < 3:
                inner += 1
                swapped_now = False
                for donor_a in sorted((idx for idx in groups[group_a] if idx not in swapped),
                                      reverse=True):
                    a_count = len(donors[donor_a].decoys)
                    if decoy_surplus > 0:
                        if a_count < 1:
                            continue
                        for donor_b in sorted((idx for idx in groups[group_b] if idx not in swapped),
                                              reverse=True):
                            if len(donors[donor_b].decoys) < a_count:
                                target_surplus, decoy_surplus = _group_swap(
                                    donors, groups, donor_a, donor_b, group_a, group_b,
                                    cutoff, target_surplus, decoy_surplus)
                                swapped_now = True
                                break
                    else:
                        for donor_b in sorted((idx for idx in groups[group_b] if idx not in swapped),
                                              reverse=True):
                            if len(donors[donor_b].decoys) > a_count:
                                target_surplus, decoy_surplus = _group_swap(
                                    donors, groups, donor_a, donor_b, group_a, group_b,
                                    cutoff, target_surplus, decoy_surplus)
                                swapped_now = True
                                break
                    if donor_a == min_index:
                        stuck = True
                        break
                    if swapped_now:
                        break


def _floor_cutoff(values_sorted):
    """C#: ``list[(int)Math.Floor(list.Count * 0.25)]`` on an already-sorted list."""
    return values_sorted[int(math.floor(len(values_sorted) * _TRAINING_FRACTION))]


def _create_peak_data(donors, donor_indices, use_pep):
    """Builds (X, y) for one CV partition. First pass selects positive examples by MBR
    score; iterative passes select by the previous round's PEP (C# CreateChromatographicPeakData
    / CreateChromatographicPeakDataIteration)."""
    if not use_pep:
        scores = sorted((p.mbr_score for i in donor_indices for p in donors[i].all), reverse=True)
        cutoff = _floor_cutoff(scores)
    else:
        # cutoff list maps null PEP -> 1 (C# `?? 1`); ascending sort
        peps = sorted((p.pep if p.pep is not None else 1.0
                       for i in donor_indices for p in donors[i].all))
        cutoff = _floor_cutoff(peps)

    X, y = [], []
    for i in donor_indices:
        for p in donors[i].all:
            if p.random_rt:
                X.append(p.features)
                y.append(False)
            elif not use_pep:
                if p.mbr_score >= cutoff:
                    X.append(p.features)
                    y.append(True)
            else:
                # C#: `peak.MbrPep <= cutoff`; a null PEP compares false -> not selected.
                if p.pep is not None and p.pep <= cutoff:
                    X.append(p.features)
                    y.append(True)
    if not X:
        return np.empty((0, 10)), np.empty((0,), dtype=bool)
    return np.asarray(X, dtype=np.float64), np.asarray(y, dtype=bool)


def _make_classifier():
    # FastTreeBinaryTrainer.Options -> HistGradientBoostingClassifier.
    return HistGradientBoostingClassifier(
        max_iter=100,             # NumberOfTrees
        max_leaf_nodes=20,        # NumberOfLeaves
        min_samples_leaf=10,      # MinimumExampleCountPerLeaf
        learning_rate=0.2,        # LearningRate
        l2_regularization=0.0,
        early_stopping=False,     # fixed number of trees, like FastTree
        class_weight="balanced",  # UnbalancedSets = true
        random_state=_RANDOM_SEED,
    )


def _both_classes(y):
    return y.any() and not y.all()


def compute_pep(table, training_fraction=_TRAINING_FRACTION, num_groups=_NUM_GROUPS,
                num_iterations=_NUM_ITERATIONS):
    """Trains the cross-validated PEP model on the MBR feature `table` and returns one
    PEP (``1 - P(target)``) per row, in the table's row order, as ``list[float]``."""
    table = _to_table(table)
    peaks = _extract_peaks(table)
    n = len(peaks)
    if n == 0:
        return []

    donors = _order_donor_groups(_build_donor_groups(peaks))
    all_scores = sorted((p.mbr_score for d in donors for p in d.all), reverse=True)
    pip_cutoff = all_scores[int(math.floor(len(all_scores) * training_fraction))]
    group_indices = _get_donor_group_indices(donors, num_groups, pip_cutoff)

    for iteration in range(num_iterations):
        fold_data = [
            _create_peak_data(donors, group_indices[g], use_pep=(iteration > 0))
            for g in range(num_groups)
        ]
        # C# returns a failure string (and stops) if any partition lacks a class.
        if not all(_both_classes(y) for _, y in fold_data):
            break

        for g in range(num_groups):
            others = [gg for gg in range(num_groups) if gg != g]
            X = np.vstack([fold_data[others[0]][0], fold_data[others[1]][0]])
            y = np.concatenate([fold_data[others[0]][1], fold_data[others[1]][1]])
            if not _both_classes(y):
                continue
            clf = _make_classifier()
            clf.fit(X, y)
            positive_col = list(clf.classes_).index(True)

            held = [p for idx in group_indices[g] for p in donors[idx].all]
            if not held:
                continue
            proba = clf.predict_proba(np.vstack([p.features for p in held]))[:, positive_col]
            for p, prob in zip(held, proba):
                p.pep = 1.0 - float(prob)

    # Map back to table row order; an unscored peak (training failed) defaults to 1.0 (worst).
    peps = [1.0] * n
    for p in peaks:
        peps[p.row] = p.pep if p.pep is not None else 1.0
    return peps
