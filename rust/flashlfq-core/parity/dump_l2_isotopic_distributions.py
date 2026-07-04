#!/usr/bin/env python3
"""L2 golden dumper — chemical formula -> isotopic distribution (masses[], intensities[]).

Produces the L2 parity artifact (`golden/L2_isotopic_distributions.tsv`) for the P1.6 gate:
for every distinct peptide in `mzLib/Test/FlashLFQ/TestData/AllPSMs.psmtsv`, the fine-grained
isotopic envelope `(masses[], intensities[])` that mzLib's
`Chemistry/IsotopicDistribution.GetDistribution(ChemicalFormula)` produces.

This is a faithful Python replica of `IsotopicDistribution.cs` (the Kubinyi/MIDAs fine-grained
polynomial algorithm), the same algorithm the Rust `isotopic_distribution` module ports. As for
L0/L1, the .NET SDK is unavailable in this environment, so the C# dumper the plan envisioned is
replaced by this replica; the Rust L2 test diffs its own `IsotopicDistribution::get_distribution`
against this golden (intensities rel-1e-6, masses abs-tol).

**Element ordering — the L2 caveat resolved here.** The convolution `MultiplyFineFinalPolynomial`
accumulates products onto a floating-point grid, so its result is element-iteration-order
sensitive. C# iterates `formula.Elements` (a `Dictionary`) in insertion order; the Rust port
iterates `ChemicalFormula::elements()` (a `BTreeMap`) in **atomic-number order**. To make an
index-aligned, near-bit-identical comparison meaningful, this replica iterates elements in the
SAME atomic-number order the Rust side uses (`for an in sorted(counts)`). Both sides are then two
independent re-implementations of the identical arithmetic in the identical order; the gate is
their agreement plus the chemistry sanity checks baked into the Rust test.

Run from anywhere:  python rust/flashlfq-core/parity/dump_l2_isotopic_distributions.py
"""

import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CORE = os.path.dirname(HERE)                       # rust/flashlfq-core
REPO = os.environ.get("MZLIB_DIR", r"F:\mzLib")    # external mzLib checkout (see README)
PSMTSV = os.path.join(REPO, "mzLib", "Test", "FlashLFQ", "TestData", "AllPSMs.psmtsv")
GOLDEN_TSV = os.path.join(HERE, "golden", "L2_isotopic_distributions.tsv")

# Reuse the L0 loader (element data) and the L1 formula builder (residue table + mod column),
# so all three goldens share one source of truth.
sys.path.insert(0, HERE)
from dump_periodic_table import extract_nist_literal, load_periodic_table  # noqa: E402
from dump_l1_peptide_formulas import (  # noqa: E402
    CS_PATH,
    build_element_index,
    full_formula_counts,
)

# Constants — verbatim from IsotopicDistribution.cs.
DEFAULT_FINE_RESOLUTION = 0.01
DEFAULT_MIN_PROBABILITY = 1e-200
DEFAULT_MOLECULAR_WEIGHT_RESOLUTION = 1e-12

NAN = float("nan")


class Poly:
    """One term of the Kubinyi polynomial (C# struct `Polynomial`): centroid mass + weight."""

    __slots__ = ("power", "probability")

    def __init__(self, power, probability):
        self.power = power
        self.probability = probability


# ---------------------------------------------------------------------------
# FactorLn — memoised cumulative log-sum, bit-identical summation order to the C#
# factorLnArray cache: factorLnArray[j+1] = factorLnArray[j] + log(j+1).
# ---------------------------------------------------------------------------
_factor_ln_cache = [0.0, 0.0]  # ln(0!) = ln(1!) = 0


def factor_ln(n):
    if n <= 1:
        return 0.0
    while len(_factor_ln_cache) <= n:
        k = len(_factor_ln_cache)  # next index == the factor being multiplied in
        _factor_ln_cache.append(_factor_ln_cache[k - 1] + math.log(k))
    return _factor_ln_cache[n]


def multiple_fine_polynomial_recursive_helper(
    mins, maxs, indices, index, f_polynomial, composition, atoms, min_prob, max_value
):
    """Port of MultipleFinePolynomialRecursiveHelper (multinomial index enumeration)."""
    indices[index] = mins[index]
    while indices[index] <= maxs[index]:
        if index < len(mins) - 1:
            multiple_fine_polynomial_recursive_helper(
                mins, maxs, indices, index + 1, f_polynomial, composition, atoms, min_prob, max_value
            )
        else:
            l = atoms - sum(indices)
            if not (l < 0 or l > max_value):
                last = len(composition) - 1
                prob = factor_ln(atoms) - factor_ln(l) + l * composition[last]["log_probability"]
                power = l * composition[last]["power"]
                for i in range(0, last):  # 0 .. last-1 inclusive
                    index_value = indices[i]
                    t_comp = composition[i]
                    prob -= factor_ln(index_value)
                    prob += index_value * t_comp["log_probability"]
                    power += index_value * t_comp["power"]
                prob = math.exp(prob)
                if prob >= min_prob:
                    f_polynomial.append(Poly(power, prob))
        indices[index] += 1


def multiply_fine_polynomial(elemental_composition, fine_resolution, mw_resolution, fine_min_prob):
    """Port of MultiplyFinePolynomial: per-element expansion, then cross-element convolution."""
    NC = 10
    NC_ADD_VALUE = 1
    N_ATOMS = 200

    n = sum(1 for c in elemental_composition if len(c) > 0)
    f_polynomial = [[] for _ in range(n)]

    for k in range(n):
        composition = elemental_composition[k]
        size = len(composition)
        atoms = composition[0]["atoms"]
        nc_add = 10 if atoms < N_ATOMS else NC_ADD_VALUE

        if size == 1:
            probability = composition[0]["probability"]
            n1 = int(atoms * probability)
            prob = factor_ln(atoms) - factor_ln(n1) + n1 * composition[0]["log_probability"]
            prob = math.exp(prob)
            f_polynomial[k].append(Poly(n1 * composition[0]["power"], prob))
        else:
            means = [0] * size
            stds = [0] * size
            for i in range(size):
                n1 = int(composition[0]["atoms"] * composition[i]["probability"])
                s1 = int(
                    math.ceil(
                        nc_add
                        + NC
                        * math.sqrt(
                            composition[i]["atoms"]
                            * composition[i]["probability"]
                            * (1.0 - composition[i]["probability"])
                        )
                    )
                )
                means[i] = n1
                stds[i] = s1

            mins = [0] * (size - 1)
            maxs = [0] * (size - 1)
            indices = [0] * (size - 1)
            for i in range(size - 1):
                mx = max(0, means[i] - stds[i])
                indices[i] = mx
                mins[i] = mx
                maxs[i] = means[i] + stds[i]
            max_value = means[size - 1] + stds[size - 1]

            multiple_fine_polynomial_recursive_helper(
                mins, maxs, indices, 0, f_polynomial[k], composition, atoms, fine_min_prob, max_value
            )

    t_polynomial = f_polynomial[0]  # C#: reference assignment (fPolynomial[0] not reused after)

    if n <= 1:
        return t_polynomial

    fgid_polynomial = []
    for k in range(1, n):
        multiply_fine_final_polynomial(
            t_polynomial, f_polynomial[k], fgid_polynomial, fine_resolution, mw_resolution, fine_min_prob
        )

    return t_polynomial


def multiply_fine_final_polynomial(
    t_polynomial, f_polynomial, fgid_polynomial, fine_resolution, mw_resolution, fine_min_prob
):
    """Port of MultiplyFineFinalPolynomial: convolve t with f onto the shared fgid grid."""
    i_count = len(t_polynomial)
    j_count = len(f_polynomial)
    if i_count == 0 or j_count == 0:
        return

    delta_mass = fine_resolution / mw_resolution
    min_probability = fine_min_prob

    min_mass = t_polynomial[0].power + f_polynomial[0].power
    max_mass = t_polynomial[i_count - 1].power + f_polynomial[j_count - 1].power

    max_index = int(abs(max_mass - min_mass) / delta_mass + 0.5)
    if max_index >= len(fgid_polynomial):
        extra = max_index - len(fgid_polynomial)
        for _ in range(extra + 1):
            fgid_polynomial.append(Poly(NAN, NAN))

    for t in range(len(t_polynomial)):
        for f in range(len(f_polynomial)):
            prob = t_polynomial[t].probability * f_polynomial[f].probability
            if prob <= min_probability:
                continue
            power = t_polynomial[t].power + f_polynomial[f].power
            indext = int(abs(power - min_mass) / delta_mass + 0.5)
            temp = fgid_polynomial[indext]
            if math.isnan(temp.power) or math.isnan(prob):
                fgid_polynomial[indext] = Poly(power * prob, prob)
            else:
                fgid_polynomial[indext] = Poly(temp.power + power * prob, temp.probability + prob)

    index = len(t_polynomial)
    j = 0
    for i in range(len(fgid_polynomial)):
        if not math.isnan(fgid_polynomial[i].probability):
            term = Poly(
                fgid_polynomial[i].power / fgid_polynomial[i].probability,
                fgid_polynomial[i].probability,
            )
            if j < index:
                t_polynomial[j] = term
                j += 1
            else:
                t_polynomial.append(term)
        fgid_polynomial[i] = Poly(NAN, NAN)

    if j < index:
        del t_polynomial[j:]


def merge_fine_polynomial(t_polynomial, mw_resolution, merge_fine_resolution):
    """Port of MergeFinePolynomial: 9 passes coalescing adjacent terms by a k-scaled threshold."""
    t_polynomial.sort(key=lambda p: p.power)
    count = len(t_polynomial)

    for k in range(1, 10):
        for i in range(count):
            power = t_polynomial[i].power
            if math.isnan(power):
                continue
            probability = t_polynomial[i].probability
            temp_power = power * probability
            temp_prob = probability

            for j in range(i + 1, count):
                # Reads t_polynomial[i].power fresh each iteration: it is the running merged
                # centroid (updated in place below) — this ordering matters for parity.
                value = abs(
                    t_polynomial[i].power * mw_resolution - t_polynomial[j].power * mw_resolution
                )
                threshold = (
                    k * merge_fine_resolution / 8
                    if k <= 8
                    else merge_fine_resolution + merge_fine_resolution / 100
                )
                if value <= threshold:
                    temp_power += t_polynomial[j].power * t_polynomial[j].probability
                    temp_prob += t_polynomial[j].probability
                    t_polynomial[i] = Poly(temp_power / temp_prob, temp_prob)
                    t_polynomial[j] = Poly(NAN, NAN)
                else:
                    break

            t_polynomial[i] = Poly(temp_power / temp_prob, temp_prob)

    return [p for p in t_polynomial if not math.isnan(p.power)]


def calculate_fine_grain(
    elemental_composition, mw_resolution, merge_fine_resolution, fine_resolution, fine_min_prob
):
    f_polynomial = multiply_fine_polynomial(
        elemental_composition, fine_resolution, mw_resolution, fine_min_prob
    )
    f_polynomial = merge_fine_polynomial(f_polynomial, mw_resolution, merge_fine_resolution)
    masses = [p.power * mw_resolution for p in f_polynomial]
    intensities = [p.probability for p in f_polynomial]
    return masses, intensities


def get_distribution(
    counts,
    elements_by_number,
    fine_resolution=DEFAULT_FINE_RESOLUTION,
    min_probability=DEFAULT_MIN_PROBABILITY,
    mw_resolution=DEFAULT_MOLECULAR_WEIGHT_RESOLUTION,
):
    """Port of IsotopicDistribution.GetDistribution for plain-element formulas.

    `counts` is {atomic_number: count}. Elements are iterated in atomic-number order to match
    the Rust BTreeMap (see module docstring). No isotope-specified atoms occur in this corpus,
    so the C# `additionalMass` path contributes 0 and is omitted.
    """
    new_fine_resolution = fine_resolution / 2.0
    merge_fine_resolution = fine_resolution

    elemental_composition = []
    for an in sorted(counts):
        count = counts[an]
        el = elements_by_number[an]
        # OrderBy(iso => iso.AtomicMass)
        isotopes = sorted(el["isotopes"], key=lambda iso: iso[1])
        comp = []
        for (_mn, am, ab) in isotopes:
            comp.append(
                {
                    "atoms": count,
                    "molecular_weight": am,
                    "power": am,
                    "probability": ab,
                    "log_probability": 0.0,
                }
            )
        elemental_composition.append(comp)

    for comp in elemental_composition:
        sum_prob = sum(c["probability"] for c in comp)
        for c in comp:
            c["probability"] /= sum_prob
            c["log_probability"] = math.log(c["probability"])
            c["power"] = math.floor(c["molecular_weight"] / mw_resolution + 0.5)

    return calculate_fine_grain(
        elemental_composition, mw_resolution, merge_fine_resolution, new_fine_resolution, min_probability
    )


def main():
    raw = extract_nist_literal(CS_PATH)
    elements = load_periodic_table(raw)
    element_index = build_element_index(elements)  # symbol -> (atomic_number, principal mass)
    elements_by_number = {el["atomic_number"]: el for el in elements}

    with open(PSMTSV, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    # One golden entry per distinct (Base Sequence, Mods Combined Chemical Formula) — same set
    # and first-seen order as the L1 gate.
    distinct = {}
    order = []
    for r in rows:
        base = r["Base Sequence"]
        mod = r["Mods Combined Chemical Formula"]
        key = (base, mod)
        if key not in distinct:
            distinct[key] = r
            order.append(key)

    lines = []
    lines.append("# L2 golden — chemical formula -> isotopic distribution (masses[], intensities[])")
    lines.append("# replicated from mzLib Chemistry/IsotopicDistribution.cs (Kubinyi fine-grained")
    lines.append("# polynomial), elements iterated in atomic-number order to match the Rust BTreeMap.")
    lines.append("# Records:")
    lines.append("#   D<TAB>baseSeq<TAB>combinedMod<TAB>nPeaks")
    lines.append("#   M<TAB>mass(repr) ... one per peak, ascending mass (merge order)")
    lines.append("#   I<TAB>intensity(repr) ... one per peak, parallel to the M record")
    lines.append("# floats are Python repr() (shortest round-trip); Rust parses them back bit-identically")

    max_peaks = 0
    for key in order:
        base, mod = key
        counts = full_formula_counts(base, mod, element_index)
        masses, intensities = get_distribution(counts, elements_by_number)
        assert len(masses) == len(intensities)
        max_peaks = max(max_peaks, len(masses))
        mod_field = mod if mod.strip() else "-"
        lines.append("D\t%s\t%s\t%d" % (base, mod_field, len(masses)))
        lines.append("M\t" + " ".join(repr(m) for m in masses))
        lines.append("I\t" + " ".join(repr(it) for it in intensities))

    os.makedirs(os.path.dirname(GOLDEN_TSV), exist_ok=True)
    with open(GOLDEN_TSV, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")

    print("psmtsv rows: %d   distinct peptides: %d" % (len(rows), len(order)))
    print("max peaks in any envelope: %d" % max_peaks)
    print("wrote golden: %s" % os.path.relpath(GOLDEN_TSV, REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
