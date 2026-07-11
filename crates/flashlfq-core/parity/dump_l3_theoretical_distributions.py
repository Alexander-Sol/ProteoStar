#!/usr/bin/env python3
"""L3 golden dumper — modified sequence -> final (massShift, normAbundance)[] + PeakfindingMass.

Produces the L3 parity artifact (`golden/L3_theoretical_distributions.tsv`) for the P1.8 gate:
for every distinct *modified* peptide (Full Sequence) in
`mzLib/Test/FlashLFQ/TestData/AllPSMs.psmtsv`, the list of expected isotope peaks
`(massShift, normalizedAbundance)` and the peakfinding mass that mzLib's
`FlashLfqEngine.CalculateTheoreticalIsotopeDistributions()` produces.

This is a faithful Python replica of that method, layered on the L0 element loader, the L1
residue table, and the L2 Kubinyi `IsotopicDistribution.GetDistribution` replica. As for
L0/L1/L2, the .NET SDK is unavailable in this environment, so the C# dumper the plan
envisioned is replaced by this replica; the Rust L3 test diffs its own
`theoretical_isotope_distribution::expected_isotope_peaks` + `peakfinding_mass` against this
golden (massShift abs-1e-6, abundance rel-1e-6, array lengths exact).

**Important fidelity point — the formula in the `None` branch is built from the BASE SEQUENCE
ONLY (no modifications).** FlashLFQ's psmtsv path leaves `Identification.OptionalChemicalFormula`
null, so the method builds `new Peptide(BaseSequence).GetChemicalFormula()` (backbone + terminal
water, *no mods*) and then, if the known mass and that bare formula disagree by more than 20 Da,
tops the formula up with **averagine** scaled to cover the gap. So a modified peptide's mod mass
is approximated by averagine here, NOT taken from the `Mods Combined Chemical Formula` column
(that column drives L1/L2, but is irrelevant to this L3 path).

The known monoisotopic mass used for re-centering and for the 20 Da test is the C#
**`Peptide Monoisotopic Mass`** column (= `Identification.MonoisotopicMass`), NOT the recomputed
formula mass — using the column makes the averagine top-up trigger identically to C#.

Distinct key is the **Full Sequence** (modified sequence), matching the C#
`ModifiedSequenceToIsotopicDistribution` dictionary keyed on `id.ModifiedSequence`.

Run from anywhere:  python rust/flashlfq-core/parity/dump_l3_theoretical_distributions.py
"""

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CORE = os.path.dirname(HERE)                       # rust/flashlfq-core
REPO = os.environ.get("MZLIB_DIR", r"F:\mzLib")    # external mzLib checkout (see README)
PSMTSV = os.path.join(REPO, "mzLib", "Test", "FlashLFQ", "TestData", "AllPSMs.psmtsv")
GOLDEN_TSV = os.path.join(HERE, "golden", "L3_theoretical_distributions.tsv")

# Reuse the L0 loader, the L1 formula builder, and the L2 distribution replica so all four
# goldens share one source of truth.
sys.path.insert(0, HERE)
from dump_periodic_table import extract_nist_literal, load_periodic_table  # noqa: E402
from dump_l1_peptide_formulas import (  # noqa: E402
    CS_PATH,
    RESIDUES,
    add_into,
    build_element_index,
    parse_formula,
)
from dump_l2_isotopic_distributions import get_distribution  # noqa: E402

# Averagine composition + distribution params — verbatim from
# FlashLfqEngine.CalculateTheoreticalIsotopeDistributions.
AVERAGE_C = 4.9384
AVERAGE_H = 7.7583
AVERAGE_O = 1.4773
AVERAGE_N = 1.3577
AVERAGE_S = 0.0417

FINE_RESOLUTION = 0.125
MIN_PROBABILITY = 1e-8
NUM_ISOTOPES_REQUIRED = 2
AVERAGINE_MASS_DIFF_THRESHOLD = 20.0


def base_formula_counts(base_seq, element_index):
    """{atomic_number: count} for the bare peptide backbone (residues + H + OH), NO mods.

    Mirrors `new Peptide(baseSequence).GetChemicalFormula()`: N-terminus H + C-terminus OH
    (= one water) + every residue's condensed formula.
    """
    counts = {}
    add_into(counts, parse_formula("H", element_index))
    add_into(counts, parse_formula("OH", element_index))
    for letter in base_seq:
        if letter not in RESIDUES:
            raise ValueError("unknown residue %r in %r" % (letter, base_seq))
        add_into(counts, parse_formula(RESIDUES[letter], element_index))
    return {an: c for an, c in counts.items() if c != 0}


def mono_mass(counts, an_to_mass):
    """Monoisotopic mass, summed in atomic-number order to match the Rust BTreeMap."""
    total = 0.0
    for an in sorted(counts):
        total += an_to_mass[an] * counts[an]
    return total


def all_residues_valid(base_seq):
    return len(base_seq) > 0 and all(c in RESIDUES for c in base_seq)


def add_averagine(counts, averagines, element_index):
    """Add averagine atoms with C# banker's-rounding integer counts.

    Mirrors the five `formula.Add("X", (int)Math.Round(averagines * averageX, 0))` calls;
    a zero count is a no-op (ChemicalFormula.Add returns early), matching add_element.
    """
    for symbol, avg in (("C", AVERAGE_C), ("H", AVERAGE_H), ("O", AVERAGE_O),
                        ("N", AVERAGE_N), ("S", AVERAGE_S)):
        # Python round() is round-half-to-even (banker's), like C# Math.Round(x, 0) and
        # Rust f64::round_ties_even. round() on a float returns an int in Python 3.
        count = round(averagines * avg)
        if count == 0:
            continue
        an = element_index[symbol][0]
        counts[an] = counts.get(an, 0) + count
    return {an: c for an, c in counts.items() if c != 0}


def resolve_formula_counts(base_seq, mono, averagine_mass, element_index, an_to_mass):
    """Resolve the {atomic_number: count} used for the isotope distribution (None-formula path)."""
    if all_residues_valid(base_seq):
        counts = base_formula_counts(base_seq, element_index)
        mass_diff = mono - mono_mass(counts, an_to_mass)
        if abs(mass_diff) > AVERAGINE_MASS_DIFF_THRESHOLD:
            counts = add_averagine(counts, mass_diff / averagine_mass, element_index)
        return counts
    else:
        counts = {}
        counts = add_averagine(counts, mono / averagine_mass, element_index)
        return counts


def expected_isotope_peaks(base_seq, mono, averagine_mass, element_index,
                           elements_by_number, an_to_mass):
    """Replica of the per-id body of CalculateTheoreticalIsotopeDistributions.

    Returns (peaks, peakfinding_mass) where peaks is a list of (massShift, normAbundance).
    """
    counts = resolve_formula_counts(base_seq, mono, averagine_mass, element_index, an_to_mass)
    formula_mono = mono_mass(counts, an_to_mass)

    masses, abundances = get_distribution(
        counts, elements_by_number, FINE_RESOLUTION, MIN_PROBABILITY
    )

    # Re-center distribution masses onto the known monoisotopic mass (kept as the two separate
    # C# steps `+= (mono - formulaMono)` then `-= mono` so the FP matches bit-for-bit).
    masses = [m + (mono - formula_mono) for m in masses]

    highest = max(abundances)

    peaks = []
    for i in range(len(masses)):
        mass_shift = masses[i] - mono
        abundance = abundances[i] / highest
        if len(peaks) < NUM_ISOTOPES_REQUIRED or abundance > 0.1:
            peaks.append((mass_shift, abundance))

    # PeakfindingMass = MonoisotopicMass + shift of the peak whose normalized abundance == 1.0.
    most_abundant_shift = next(ms for (ms, ab) in peaks if ab == 1.0)
    peakfinding_mass = mono + most_abundant_shift
    return peaks, peakfinding_mass


def main():
    raw = extract_nist_literal(CS_PATH)
    elements = load_periodic_table(raw)
    element_index = build_element_index(elements)  # symbol -> (atomic_number, principal mass)
    elements_by_number = {el["atomic_number"]: el for el in elements}
    an_to_mass = {an: m for sym, (an, m) in element_index.items()}

    avg = {el["symbol"]: el["average_mass"] for el in elements}
    averagine_mass = (avg["C"] * AVERAGE_C + avg["H"] * AVERAGE_H + avg["O"] * AVERAGE_O
                      + avg["N"] * AVERAGE_N + avg["S"] * AVERAGE_S)

    with open(PSMTSV, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    # One golden entry per distinct Full Sequence (modified sequence), first-seen order —
    # matching the C# dictionary keyed on id.ModifiedSequence.
    distinct = {}
    order = []
    for r in rows:
        key = r["Full Sequence"]
        if key not in distinct:
            distinct[key] = r
            order.append(key)

    lines = []
    lines.append("# L3 golden — modified sequence -> (massShift, normAbundance)[] + PeakfindingMass")
    lines.append("# replicated from mzLib FlashLfqEngine.CalculateTheoreticalIsotopeDistributions.")
    lines.append("# Formula is built from the BASE SEQUENCE ONLY (no mods) + averagine top-up when")
    lines.append("# |mono - baseFormulaMono| > 20; mono is the psmtsv 'Peptide Monoisotopic Mass'.")
    lines.append("# Records:")
    lines.append("#   D<TAB>fullSeq<TAB>baseSeq<TAB>mono(repr)<TAB>nPeaks<TAB>peakfindingMass(repr)")
    lines.append("#   S<TAB>massShift(repr) ... one per peak (distribution/mass-ascending order)")
    lines.append("#   A<TAB>normAbundance(repr) ... one per peak, parallel to the S record")
    lines.append("# floats are Python repr() (shortest round-trip); Rust parses them back bit-identically")

    max_peaks = 0
    n_averagine = 0
    for key in order:
        r = distinct[key]
        base = r["Base Sequence"]
        mono = float(r["Peptide Monoisotopic Mass"])

        # Track whether averagine fired (informational).
        if all_residues_valid(base):
            md = mono - mono_mass(base_formula_counts(base, element_index), an_to_mass)
            if abs(md) > AVERAGINE_MASS_DIFF_THRESHOLD:
                n_averagine += 1

        peaks, pf = expected_isotope_peaks(
            base, mono, averagine_mass, element_index, elements_by_number, an_to_mass
        )
        max_peaks = max(max_peaks, len(peaks))
        lines.append("D\t%s\t%s\t%r\t%d\t%r" % (key, base, mono, len(peaks), pf))
        lines.append("S\t" + " ".join(repr(ms) for (ms, _ab) in peaks))
        lines.append("A\t" + " ".join(repr(ab) for (_ms, ab) in peaks))

    os.makedirs(os.path.dirname(GOLDEN_TSV), exist_ok=True)
    with open(GOLDEN_TSV, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")

    print("psmtsv rows: %d   distinct modified sequences: %d" % (len(rows), len(order)))
    print("peptides where averagine top-up fired: %d" % n_averagine)
    print("max peaks in any envelope: %d" % max_peaks)
    print("averagine unit mass: %r" % averagine_mass)
    print("wrote golden: %s" % os.path.relpath(GOLDEN_TSV, REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
