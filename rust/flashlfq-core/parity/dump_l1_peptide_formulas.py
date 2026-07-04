#!/usr/bin/env python3
"""L1 golden dumper — modified sequence -> element counts + monoisotopic mass.

Produces the L1 parity artifact (`golden/L1_peptide_formulas.tsv`) for the P1.4 gate:
for every distinct peptide in `mzLib/Test/FlashLFQ/TestData/AllPSMs.psmtsv`, the full
chemical formula (backbone residues + terminal water + modifications) as exact element
counts, plus the monoisotopic mass.

Why this is tractable without a modified-sequence parser: the psmtsv pre-resolves every
modification to a chemical formula. The **Mods Combined Chemical Formula** column is the
sum of all modification formulas on the peptide (mzLib writes it from the same
`Modification.ChemicalFormula` objects that `PeptideWithSetModifications.FullChemicalFormula`
adds). So the full formula is exactly

    peptide_base_formula(Base Sequence)  +  parse(Mods Combined Chemical Formula)

which mirrors `new Peptide(baseSequence).GetChemicalFormula()` then `Add(modFormula)` — the
same arithmetic the Rust `peptide::peptide_formula_with_mods` performs. P1.3 deliberately
left mod-string resolution out; here the resolution is "read the column mzLib already wrote".

This script is the independent reference implementation (a faithful Python replica of the
residue table + chemical-formula mono-mass sum over the *same* embedded periodic table the
Rust side uses, `data/periodic_table.nist.txt`). The Rust L1 test diffs its own computation
against this golden (counts exact, mass rel-1e-9) AND cross-checks the mono mass against the
psmtsv's C# `Peptide Monoisotopic Mass` column (the real ground truth, at its 5-decimal
rounding tolerance). Agreement of the two independent replicas plus the C# cross-check is
the L1 gate.

No .NET SDK in this environment, so (as for L0) the C# dumper the plan envisioned is replaced
by this replica; the psmtsv column keeps a real tie to C# output.

Run from anywhere:  python rust/flashlfq-core/parity/dump_l1_peptide_formulas.py
"""

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CORE = os.path.dirname(HERE)                       # rust/flashlfq-core
REPO = os.environ.get("MZLIB_DIR", r"F:\mzLib")    # external mzLib checkout (see README)
PSMTSV = os.path.join(REPO, "mzLib", "Test", "FlashLFQ", "TestData", "AllPSMs.psmtsv")
GOLDEN_TSV = os.path.join(HERE, "golden", "L1_peptide_formulas.tsv")

# Reuse the L0 loader so both goldens share one source of truth for element data.
sys.path.insert(0, HERE)
from dump_periodic_table import extract_nist_literal, load_periodic_table  # noqa: E402

CS_PATH = os.path.join(REPO, "mzLib", "Chemistry", "PeriodicTable.cs")

# Residue one-letter -> condensed (water-loss) formula. Verbatim from Residue.cs's
# ResiduesDictionary; identical to peptide::residue_formula on the Rust side. B/J/X/Z are
# null in mzLib and intentionally absent.
RESIDUES = {
    "A": "C3H5NO", "R": "C6H12N4O", "N": "C4H6N2O2", "D": "C4H5NO3", "C": "C3H5NOS",
    "E": "C5H7NO3", "Q": "C5H8N2O2", "G": "C2H3NO", "H": "C6H7N3O", "I": "C6H11NO",
    "L": "C6H11NO", "K": "C6H12N2O", "M": "C5H9NOS", "F": "C9H9NO", "P": "C5H7NO",
    "O": "C12H19N3O2", "U": "C3H5NOSe", "S": "C3H5NO2", "T": "C4H7NO2", "W": "C11H10N2O",
    "Y": "C9H9NO2", "V": "C5H9NO",
}


def build_element_index(elements):
    """symbol -> (atomic_number, principal_isotope_atomic_mass)."""
    idx = {}
    for el in elements:
        pmn = el["principal_mass_number"]
        pmass = next(am for (mn, am, _ab) in el["isotopes"] if mn == pmn)
        idx[el["symbol"]] = (el["atomic_number"], pmass)
    return idx


def parse_formula(formula, element_index):
    """Parse a plain chemical formula -> {atomic_number: count}.

    Faithful to ChemicalFormula.ParseFormula for the subset that appears here: symbol then
    optional leading '-' then optional integer count (default 1). No isotope braces occur in
    the residue table or in any AllPSMs mod formula, so they are unnecessary; if one ever
    appears the tokenizer raises rather than silently mis-parsing.
    """
    counts = {}
    i, n = 0, len(formula)
    while i < n:
        c = formula[i]
        if c.isspace():
            i += 1
            continue
        if not ("A" <= c <= "Z"):
            raise ValueError("bad formula token at %d in %r" % (i, formula))
        j = i + 1
        while j < n and "a" <= formula[j] <= "z":
            j += 1
        symbol = formula[i:j]
        i = j
        if i < n and formula[i] == "{":
            raise ValueError("isotope notation unsupported in L1 corpus: %r" % formula)
        sign = 1
        if i < n and formula[i] == "-":
            sign = -1
            i += 1
        k = i
        while i < n and formula[i].isdigit():
            i += 1
        count = int(formula[k:i]) if i > k else 1
        if symbol not in element_index:
            raise ValueError("unknown element %r in %r" % (symbol, formula))
        an = element_index[symbol][0]
        counts[an] = counts.get(an, 0) + sign * count
    return counts


def add_into(acc, other):
    for an, c in other.items():
        acc[an] = acc.get(an, 0) + c


def full_formula_counts(base_seq, combined_mod, element_index):
    """{atomic_number: count} for backbone (residues + H + OH) + combined mod formula."""
    counts = {}
    # N-terminus H + C-terminus OH = one water.
    add_into(counts, parse_formula("H", element_index))
    add_into(counts, parse_formula("OH", element_index))
    for letter in base_seq:
        if letter not in RESIDUES:
            raise ValueError("unknown residue %r in %r" % (letter, base_seq))
        add_into(counts, parse_formula(RESIDUES[letter], element_index))
    if combined_mod.strip():
        add_into(counts, parse_formula(combined_mod, element_index))
    # Mirror the C# removal threshold: drop elements whose running count is exactly 0.
    return {an: c for an, c in counts.items() if c != 0}


def mono_mass(counts, an_to_mass):
    total = 0.0
    for an, c in counts.items():
        total += an_to_mass[an] * c
    return total


def main():
    raw = extract_nist_literal(CS_PATH)
    elements = load_periodic_table(raw)
    element_index = build_element_index(elements)
    an_to_symbol = {an: sym for sym, (an, _m) in element_index.items()}
    an_to_mass = {an: m for sym, (an, m) in element_index.items()}

    with open(PSMTSV, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    # One golden entry per distinct (Base Sequence, Mods Combined Chemical Formula). The C#
    # mono mass cross-check column is identical across duplicate rows, so keep the first.
    distinct = {}
    order = []
    for r in rows:
        base = r["Base Sequence"]
        mod = r["Mods Combined Chemical Formula"]
        key = (base, mod)
        if key not in distinct:
            distinct[key] = r
            order.append(key)

    worst_xcheck = (None, 0.0)  # (key, abs error vs C# column)
    lines = []
    lines.append("# L1 golden — modified peptide -> element counts + monoisotopic mass")
    lines.append("# replicated from mzLib residue table + ChemicalFormula over the embedded")
    lines.append("# periodic table; mod formulas read from AllPSMs.psmtsv 'Mods Combined")
    lines.append("# Chemical Formula'. csharpMono is the psmtsv 'Peptide Monoisotopic Mass'.")
    lines.append("# P<TAB>baseSeq<TAB>combinedMod<TAB>monoMass(repr)<TAB>csharpMono<TAB>counts")
    lines.append("#   counts = space-joined symbol:count, sorted by atomic number")
    for key in order:
        base, mod = key
        r = distinct[key]
        counts = full_formula_counts(base, mod, element_index)
        mm = mono_mass(counts, an_to_mass)
        csharp = float(r["Peptide Monoisotopic Mass"])
        err = abs(mm - csharp)
        if err > worst_xcheck[1]:
            worst_xcheck = (key, err)
        counts_str = " ".join(
            "%s:%d" % (an_to_symbol[an], counts[an]) for an in sorted(counts)
        )
        mod_field = mod if mod.strip() else "-"
        lines.append("P\t%s\t%s\t%r\t%s\t%s" % (base, mod_field, mm, repr(csharp), counts_str))

    os.makedirs(os.path.dirname(GOLDEN_TSV), exist_ok=True)
    with open(GOLDEN_TSV, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")

    print("psmtsv rows: %d   distinct peptides: %d" % (len(rows), len(order)))
    wk, we = worst_xcheck
    print("worst mono-mass error vs C# 'Peptide Monoisotopic Mass' column: %.3e Da" % we)
    print("  at %r" % (wk,))
    print("wrote golden: %s" % os.path.relpath(GOLDEN_TSV, REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
