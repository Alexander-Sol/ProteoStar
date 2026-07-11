#!/usr/bin/env python3
"""L0 golden dumper — periodic-table isotope data.

Produces the L0 parity artifact (`golden/L0_periodic_table.tsv`) by replicating
mzLib's `Chemistry/PeriodicTable.cs` static loader *exactly*. The plan (PLAN.md P1.1)
calls for a C# dumper, but the .NET SDK is not available in this environment, so we
re-run the loader's deterministic parsing logic here instead. The parsing is a faithful
line-by-line transcription of the C# `static PeriodicTable()` constructor; see the
inline comments mapping each branch to the C# source.

The source decimal strings come from the *verbatim* `thePeriodicTable` literal embedded
in PeriodicTable.cs (extracted to `data/periodic_table.nist.txt`), NOT re-sourced from
NIST — parity is against what mzLib loads, per the feasibility doc.

Two consumers read the same `data/periodic_table.nist.txt`:
  * this script (golden / reference loader), and
  * the Rust `periodic_table` module (the embedded table under test).
The Rust L0 test diffs its loaded table against the golden tsv this emits; agreement of
the two independent re-implementations + `validate_abundances` is the L0 gate.

Run from anywhere:  python rust/flashlfq-core/parity/dump_periodic_table.py
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CORE = os.path.dirname(HERE)                       # rust/flashlfq-core
REPO = os.environ.get("MZLIB_DIR", r"F:\mzLib")    # external mzLib checkout (see README)
CS_PATH = os.path.join(REPO, "mzLib", "Chemistry", "PeriodicTable.cs")
NIST_TXT = os.path.join(CORE, "data", "periodic_table.nist.txt")
GOLDEN_TSV = os.path.join(HERE, "golden", "L0_periodic_table.tsv")

_NUM = re.compile(r"[\d\.]+")
_AFTER_BRACKET = re.compile(r"(?<=\[)[\d\.]+")
_AFTER_COMMA = re.compile(r"(?<=,)[\d\.]+")


def extract_nist_literal(cs_path):
    """Pull the verbatim `thePeriodicTable` string literal out of PeriodicTable.cs.

    The literal is a C# verbatim string `@"..."`; the data contains no embedded quotes,
    so the first `"` after the opening `@"` terminates it.
    """
    with open(cs_path, "r", encoding="utf-8") as f:
        text = f.read()
    # Anchor on the const *declaration*, not its earlier use in the constructor — the
    # constructor body contains `@"..."` regex literals that would otherwise match first.
    anchor = text.index("private const string thePeriodicTable")
    start = text.index('@"', anchor) + 2
    end = text.index('"', start)
    return text[start:end]


def load_periodic_table(raw):
    """Replicate the C# `static PeriodicTable()` constructor.

    Returns a list of element dicts in insertion order (atomic-number order as they
    first appear in the data), each:
        {symbol, atomic_number, average_mass, principal_mass_number,
         isotopes: [(mass_number, atomic_mass, abundance), ...]}  # added-order
    """
    elements_by_number = {}   # atomic_number -> element dict  (mirrors _elementsArray)
    order = []                # insertion order of elements

    atomic_number = symbol = isotope_number = None
    atomic_mass = abundance = average_mass = None

    for line in raw.split("\n"):
        if not line.strip():
            continue
        key = line[: line.index("=")].strip()
        value = line[line.index("=") + 1 :].strip()

        if key == "Atomic Number":
            atomic_number = int(value)
        elif key == "Atomic Symbol":
            symbol = value
        elif key == "Mass Number":
            isotope_number = int(value)
        elif key == "Relative Atomic Mass":
            atomic_mass = float(_NUM.search(line).group())
        elif key == "Isotopic Composition":
            m = _NUM.search(line)
            if m:
                abundance = float(m.group())
        elif key == "Standard Atomic Weight":
            if "[" in line:
                avg1 = float(_AFTER_BRACKET.search(line).group())
                m2 = _AFTER_COMMA.search(line)
                if m2:  # C#: double.TryParse(...) succeeded
                    avg2 = float(m2.group())
                    average_mass = (avg1 + avg2) / 2
                else:
                    average_mass = avg1
            else:
                m = _NUM.search(line)
                if m:
                    average_mass = float(m.group())
        elif key == "Notes":
            existing = elements_by_number.get(atomic_number)
            if existing is None:
                # C#: create the element only if every field is present (abundance incl.)
                if (
                    symbol is not None
                    and atomic_number is not None
                    and average_mass is not None
                    and isotope_number is not None
                    and atomic_mass is not None
                    and abundance is not None
                ):
                    el = {
                        "symbol": symbol,
                        "atomic_number": atomic_number,
                        "average_mass": average_mass,
                        "isotopes": [],
                        "principal_mass_number": None,
                        "_principal_abundance": None,
                    }
                    _add_isotope(el, isotope_number, atomic_mass, abundance)
                    elements_by_number[atomic_number] = el
                    order.append(el)
            elif abundance is not None:
                _add_isotope(existing, isotope_number, atomic_mass, abundance)

            atomic_number = symbol = isotope_number = None
            atomic_mass = abundance = average_mass = None
        else:
            raise ValueError("Could not parse line from periodic table: " + line)

    order.sort(key=lambda e: e["atomic_number"])
    return order


def _add_isotope(el, mass_number, atomic_mass, abundance):
    """Mirror Element.AddIsotope: append in added-order, track strict-greater principal."""
    el["isotopes"].append((mass_number, atomic_mass, abundance))
    if el["_principal_abundance"] is None or abundance > el["_principal_abundance"]:
        el["_principal_abundance"] = abundance
        el["principal_mass_number"] = mass_number


def validate_abundances(elements, epsilon):
    """Mirror PeriodicTable.ValidateAbundances: per-element sum-of-abundances ~= 1."""
    worst_sym, worst_err = None, 0.0
    for el in elements:
        total = 0.0  # sum in added-order, as C# Isotopes enumerates
        for (_mn, _am, ab) in el["isotopes"]:
            total += ab
        err = abs(total - 1.0)
        if err > worst_err:
            worst_err, worst_sym = err, el["symbol"]
        if err > epsilon:
            return False, worst_sym, worst_err
    return True, worst_sym, worst_err


def write_golden(elements, path):
    lines = []
    lines.append("# L0 periodic-table golden dump — replicated from mzLib Chemistry/PeriodicTable.cs")
    lines.append("# E<TAB>symbol<TAB>atomicNumber<TAB>averageMass<TAB>principalMassNumber<TAB>isotopeCount")
    lines.append("# I<TAB>massNumber<TAB>atomicMass<TAB>relativeAbundance   (C# added-order)")
    lines.append("# floats are Python repr() (shortest round-trip); Rust parses them back bit-identically")
    for el in elements:
        lines.append(
            "E\t{symbol}\t{atomic_number}\t{average_mass!r}\t{principal_mass_number}\t{n}".format(
                n=len(el["isotopes"]), **el
            )
        )
        for (mn, am, ab) in el["isotopes"]:
            lines.append("I\t{}\t{!r}\t{!r}".format(mn, am, ab))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # Force \n endings so the byte-for-byte Rust diff is platform-independent.
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


def main():
    raw = extract_nist_literal(CS_PATH)
    os.makedirs(os.path.dirname(NIST_TXT), exist_ok=True)
    with open(NIST_TXT, "w", encoding="utf-8", newline="\n") as f:
        f.write(raw)

    elements = load_periodic_table(raw)
    ok, worst_sym, worst_err = validate_abundances(elements, 1e-15)
    n_iso = sum(len(e["isotopes"]) for e in elements)
    print("loaded elements: {}  isotopes: {}".format(len(elements), n_iso))
    print("validate_abundances(1e-15): {}  (worst: {} err={:.3e})".format(ok, worst_sym, worst_err))
    if not ok:
        print("ERROR: abundance validation failed", file=sys.stderr)
        return 1

    write_golden(elements, GOLDEN_TSV)
    print("wrote golden: {}".format(os.path.relpath(GOLDEN_TSV, REPO)))
    print("wrote nist  : {}".format(os.path.relpath(NIST_TXT, REPO)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
