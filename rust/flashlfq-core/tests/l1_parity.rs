//! L1 parity gate (PLAN.md P1.4): modified sequence -> element counts + monoisotopic mass.
//!
//! For every distinct peptide in `mzLib/Test/FlashLFQ/TestData/AllPSMs.psmtsv`, the Rust
//! chemistry port (`peptide::peptide_formula_with_mods` over `chemical_formula` + the L0
//! `periodic_table`) must reproduce:
//!   * **element counts** — exact, and
//!   * **monoisotopic mass** — relative `1e-9`,
//! versus the independent Python replica golden (`parity/golden/L1_peptide_formulas.tsv`,
//! emitted by `parity/dump_l1_peptide_formulas.py`).
//!
//! The full formula is the backbone (`new Peptide(baseSeq).GetChemicalFormula()`) plus the
//! modifications, whose combined formula mzLib already wrote into the psmtsv column
//! **Mods Combined Chemical Formula** — so no modified-sequence string parser is needed
//! (P1.3 deferred that; here we read the column mzLib resolved). This mirrors C#
//! `PeptideWithSetModifications.FullChemicalFormula`.
//!
//! Two independent checks make this a real parity gate despite the absent .NET SDK:
//!   1. Rust vs the Python replica golden (counts exact, mass rel-1e-9) — two independent
//!      re-implementations of the same arithmetic over the same embedded periodic table.
//!   2. Rust vs the psmtsv **Peptide Monoisotopic Mass** column — the genuine C# output,
//!      compared at its 5-decimal print rounding (abs `< 1e-5`; observed worst 4.98e-6).
//!
//! The test also drives the *distinct peptide set straight from the live psmtsv* and asserts
//! it matches the golden set, so the gate tracks the corpus file rather than a frozen copy.

use std::collections::BTreeMap;

use flashlfq_core::chemical_formula::ChemicalFormula;
use flashlfq_core::peptide::peptide_formula_with_mods;
use flashlfq_core::periodic_table::periodic_table;

/// The golden L1 artifact emitted by `parity/dump_l1_peptide_formulas.py`.
const GOLDEN: &str = include_str!("../parity/golden/L1_peptide_formulas.tsv");

/// Counts exact; the L1 mass tolerance is relative `1e-9`.
const MASS_REL_TOL: f64 = 1e-9;
/// Cross-check vs the C# psmtsv column, bounded by its 5-decimal print rounding.
const CSHARP_ABS_TOL: f64 = 1e-5;

struct GoldenEntry {
    base_seq: String,
    combined_mod: String,
    mono_mass: f64,
    csharp_mono: f64,
    /// (atomic_number, count), sorted by atomic number — exact element counts.
    counts: Vec<(u16, i32)>,
}

fn parse_golden() -> Vec<GoldenEntry> {
    let mut out = Vec::new();
    for line in GOLDEN.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[0], "P", "unexpected golden record tag {:?}", f[0]);
        let combined_mod = if f[2] == "-" { String::new() } else { f[2].to_string() };
        let counts = f[5]
            .split(' ')
            .filter(|s| !s.is_empty())
            .map(|tok| {
                let (sym, cnt) = tok.split_once(':').expect("count token is symbol:count");
                let an = periodic_table()
                    .element_by_symbol(sym)
                    .unwrap_or_else(|| panic!("golden symbol {sym} not in periodic table"))
                    .atomic_number;
                (an, cnt.parse::<i32>().expect("integer count"))
            })
            .collect();
        out.push(GoldenEntry {
            base_seq: f[1].to_string(),
            combined_mod,
            mono_mass: f[3].parse().expect("mono mass repr parses"),
            csharp_mono: f[4].parse().expect("csharp mono parses"),
            counts,
        });
    }
    out
}

/// Distinct `(base_sequence, combined_mod) -> csharp_mono`, first-seen order, read straight
/// from the live psmtsv so the gate tracks the corpus file, not a frozen copy.
fn read_psmtsv_distinct() -> Vec<(String, String, f64)> {
    let path = flashlfq_core::mzlib_test_data("AllPSMs.psmtsv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));

    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("header line").split('\t').collect();
    let col = |name: &str| {
        header
            .iter()
            .position(|h| *h == name)
            .unwrap_or_else(|| panic!("psmtsv missing column {name:?}"))
    };
    let c_base = col("Base Sequence");
    let c_mod = col("Mods Combined Chemical Formula");
    let c_mass = col("Peptide Monoisotopic Mass");

    let mut seen: BTreeMap<(String, String), ()> = BTreeMap::new();
    let mut order = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let base = fields[c_base].to_string();
        let combined_mod = fields[c_mod].trim().to_string();
        let csharp: f64 = fields[c_mass].parse().expect("Peptide Monoisotopic Mass parses");
        let key = (base.clone(), combined_mod.clone());
        if seen.insert(key, ()).is_none() {
            order.push((base, combined_mod, csharp));
        }
    }
    order
}

#[test]
fn l1_modified_peptide_formulas_match_golden_and_csharp() {
    let golden = parse_golden();
    let golden_by_key: BTreeMap<(&str, &str), &GoldenEntry> = golden
        .iter()
        .map(|g| ((g.base_seq.as_str(), g.combined_mod.as_str()), g))
        .collect();

    let distinct = read_psmtsv_distinct();
    assert_eq!(
        distinct.len(),
        golden.len(),
        "distinct peptide count: psmtsv {} vs golden {} (regenerate the golden if the corpus changed)",
        distinct.len(),
        golden.len(),
    );

    let mut worst_mass_rel = 0.0_f64;
    let mut worst_csharp_abs = 0.0_f64;

    for (base, combined_mod, csharp_mono) in &distinct {
        let g = golden_by_key
            .get(&(base.as_str(), combined_mod.as_str()))
            .unwrap_or_else(|| panic!("no golden entry for ({base:?}, {combined_mod:?})"));

        // The golden must have been generated from this same psmtsv row.
        assert_eq!(
            g.csharp_mono, *csharp_mono,
            "golden C# mono disagrees with live psmtsv for {base:?}"
        );

        // Build the Rust formula: backbone + the combined modification formula (if any).
        let mods: Vec<ChemicalFormula> = if combined_mod.is_empty() {
            Vec::new()
        } else {
            vec![ChemicalFormula::parse_formula(combined_mod)
                .unwrap_or_else(|e| panic!("parsing mod formula {combined_mod:?}: {e}"))]
        };
        let formula = peptide_formula_with_mods(base, &mods)
            .unwrap_or_else(|e| panic!("building formula for {base:?}: {e}"));

        // Element counts: exact, every golden element, with no extras. The corpus carries no
        // isotope-specified atoms, so atom_count equals the sum of element counts — comparing
        // it catches any element the Rust side has that the golden lacks.
        let mut expected_atoms = 0_i32;
        for (an, count) in &g.counts {
            assert_eq!(
                formula.count_of_element(*an),
                *count,
                "element count (Z={an}) for {base:?} + {combined_mod:?}"
            );
            expected_atoms += *count;
        }
        assert_eq!(
            formula.atom_count(),
            expected_atoms,
            "total atom count for {base:?} + {combined_mod:?} (extra element present?)"
        );

        // Monoisotopic mass: rel-1e-9 vs the Python replica.
        let mono = formula.monoisotopic_mass();
        let rel = (mono - g.mono_mass).abs() / g.mono_mass.abs().max(1.0);
        worst_mass_rel = worst_mass_rel.max(rel);
        assert!(
            rel < MASS_REL_TOL,
            "mono mass vs golden for {base:?} + {combined_mod:?}: rust={mono}, golden={}, rel={rel:e}",
            g.mono_mass
        );

        // Cross-check vs the genuine C# output, bounded by its print rounding.
        let abs = (mono - csharp_mono).abs();
        worst_csharp_abs = worst_csharp_abs.max(abs);
        assert!(
            abs < CSHARP_ABS_TOL,
            "mono mass vs C# psmtsv column for {base:?} + {combined_mod:?}: rust={mono}, csharp={csharp_mono}, abs={abs:e}"
        );
    }

    eprintln!(
        "L1 parity: {} distinct peptides — worst mass rel vs golden {:e}, worst abs vs C# {:e} Da",
        distinct.len(),
        worst_mass_rel,
        worst_csharp_abs
    );
}
