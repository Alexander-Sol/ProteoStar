//! L2 parity gate (PLAN.md P1.6): chemical formula -> isotopic distribution.
//!
//! For every distinct peptide in `mzLib/Test/FlashLFQ/TestData/AllPSMs.psmtsv`, the Rust
//! `IsotopicDistribution::get_distribution` (the Kubinyi fine-grained polynomial port) must
//! reproduce the `(masses[], intensities[])` envelope of the independent Python replica golden
//! (`parity/golden/L2_isotopic_distributions.tsv`, emitted by
//! `parity/dump_l2_isotopic_distributions.py`):
//!   * **intensities** — relative `1e-6`, and
//!   * **masses** — absolute `1e-6` Da,
//! with **array lengths exact**.
//!
//! Both sides are faithful re-implementations of `Chemistry/IsotopicDistribution.cs`. The
//! envelope convolution is floating-point order-sensitive, so the replica iterates elements in
//! the same **atomic-number order** the Rust `ChemicalFormula::elements()` BTreeMap uses,
//! making this an index-aligned, near-bit-identical comparison (see the dumper's docstring).
//!
//! As for L0/L1 the .NET SDK is absent, so the golden is the Python replica rather than a live
//! C# dump; agreement of the two independent re-implementations is the gate. The full formula
//! is the backbone + the psmtsv `Mods Combined Chemical Formula` column, exactly as in L1.
//!
//! The distinct peptide set is driven straight from the live psmtsv and asserted to match the
//! golden set, so the gate tracks the corpus file rather than a frozen copy.

use std::collections::BTreeMap;

use flashlfq_core::chemical_formula::ChemicalFormula;
use flashlfq_core::isotopic_distribution::IsotopicDistribution;
use flashlfq_core::peptide::peptide_formula_with_mods;

/// The golden L2 artifact emitted by `parity/dump_l2_isotopic_distributions.py`.
const GOLDEN: &str = include_str!("../parity/golden/L2_isotopic_distributions.tsv");

/// Intensities are compared at relative `1e-6`.
const INTENSITY_REL_TOL: f64 = 1e-6;
/// Masses are compared at absolute `1e-6` Da.
const MASS_ABS_TOL: f64 = 1e-6;

struct GoldenEntry {
    base_seq: String,
    combined_mod: String,
    masses: Vec<f64>,
    intensities: Vec<f64>,
}

/// Parse the D/M/I record stream into one entry per peptide.
fn parse_golden() -> Vec<GoldenEntry> {
    let mut out: Vec<GoldenEntry> = Vec::new();
    let mut lines = GOLDEN.lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty());

    while let Some(line) = lines.next() {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[0], "D", "expected a D (descriptor) record, got {:?}", line);
        let combined_mod = if f[2] == "-" { String::new() } else { f[2].to_string() };
        let n_peaks: usize = f[3].parse().expect("nPeaks parses");

        let m_line = lines.next().expect("M record follows D record");
        let masses = parse_floats(m_line, 'M');
        let i_line = lines.next().expect("I record follows M record");
        let intensities = parse_floats(i_line, 'I');

        assert_eq!(masses.len(), n_peaks, "M record length matches nPeaks for {:?}", f[1]);
        assert_eq!(intensities.len(), n_peaks, "I record length matches nPeaks for {:?}", f[1]);

        out.push(GoldenEntry {
            base_seq: f[1].to_string(),
            combined_mod,
            masses,
            intensities,
        });
    }
    out
}

/// Parse a `TAG<TAB>v v v ...` record's space-joined repr() floats.
fn parse_floats(line: &str, tag: char) -> Vec<f64> {
    let (t, rest) = line.split_once('\t').expect("record has a tab");
    assert_eq!(t.chars().next().unwrap(), tag, "record tag {tag}");
    rest.split(' ')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<f64>().expect("float repr parses"))
        .collect()
}

/// Distinct `(base_sequence, combined_mod)`, first-seen order, read straight from the live
/// psmtsv so the gate tracks the corpus file, not a frozen copy. (Same logic as the L1 gate.)
fn read_psmtsv_distinct() -> Vec<(String, String)> {
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

    let mut seen: BTreeMap<(String, String), ()> = BTreeMap::new();
    let mut order = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let base = fields[c_base].to_string();
        let combined_mod = fields[c_mod].trim().to_string();
        let key = (base.clone(), combined_mod.clone());
        if seen.insert(key, ()).is_none() {
            order.push((base, combined_mod));
        }
    }
    order
}

#[test]
fn l2_isotopic_distributions_match_golden() {
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

    let mut worst_mass_abs = 0.0_f64;
    let mut worst_int_rel = 0.0_f64;
    let mut total_peaks = 0usize;

    for (base, combined_mod) in &distinct {
        let g = golden_by_key
            .get(&(base.as_str(), combined_mod.as_str()))
            .unwrap_or_else(|| panic!("no golden entry for ({base:?}, {combined_mod:?})"));

        let mods: Vec<ChemicalFormula> = if combined_mod.is_empty() {
            Vec::new()
        } else {
            vec![ChemicalFormula::parse_formula(combined_mod)
                .unwrap_or_else(|e| panic!("parsing mod formula {combined_mod:?}: {e}"))]
        };
        let formula = peptide_formula_with_mods(base, &mods)
            .unwrap_or_else(|e| panic!("building formula for {base:?}: {e}"));

        let dist = IsotopicDistribution::get_distribution(&formula);

        // Array lengths exact: same algorithm, same element order → same peak count.
        assert_eq!(
            dist.masses.len(),
            g.masses.len(),
            "peak count for {base:?} + {combined_mod:?}: rust {} vs golden {}",
            dist.masses.len(),
            g.masses.len()
        );
        assert_eq!(dist.masses.len(), dist.intensities.len());

        for i in 0..g.masses.len() {
            let mass_abs = (dist.masses[i] - g.masses[i]).abs();
            worst_mass_abs = worst_mass_abs.max(mass_abs);
            assert!(
                mass_abs < MASS_ABS_TOL,
                "mass[{i}] for {base:?} + {combined_mod:?}: rust={}, golden={}, abs={mass_abs:e}",
                dist.masses[i],
                g.masses[i]
            );

            let denom = g.intensities[i].abs().max(dist.intensities[i].abs());
            let int_rel = if denom > 0.0 {
                (dist.intensities[i] - g.intensities[i]).abs() / denom
            } else {
                0.0
            };
            worst_int_rel = worst_int_rel.max(int_rel);
            assert!(
                int_rel < INTENSITY_REL_TOL,
                "intensity[{i}] for {base:?} + {combined_mod:?}: rust={}, golden={}, rel={int_rel:e}",
                dist.intensities[i],
                g.intensities[i]
            );
        }
        total_peaks += g.masses.len();
    }

    eprintln!(
        "L2 parity: {} peptides, {} total peaks — worst mass abs {:e} Da, worst intensity rel {:e}",
        distinct.len(),
        total_peaks,
        worst_mass_abs,
        worst_int_rel
    );
}
