//! L3 parity gate (PLAN.md P1.8): modified sequence -> final `(massShift, normAbundance)[]`
//! and `PeakfindingMass`. **This is the chemistry gate — all later (tracing-core) work depends
//! on it passing.**
//!
//! For every distinct modified peptide (Full Sequence) in
//! `mzLib/Test/FlashLFQ/TestData/AllPSMs.psmtsv`, the Rust
//! `theoretical_isotope_distribution::expected_isotope_peaks` (the
//! `CalculateTheoreticalIsotopeDistributions` port) must reproduce the expected isotope peaks and
//! peakfinding mass of the independent Python replica golden
//! (`parity/golden/L3_theoretical_distributions.tsv`, emitted by
//! `parity/dump_l3_theoretical_distributions.py`):
//!   * **mass shifts** — absolute `1e-6` Da,
//!   * **abundances** — relative `1e-6`,
//!   * **peakfinding mass** — absolute `1e-6` Da,
//! with **array lengths exact**.
//!
//! **The formula resolution mirrors C# exactly.** FlashLFQ leaves
//! `Identification.OptionalChemicalFormula` null on the psmtsv path, so the formula is built from
//! the **base sequence only** (no mods) and topped up with averagine when the known mass differs
//! from the bare-formula mass by more than 20 Da. The known monoisotopic mass driving both the
//! re-centering and the 20 Da test is the psmtsv **`Peptide Monoisotopic Mass`** column
//! (= `Identification.MonoisotopicMass`), parsed to `f64` identically on both sides — so the
//! averagine top-up fires identically. (The `Mods Combined Chemical Formula` column that drives
//! L1/L2 is deliberately NOT used here.)
//!
//! As for L0/L1/L2 the .NET SDK is absent, so the golden is the Python replica rather than a live
//! C# dump; agreement of the two independent re-implementations is the gate. The distinct modified
//! sequence set is driven straight from the live psmtsv and asserted to match the golden set, so
//! the gate tracks the corpus file rather than a frozen copy.

use std::collections::BTreeMap;

use flashlfq_core::theoretical_isotope_distribution::{
    expected_isotope_peaks, peakfinding_mass, DEFAULT_NUM_ISOTOPES_REQUIRED,
};

/// The golden L3 artifact emitted by `parity/dump_l3_theoretical_distributions.py`.
const GOLDEN: &str = include_str!("../parity/golden/L3_theoretical_distributions.tsv");

/// Mass shifts and the peakfinding mass are compared at absolute `1e-6` Da.
const MASS_ABS_TOL: f64 = 1e-6;
/// Abundances are compared at relative `1e-6`.
const ABUNDANCE_REL_TOL: f64 = 1e-6;

struct GoldenEntry {
    full_seq: String,
    base_seq: String,
    mono: f64,
    peakfinding_mass: f64,
    mass_shifts: Vec<f64>,
    abundances: Vec<f64>,
}

/// Parse the D/S/A record stream into one entry per modified peptide.
fn parse_golden() -> Vec<GoldenEntry> {
    let mut out: Vec<GoldenEntry> = Vec::new();
    let mut lines = GOLDEN.lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty());

    while let Some(line) = lines.next() {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[0], "D", "expected a D (descriptor) record, got {:?}", line);
        let full_seq = f[1].to_string();
        let base_seq = f[2].to_string();
        let mono: f64 = f[3].parse().expect("mono parses");
        let n_peaks: usize = f[4].parse().expect("nPeaks parses");
        let peakfinding_mass: f64 = f[5].parse().expect("peakfinding mass parses");

        let s_line = lines.next().expect("S record follows D record");
        let mass_shifts = parse_floats(s_line, 'S');
        let a_line = lines.next().expect("A record follows S record");
        let abundances = parse_floats(a_line, 'A');

        assert_eq!(mass_shifts.len(), n_peaks, "S record length matches nPeaks for {full_seq:?}");
        assert_eq!(abundances.len(), n_peaks, "A record length matches nPeaks for {full_seq:?}");

        out.push(GoldenEntry {
            full_seq,
            base_seq,
            mono,
            peakfinding_mass,
            mass_shifts,
            abundances,
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

/// Distinct `(full_sequence, base_sequence, monoisotopic_mass)`, first-seen order, read straight
/// from the live psmtsv so the gate tracks the corpus file. Keyed on the Full Sequence (modified
/// sequence), matching the C# `ModifiedSequenceToIsotopicDistribution` dictionary.
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
    let c_full = col("Full Sequence");
    let c_base = col("Base Sequence");
    let c_mono = col("Peptide Monoisotopic Mass");

    let mut seen: BTreeMap<String, ()> = BTreeMap::new();
    let mut order = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let full = fields[c_full].to_string();
        if seen.insert(full.clone(), ()).is_none() {
            let base = fields[c_base].to_string();
            let mono: f64 = fields[c_mono]
                .parse()
                .unwrap_or_else(|e| panic!("parsing mono {:?}: {e}", fields[c_mono]));
            order.push((full, base, mono));
        }
    }
    order
}

#[test]
fn l3_theoretical_distributions_match_golden() {
    let golden = parse_golden();
    let golden_by_key: BTreeMap<&str, &GoldenEntry> =
        golden.iter().map(|g| (g.full_seq.as_str(), g)).collect();

    let distinct = read_psmtsv_distinct();
    assert_eq!(
        distinct.len(),
        golden.len(),
        "distinct modified-sequence count: psmtsv {} vs golden {} (regenerate the golden if the corpus changed)",
        distinct.len(),
        golden.len(),
    );

    let mut worst_shift_abs = 0.0_f64;
    let mut worst_abundance_rel = 0.0_f64;
    let mut worst_pf_abs = 0.0_f64;
    let mut total_peaks = 0usize;

    for (full, base, mono) in &distinct {
        let g = golden_by_key
            .get(full.as_str())
            .unwrap_or_else(|| panic!("no golden entry for {full:?}"));

        // The golden's base/mono should match the psmtsv (same source); guard against drift.
        assert_eq!(&g.base_seq, base, "golden base sequence for {full:?}");
        assert_eq!(g.mono.to_bits(), mono.to_bits(), "golden mono for {full:?}");

        let peaks = expected_isotope_peaks(None, base, *mono, DEFAULT_NUM_ISOTOPES_REQUIRED);
        let pf = peakfinding_mass(*mono, &peaks);

        // Array lengths exact: same algorithm, same element order, same truncation rule.
        assert_eq!(
            peaks.len(),
            g.mass_shifts.len(),
            "peak count for {full:?}: rust {} vs golden {}",
            peaks.len(),
            g.mass_shifts.len()
        );

        for i in 0..g.mass_shifts.len() {
            let shift_abs = (peaks[i].mass_shift - g.mass_shifts[i]).abs();
            worst_shift_abs = worst_shift_abs.max(shift_abs);
            assert!(
                shift_abs < MASS_ABS_TOL,
                "massShift[{i}] for {full:?}: rust={}, golden={}, abs={shift_abs:e}",
                peaks[i].mass_shift,
                g.mass_shifts[i]
            );

            let denom = g.abundances[i].abs().max(peaks[i].normalized_abundance.abs());
            let abundance_rel = if denom > 0.0 {
                (peaks[i].normalized_abundance - g.abundances[i]).abs() / denom
            } else {
                0.0
            };
            worst_abundance_rel = worst_abundance_rel.max(abundance_rel);
            assert!(
                abundance_rel < ABUNDANCE_REL_TOL,
                "abundance[{i}] for {full:?}: rust={}, golden={}, rel={abundance_rel:e}",
                peaks[i].normalized_abundance,
                g.abundances[i]
            );
        }

        let pf_abs = (pf - g.peakfinding_mass).abs();
        worst_pf_abs = worst_pf_abs.max(pf_abs);
        assert!(
            pf_abs < MASS_ABS_TOL,
            "peakfinding mass for {full:?}: rust={pf}, golden={}, abs={pf_abs:e}",
            g.peakfinding_mass
        );

        total_peaks += g.mass_shifts.len();
    }

    eprintln!(
        "L3 parity: {} modified peptides, {} total peaks — worst massShift abs {:e} Da, \
         worst abundance rel {:e}, worst peakfinding mass abs {:e} Da",
        distinct.len(),
        total_peaks,
        worst_shift_abs,
        worst_abundance_rel,
        worst_pf_abs
    );
}
