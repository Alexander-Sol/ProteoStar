//! Theoretical isotope distributions for FlashLFQ identifications (parity layer **L3**).
//!
//! A port of `FlashLfqEngine.CalculateTheoreticalIsotopeDistributions()`. For each peptide it
//! turns a base sequence (and known monoisotopic mass) into the list of *expected isotope mass
//! shifts and normalized abundances* that the tracing core later searches for in MS1 scans, plus
//! the **peakfinding mass** (the monoisotopic mass shifted onto the most abundant isotope).
//!
//! The pipeline per peptide:
//! 1. **Resolve a chemical formula.** If the identification already carries one, use it. Otherwise
//!    build it from the base sequence ([`crate::peptide::peptide_base_formula`]); if the
//!    sequence/known-mass disagree by more than 20 Da (an unlocalized/unknown modification), top
//!    the formula up with **averagine** — an averaged amino-acid unit — scaled to cover the gap.
//!    If the sequence has non-residue characters, build the whole formula from averagine.
//! 2. **Compute the fine isotopic distribution** of that formula
//!    ([`IsotopicDistribution::get_distribution_with`] at `fineResolution = 0.125`,
//!    `minProbability = 1e-8`).
//! 3. **Re-center on the known monoisotopic mass** (the formula's mono mass may differ by the
//!    averagine residual), convert each peak to a *mass shift* relative to the monoisotopic mass,
//!    and **normalize abundances to the most abundant peak**.
//! 4. **Truncate**: keep a peak iff fewer than `NumIsotopesRequired` are kept so far, or its
//!    normalized abundance exceeds `0.1`.
//!
//! The averagine coefficients, the 20 Da threshold, the `0.125 / 1e-8` distribution parameters,
//! the `Math.Round` (banker's rounding) integer counts, and the `0.1` truncation cutoff are all
//! preserved verbatim from the C# so the L3 parity gate (P1.8) can be exact.

use crate::chemical_formula::ChemicalFormula;
use crate::isotopic_distribution::{IsotopicDistribution, DEFAULT_MOLECULAR_WEIGHT_RESOLUTION};
use crate::peptide::peptide_base_formula;
use crate::periodic_table::periodic_table;

/// Averagine composition — the average number of each element per "averagine" unit, from
/// `CalculateTheoreticalIsotopeDistributions`. Used to fill in the formula of a peptide whose
/// modification chemistry is unknown.
pub const AVERAGE_C: f64 = 4.9384;
pub const AVERAGE_H: f64 = 7.7583;
pub const AVERAGE_O: f64 = 1.4773;
pub const AVERAGE_N: f64 = 1.3577;
pub const AVERAGE_S: f64 = 0.0417;

/// Distribution parameters used by `CalculateTheoreticalIsotopeDistributions` — note these are
/// *not* the `IsotopicDistribution` defaults: FlashLFQ calls `GetDistribution(formula, 0.125,
/// 1e-8)` (the molecular-weight resolution stays at the default `1e-12`).
pub const FINE_RESOLUTION: f64 = 0.125;
pub const MIN_PROBABILITY: f64 = 1e-8;

/// The default number of isotope peaks FlashLFQ requires (`FlashLfqParameters.NumIsotopesRequired`).
pub const DEFAULT_NUM_ISOTOPES_REQUIRED: usize = 2;

/// Mass-difference threshold (Da): if the known mass and the bare-sequence formula mass disagree
/// by more than this, the peptide carries an unaccounted modification and the formula is topped up
/// with averagine.
const AVERAGINE_MASS_DIFF_THRESHOLD: f64 = 20.0;

/// One expected isotope peak: a mass shift relative to the monoisotopic mass, and an abundance
/// normalized so the most abundant peak is `1.0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExpectedIsotopePeak {
    /// Mass shift from the peptide's monoisotopic mass (Da). The monoisotopic peak itself is the
    /// shift closest to 0; heavier isotopologues have positive shifts.
    pub mass_shift: f64,
    /// Abundance normalized to the most abundant isotope (which is exactly `1.0`).
    pub normalized_abundance: f64,
}

/// The mass of one averagine unit: `Σ averageElementCount * element.AverageMass` over C, H, O, N, S.
///
/// Mirrors the `averagineMass` computed at the top of `CalculateTheoreticalIsotopeDistributions`.
pub fn averagine_mass() -> f64 {
    let pt = periodic_table();
    let avg = |symbol: &str| {
        pt.element_by_symbol(symbol)
            .unwrap_or_else(|| panic!("averagine element {symbol} missing from periodic table"))
            .average_mass
    };
    avg("C") * AVERAGE_C
        + avg("H") * AVERAGE_H
        + avg("O") * AVERAGE_O
        + avg("N") * AVERAGE_N
        + avg("S") * AVERAGE_S
}

/// `true` iff every character of `base_sequence` is a known amino-acid residue and the sequence is
/// non-empty. Mirrors `AminoAcidPolymerExtensions.AllSequenceResiduesAreValid`.
pub fn all_sequence_residues_are_valid(base_sequence: &str) -> bool {
    !base_sequence.is_empty()
        && base_sequence
            .chars()
            .all(|c| crate::peptide::residue_formula(c).is_some())
}

/// Add `averagines` worth of averagine atoms to `formula`, with C# banker's-rounding integer
/// counts. Mirrors the five `formula.Add("X", (int)Math.Round(averagines * averageX, 0))` calls.
fn add_averagine(formula: &mut ChemicalFormula, averagines: f64) {
    let pt = periodic_table();
    let mut add = |symbol: &str, average_count: f64| {
        let element = pt
            .element_by_symbol(symbol)
            .unwrap_or_else(|| panic!("averagine element {symbol} missing from periodic table"));
        // C# `(int)Math.Round(value, 0)` is round-half-to-even (banker's rounding).
        let count = (averagines * average_count).round_ties_even() as i32;
        formula.add_element(element.atomic_number, count);
    };
    add("C", AVERAGE_C);
    add("H", AVERAGE_H);
    add("O", AVERAGE_O);
    add("N", AVERAGE_N);
    add("S", AVERAGE_S);
}

/// Resolve the chemical formula used to compute a peptide's isotope distribution.
///
/// Mirrors the formula-resolution block of `CalculateTheoreticalIsotopeDistributions`:
/// - an explicit `optional_formula` is used as-is;
/// - otherwise, for an all-valid sequence, the bare-sequence formula is topped up with averagine
///   when the known mass differs by more than 20 Da;
/// - otherwise (non-residue characters) the formula is built entirely from averagine.
fn resolve_formula(
    optional_formula: Option<&ChemicalFormula>,
    base_sequence: &str,
    monoisotopic_mass: f64,
) -> ChemicalFormula {
    if let Some(formula) = optional_formula {
        return formula.clone();
    }

    let avg_mass = averagine_mass();

    if all_sequence_residues_are_valid(base_sequence) {
        // A valid sequence always parses (every char is a residue).
        let mut formula = peptide_base_formula(base_sequence)
            .expect("all-valid base sequence must build a peptide formula");
        let mass_diff = monoisotopic_mass - formula.monoisotopic_mass();

        if mass_diff.abs() > AVERAGINE_MASS_DIFF_THRESHOLD {
            add_averagine(&mut formula, mass_diff / avg_mass);
        }
        formula
    } else {
        let mut formula = ChemicalFormula::new();
        add_averagine(&mut formula, monoisotopic_mass / avg_mass);
        formula
    }
}

/// Compute the expected isotope peaks for one peptide.
///
/// Faithful port of the per-identification body of `CalculateTheoreticalIsotopeDistributions`:
/// resolve the formula, compute the fine distribution at `(0.125, 1e-8)`, re-center on
/// `monoisotopic_mass`, normalize to the most abundant peak, and truncate by
/// `num_isotopes_required` / the `0.1` abundance cutoff.
///
/// `optional_formula` corresponds to `Identification.OptionalChemicalFormula` (pass `None` to
/// build from the sequence). The returned peaks are in distribution (mass-ascending) order, with
/// the most abundant peak's `normalized_abundance == 1.0`.
pub fn expected_isotope_peaks(
    optional_formula: Option<&ChemicalFormula>,
    base_sequence: &str,
    monoisotopic_mass: f64,
    num_isotopes_required: usize,
) -> Vec<ExpectedIsotopePeak> {
    let formula = resolve_formula(optional_formula, base_sequence, monoisotopic_mass);

    let distribution = IsotopicDistribution::get_distribution_with(
        &formula,
        FINE_RESOLUTION,
        MIN_PROBABILITY,
        DEFAULT_MOLECULAR_WEIGHT_RESOLUTION,
    );

    let mut masses = distribution.masses.clone();
    let mut abundances = distribution.intensities.clone();

    // Re-center the distribution masses onto the known monoisotopic mass (kept as the two
    // separate `+= (mono - formulaMono)` then `-= mono` steps the C# performs, so the floating
    // point matches bit-for-bit).
    let formula_mono = formula.monoisotopic_mass();
    for mass in masses.iter_mut() {
        *mass += monoisotopic_mass - formula_mono;
    }

    let highest_abundance = abundances
        .iter()
        .copied()
        .fold(f64::MIN, f64::max);

    let mut peaks: Vec<ExpectedIsotopePeak> = Vec::new();
    for i in 0..masses.len() {
        // expected isotopic mass shift for this peptide
        let mass_shift = masses[i] - monoisotopic_mass;
        // normalized abundance of each isotope
        abundances[i] /= highest_abundance;

        if peaks.len() < num_isotopes_required || abundances[i] > 0.1 {
            peaks.push(ExpectedIsotopePeak {
                mass_shift,
                normalized_abundance: abundances[i],
            });
        }
    }

    peaks
}

/// The mass shift of the most abundant isotope (the peak whose normalized abundance is exactly
/// `1.0`). Mirrors `.First(p => p.Item2 == 1.0).Item1`.
///
/// Returns `None` if no peak has abundance `1.0` (never the case for a distribution produced by
/// [`expected_isotope_peaks`], which always normalizes one peak to `1.0`).
pub fn most_abundant_isotope_shift(peaks: &[ExpectedIsotopePeak]) -> Option<f64> {
    peaks
        .iter()
        .find(|p| p.normalized_abundance == 1.0)
        .map(|p| p.mass_shift)
}

/// The peakfinding mass for an identification: its monoisotopic mass shifted onto the most
/// abundant isotope. Mirrors `identification.PeakfindingMass = MonoisotopicMass +
/// mostAbundantIsotopeShift`.
pub fn peakfinding_mass(monoisotopic_mass: f64, peaks: &[ExpectedIsotopePeak]) -> f64 {
    let shift = most_abundant_isotope_shift(peaks)
        .expect("expected_isotope_peaks always normalizes one peak to abundance 1.0");
    monoisotopic_mass + shift
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Relative-tolerance assert (mirrors the helper used elsewhere in the crate).
    fn assert_rel(actual: f64, expected: f64, tol: f64, what: &str) {
        let rel = (actual - expected).abs() / expected.abs().max(1.0);
        assert!(
            rel < tol,
            "{what}: actual={actual}, expected={expected}, rel={rel:e} >= {tol:e}"
        );
    }

    #[test]
    fn averagine_mass_matches_reference() {
        // Σ averageX * AverageMass(X) — the averagine unit mass is ~111.12 Da (the classic
        // Senko averagine residue). Value computed from the L0 average masses.
        let m = averagine_mass();
        assert_rel(m, 111.1235549, 1e-6, "averagine unit mass");
    }

    #[test]
    fn round_ties_even_in_averagine_counts() {
        // Sanity-check that Rust's round_ties_even matches C# Math.Round midpoint behavior:
        // 2.5 -> 2, 3.5 -> 4, 0.5 -> 0.
        assert_eq!((2.5_f64).round_ties_even() as i32, 2);
        assert_eq!((3.5_f64).round_ties_even() as i32, 4);
        assert_eq!((0.5_f64).round_ties_even() as i32, 0);
    }

    #[test]
    fn valid_sequence_no_modification_uses_bare_formula() {
        // PEPTIDE, with its own exact monoisotopic mass: mass_diff ~0, so no averagine top-up.
        let formula = peptide_base_formula("PEPTIDE").unwrap();
        let mono = formula.monoisotopic_mass();
        let resolved = resolve_formula(None, "PEPTIDE", mono);
        assert_eq!(resolved, formula, "no averagine should be added");
    }

    #[test]
    fn large_mass_gap_adds_averagine() {
        // Pretend PEPTIDE carries a big unknown modification (+500 Da): formula should grow.
        let bare = peptide_base_formula("PEPTIDE").unwrap();
        let bare_mono = bare.monoisotopic_mass();
        let resolved = resolve_formula(None, "PEPTIDE", bare_mono + 500.0);
        assert_ne!(resolved, bare, "averagine top-up should change the formula");
        // The topped-up formula's mono mass should move toward the inflated target.
        assert!(
            resolved.monoisotopic_mass() > bare_mono + 400.0,
            "averagine should account for most of the +500 Da gap"
        );
    }

    #[test]
    fn invalid_sequence_builds_pure_averagine() {
        // A sequence with a non-residue char ('B') falls back to a pure-averagine formula sized
        // by the known mass.
        let resolved = resolve_formula(None, "PEPBTIDE", 1000.0);
        assert!(resolved.atom_count() > 0, "averagine formula should have atoms");
        assert_rel(resolved.monoisotopic_mass(), 1000.0, 0.05, "averagine mass approximates target");
    }

    #[test]
    fn explicit_formula_is_used_verbatim() {
        let explicit = ChemicalFormula::parse_formula("C50H80N15O20").unwrap();
        let resolved = resolve_formula(Some(&explicit), "PEPTIDE", 9999.0);
        assert_eq!(resolved, explicit, "explicit formula overrides sequence/mass");
    }

    #[test]
    fn envelope_is_normalized_and_truncated() {
        let formula = peptide_base_formula("PEPTIDE").unwrap();
        let mono = formula.monoisotopic_mass();
        let peaks = expected_isotope_peaks(None, "PEPTIDE", mono, DEFAULT_NUM_ISOTOPES_REQUIRED);

        // At least NumIsotopesRequired peaks are always kept.
        assert!(peaks.len() >= DEFAULT_NUM_ISOTOPES_REQUIRED);

        // Exactly one peak normalizes to 1.0 (the most abundant).
        let ones = peaks.iter().filter(|p| p.normalized_abundance == 1.0).count();
        assert_eq!(ones, 1, "exactly one peak should be the normalization anchor");

        // Every retained peak is either among the first NumIsotopesRequired or above the 0.1 cut.
        for (i, p) in peaks.iter().enumerate() {
            assert!(
                i < DEFAULT_NUM_ISOTOPES_REQUIRED || p.normalized_abundance > 0.1,
                "peak {i} (abundance {}) violates the truncation rule",
                p.normalized_abundance
            );
        }

        // All abundances in (0, 1].
        for p in &peaks {
            assert!(p.normalized_abundance > 0.0 && p.normalized_abundance <= 1.0);
        }
    }

    #[test]
    fn peakfinding_mass_is_mono_plus_most_abundant_shift() {
        let formula = peptide_base_formula("PEPTIDE").unwrap();
        let mono = formula.monoisotopic_mass();
        let peaks = expected_isotope_peaks(None, "PEPTIDE", mono, DEFAULT_NUM_ISOTOPES_REQUIRED);

        let shift = most_abundant_isotope_shift(&peaks).expect("a 1.0 peak exists");
        let pf = peakfinding_mass(mono, &peaks);
        assert_eq!(pf, mono + shift);

        // For a small peptide the monoisotopic peak is usually the most abundant → shift ~0.
        // PEPTIDE's mono peak dominates, so peakfinding mass ≈ monoisotopic mass.
        assert!((pf - mono).abs() < 2.0, "peakfinding shift should be a small number of isotopes");
    }

    #[test]
    fn mass_shifts_ascend_from_near_zero() {
        let formula = peptide_base_formula("PEPTIDE").unwrap();
        let mono = formula.monoisotopic_mass();
        let peaks = expected_isotope_peaks(None, "PEPTIDE", mono, DEFAULT_NUM_ISOTOPES_REQUIRED);

        // First retained peak (the monoisotopic peak of the bare formula) has shift ~0.
        assert!(peaks[0].mass_shift.abs() < 1e-6, "first shift ~ 0 for an unmodified peptide");
        // Shifts are mass-ascending (distribution order).
        for w in peaks.windows(2) {
            assert!(w[1].mass_shift >= w[0].mass_shift, "shifts should be non-decreasing");
        }
    }
}
