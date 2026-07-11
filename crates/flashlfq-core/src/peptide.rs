//! Peptide → chemical formula (parity layer **L1**, part 2).
//!
//! A port of the formula side of mzLib's amino-acid polymer code: the residue→formula table
//! from `Proteomics/AminoAcidPolymer/Residue.cs` and the formula assembly in
//! `AminoAcidPolymer.GetChemicalFormula()`.
//!
//! `Peptide.GetChemicalFormula()` (called as
//! `new Peptide(baseSequence).GetChemicalFormula()`) sums, for an *unmodified* peptide:
//! - the N-terminus formula (`H`),
//! - the C-terminus formula (`OH`),
//! - each residue's formula.
//!
//! The two termini together contribute one water (`H2O`) — the residues are the
//! water-loss (condensed) forms, so a peptide is `Σ residues + H2O`. This matches the
//! free amino acids exactly: glycine residue `C2H3NO` + `H2O` = `C2H5NO2` (free glycine).
//!
//! Modifications are layered on top by mzLib's `PeptideWithSetModifications.FullChemicalFormula`,
//! which simply `Add`s each modification's `ChemicalFormula`. Here that is
//! [`peptide_formula_with_mods`]: build the backbone, then add modification formulas. The
//! *source* of a modification's formula (mod list / GPTMD lookup, modified-sequence parsing)
//! is deferred to the L1 parity harness (P1.4); this layer only owns the residue table and
//! the assembly arithmetic.

use crate::chemical_formula::{ChemicalFormula, FormulaParseError};

/// Error building a peptide's chemical formula.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PeptideError {
    /// A sequence character is not a known amino-acid residue (no formula in the table).
    UnknownResidue(char),
}

impl std::fmt::Display for PeptideError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PeptideError::UnknownResidue(c) => {
                write!(f, "Sequence character '{c}' is not a known amino-acid residue")
            }
        }
    }
}

impl std::error::Error for PeptideError {}

/// Residue formula string for a one-letter amino-acid code, or `None` if unknown.
///
/// Verbatim from `Residue.cs`'s `ResiduesDictionary` (the formula passed to each
/// `ChemicalFormula.ParseFormula`). These are the condensed (water-loss) residue forms.
/// Includes the two non-standard residues mzLib carries: `O` pyrrolysine and `U`
/// selenocysteine. The letters mzLib leaves `null` in `ResiduesByLetter` (B, J, X, Z) have
/// no residue formula and return `None`.
pub fn residue_formula(letter: char) -> Option<&'static str> {
    Some(match letter {
        'A' => "C3H5NO",     // Alanine
        'R' => "C6H12N4O",   // Arginine
        'N' => "C4H6N2O2",   // Asparagine
        'D' => "C4H5NO3",    // Aspartic Acid
        'C' => "C3H5NOS",    // Cysteine
        'E' => "C5H7NO3",    // Glutamic Acid
        'Q' => "C5H8N2O2",   // Glutamine
        'G' => "C2H3NO",     // Glycine
        'H' => "C6H7N3O",    // Histidine
        'I' => "C6H11NO",    // Isoleucine
        'L' => "C6H11NO",    // Leucine
        'K' => "C6H12N2O",   // Lysine
        'M' => "C5H9NOS",    // Methionine
        'F' => "C9H9NO",     // Phenylalanine
        'P' => "C5H7NO",     // Proline
        'O' => "C12H19N3O2", // Pyrrolysine
        'U' => "C3H5NOSe",   // Selenocysteine
        'S' => "C3H5NO2",    // Serine
        'T' => "C4H7NO2",    // Threonine
        'W' => "C11H10N2O",  // Tryptophan
        'Y' => "C9H9NO2",    // Tyrosine
        'V' => "C5H9NO",     // Valine
        _ => return None,
    })
}

/// Chemical formula of an unmodified peptide from its base (one-letter) sequence.
///
/// Mirrors `new Peptide(baseSequence).GetChemicalFormula()`: N-terminus `H` + C-terminus
/// `OH` + the formula of every residue. Returns [`PeptideError::UnknownResidue`] for any
/// character not in the residue table.
pub fn peptide_base_formula(sequence: &str) -> Result<ChemicalFormula, PeptideError> {
    let mut formula = ChemicalFormula::new();

    // N-terminus (H) and C-terminus (OH) — together one water. These literals parse
    // cleanly, so the parse cannot fail.
    formula.add_formula(&parse_known("H"));
    formula.add_formula(&parse_known("OH"));

    for letter in sequence.chars() {
        let residue = residue_formula(letter).ok_or(PeptideError::UnknownResidue(letter))?;
        formula.add_formula(&parse_known(residue));
    }

    Ok(formula)
}

/// Peptide formula with modifications added (mirrors `FullChemicalFormula`).
///
/// Builds the backbone via [`peptide_base_formula`], then `Add`s each modification formula.
/// The modification formulas are taken as already-resolved [`ChemicalFormula`]s — resolving
/// a modified-sequence string to its mod formulas is the L1 harness's job (P1.4).
pub fn peptide_formula_with_mods(
    sequence: &str,
    modifications: &[ChemicalFormula],
) -> Result<ChemicalFormula, PeptideError> {
    let mut formula = peptide_base_formula(sequence)?;
    for modification in modifications {
        formula.add_formula(modification);
    }
    Ok(formula)
}

/// Parse a formula literal that is statically known to be well-formed (residue/terminus
/// strings from the embedded table). A parse failure here is a bug in the table, not input.
fn parse_known(formula: &str) -> ChemicalFormula {
    ChemicalFormula::parse_formula(formula).unwrap_or_else(|e: FormulaParseError| {
        panic!("residue/terminus formula {formula:?} failed to parse: {e}")
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::periodic_table::periodic_table;

    /// Relative tolerance helper (the L1 gate uses rel-1e-9; literature masses are quoted a
    /// touch coarser, so we assert against rel-1e-7 for the known-mass checks).
    fn assert_rel(actual: f64, expected: f64, tol: f64, what: &str) {
        let rel = (actual - expected).abs() / expected.abs().max(1.0);
        assert!(
            rel < tol,
            "{what}: actual={actual}, expected={expected}, rel={rel:e} >= {tol:e}"
        );
    }

    /// Plain-element count for a symbol, via the periodic table → atomic number.
    fn count(f: &ChemicalFormula, symbol: &str) -> i32 {
        let an = periodic_table()
            .element_by_symbol(symbol)
            .unwrap_or_else(|| panic!("unknown symbol {symbol}"))
            .atomic_number;
        f.count_of_element(an)
    }

    #[test]
    fn single_residues_equal_free_amino_acids() {
        // A length-1 "peptide" is the free amino acid: residue + H2O.
        // Glycine C2H5NO2, mono 75.03203.
        let gly = peptide_base_formula("G").unwrap();
        assert_eq!(count(&gly, "C"), 2);
        assert_eq!(count(&gly, "H"), 5);
        assert_eq!(count(&gly, "N"), 1);
        assert_eq!(count(&gly, "O"), 2);
        assert_rel(gly.monoisotopic_mass(), 75.03202840, 1e-7, "Gly mono");

        // Alanine C3H7NO2, mono 89.04768.
        let ala = peptide_base_formula("A").unwrap();
        assert_eq!(count(&ala, "C"), 3);
        assert_eq!(count(&ala, "H"), 7);
        assert_eq!(count(&ala, "N"), 1);
        assert_eq!(count(&ala, "O"), 2);
        assert_rel(ala.monoisotopic_mass(), 89.04767847, 1e-7, "Ala mono");
    }

    #[test]
    fn peptide_element_counts_and_mass() {
        // The canonical example: "PEPTIDE" = C34H53N7O15, monoisotopic 799.35996 Da.
        // Hand sum of residue formulas + H2O:
        //   P C5H7NO  E C5H7NO3  P C5H7NO  T C4H7NO2  I C6H11NO  D C4H5NO3  E C5H7NO3
        //   ΣC=34 ΣH=51 ΣN=7 ΣO=14, + H2O → H53 O15.
        let f = peptide_base_formula("PEPTIDE").unwrap();
        assert_eq!(count(&f, "C"), 34, "PEPTIDE carbons");
        assert_eq!(count(&f, "H"), 53, "PEPTIDE hydrogens");
        assert_eq!(count(&f, "N"), 7, "PEPTIDE nitrogens");
        assert_eq!(count(&f, "O"), 15, "PEPTIDE oxygens");
        assert_eq!(f.atom_count(), 34 + 53 + 7 + 15);
        assert_rel(f.monoisotopic_mass(), 799.35996, 1e-7, "PEPTIDE mono");
    }

    #[test]
    fn cysteine_and_methionine_carry_sulfur() {
        // "MCK": M C5H9NOS, C C3H5NOS, K C6H12N2O → ΣC14 H26 N4 O3 S2, +H2O → H28 O4.
        let f = peptide_base_formula("MCK").unwrap();
        assert_eq!(count(&f, "C"), 14);
        assert_eq!(count(&f, "H"), 28);
        assert_eq!(count(&f, "N"), 4);
        assert_eq!(count(&f, "O"), 4);
        assert_eq!(count(&f, "S"), 2);
    }

    #[test]
    fn selenocysteine_residue_carries_selenium() {
        // Free selenocysteine: C3H5NOSe + H2O = C3H7NO2Se.
        let f = peptide_base_formula("U").unwrap();
        assert_eq!(count(&f, "C"), 3);
        assert_eq!(count(&f, "H"), 7);
        assert_eq!(count(&f, "N"), 1);
        assert_eq!(count(&f, "O"), 2);
        assert_eq!(count(&f, "Se"), 1);
    }

    #[test]
    fn modifications_are_added_on_top() {
        // PEPTIDE with a phosphorylation (HPO3, mod mass +79.96633): C34H54N7O18P.
        let phospho = ChemicalFormula::parse_formula("HO3P").unwrap();
        let base = peptide_base_formula("PEPTIDE").unwrap();
        let modded = peptide_formula_with_mods("PEPTIDE", &[phospho]).unwrap();

        assert_eq!(count(&modded, "C"), 34, "phospho keeps carbons");
        assert_eq!(count(&modded, "H"), 54, "phospho adds 1 H");
        assert_eq!(count(&modded, "O"), 18, "phospho adds 3 O");
        assert_eq!(count(&modded, "P"), 1, "phospho adds 1 P");

        let delta = modded.monoisotopic_mass() - base.monoisotopic_mass();
        assert_rel(delta, 79.96633, 1e-6, "phospho mass shift");
    }

    #[test]
    fn unknown_residue_is_rejected() {
        // 'B' (Asx) has no residue formula in mzLib's table.
        assert_eq!(
            peptide_base_formula("PEPBTIDE"),
            Err(PeptideError::UnknownResidue('B'))
        );
        // 'J' is a mass-only special case in C# (Ile mass) but null in ResiduesByLetter,
        // so it has no formula either.
        assert_eq!(
            peptide_base_formula("J"),
            Err(PeptideError::UnknownResidue('J'))
        );
    }
}
