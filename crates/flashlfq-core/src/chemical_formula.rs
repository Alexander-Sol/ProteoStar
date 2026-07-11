//! Chemical formulas (parity layer **L1**, part 1).
//!
//! A port of mzLib's `Chemistry/ChemicalFormula.cs`. A `ChemicalFormula` is a multiset of
//! chemical elements and *specific* isotopes (e.g. carbon-13), each with an integer count.
//! It carries no physical identity — peptides, modifications, etc. *have* a formula.
//!
//! Two stores, exactly as in C#:
//! - `elements`: plain elements with an unspecified isotope (their mass uses the element's
//!   average / principal-isotope mass), keyed by atomic number.
//! - `isotopes`: isotope-specified entries (mass uses that isotope's exact atomic mass),
//!   keyed by `(atomic_number, mass_number)`.
//!
//! The two masses that matter downstream:
//! - [`ChemicalFormula::monoisotopic_mass`] — elements contribute their *principal* isotope
//!   mass; isotope entries contribute their exact atomic mass.
//! - [`ChemicalFormula::average_mass`] — elements contribute their abundance-weighted average
//!   mass; isotope entries (already mass-specific) contribute their exact atomic mass.
//!
//! Element/isotope data comes from the L0 [`crate::periodic_table`]. Masses are summed in a
//! deterministic (sorted) order; the L1 parity tolerance is relative `1e-9`, so summation
//! order vs the C# `Dictionary` enumeration order is immaterial.

use std::collections::BTreeMap;

use crate::periodic_table::periodic_table;

/// A chemical formula: a multiset of elements and isotope-specified atoms with integer counts.
///
/// Counts may transiently go negative while building (subtracting more than present); an entry
/// is dropped when it reaches the C# removal threshold (`== 0` for elements, `<= 0` for
/// isotopes), so the stores never retain a zero count.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ChemicalFormula {
    /// Plain elements: `atomic_number -> count`.
    elements: BTreeMap<u16, i32>,
    /// Isotope-specified atoms: `(atomic_number, mass_number) -> count`.
    isotopes: BTreeMap<(u16, u16), i32>,
}

/// Error parsing a chemical-formula string (mirrors the C# `MzLibException` cases).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FormulaParseError {
    /// The string is not in the chemical-formula grammar (C# `ValidateFormulaRegex` fails).
    BadFormat(String),
    /// A chemical symbol was not found in the periodic table.
    UnknownElement(String),
    /// An isotope `{n}` was requested for an element that has no loaded isotope of that mass.
    UnknownIsotope { symbol: String, mass_number: u16 },
}

impl std::fmt::Display for FormulaParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FormulaParseError::BadFormat(s) => write!(
                f,
                "Input string for chemical formula was in an incorrect format: {s}"
            ),
            FormulaParseError::UnknownElement(s) => {
                write!(f, "Element with atomic symbol {s} not found")
            }
            FormulaParseError::UnknownIsotope {
                symbol,
                mass_number,
            } => write!(f, "Element {symbol} has no isotope with mass number {mass_number}"),
        }
    }
}

impl std::error::Error for FormulaParseError {}

impl ChemicalFormula {
    /// An empty formula.
    pub fn new() -> Self {
        ChemicalFormula::default()
    }

    /// Add `count` plain atoms of the element with the given atomic number.
    ///
    /// Mirrors `ChemicalFormula.Add(Element, int)`: a zero `count` is a no-op, and an entry
    /// whose running total reaches exactly `0` is removed.
    pub fn add_element(&mut self, atomic_number: u16, count: i32) {
        if count == 0 {
            return;
        }
        let entry = self.elements.entry(atomic_number).or_insert(0);
        *entry += count;
        if *entry == 0 {
            self.elements.remove(&atomic_number);
        }
    }

    /// Add `count` atoms of a specific isotope `(atomic_number, mass_number)`.
    ///
    /// Mirrors `ChemicalFormula.Add(Isotope, int)`: a zero `count` is a no-op, and an entry
    /// whose running total reaches `<= 0` is removed.
    pub fn add_isotope(&mut self, atomic_number: u16, mass_number: u16, count: i32) {
        if count == 0 {
            return;
        }
        let key = (atomic_number, mass_number);
        let entry = self.isotopes.entry(key).or_insert(0);
        *entry += count;
        if *entry <= 0 {
            self.isotopes.remove(&key);
        }
    }

    /// Merge another formula into this one (mirrors `ChemicalFormula.Add(ChemicalFormula)`).
    pub fn add_formula(&mut self, other: &ChemicalFormula) {
        for (&an, &count) in &other.elements {
            self.add_element(an, count);
        }
        for (&(an, mn), &count) in &other.isotopes {
            self.add_isotope(an, mn, count);
        }
    }

    /// Total number of atoms (elements + isotopes), counts summed (`ChemicalFormula.AtomCount`).
    pub fn atom_count(&self) -> i32 {
        self.elements.values().sum::<i32>() + self.isotopes.values().sum::<i32>()
    }

    /// Count of plain (non-isotope-specified) atoms of the element with this atomic number.
    /// Returns `0` if the element is absent. Used by the parity harness to assert *exact*
    /// element counts (mirrors reading a count out of the C# formula's element dictionary).
    pub fn count_of_element(&self, atomic_number: u16) -> i32 {
        self.elements.get(&atomic_number).copied().unwrap_or(0)
    }

    /// Count of a specific isotope `(atomic_number, mass_number)`. Returns `0` if absent.
    pub fn count_of_isotope(&self, atomic_number: u16, mass_number: u16) -> i32 {
        self.isotopes
            .get(&(atomic_number, mass_number))
            .copied()
            .unwrap_or(0)
    }

    /// Iterate plain (non-isotope-specified) elements as `(atomic_number, count)`, in
    /// atomic-number order. Mirrors enumerating the C# `ChemicalFormula.Elements` dictionary
    /// (order differs from C#'s `Dictionary`, but downstream consumers — e.g. the isotopic
    /// distribution — are order-tolerant within the parity tolerance).
    pub fn elements(&self) -> impl Iterator<Item = (u16, i32)> + '_ {
        self.elements.iter().map(|(&an, &count)| (an, count))
    }

    /// Iterate isotope-specified atoms as `((atomic_number, mass_number), count)`.
    /// Mirrors enumerating the C# `ChemicalFormula.Isotopes` dictionary.
    pub fn isotopes(&self) -> impl Iterator<Item = ((u16, u16), i32)> + '_ {
        self.isotopes.iter().map(|(&key, &count)| (key, count))
    }

    /// The monoisotopic mass: elements use their principal-isotope mass, isotope entries use
    /// their exact atomic mass (`ChemicalFormula.MonoisotopicMass`).
    pub fn monoisotopic_mass(&self) -> f64 {
        let pt = periodic_table();
        let mut sum = 0.0_f64;
        for (&(an, mn), &count) in &self.isotopes {
            let element = pt
                .element_by_number(an)
                .expect("isotope entry references a loaded element");
            let isotope = element
                .isotope_by_mass_number(mn)
                .expect("isotope entry references a loaded isotope");
            sum += isotope.atomic_mass * count as f64;
        }
        for (&an, &count) in &self.elements {
            let element = pt
                .element_by_number(an)
                .expect("element entry references a loaded element");
            sum += element.principal_isotope().atomic_mass * count as f64;
        }
        sum
    }

    /// The average mass: elements use their abundance-weighted average mass, isotope entries
    /// use their exact atomic mass (`ChemicalFormula.AverageMass`).
    pub fn average_mass(&self) -> f64 {
        let pt = periodic_table();
        let mut sum = 0.0_f64;
        for (&(an, mn), &count) in &self.isotopes {
            let element = pt
                .element_by_number(an)
                .expect("isotope entry references a loaded element");
            let isotope = element
                .isotope_by_mass_number(mn)
                .expect("isotope entry references a loaded isotope");
            sum += isotope.atomic_mass * count as f64;
        }
        for (&an, &count) in &self.elements {
            let element = pt
                .element_by_number(an)
                .expect("element entry references a loaded element");
            sum += element.average_mass * count as f64;
        }
        sum
    }

    /// Parse a formula string such as `"C6H12O6"` or `"C{13}6H12O6"`.
    ///
    /// Faithful to the C# `ParseFormula` grammar
    /// (`\s*([A-Z][a-z]*)(?:\{([0-9]+)\})?(-)?([0-9]+)?\s*`, repeated to cover the whole
    /// string): an uppercase letter then optional lowercase letters give the symbol; an
    /// optional `{n}` selects a specific isotope; an optional leading `-` subtracts; an
    /// optional integer gives the count (default `1`). `"C-2C{13}5"` subtracts 2 carbons
    /// then adds 5 carbon-13.
    pub fn parse_formula(formula: &str) -> Result<ChemicalFormula, FormulaParseError> {
        let pt = periodic_table();
        let mut f = ChemicalFormula::new();
        let bytes = formula.as_bytes();
        let mut i = 0;
        let bad = || FormulaParseError::BadFormat(formula.to_string());

        while i < bytes.len() {
            // \s*
            while i < bytes.len() && bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            if i >= bytes.len() {
                break;
            }

            // ([A-Z][a-z]*)  — required chemical symbol
            if !bytes[i].is_ascii_uppercase() {
                return Err(bad());
            }
            let sym_start = i;
            i += 1;
            while i < bytes.len() && bytes[i].is_ascii_lowercase() {
                i += 1;
            }
            let symbol = &formula[sym_start..i];

            // (?:\{([0-9]+)\})?  — optional isotope mass number
            let mut mass_number: Option<u16> = None;
            if i < bytes.len() && bytes[i] == b'{' {
                i += 1;
                let num_start = i;
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                if i == num_start || i >= bytes.len() || bytes[i] != b'}' {
                    return Err(bad());
                }
                mass_number = Some(formula[num_start..i].parse().map_err(|_| bad())?);
                i += 1; // consume '}'
            }

            // (-)?  — optional sign
            let sign: i32 = if i < bytes.len() && bytes[i] == b'-' {
                i += 1;
                -1
            } else {
                1
            };

            // ([0-9]+)?  — optional count (defaults to 1)
            let num_start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            let count: i32 = if i == num_start {
                1
            } else {
                formula[num_start..i].parse().map_err(|_| bad())?
            };

            let element = pt
                .element_by_symbol(symbol)
                .ok_or_else(|| FormulaParseError::UnknownElement(symbol.to_string()))?;

            match mass_number {
                Some(mn) => {
                    if element.isotope_by_mass_number(mn).is_none() {
                        return Err(FormulaParseError::UnknownIsotope {
                            symbol: symbol.to_string(),
                            mass_number: mn,
                        });
                    }
                    f.add_isotope(element.atomic_number, mn, sign * count);
                }
                None => f.add_element(element.atomic_number, sign * count),
            }
        }

        Ok(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Relative tolerance for mass parity (the L1 gate uses rel-1e-9; we assert against
    /// literature values a touch looser to absorb reference rounding).
    fn assert_rel(actual: f64, expected: f64, tol: f64, what: &str) {
        let rel = (actual - expected).abs() / expected.abs().max(1.0);
        assert!(
            rel < tol,
            "{what}: actual={actual}, expected={expected}, rel={rel:e} >= {tol:e}"
        );
    }

    #[test]
    fn water_masses() {
        let water = ChemicalFormula::parse_formula("H2O").unwrap();
        assert_eq!(water.atom_count(), 3);
        // Monoisotopic H2O = 2*1.00782503223 + 15.99491461957 = 18.0105646840
        assert_rel(water.monoisotopic_mass(), 18.0105646863, 1e-9, "H2O mono");
        // Average ~18.0153 g/mol.
        assert_rel(water.average_mass(), 18.01528, 1e-5, "H2O average");
    }

    #[test]
    fn glucose_masses() {
        let glucose = ChemicalFormula::parse_formula("C6H12O6").unwrap();
        assert_eq!(glucose.atom_count(), 24);
        // Monoisotopic C6H12O6 = 6*12 + 12*1.00782503223 + 6*15.99491461957 = 180.0633881
        assert_rel(glucose.monoisotopic_mass(), 180.0633881, 1e-9, "glucose mono");
        // Average ~180.156 g/mol.
        assert_rel(glucose.average_mass(), 180.15588, 1e-5, "glucose average");
    }

    #[test]
    fn isotope_specified_carbon13() {
        // Fully 13C-labelled glucose: each carbon is carbon-13.
        let heavy = ChemicalFormula::parse_formula("C{13}6H12O6").unwrap();
        assert_eq!(heavy.atom_count(), 24);
        // 6*13.00335483507 + 12*1.00782503223 + 6*15.99491461957 = 186.0835171
        assert_rel(heavy.monoisotopic_mass(), 186.0835171, 1e-9, "13C6 glucose mono");

        // The light/heavy mono-mass delta is 6 * (13C - 12C) neutron-ish shift.
        let light = ChemicalFormula::parse_formula("C6H12O6").unwrap();
        let delta = heavy.monoisotopic_mass() - light.monoisotopic_mass();
        assert_rel(delta, 6.0 * (13.00335483507 - 12.0), 1e-9, "6x 13C shift");
    }

    #[test]
    fn parsing_handles_signs_and_multiletter_symbols() {
        // Sodium (two-letter symbol) and a subtraction that cancels carbons to zero.
        let f = ChemicalFormula::parse_formula("C2C-2Na2").unwrap();
        // Both carbon entries cancel → element removed; only Na2 remains.
        assert_eq!(f.atom_count(), 2);
        let na = periodic_table().element_by_symbol("Na").unwrap();
        assert_rel(
            f.monoisotopic_mass(),
            2.0 * na.principal_isotope().atomic_mass,
            1e-12,
            "Na2 mono",
        );
    }

    #[test]
    fn add_formula_combines_counts() {
        let mut total = ChemicalFormula::parse_formula("H2O").unwrap();
        total.add_formula(&ChemicalFormula::parse_formula("H2O").unwrap());
        // 2 waters = H4O2.
        assert_eq!(total.atom_count(), 6);
        assert_rel(
            total.monoisotopic_mass(),
            2.0 * 18.0105646863,
            1e-9,
            "2x water mono",
        );
    }

    #[test]
    fn rejects_malformed_and_unknown() {
        assert!(matches!(
            ChemicalFormula::parse_formula("123"),
            Err(FormulaParseError::BadFormat(_))
        ));
        assert!(matches!(
            ChemicalFormula::parse_formula("Xx2"),
            Err(FormulaParseError::UnknownElement(_))
        ));
        assert!(matches!(
            ChemicalFormula::parse_formula("C{99}2"),
            Err(FormulaParseError::UnknownIsotope { .. })
        ));
    }
}
