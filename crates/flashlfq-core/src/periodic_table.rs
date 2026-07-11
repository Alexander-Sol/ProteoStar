//! Periodic-table isotope data (parity layer **L0**).
//!
//! A verbatim port of mzLib's `Chemistry/PeriodicTable.cs` static loader. The decimal
//! source data is embedded from `data/periodic_table.nist.txt` — the *verbatim*
//! `thePeriodicTable` literal extracted from the C# source — and parsed here with the
//! exact same line-by-line logic as the C# `static PeriodicTable()` constructor. We do
//! **not** re-source abundances from NIST: parity is against what mzLib loads, because
//! sub-ppm abundance drift propagates into every theoretical isotope envelope (see
//! `agent_info/FlashLFQ-Rust-Rewrite-Feasibility.md`, "Isotope-distribution port plan").
//!
//! The L0 gate: this loaded table is diffed (exact) against the golden artifact
//! `parity/golden/L0_periodic_table.tsv`, produced by the independent Python replica in
//! `parity/dump_periodic_table.py`. Two independent re-implementations of the loader
//! agreeing exactly — plus `validate_abundances` — is the L0 acceptance criterion.

use std::collections::HashMap;
use std::sync::OnceLock;

/// The verbatim mzLib periodic-table source data (the `thePeriodicTable` C# literal).
const NIST_DATA: &str = include_str!("../data/periodic_table.nist.txt");

/// A single isotope of an element: `(massNumber, atomicMass, relativeAbundance)`.
#[derive(Debug, Clone, PartialEq)]
pub struct Isotope {
    pub mass_number: u16,
    pub atomic_mass: f64,
    pub relative_abundance: f64,
}

/// A chemical element with the isotopes mzLib actually loads (those with a measured
/// abundance), in the order they were added — mirroring `Element.IsotopesInOrderTheyWereAdded`.
#[derive(Debug, Clone, PartialEq)]
pub struct Element {
    pub symbol: String,
    pub atomic_number: u16,
    pub average_mass: f64,
    /// Mass number of the most abundant (principal) isotope; first-added wins ties
    /// (mirrors the strict `>` in `Element.AddIsotope`).
    pub principal_mass_number: u16,
    pub isotopes: Vec<Isotope>,
}

impl Element {
    /// The most abundant (principal) isotope.
    pub fn principal_isotope(&self) -> &Isotope {
        self.isotopes
            .iter()
            .find(|i| i.mass_number == self.principal_mass_number)
            .expect("principal_mass_number always references a loaded isotope")
    }

    /// Look up a specific isotope by its mass number (mirrors the C# `Element this[int]`
    /// indexer). Returns `None` if this element carries no loaded isotope of that mass.
    pub fn isotope_by_mass_number(&self, mass_number: u16) -> Option<&Isotope> {
        self.isotopes.iter().find(|i| i.mass_number == mass_number)
    }
}

/// The loaded periodic table: elements in atomic-number order, indexed by symbol and number.
#[derive(Debug, Clone)]
pub struct PeriodicTable {
    elements: Vec<Element>,
    by_symbol: HashMap<String, usize>,
    by_number: HashMap<u16, usize>,
}

impl PeriodicTable {
    /// All loaded elements, sorted by atomic number.
    pub fn elements(&self) -> &[Element] {
        &self.elements
    }

    /// Look up an element by its atomic symbol (e.g. `"C"`).
    pub fn element_by_symbol(&self, symbol: &str) -> Option<&Element> {
        self.by_symbol.get(symbol).map(|&i| &self.elements[i])
    }

    /// Look up an element by its atomic number (e.g. `6` for carbon).
    pub fn element_by_number(&self, atomic_number: u16) -> Option<&Element> {
        self.by_number.get(&atomic_number).map(|&i| &self.elements[i])
    }

    /// Mirror of `PeriodicTable.ValidateAbundances`: every element's isotope abundances
    /// (summed in added-order) lie within `epsilon` of 1.0. C# asserts this true at 1e-15.
    pub fn validate_abundances(&self, epsilon: f64) -> bool {
        self.elements.iter().all(|e| {
            let mut total = 0.0_f64;
            for iso in &e.isotopes {
                total += iso.relative_abundance;
            }
            (total - 1.0).abs() <= epsilon
        })
    }
}

/// The process-wide loaded periodic table (parsed once).
pub fn periodic_table() -> &'static PeriodicTable {
    static TABLE: OnceLock<PeriodicTable> = OnceLock::new();
    TABLE.get_or_init(|| load(NIST_DATA))
}

/// First maximal run of `[0-9.]` in `s`, mirroring the C# regex `[\d\.]+`.
fn first_number(s: &str) -> Option<&str> {
    let bytes = s.as_bytes();
    let start = bytes
        .iter()
        .position(|&b| b.is_ascii_digit() || b == b'.')?;
    let mut end = start;
    while end < bytes.len() && (bytes[end].is_ascii_digit() || bytes[end] == b'.') {
        end += 1;
    }
    Some(&s[start..end])
}

/// Parse a verbatim periodic-table data block, replicating the C# `static PeriodicTable()`
/// constructor exactly (field accumulation, abundance-gated insertion, added-order isotopes).
fn load(raw: &str) -> PeriodicTable {
    let mut elements: Vec<Element> = Vec::new();
    let mut index_by_number: HashMap<u16, usize> = HashMap::new();

    let mut atomic_number: Option<u16> = None;
    let mut symbol: Option<String> = None;
    let mut isotope_number: Option<u16> = None;
    let mut atomic_mass: Option<f64> = None;
    let mut abundance: Option<f64> = None;
    let mut average_mass: Option<f64> = None;

    for line in raw.split('\n') {
        if line.trim().is_empty() {
            continue;
        }
        let eq = line
            .find('=')
            .expect("every non-blank periodic-table line has '='");
        let key = line[..eq].trim();
        let value = line[eq + 1..].trim();

        match key {
            "Atomic Number" => atomic_number = Some(value.parse().expect("atomic number")),
            "Atomic Symbol" => symbol = Some(value.to_string()),
            "Mass Number" => isotope_number = Some(value.parse().expect("mass number")),
            "Relative Atomic Mass" => {
                atomic_mass = Some(
                    first_number(line)
                        .expect("relative atomic mass has a number")
                        .parse()
                        .expect("atomic mass parses"),
                );
            }
            "Isotopic Composition" => {
                if let Some(num) = first_number(line) {
                    abundance = Some(num.parse().expect("abundance parses"));
                }
            }
            "Standard Atomic Weight" => {
                if line.contains('[') {
                    let after_bracket = &line[line.find('[').unwrap() + 1..];
                    let avg1: f64 = first_number(after_bracket)
                        .expect("bracketed weight has a number")
                        .parse()
                        .expect("avg1 parses");
                    let avg2 = line
                        .find(',')
                        .and_then(|c| first_number(&line[c + 1..]))
                        .and_then(|n| n.parse::<f64>().ok());
                    average_mass = Some(match avg2 {
                        Some(a2) => (avg1 + a2) / 2.0,
                        None => avg1,
                    });
                } else if let Some(num) = first_number(line) {
                    average_mass = Some(num.parse().expect("standard atomic weight parses"));
                }
            }
            "Notes" => {
                match atomic_number.and_then(|n| index_by_number.get(&n).copied()) {
                    None => {
                        // C#: create the element only when every field (abundance incl.) is present.
                        if let (Some(sym), Some(num), Some(avg), Some(iso), Some(am), Some(ab)) = (
                            symbol.as_ref(),
                            atomic_number,
                            average_mass,
                            isotope_number,
                            atomic_mass,
                            abundance,
                        ) {
                            let mut element = Element {
                                symbol: sym.clone(),
                                atomic_number: num,
                                average_mass: avg,
                                principal_mass_number: iso,
                                isotopes: Vec::new(),
                            };
                            add_isotope(&mut element, iso, am, ab);
                            index_by_number.insert(num, elements.len());
                            elements.push(element);
                        }
                    }
                    Some(idx) => {
                        // Element exists: append the isotope only if this block has an abundance.
                        if let (Some(iso), Some(am), Some(ab)) =
                            (isotope_number, atomic_mass, abundance)
                        {
                            add_isotope(&mut elements[idx], iso, am, ab);
                        }
                    }
                }

                atomic_number = None;
                symbol = None;
                isotope_number = None;
                atomic_mass = None;
                abundance = None;
                average_mass = None;
            }
            other => panic!("Could not parse line from periodic table: key={other:?}"),
        }
    }

    elements.sort_by_key(|e| e.atomic_number);
    let by_symbol = elements
        .iter()
        .enumerate()
        .map(|(i, e)| (e.symbol.clone(), i))
        .collect();
    let by_number = elements
        .iter()
        .enumerate()
        .map(|(i, e)| (e.atomic_number, i))
        .collect();
    PeriodicTable {
        elements,
        by_symbol,
        by_number,
    }
}

/// Mirror of `Element.AddIsotope`: append in added-order, track strict-greater principal.
fn add_isotope(element: &mut Element, mass_number: u16, atomic_mass: f64, abundance: f64) {
    let new_principal = element
        .isotopes
        .iter()
        .find(|i| i.mass_number == element.principal_mass_number)
        .map_or(true, |p| abundance > p.relative_abundance);
    if new_principal {
        element.principal_mass_number = mass_number;
    }
    element.isotopes.push(Isotope {
        mass_number,
        atomic_mass,
        relative_abundance: abundance,
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The golden L0 artifact emitted by `parity/dump_periodic_table.py`.
    const GOLDEN: &str = include_str!("../parity/golden/L0_periodic_table.tsv");

    /// Parsed golden: `(symbol, atomic_number, average_mass, principal, isotopes)`.
    struct GoldenElement {
        symbol: String,
        atomic_number: u16,
        average_mass: f64,
        principal_mass_number: u16,
        isotopes: Vec<Isotope>,
    }

    fn parse_golden() -> Vec<GoldenElement> {
        let mut out: Vec<GoldenElement> = Vec::new();
        for line in GOLDEN.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            match f[0] {
                "E" => out.push(GoldenElement {
                    symbol: f[1].to_string(),
                    atomic_number: f[2].parse().unwrap(),
                    average_mass: f[3].parse().unwrap(),
                    principal_mass_number: f[4].parse().unwrap(),
                    isotopes: Vec::new(),
                }),
                "I" => out.last_mut().unwrap().isotopes.push(Isotope {
                    mass_number: f[1].parse().unwrap(),
                    atomic_mass: f[2].parse().unwrap(),
                    relative_abundance: f[3].parse().unwrap(),
                }),
                other => panic!("unexpected golden record tag {other:?}"),
            }
        }
        out
    }

    #[test]
    fn loads_expected_element_and_isotope_counts() {
        let pt = periodic_table();
        // Parity numbers from the Python replica: mzLib loads the 84 naturally-occurring
        // elements (those with a measured isotopic abundance) and 288 isotopes total.
        assert_eq!(pt.elements().len(), 84, "loaded element count");
        let isotopes: usize = pt.elements().iter().map(|e| e.isotopes.len()).sum();
        assert_eq!(isotopes, 288, "loaded isotope count");
    }

    #[test]
    fn hydrogen_carbon_known_values() {
        let pt = periodic_table();

        // Hydrogen carries protium + deuterium (D is folded into element 1 by the loader);
        // tritium has no abundance and is dropped.
        let h = pt.element_by_symbol("H").expect("H loaded");
        assert_eq!(h.atomic_number, 1);
        assert_eq!(h.isotopes.len(), 2);
        assert_eq!(h.principal_mass_number, 1);
        assert_eq!(h.isotopes[0].mass_number, 1);
        assert_eq!(h.isotopes[0].atomic_mass, 1.00782503223);
        assert_eq!(h.isotopes[0].relative_abundance, 0.999885);
        assert_eq!(h.isotopes[1].mass_number, 2); // deuterium

        // Carbon-12 is exactly 12 by definition and is the principal isotope.
        let c = pt.element_by_number(6).expect("C loaded");
        assert_eq!(c.symbol, "C");
        assert_eq!(c.principal_mass_number, 12);
        assert_eq!(c.principal_isotope().atomic_mass, 12.0);
    }

    #[test]
    fn validate_abundances_passes_at_1e15() {
        // Mirrors the C# test ChemicalFormulaPropertiesTests.ValidatePeriodicTable:
        // ValidateAbundances(1e-15) is true; at epsilon 0 it is false.
        let pt = periodic_table();
        assert!(pt.validate_abundances(1e-15), "abundances must sum to 1 within 1e-15");
        assert!(!pt.validate_abundances(0.0), "not all sums are bit-exactly 1.0");
    }

    /// L0 parity gate: the embedded/loaded table must match the golden dump **exactly**.
    /// Integer fields compare exact; floats compare with `==` (bit-identical, because both
    /// sides parse the same decimal tokens via correctly-rounded f64 parsing).
    #[test]
    fn l0_diff_vs_golden_is_exact() {
        let pt = periodic_table();
        let golden = parse_golden();
        assert_eq!(pt.elements().len(), golden.len(), "element count vs golden");

        for (loaded, gold) in pt.elements().iter().zip(golden.iter()) {
            assert_eq!(loaded.symbol, gold.symbol, "symbol");
            assert_eq!(loaded.atomic_number, gold.atomic_number, "atomic number ({})", gold.symbol);
            assert_eq!(
                loaded.average_mass, gold.average_mass,
                "average mass ({})",
                gold.symbol
            );
            assert_eq!(
                loaded.principal_mass_number, gold.principal_mass_number,
                "principal isotope ({})",
                gold.symbol
            );
            assert_eq!(
                loaded.isotopes.len(),
                gold.isotopes.len(),
                "isotope count ({})",
                gold.symbol
            );
            for (li, gi) in loaded.isotopes.iter().zip(gold.isotopes.iter()) {
                assert_eq!(li, gi, "isotope {} of {}", gi.mass_number, gold.symbol);
            }
        }
    }
}
