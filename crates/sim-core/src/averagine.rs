//! mzLib's `Averagine` static tables, embedded verbatim (contract §4).
//!
//! The envelope shapes must come from the identical table C# used, so rather than regenerate it
//! from the averagine composition (as `flashlfq_core::deconvolution` does for its own purposes)
//! the table exported by the Phase-0 golden generator is embedded and parsed once.

use std::sync::OnceLock;

use serde::Deserialize;

static AVERAGINE_JSON: &str = include_str!("../data/averagine.v1.json");

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AveragineTable {
    /// Per row, isotopologue neutral masses, most-intense-first.
    pub all_masses: Vec<Vec<f64>>,
    /// Per row, isotopologue intensities, most-intense-first.
    pub all_intensities: Vec<Vec<f64>>,
    /// Ascending; the lookup key for a target mass.
    pub most_intense_masses: Vec<f64>,
    /// Per row, most-intense mass − monoisotopic mass.
    pub diff_to_monoisotopic: Vec<f64>,
}

impl AveragineTable {
    pub fn get() -> &'static AveragineTable {
        static TABLE: OnceLock<AveragineTable> = OnceLock::new();
        TABLE.get_or_init(|| {
            serde_json::from_str(AVERAGINE_JSON).expect("embedded averagine.v1.json is valid")
        })
    }

    /// C# `Averagine.GetMostIntenseMassIndex` → `MostIntenseMasses.GetClosestIndex(mass)`,
    /// including its tie-break: equidistant neighbours resolve to the upper one.
    pub fn most_intense_mass_index(&self, test_mass: f64) -> usize {
        let a = &self.most_intense_masses;
        // Array.BinarySearch: exact match wins; otherwise the insertion point.
        let i = match a.binary_search_by(|m| m.total_cmp(&test_mass)) {
            Ok(i) => return i,
            Err(i) => i,
        };
        if i == a.len() {
            return i - 1;
        }
        if i == 0 {
            return 0;
        }
        let backwards = test_mass - a[i - 1];
        let forwards = a[i] - test_mass;
        if backwards < forwards { i - 1 } else { i }
    }
}
