//! Mass tolerance — port of the slice of mzLib `MzLibUtil/PpmTolerance.cs` that the
//! peak-indexing engine needs (`GetMinimumValue`, `GetMaximumValue`, `Within`).
//!
//! Only the parts-per-million flavour is ported here; FlashLFQ's peak finding is always
//! ppm-based. The arithmetic mirrors the C# `PpmTolerance` exactly so bin ranges and
//! `Within` decisions line up bit-for-bit during parity gating.

/// A symmetric parts-per-million tolerance around a mass.
///
/// Mirrors `MzLibUtil.PpmTolerance`: a tolerance of `value` ppm spans
/// `mean * (1 ± value / 1e6)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PpmTolerance {
    /// The tolerance width in parts per million.
    pub value: f64,
}

impl PpmTolerance {
    /// Creates a ppm tolerance with the given width (in ppm).
    pub fn new(value: f64) -> Self {
        PpmTolerance { value }
    }

    /// The smallest mass within tolerance of `mean`: `mean * (1 - value/1e6)`.
    ///
    /// Faithful to C# `PpmTolerance.GetMinimumValue`.
    pub fn get_minimum_value(&self, mean: f64) -> f64 {
        mean * (1.0 - (self.value / 1e6))
    }

    /// The largest mass within tolerance of `mean`: `mean * (1 + value/1e6)`.
    ///
    /// Faithful to C# `PpmTolerance.GetMaximumValue`.
    pub fn get_maximum_value(&self, mean: f64) -> f64 {
        mean * (1.0 + (self.value / 1e6))
    }

    /// Whether `experimental` is within tolerance of `theoretical`.
    ///
    /// Faithful to C# `PpmTolerance.Within`:
    /// `|((experimental - theoretical) / theoretical) * 1e6| <= value`.
    pub fn within(&self, experimental: f64, theoretical: f64) -> bool {
        ((experimental - theoretical) / theoretical * 1e6).abs() <= self.value
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn min_max_bracket_the_mean() {
        let tol = PpmTolerance::new(10.0);
        let mean = 1000.0;
        assert!((tol.get_minimum_value(mean) - 999.99).abs() < 1e-9);
        assert!((tol.get_maximum_value(mean) - 1000.01).abs() < 1e-9);
    }

    #[test]
    fn within_matches_ppm_definition() {
        let tol = PpmTolerance::new(20.0);
        // 500 vs 500.005 → 10 ppm, inside 20 ppm
        assert!(tol.within(500.005, 500.0));
        // 500 vs 500.02 → 40 ppm, outside 20 ppm
        assert!(!tol.within(500.02, 500.0));
        // exactly on the boundary is inclusive (<=)
        let theo = 500.0;
        let exp = theo * (1.0 + 20.0 / 1e6);
        assert!(tol.within(exp, theo));
    }
}
