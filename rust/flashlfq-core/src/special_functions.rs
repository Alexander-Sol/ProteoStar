//! Special functions + probability distributions used by the MBR scorer (PLAN.md P3.2c).
//!
//! FlashLFQ's `MbrScorer` scores match-between-runs transfers with three MathNet objects:
//! `MathNet.Numerics.Distributions.Normal` (intensity / ppm / scan-count / RT distributions) and
//! `MathNet.Numerics.Distributions.Gamma` (the isotopic-correlation distribution). Their CDFs reduce
//! to two special functions — the complementary error function `erfc` (Normal) and the regularized
//! lower incomplete gamma `P(a, x)` (Gamma). This module ports both, plus the log-gamma they need,
//! and wraps them in [`Normal`] / [`Gamma`] value types that mirror the MathNet API the scorer uses.
//!
//! ## Faithfulness
//! - [`gamma_ln`] is MathNet's Lanczos approximation (`SpecialFunctions.GammaLn`), constants verbatim.
//! - [`gamma_lower_regularized`] / [`gamma_upper_regularized`] are MathNet's Cephes-derived
//!   `GammaLowerRegularized` / `GammaUpperRegularized` (series for `x ≤ 1 || x ≤ a`, Lentz-style
//!   continued fraction otherwise), ported line-for-line including the `big`/`bigInv` rescaling.
//! - `erf`/`erfc` are expressed through the regularized incomplete gamma via the exact identities
//!   `erf(x) = P(½, x²)` and `erfc(x) = Q(½, x²)` for `x ≥ 0` (rather than MathNet's Boost rational
//!   polynomial). Both routes are accurate to ≈1e-13 — far inside the FlashLFQ rel-1e-6 parity
//!   budget — and the `Q` route avoids cancellation for large `x`, so the Normal CDF (which only ever
//!   evaluates `erfc` of a non-negative argument inside the scorer) matches C# to full working
//!   precision.
//! - [`Normal::cumulative_distribution`] is `0.5·erfc((mean − x)/(stddev·√2))` (MathNet `Normal.CDF`),
//!   and [`Gamma::cumulative_distribution`] is `P(shape, x·rate)` (MathNet `Gamma.CDF`, rate-parameterized).

use std::f64::consts::{E, LN_2, PI, SQRT_2};

// MathNet Lanczos coefficients (g = 10.900511, n = 10) — `SpecialFunctions.GammaLn`.
const GAMMA_R: f64 = 10.900511;
const GAMMA_DK: [f64; 11] = [
    2.48574089138753565546e-5,
    1.05142378581721974210,
    -3.45687097222016235469,
    4.51227709466894823700,
    -2.98285225323576655721,
    1.05639711577126713077,
    -1.95428773191645869583e-1,
    1.70970543404441224307e-2,
    -5.71926117404305781283e-4,
    4.63399473359905636708e-6,
    -2.71994908488607703910e-9,
];
const LOG_TWO_SQRT_E_OVER_PI: f64 = 0.6207822376352452223455184457816472122518527279025978;
const LN_PI: f64 = 1.1447298858494001741434273513530587116472948129153;

/// Natural log of the Gamma function, `ln Γ(z)`. Port of MathNet `SpecialFunctions.GammaLn`
/// (Lanczos approximation with reflection for `z < 0.5`).
pub fn gamma_ln(z: f64) -> f64 {
    if z < 0.5 {
        let mut s = GAMMA_DK[0];
        for i in 1..=10 {
            s += GAMMA_DK[i] / (i as f64 - z);
        }
        LN_PI
            - (PI * z).sin().ln()
            - s.ln()
            - LOG_TWO_SQRT_E_OVER_PI
            - ((0.5 - z) * ((0.5 - z + GAMMA_R) / E).ln())
    } else {
        let mut s = GAMMA_DK[0];
        for i in 1..=10 {
            s += GAMMA_DK[i] / (z + i as f64 - 1.0);
        }
        s.ln() + LOG_TWO_SQRT_E_OVER_PI + ((z - 0.5) * ((z - 0.5 + GAMMA_R) / E).ln())
    }
}

// Convergence epsilon and overflow guards from MathNet's incomplete-gamma routines.
const GAMMA_EPSILON: f64 = 0.000000000000001; // 1e-15
const GAMMA_BIG: f64 = 4503599627370496.0;
const GAMMA_BIG_INV: f64 = 2.22044604925031308085e-16;
// `Math.Log(double.MinValue)` ≈ -709.7827…; below this `exp(ax)` underflows to 0.
const LOG_MIN: f64 = -709.78271289338399;

/// Regularized lower incomplete gamma `P(a, x) = γ(a, x) / Γ(a)`. Port of MathNet
/// `SpecialFunctions.GammaLowerRegularized`.
pub fn gamma_lower_regularized(a: f64, x: f64) -> f64 {
    if a < 0.0 || x < 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return if x == 0.0 { f64::NAN } else { 1.0 };
    }
    if x == 0.0 {
        return 0.0;
    }

    let ax = a * x.ln() - x - gamma_ln(a);
    if ax < LOG_MIN {
        return if a < x { 1.0 } else { 0.0 };
    }

    if x <= 1.0 || x <= a {
        // Series expansion.
        let mut r2 = a;
        let mut c2 = 1.0;
        let mut ans2 = 1.0;
        loop {
            r2 += 1.0;
            c2 *= x / r2;
            ans2 += c2;
            if (c2 / ans2) <= GAMMA_EPSILON {
                break;
            }
        }
        return ax.exp() * ans2 / a;
    }

    // Continued fraction (Lentz), returning 1 - Q.
    let mut c = 0.0_f64;
    let mut y = 1.0 - a;
    let mut z = x + y + 1.0;
    let mut p3 = 1.0;
    let mut q3 = x;
    let mut p2 = x + 1.0;
    let mut q2 = z * x;
    let mut ans = p2 / q2;
    loop {
        c += 1.0;
        y += 1.0;
        z += 2.0;
        let yc = y * c;
        let p = p2 * z - p3 * yc;
        let q = q2 * z - q3 * yc;
        let error;
        if q != 0.0 {
            let nextans = p / q;
            error = ((ans - nextans) / nextans).abs();
            ans = nextans;
        } else {
            error = 1.0;
        }
        p3 = p2;
        p2 = p;
        q3 = q2;
        q2 = q;
        if p.abs() > GAMMA_BIG {
            p3 *= GAMMA_BIG_INV;
            p2 *= GAMMA_BIG_INV;
            q3 *= GAMMA_BIG_INV;
            q2 *= GAMMA_BIG_INV;
        }
        if error <= GAMMA_EPSILON {
            break;
        }
    }
    1.0 - ax.exp() * ans
}

/// Regularized upper incomplete gamma `Q(a, x) = Γ(a, x) / Γ(a) = 1 − P(a, x)`. Port of MathNet
/// `SpecialFunctions.GammaUpperRegularized` (continued fraction for `x ≥ 1 && x > a`, else `1 − P`).
pub fn gamma_upper_regularized(a: f64, x: f64) -> f64 {
    if x < 1.0 || x <= a {
        return 1.0 - gamma_lower_regularized(a, x);
    }

    let ax0 = a * x.ln() - x - gamma_ln(a);
    if ax0 < LOG_MIN {
        return if a < x { 0.0 } else { 1.0 };
    }
    let ax = ax0.exp();

    let mut y = 1.0 - a;
    let mut z = x + y + 1.0;
    let mut c = 0.0_f64;
    let mut pkm2 = 1.0;
    let mut qkm2 = x;
    let mut pkm1 = x + 1.0;
    let mut qkm1 = z * x;
    let mut ans = pkm1 / qkm1;
    loop {
        c += 1.0;
        y += 1.0;
        z += 2.0;
        let yc = y * c;
        let pk = pkm1 * z - pkm2 * yc;
        let qk = qkm1 * z - qkm2 * yc;
        let t;
        if qk != 0.0 {
            let r = pk / qk;
            t = ((ans - r) / r).abs();
            ans = r;
        } else {
            t = 1.0;
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if pk.abs() > GAMMA_BIG {
            pkm2 *= GAMMA_BIG_INV;
            pkm1 *= GAMMA_BIG_INV;
            qkm2 *= GAMMA_BIG_INV;
            qkm1 *= GAMMA_BIG_INV;
        }
        if t <= GAMMA_EPSILON {
            break;
        }
    }
    ax * ans
}

/// The error function `erf(x)`. Computed via `erf(x) = P(½, x²)` for `x ≥ 0` and odd symmetry
/// `erf(−x) = −erf(x)` — see the module note on the gamma route.
pub fn erf(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == 0.0 {
        return 0.0;
    }
    if x < 0.0 {
        return -erf(-x);
    }
    if x.is_infinite() {
        return 1.0;
    }
    gamma_lower_regularized(0.5, x * x)
}

/// The complementary error function `erfc(x) = 1 − erf(x)`. Computed via `erfc(x) = Q(½, x²)` for
/// `x ≥ 0` (no cancellation for large `x`) and `erfc(−x) = 2 − erfc(x)`.
pub fn erfc(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == 0.0 {
        return 1.0;
    }
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    if x.is_infinite() {
        return 0.0;
    }
    gamma_upper_regularized(0.5, x * x)
}

/// Base-2 logarithm computed the way C# `Math.Log(x, 2)` does — `ln(x) / ln(2)` — so the scorer's
/// log-intensity / log-fold-change values match bit-for-bit (distinct from the fused `f64::log2`).
pub fn log2(x: f64) -> f64 {
    x.ln() / LN_2
}

/// A normal (Gaussian) distribution, mirroring the slice of `MathNet.Numerics.Distributions.Normal`
/// the MBR scorer uses: its [`mean`](Self::mean), [`std_dev`](Self::std_dev), and CDF.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Normal {
    /// Distribution mean (C# `Mean`).
    pub mean: f64,
    /// Distribution standard deviation (C# `StdDev`).
    pub std_dev: f64,
}

impl Normal {
    /// Constructs a normal distribution (C# `new Normal(mean, stddev)`). No validation; callers
    /// guard parameters via [`Normal::is_valid_parameter_set`] exactly as the C# scorer does.
    pub fn new(mean: f64, std_dev: f64) -> Normal {
        Normal { mean, std_dev }
    }

    /// CDF at `x`: `0.5·erfc((mean − x)/(stddev·√2))` (MathNet `Normal.CDF`).
    pub fn cumulative_distribution(&self, x: f64) -> f64 {
        0.5 * erfc((self.mean - x) / (self.std_dev * SQRT_2))
    }

    /// Whether `(mean, stddev)` is a valid Normal parameterization — `stddev ≥ 0 && !mean.is_nan()`
    /// (MathNet `Normal.IsValidParameterSet`).
    pub fn is_valid_parameter_set(mean: f64, std_dev: f64) -> bool {
        std_dev >= 0.0 && !mean.is_nan()
    }
}

/// A gamma distribution in MathNet's **rate** parameterization (`shape` = α, `rate` = β), mirroring
/// the slice of `MathNet.Numerics.Distributions.Gamma` the MBR scorer uses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Gamma {
    /// Shape parameter α (C# `Shape`).
    pub shape: f64,
    /// Rate parameter β (C# `Rate`).
    pub rate: f64,
}

impl Gamma {
    /// Constructs a gamma distribution (C# `new Gamma(shape, rate)`).
    pub fn new(shape: f64, rate: f64) -> Gamma {
        Gamma { shape, rate }
    }

    /// CDF at `x`: the regularized lower incomplete gamma `P(shape, x·rate)` (MathNet `Gamma.CDF`).
    pub fn cumulative_distribution(&self, x: f64) -> f64 {
        gamma_lower_regularized(self.shape, x * self.rate)
    }

    /// Whether `(shape, rate)` is a valid Gamma parameterization — both `≥ 0` and non-NaN
    /// (MathNet `Gamma.IsValidParameterSet`).
    pub fn is_valid_parameter_set(shape: f64, rate: f64) -> bool {
        shape >= 0.0 && rate >= 0.0 && !shape.is_nan() && !rate.is_nan()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gamma_ln_matches_known_values() {
        // ln Γ(0.5) = ln(√π) = 0.5·ln(π).
        assert!((gamma_ln(0.5) - 0.5723649429247001).abs() < 1e-13);
        // Γ(5) = 24 -> ln 24.
        assert!((gamma_ln(5.0) - 24.0_f64.ln()).abs() < 1e-12);
        // Γ(1) = 1 -> 0; Γ(2) = 1 -> 0.
        assert!(gamma_ln(1.0).abs() < 1e-12);
        assert!(gamma_ln(2.0).abs() < 1e-12);
    }

    #[test]
    fn erf_matches_known_values() {
        assert_eq!(erf(0.0), 0.0);
        assert!((erf(1.0) - 0.8427007929497149).abs() < 1e-12);
        assert!((erf(0.5) - 0.5204998778130465).abs() < 1e-12);
        assert!((erf(2.0) - 0.9953222650189527).abs() < 1e-12);
        // Odd symmetry.
        assert!((erf(-1.3) + erf(1.3)).abs() < 1e-15);
    }

    #[test]
    fn erfc_matches_known_values_and_complements_erf() {
        assert_eq!(erfc(0.0), 1.0);
        assert!((erfc(1.0) - 0.15729920705028513).abs() < 1e-12);
        // erfc(x) = 1 - erf(x) across a range, including the large-x tail.
        for &x in &[0.1, 0.7, 1.5, 3.0, 5.0] {
            assert!((erfc(x) - (1.0 - erf(x))).abs() < 1e-12, "x = {x}");
        }
        // Negative-argument reflection.
        assert!((erfc(-1.0) - (2.0 - erfc(1.0))).abs() < 1e-15);
    }

    #[test]
    fn normal_cdf_matches_known_values() {
        let std_normal = Normal::new(0.0, 1.0);
        assert!((std_normal.cumulative_distribution(0.0) - 0.5).abs() < 1e-13);
        // Standard normal CDF at 1 / -1 (≈ 0.8413 / 0.1587).
        assert!((std_normal.cumulative_distribution(1.0) - 0.8413447460685429).abs() < 1e-12);
        assert!((std_normal.cumulative_distribution(-1.0) - 0.15865525393145707).abs() < 1e-12);
        // Shifted/scaled normal: CDF at the mean is 0.5.
        let shifted = Normal::new(3.0, 2.0);
        assert!((shifted.cumulative_distribution(3.0) - 0.5).abs() < 1e-13);
    }

    #[test]
    fn gamma_cdf_matches_known_values() {
        // Exponential as Gamma(shape=1, rate=1): CDF = 1 - e^{-x}.
        let exp = Gamma::new(1.0, 1.0);
        for &x in &[0.0, 0.5, 1.0, 2.5] {
            assert!(
                (exp.cumulative_distribution(x) - (1.0 - (-x).exp())).abs() < 1e-12,
                "x = {x}"
            );
        }
        // Gamma(2, 1) CDF = 1 - e^{-x}(1 + x).
        let g2 = Gamma::new(2.0, 1.0);
        let x = 1.5;
        let expected = 1.0 - (-x as f64).exp() * (1.0 + x);
        assert!((g2.cumulative_distribution(x) - expected).abs() < 1e-12);
    }

    #[test]
    fn parameter_validity_checks() {
        assert!(Normal::is_valid_parameter_set(0.0, 1.0));
        assert!(Normal::is_valid_parameter_set(0.0, 0.0)); // stddev >= 0
        assert!(!Normal::is_valid_parameter_set(0.0, -1.0));
        assert!(!Normal::is_valid_parameter_set(f64::NAN, 1.0));

        assert!(Gamma::is_valid_parameter_set(2.0, 3.0));
        assert!(!Gamma::is_valid_parameter_set(-1.0, 3.0));
        assert!(!Gamma::is_valid_parameter_set(2.0, f64::NAN));
    }
}
