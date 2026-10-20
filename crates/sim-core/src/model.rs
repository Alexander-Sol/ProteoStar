//! The forward-model pieces: `g_p(t)`, `f_p(z)`, `σ_m(mz)` and the isotope envelope `φ`.
//!
//! Ports of `mzLib/TopDownSimulator/Model/*`. Arithmetic is written in the same operation order
//! as the C# so that results agree to the bit, not just within tolerance.

use std::f64::consts::PI;

use flashlfq_core::isotopic_envelope::PROTON_MASS;

use crate::averagine::AveragineTable;
use crate::spec::{SimChargeDistribution, SimPeakWidthModel, SimProteoform, SimRtProfile};

/// mzLib `ToMz` for a positive charge.
#[inline]
pub fn to_mz(mass: f64, charge: i32) -> f64 {
    mass / charge.abs() as f64 + charge.signum() as f64 * PROTON_MASS
}

/// C# `Math.Sign` as a double: 0 at 0, unlike `f64::signum`.
#[inline]
fn sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// Exponentially-modified Gaussian elution profile, unit area. Port of `EmgProfile`.
#[derive(Debug, Clone, Copy)]
pub struct EmgProfile {
    pub mu: f64,
    pub sigma: f64,
    pub tau: f64,
}

impl From<SimRtProfile> for EmgProfile {
    fn from(r: SimRtProfile) -> Self {
        Self { mu: r.mu, sigma: r.sigma, tau: r.tau }
    }
}

impl EmgProfile {
    pub fn evaluate(&self, t: f64) -> f64 {
        let (mu, sigma, tau) = (self.mu, self.sigma, self.tau);
        if tau <= 0.0 {
            let z0 = (t - mu) / sigma;
            return (-0.5 * z0 * z0).exp() / (sigma * (2.0 * PI).sqrt());
        }

        let arg = (mu - t) / tau + sigma * sigma / (2.0 * tau * tau);
        let z = (t - mu - sigma * sigma / tau) / (sigma * 2f64.sqrt());
        arg.exp() / (2.0 * tau) * erfc(-z)
    }
}

fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}

/// Abramowitz & Stegun 7.1.26 — the C# approximation, kept rather than a better `erf` so the
/// RT profile matches it exactly.
fn erf(x: f64) -> f64 {
    let sign = sign(x);
    let x = x.abs();
    const A1: f64 = 0.254829592;
    const A2: f64 = -0.284496736;
    const A3: f64 = 1.421413741;
    const A4: f64 = -1.453152027;
    const A5: f64 = 1.061405429;
    const P: f64 = 0.3275911;
    let t = 1.0 / (1.0 + P * x);
    let y = 1.0 - (((((A5 * t + A4) * t) + A3) * t + A2) * t + A1) * t * (-x * x).exp();
    sign * y
}

/// `f_p(z)`. Port of `GaussianChargeDistribution`.
///
/// Note this is the Gaussian *density* in z, not renormalized over the charge range — the C#
/// does not normalize either, and abundances are fitted against this exact form.
#[derive(Debug, Clone, Copy)]
pub struct ChargeDistribution {
    pub mu_z: f64,
    pub sigma_z: f64,
}

impl From<SimChargeDistribution> for ChargeDistribution {
    fn from(c: SimChargeDistribution) -> Self {
        let SimChargeDistribution::Gaussian { mu_z, sigma_z } = c;
        Self { mu_z, sigma_z }
    }
}

impl ChargeDistribution {
    pub fn evaluate(&self, z: i32) -> f64 {
        let d = (z as f64 - self.mu_z) / self.sigma_z;
        (-0.5 * d * d).exp() / (self.sigma_z * (2.0 * PI).sqrt())
    }
}

/// `σ_m` as a function of m/z. Port of `IPeakWidthModel`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PeakWidthModel {
    /// `σ = k·(m/z)^1.5`.
    Orbitrap { k: f64 },
    Constant { sigma: f64 },
}

impl From<SimPeakWidthModel> for PeakWidthModel {
    fn from(m: SimPeakWidthModel) -> Self {
        match m {
            SimPeakWidthModel::Orbitrap { k } => Self::Orbitrap { k },
            SimPeakWidthModel::Constant { sigma } => Self::Constant { sigma },
        }
    }
}

impl PeakWidthModel {
    #[inline]
    pub fn sigma_at(&self, mz: f64) -> f64 {
        match *self {
            Self::Orbitrap { k } => {
                if mz <= 0.0 {
                    0.0
                } else {
                    k * mz.powf(1.5)
                }
            }
            Self::Constant { sigma } => sigma,
        }
    }
}

/// One proteoform's parameter vector, decoded from the spec.
#[derive(Debug, Clone)]
pub struct ProteoformModel {
    pub monoisotopic_mass: f64,
    pub abundance: f64,
    pub rt: EmgProfile,
    pub charge: ChargeDistribution,
}

impl From<&SimProteoform> for ProteoformModel {
    fn from(p: &SimProteoform) -> Self {
        Self {
            monoisotopic_mass: p.monoisotopic_mass,
            abundance: p.abundance,
            rt: p.rt.into(),
            charge: p.charge.into(),
        }
    }
}

/// The averagine isotope envelope anchored on a monoisotopic mass. Port of the constructor of
/// `IsotopeEnvelopeKernel` — the evaluation half lives in the centroider, which is the only caller.
///
/// Isotopologues are in **ascending** neutral mass. Index 0 is not necessarily the monoisotopic
/// peak: averagine's 1e-8 cutoff drops it above ~31 kDa, so compute the monoisotopic m/z from the
/// mass instead of reading it from here.
#[derive(Debug, Clone)]
pub struct IsotopeEnvelope {
    pub neutral_masses: Vec<f64>,
    /// Sums to 1.
    pub intensities: Vec<f64>,
}

impl IsotopeEnvelope {
    pub fn new(monoisotopic_mass: f64) -> Self {
        let table = AveragineTable::get();
        let idx = table.most_intense_mass_index(monoisotopic_mass);
        // Both rows are most-intense-first, so [0] is the most intense isotopologue, which is
        // what DiffToMonoisotopic is measured from.
        let theor_masses = &table.all_masses[idx];
        let theor_intensities = &table.all_intensities[idx];

        let averagine_monoisotopic = theor_masses[0] - table.diff_to_monoisotopic[idx];
        let shift = monoisotopic_mass - averagine_monoisotopic;

        // Sum in the table's own order so the normalization matches C# to the bit.
        let sum: f64 = theor_intensities.iter().fold(0.0, |acc, &v| acc + v);
        let mut pairs: Vec<(f64, f64)> = theor_masses
            .iter()
            .zip(theor_intensities)
            .map(|(&m, &i)| (m + shift, i / sum))
            .collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));

        Self {
            neutral_masses: pairs.iter().map(|p| p.0).collect(),
            intensities: pairs.iter().map(|p| p.1).collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.neutral_masses.len()
    }

    pub fn is_empty(&self) -> bool {
        self.neutral_masses.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erf_is_zero_at_zero_and_odd() {
        assert_eq!(erf(0.0), 0.0);
        assert_eq!(erf(0.7), -erf(-0.7));
    }

    #[test]
    fn emg_integrates_to_one() {
        let g = EmgProfile { mu: 20.0, sigma: 0.22, tau: 0.08 };
        let dt = 1e-4;
        let area: f64 = (0..100_000).map(|i| g.evaluate(15.0 + i as f64 * dt) * dt).sum();
        // A&S erf's ~1e-7 absolute error is amplified in the leading tail, where the EMG is a
        // huge exp() times a tiny erfc(); the C# profile carries the same ~1e-4 excess area.
        assert!((area - 1.0).abs() < 1e-3, "area {area}");
    }

    #[test]
    fn envelope_is_ascending_normalized_and_anchored() {
        for mass in [5000.25, 35000.75] {
            let env = IsotopeEnvelope::new(mass);
            assert!(env.neutral_masses.windows(2).all(|w| w[1] > w[0]));
            let total: f64 = env.intensities.iter().sum();
            assert!((total - 1.0).abs() < 1e-12);
            // Anchored: the isotopologues sit on the 1.00235 Da ladder from the monoisotopic mass.
            let offset = (env.neutral_masses[0] - mass) / 1.002_35;
            assert!((offset - offset.round()).abs() < 0.05, "mass {mass}: offset {offset}");
        }
        // Above ~31 kDa averagine drops the monoisotopic peak.
        assert!(IsotopeEnvelope::new(35000.75).neutral_masses[0] > 35000.75 + 0.5);
    }
}
