//! Isotopic-envelope search — port of `FlashLfqEngine.GetIsotopicEnvelopes` and
//! `FlashLfqEngine.CheckIsotopicEnvelopeCorrelation`.
//!
//! Given an XIC (a list of "most abundant mass" peaks traced across retention time, one per
//! scan) and the theoretical isotope distribution of an identification, this finds, in each
//! scan, the set of isotopic peaks that match the theoretical envelope and sums their
//! intensities into an [`IsotopicEnvelope`].
//!
//! ## Fidelity notes
//! - The ±1 mass-shift off-by-one logic is ported verbatim: a `massShiftToIsotopePeaks`
//!   dictionary keyed `{-1, 0, 1}` records, for each hypothesis (the observed peak is the
//!   "real" monoisotopic species, or is off by one ¹³C either way), the experimental and
//!   theoretical isotopic peaks found. Only the `0` (accurate) hypothesis contributes to the
//!   returned envelope; the `-1`/`+1` hypotheses exist purely so the correlation gate can
//!   reject envelopes that are better explained by a neighbouring species.
//! - The intensity-ratio gate is `/4 .. *4`: an experimental isotope peak is accepted only if
//!   its intensity is within a factor of four of the theoretical intensity
//!   (`theoreticalAbundance * observedPeakIntensity`).
//! - The peak walk starts at `peakfindingMassIndex` (the most abundant isotope) and works
//!   **down** to the monoisotopic peak (`direction = -1`, starting one below) and **up** through
//!   the heavier isotopes (`direction = +1`), breaking each direction as soon as an expected
//!   peak is missing or out of the intensity window.
//! - `peak.M.ToMass` / `mass.ToMz` use the proton-mass conversions from mzLib
//!   `Chemistry/ClassExtensions.cs`. `peak.M` is `f32`, so `ToMass` is computed in `f32` for
//!   the `|charge| * m/z` term then widened (matching the C# `float` overload); `ToMz` is all
//!   `f64` (the C# `double` overload), since the isotope mass it converts is `double`.
//! - The Pearson correlation is the MathNet `Correlation.Pearson` Welford-style online
//!   algorithm, ported in [`pearson`] so the running means/variances accumulate in the same
//!   order C# uses. An empty or constant input yields `NaN`, which the gate maps to `-1`.
//!
//! Unlike C#, which threads a `SpectraFileInfo` through an `IndexingEngineDictionary`, this
//! takes a single [`PeakIndexingEngine`] directly (the m/z index for the one file in play).
//! The theoretical distribution is the [`ExpectedIsotopePeak`] slice produced by
//! [`crate::theoretical_isotope_distribution::expected_isotope_peaks`] (mzLib's
//! `ModifiedSequenceToIsotopicDistribution[modifiedSequence]`).

use crate::peak_indexing::{IndexedMassSpectralPeak, PeakIndexingEngine};
use crate::theoretical_isotope_distribution::ExpectedIsotopePeak;
use crate::tolerance::PpmTolerance;

/// Mass of a proton, in daltons. mzLib `Chemistry/Constants.cs` `ProtonMass`.
pub const PROTON_MASS: f64 = 1.007276466879;

/// Mass difference between ¹³C and ¹²C, in daltons. mzLib `Chemistry/Constants.cs` `C13MinusC12`.
pub const C13_MINUS_C12: f64 = 1.00335483810;

/// The summed intensities of all isotope peaks detected in a single MS1 scan for one species.
///
/// Port of `FlashLFQ.IsotopicEnvelope`. The reported [`intensity`](Self::intensity) is the
/// summed isotopic-peak intensity divided by the charge state (C# `Intensity = intensity /
/// chargeState`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IsotopicEnvelope {
    /// The most abundant isotopic peak used for peak finding (C# `IndexedPeak`).
    pub indexed_peak: IndexedMassSpectralPeak,
    /// The charge state this envelope was found at.
    pub charge_state: i32,
    /// The summed isotopic-peak intensity for this scan, divided by the charge state. May
    /// include imputed intensity for expected isotopes that were not observed.
    pub intensity: f64,
    /// Pearson correlation of the observed isotope intensities to the theoretical abundances.
    pub pearson_correlation: f64,
}

impl IsotopicEnvelope {
    /// Builds an envelope, dividing the summed `intensity` by `charge_state` (C# constructor).
    pub fn new(
        monoisotopic_peak: IndexedMassSpectralPeak,
        charge_state: i32,
        intensity: f64,
        pearson_correlation: f64,
    ) -> Self {
        IsotopicEnvelope {
            indexed_peak: monoisotopic_peak,
            charge_state,
            intensity: intensity / charge_state as f64,
            pearson_correlation,
        }
    }
}

/// One recorded isotopic peak hypothesis: the experimental intensity found, the theoretical
/// intensity expected, and the theoretical (expected) mass. Mirrors the C# value tuple
/// `(double expIntensity, double theorIntensity, double theorMass)`.
#[derive(Debug, Clone, Copy)]
struct IsotopePeakRecord {
    exp_intensity: f64,
    theor_intensity: f64,
    theor_mass: f64,
}

/// Converts an m/z (here an `f32` peak m/z) to a neutral mass at the given charge, assuming the
/// charge comes from protons. Port of `ClassExtensions.ToMass(this float, int)`: the
/// `|charge| * m/z` term is computed in `f32` (the C# `int * float` promotion) then widened.
pub fn mz_to_mass_f32(mz: f32, charge: i32) -> f64 {
    (charge.abs() as f32 * mz) as f64 - charge as f64 * PROTON_MASS
}

/// Converts a neutral mass to m/z at the given charge, assuming the charge comes from protons.
/// Port of `ClassExtensions.ToMz(this double, int)` (all-`f64`).
pub fn mass_to_mz_f64(mass: f64, charge: i32) -> f64 {
    mass / charge.abs() as f64 + (charge.signum() as f64) * PROTON_MASS
}

/// Finds, for each peak in `xic`, the isotopic envelope around it that matches the theoretical
/// distribution, returning one [`IsotopicEnvelope`] per scan whose envelope passes the
/// correlation gate. Faithful port of `FlashLfqEngine.GetIsotopicEnvelopes`.
///
/// - `engine`: the m/z index for the spectra file (C# `IndexingEngineDictionary[spectraFile]`).
/// - `xic`: the most-abundant-mass peaks traced across scans (one candidate per scan).
/// - `isotope_mass_shifts`: the theoretical `(massShift, normalizedAbundance)` distribution for
///   the identification's modified sequence (C#
///   `ModifiedSequenceToIsotopicDistribution[modifiedSequence]`).
/// - `monoisotopic_mass`, `peakfinding_mass`: the identification's masses.
/// - `charge_state`: the charge state being searched.
/// - `num_isotopes_required`: minimum number of accurate (shift-0) isotope peaks required
///   (C# `FlashParams.NumIsotopesRequired`).
/// - `isotope_ppm_tolerance`: the ppm tolerance for matching isotope peaks
///   (C# `FlashParams.IsotopePpmTolerance`).
#[allow(clippy::too_many_arguments)]
pub fn get_isotopic_envelopes(
    engine: &PeakIndexingEngine,
    xic: &[IndexedMassSpectralPeak],
    isotope_mass_shifts: &[ExpectedIsotopePeak],
    monoisotopic_mass: f64,
    peakfinding_mass: f64,
    charge_state: i32,
    num_isotopes_required: usize,
    isotope_ppm_tolerance: f64,
) -> Vec<IsotopicEnvelope> {
    let mut isotopic_envelopes: Vec<IsotopicEnvelope> = Vec::new();

    if isotope_mass_shifts.len() < num_isotopes_required {
        return isotopic_envelopes;
    }

    let isotope_tolerance = PpmTolerance::new(isotope_ppm_tolerance);

    let theoretical_isotope_mass_shifts: Vec<f64> =
        isotope_mass_shifts.iter().map(|p| p.mass_shift).collect();
    let theoretical_isotope_abundances: Vec<f64> = isotope_mass_shifts
        .iter()
        .map(|p| p.normalized_abundance)
        .collect();

    // (int)Math.Round(PeakfindingMass - MonoisotopicMass, 0) — banker's rounding.
    let peakfinding_mass_index =
        (peakfinding_mass - monoisotopic_mass).round_ties_even() as i32;

    let mut experimental_isotope_intensities = vec![0.0f64; isotope_mass_shifts.len()];

    // The -1 / 0 / +1 mass-shift hypotheses (negative off-by-one, accurate, positive off-by-one).
    // Index 0 -> key -1, index 1 -> key 0, index 2 -> key +1 (the C# dict insertion order).
    let shift_keys: [i32; 3] = [-1, 0, 1];
    let directions: [i32; 2] = [-1, 1];

    for peak in xic {
        for v in experimental_isotope_intensities.iter_mut() {
            *v = 0.0;
        }
        let mut mass_shift_to_isotope_peaks: [Vec<IsotopePeakRecord>; 3] =
            [Vec::new(), Vec::new(), Vec::new()];

        // isotope masses are calculated relative to the observed peak
        let observed_mass = mz_to_mass_f32(peak.m(), charge_state);
        let observed_mass_error = observed_mass - peakfinding_mass;

        for (shift_idx, &shift_key) in shift_keys.iter().enumerate() {
            for &direction in &directions {
                let start = if direction == -1 {
                    peakfinding_mass_index - 1
                } else {
                    peakfinding_mass_index
                };

                let mut i = start;
                while i < theoretical_isotope_abundances.len() as i32 && i >= 0 {
                    let idx = i as usize;
                    let isotope_mass = monoisotopic_mass
                        + observed_mass_error
                        + theoretical_isotope_mass_shifts[idx]
                        + shift_key as f64 * C13_MINUS_C12;
                    let theoretical_isotope_intensity =
                        theoretical_isotope_abundances[idx] * peak.intensity as f64;

                    let isotope_peak = engine.get_indexed_peak(
                        mass_to_mz_f64(isotope_mass, charge_state),
                        peak.zero_based_scan_index,
                        &isotope_tolerance,
                    );

                    let isotope_peak = match isotope_peak {
                        Some(p) => p,
                        None => break,
                    };
                    let exp_int = isotope_peak.intensity as f64;
                    if exp_int < theoretical_isotope_intensity / 4.0
                        || exp_int > theoretical_isotope_intensity * 4.0
                    {
                        break;
                    }

                    mass_shift_to_isotope_peaks[shift_idx].push(IsotopePeakRecord {
                        exp_intensity: exp_int,
                        theor_intensity: theoretical_isotope_intensity,
                        theor_mass: isotope_mass,
                    });
                    if shift_key == 0 {
                        experimental_isotope_intensities[idx] = exp_int;
                    }

                    i += direction;
                }
            }
        }

        // check number of accurate (shift-0) isotope peaks observed
        if mass_shift_to_isotope_peaks[1].len() < num_isotopes_required {
            continue;
        }

        if let Some(pearson_corr) = check_isotopic_envelope_correlation(
            &mut mass_shift_to_isotope_peaks,
            peak,
            charge_state,
            &isotope_tolerance,
            engine,
        ) {
            // impute unobserved isotope peak intensities
            let pf_idx = peakfinding_mass_index as usize;
            for i in 0..experimental_isotope_intensities.len() {
                if experimental_isotope_intensities[i] == 0.0 {
                    experimental_isotope_intensities[i] =
                        theoretical_isotope_abundances[i] * experimental_isotope_intensities[pf_idx];
                }
            }

            let summed: f64 = experimental_isotope_intensities.iter().sum();
            isotopic_envelopes.push(IsotopicEnvelope::new(
                *peak,
                charge_state,
                summed,
                pearson_corr,
            ));
        }
    }

    isotopic_envelopes
}

/// Checks whether the experimental isotope intensities are best described by the theoretical
/// abundances rather than by a ±1 ¹³C-shifted neighbour. Faithful port of
/// `FlashLfqEngine.CheckIsotopicEnvelopeCorrelation`.
///
/// Returns `Some(pearsonCorrelation)` (the pre-padding shift-0 correlation) when the envelope
/// passes the gate (`pearson > 0.7`, and neither shifted alternative beats the padded shift-0
/// correlation by `0.1` or more), else `None`. The `mass_shift_to_isotope_peaks` lists are
/// mutated in place to append the "unexpected" padding peak, exactly as C# does.
fn check_isotopic_envelope_correlation(
    mass_shift_to_isotope_peaks: &mut [Vec<IsotopePeakRecord>; 3],
    peak: &IndexedMassSpectralPeak,
    charge_state: i32,
    isotope_tolerance: &PpmTolerance,
    engine: &PeakIndexingEngine,
) -> Option<f64> {
    // correlation of experimental isotope intensities to theoretical abundances (pre-padding)
    let pearson_correlation = pearson_of(&mass_shift_to_isotope_peaks[1]);

    // check for unexpected peaks: for each non-empty hypothesis, look one ¹³C below its lowest
    // expected mass and append that peak (or a zero) as padding.
    for records in mass_shift_to_isotope_peaks.iter_mut() {
        if records.is_empty() {
            continue;
        }
        let unexpected_mass = records
            .iter()
            .map(|r| r.theor_mass)
            .fold(f64::INFINITY, f64::min)
            - C13_MINUS_C12;
        let unexpected_peak = engine.get_indexed_peak(
            mass_to_mz_f64(unexpected_mass, charge_state),
            peak.zero_based_scan_index,
            isotope_tolerance,
        );
        match unexpected_peak {
            None => records.push(IsotopePeakRecord {
                exp_intensity: 0.0,
                theor_intensity: 0.0,
                theor_mass: unexpected_mass,
            }),
            Some(p) => records.push(IsotopePeakRecord {
                exp_intensity: p.intensity as f64,
                theor_intensity: 0.0,
                theor_mass: unexpected_mass,
            }),
        }
    }

    let corr_with_padding = pearson_of(&mass_shift_to_isotope_peaks[1]);
    let mut corr_shifted_left = pearson_of(&mass_shift_to_isotope_peaks[0]);
    let mut corr_shifted_right = pearson_of(&mass_shift_to_isotope_peaks[2]);

    if corr_shifted_left.is_nan() {
        corr_shifted_left = -1.0;
    }
    if corr_shifted_right.is_nan() {
        corr_shifted_right = -1.0;
    }

    let pass = pearson_correlation > 0.7
        && corr_shifted_left - corr_with_padding < 0.1
        && corr_shifted_right - corr_with_padding < 0.1;

    if pass {
        Some(pearson_correlation)
    } else {
        None
    }
}

/// Pearson correlation of the experimental vs theoretical intensities in a record list.
fn pearson_of(records: &[IsotopePeakRecord]) -> f64 {
    let exp: Vec<f64> = records.iter().map(|r| r.exp_intensity).collect();
    let theor: Vec<f64> = records.iter().map(|r| r.theor_intensity).collect();
    pearson(&exp, &theor)
}

/// Pearson product-moment correlation coefficient of two equal-length samples.
///
/// Port of MathNet.Numerics `Correlation.Pearson`: a Welford-style single-pass accumulation of
/// the running means, variances, and covariance (deliberately using differences from the running
/// mean rather than summing products, for numerical accuracy). Returns `NaN` for empty input or
/// when either sample has zero variance (a constant array), matching MathNet (`r / sqrt(0)`).
///
/// # Panics
/// Panics if `data_a` and `data_b` differ in length (the C# code throws
/// `ArgumentOutOfRangeException`).
pub fn pearson(data_a: &[f64], data_b: &[f64]) -> f64 {
    assert_eq!(
        data_a.len(),
        data_b.len(),
        "pearson: input samples differ in length"
    );

    let mut n: i64 = 0;
    let mut r = 0.0f64;
    let mut mean_a = 0.0f64;
    let mut mean_b = 0.0f64;
    let mut var_a = 0.0f64;
    let mut var_b = 0.0f64;

    for (&current_a, &current_b) in data_a.iter().zip(data_b.iter()) {
        let delta_a = current_a - mean_a;
        n += 1;
        let scale_delta_a = delta_a / n as f64;

        let delta_b = current_b - mean_b;
        let scale_delta_b = delta_b / n as f64;

        mean_a += scale_delta_a;
        mean_b += scale_delta_b;

        var_a += scale_delta_a * delta_a * (n - 1) as f64;
        var_b += scale_delta_b * delta_b * (n - 1) as f64;
        r += delta_a * delta_b * (n - 1) as f64 / n as f64;
    }

    r / (var_a * var_b).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::peak_indexing::{PeakIndexingEngine, Scan};
    use crate::theoretical_isotope_distribution::ExpectedIsotopePeak;

    /// A clean three-isotope theoretical distribution: monoisotopic peak most abundant (1.0),
    /// then 0.5 and 0.2, spaced one ¹³C apart. The most abundant isotope is index 0, so the
    /// peakfinding mass equals the monoisotopic mass.
    fn three_isotope_distribution() -> Vec<ExpectedIsotopePeak> {
        vec![
            ExpectedIsotopePeak {
                mass_shift: 0.0,
                normalized_abundance: 1.0,
            },
            ExpectedIsotopePeak {
                mass_shift: C13_MINUS_C12,
                normalized_abundance: 0.5,
            },
            ExpectedIsotopePeak {
                mass_shift: 2.0 * C13_MINUS_C12,
                normalized_abundance: 0.2,
            },
        ]
    }

    /// Builds scans each carrying the three isotope peaks of a species at `mono` and charge
    /// `z`, with per-scan intensity multipliers forming an elution profile.
    fn envelope_scans(mono: f64, z: i32, base_intensity: f64, multipliers: &[f64]) -> Vec<Scan> {
        let dist = three_isotope_distribution();
        let mut scans = Vec::new();
        for (s, mult) in multipliers.iter().enumerate() {
            let mut mz = Vec::new();
            let mut intensity = Vec::new();
            for peak in &dist {
                // neutral mass of this isotope -> m/z at charge z
                let isotope_mass = mono + peak.mass_shift;
                mz.push(mass_to_mz_f64(isotope_mass, z));
                intensity.push(peak.normalized_abundance * base_intensity * mult);
            }
            scans.push(Scan {
                mz,
                intensity,
                one_based_scan_number: (s + 1) as i32,
                retention_time: 1.0 + s as f64 / 10.0,
                msn_order: 1,
            });
        }
        scans
    }

    /// Pulls the most-abundant-mass (mono, here) peak out of each scan to form the XIC the
    /// envelope search consumes.
    fn mono_xic(engine: &PeakIndexingEngine, mono: f64, z: i32, n: i32) -> Vec<IndexedMassSpectralPeak> {
        let mono_mz = mass_to_mz_f64(mono, z);
        let tol = PpmTolerance::new(10.0);
        let mut xic = Vec::new();
        for scan in 0..n {
            if let Some(p) = engine.get_indexed_peak(mono_mz, scan, &tol) {
                xic.push(*p);
            }
        }
        xic
    }

    #[test]
    fn pearson_of_identical_arrays_is_one() {
        let a = [1.0, 2.0, 3.0, 4.0];
        assert!((pearson(&a, &a) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pearson_of_anticorrelated_arrays_is_minus_one() {
        let a = [1.0, 2.0, 3.0, 4.0];
        let b = [4.0, 3.0, 2.0, 1.0];
        assert!((pearson(&a, &b) + 1.0).abs() < 1e-12);
    }

    #[test]
    fn pearson_of_constant_array_is_nan() {
        let a = [5.0, 5.0, 5.0];
        let b = [1.0, 2.0, 3.0];
        assert!(pearson(&a, &b).is_nan());
        assert!(pearson(&[], &[]).is_nan());
    }

    #[test]
    fn mz_mass_round_trip() {
        // ToMz then ToMass should recover the neutral mass (to f32 widening precision).
        let mass = 1000.0;
        let z = 2;
        let mz = mass_to_mz_f64(mass, z);
        let recovered = mz_to_mass_f32(mz as f32, z);
        assert!((recovered - mass).abs() < 1e-2);
    }

    #[test]
    fn finds_one_envelope_per_scan_with_high_correlation() {
        let mono = 1000.0;
        let z = 1;
        let base_intensity = 1.0e6;
        let multipliers = [1.0, 2.0, 4.0, 6.0, 4.0, 2.0, 1.0];
        let scans = envelope_scans(mono, z, base_intensity, &multipliers);
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let dist = three_isotope_distribution();
        let xic = mono_xic(&engine, mono, z, multipliers.len() as i32);
        assert_eq!(xic.len(), multipliers.len());

        // peakfinding mass == mono (most abundant isotope is index 0).
        let envelopes = get_isotopic_envelopes(&engine, &xic, &dist, mono, mono, z, 2, 5.0);

        // One envelope per scan: every scan carries the full, clean envelope.
        assert_eq!(envelopes.len(), multipliers.len());
        for env in &envelopes {
            assert!(
                env.pearson_correlation > 0.7,
                "expected high correlation, got {}",
                env.pearson_correlation
            );
            assert_eq!(env.charge_state, z);
        }

        // The apex scan (multiplier 6, index 3) should have the largest summed intensity.
        let apex = envelopes
            .iter()
            .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
            .unwrap();
        assert_eq!(apex.indexed_peak.zero_based_scan_index, 3);
        // Intensity at apex = (1.0 + 0.5 + 0.2) * base * 6 / charge.
        let expected_apex = (1.0 + 0.5 + 0.2) * base_intensity * 6.0 / z as f64;
        assert!((apex.intensity - expected_apex).abs() / expected_apex < 1e-4);
    }

    #[test]
    fn returns_empty_when_fewer_isotopes_than_required() {
        // A distribution with one peak cannot meet num_isotopes_required = 2.
        let dist = vec![ExpectedIsotopePeak {
            mass_shift: 0.0,
            normalized_abundance: 1.0,
        }];
        let scans = envelope_scans(1000.0, 1, 1.0e6, &[1.0, 2.0, 1.0]);
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let xic = mono_xic(&engine, 1000.0, 1, 3);
        let envelopes = get_isotopic_envelopes(&engine, &xic, &dist, 1000.0, 1000.0, 1, 2, 5.0);
        assert!(envelopes.is_empty());
    }

    #[test]
    fn finds_envelopes_at_charge_two() {
        let mono = 1600.0;
        let z = 2;
        let base_intensity = 5.0e5;
        let multipliers = [1.0, 3.0, 5.0, 3.0, 1.0];
        let scans = envelope_scans(mono, z, base_intensity, &multipliers);
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let dist = three_isotope_distribution();
        let xic = mono_xic(&engine, mono, z, multipliers.len() as i32);
        let envelopes = get_isotopic_envelopes(&engine, &xic, &dist, mono, mono, z, 2, 5.0);
        assert_eq!(envelopes.len(), multipliers.len());
        for env in &envelopes {
            assert!(env.pearson_correlation > 0.7);
            assert_eq!(env.charge_state, 2);
        }
    }
}
