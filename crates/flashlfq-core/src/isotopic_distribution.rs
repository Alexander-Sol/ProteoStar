//! Isotopic distribution (parity layer **L2**, part 1).
//!
//! A faithful port of mzLib's `Chemistry/IsotopicDistribution.cs` — the *fine-grained*
//! Kubinyi/MIDAs polynomial algorithm that turns a [`ChemicalFormula`] into a discrete
//! `(masses[], intensities[])` isotopic envelope. The reference is:
//!
//! > Molecular Isotopic Distribution Analysis (MIDAs) with Adjustable Mass Accuracy.
//! > Gelio Alves, Aleksy Y. Ogurtsov, and Yi-Kuo Yu. J. Am. Soc. Mass Spectrom. (2014).
//!
//! Everything here mirrors the C# line-for-line: the same constants
//! (`defaultFineResolution = 0.01`, `defaultMinProbability = 1e-200`,
//! `defaultMolecularWeightResolution = 1e-12`), the `fineResolution/2` split, the per-element
//! multinomial expansion (`MultiplyFinePolynomial` / `MultipleFinePolynomialRecursiveHelper`),
//! the cross-element convolution onto a binned grid (`MultiplyFineFinalPolynomial`), and the
//! **merge order** of `MergeFinePolynomial` (its 9 passes with the `k`-scaled threshold, where
//! the comparison reads the *running merged centroid* of bin `i`, not its original power).
//!
//! Why bit-for-bit care: this envelope feeds averagine fill-in (P1.7) and ultimately every
//! peak-tracing decision, so abundance drift compounds. The L2 parity gate is **P1.6**; this
//! task (P1.5) only builds the algorithm and produces the two arrays.

use std::cell::RefCell;

use crate::chemical_formula::ChemicalFormula;
use crate::periodic_table::periodic_table;

/// Default fine resolution (`IsotopicDistribution.defaultFineResolution`).
pub const DEFAULT_FINE_RESOLUTION: f64 = 0.01;
/// Default minimum probability cutoff (`IsotopicDistribution.defaultMinProbability`).
pub const DEFAULT_MIN_PROBABILITY: f64 = 1e-200;
/// Default molecular-weight resolution (`IsotopicDistribution.defaultMolecularWeightResolution`).
pub const DEFAULT_MOLECULAR_WEIGHT_RESOLUTION: f64 = 1e-12;

/// A discrete isotopic distribution: parallel `masses` / `intensities` arrays.
///
/// Mirrors the C# `IsotopicDistribution` (which exposes `Masses` and `Intensities`). The
/// arrays are in the merged-polynomial order produced by `CalculateFineGrain`
/// (sorted by mass during the merge, then filtered to non-empty bins).
#[derive(Debug, Clone, PartialEq)]
pub struct IsotopicDistribution {
    /// Isotopologue masses (Da).
    pub masses: Vec<f64>,
    /// Relative (un-normalised) intensities, one per mass.
    pub intensities: Vec<f64>,
}

/// One term of the Kubinyi polynomial: a centroid mass (`power`) and weight (`probability`).
///
/// Empty grid bins are marked with `NaN` in both fields, exactly as the C# struct uses
/// `double.NaN` as its sentinel. (The C# field is misspelled `Probablity`; we keep the
/// corrected name `probability`.)
#[derive(Debug, Clone, Copy, PartialEq)]
struct Polynomial {
    power: f64,
    probability: f64,
}

impl Polynomial {
    /// The empty/sentinel term (`{ Power = NaN, Probablity = NaN }`).
    const NAN: Polynomial = Polynomial {
        power: f64::NAN,
        probability: f64::NAN,
    };
}

/// One isotope's contribution within a single element's expansion (the C# `Composition` class).
#[derive(Debug, Clone, Copy)]
struct Composition {
    power: f64,
    probability: f64,
    log_probability: f64,
    molecular_weight: f64,
    atoms: i32,
}

thread_local! {
    /// Memoised `ln(n!)` table, built cumulatively in increasing `n` so the summation order
    /// is bit-identical to the C# `FactorLn` cache (`factorLnArray[j+1] = factorLnArray[j] +
    /// log(j+1)`). Indices 0 and 1 are `0.0`. The cache is process-local state in C#; a
    /// thread-local here keeps it deterministic without locking.
    static FACTOR_LN: RefCell<Vec<f64>> = const { RefCell::new(Vec::new()) };
}

/// `ln(n!)`, mirroring `IsotopicDistribution.FactorLn` (cumulative log-sum, memoised).
fn factor_ln(n: i32) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    let n = n as usize;
    FACTOR_LN.with(|cache| {
        let mut c = cache.borrow_mut();
        if c.is_empty() {
            c.push(0.0); // ln(0!) = 0
            c.push(0.0); // ln(1!) = 0
        }
        while c.len() <= n {
            let k = c.len(); // next index to fill == the factor being multiplied in
            let prev = c[k - 1];
            c.push(prev + (k as f64).ln());
        }
        c[n]
    })
}

/// `(fineResolution / 2, fineResolution)` — `GetNewFineAndMergeResolutions`.
fn get_new_fine_and_merge_resolutions(fine_resolution: f64) -> (f64, f64) {
    (fine_resolution / 2.0, fine_resolution)
}

impl IsotopicDistribution {
    /// Compute the fine-grained isotopic distribution of `formula` with mzLib's defaults.
    ///
    /// Equivalent to the C# `IsotopicDistribution.GetDistribution(ChemicalFormula)`.
    pub fn get_distribution(formula: &ChemicalFormula) -> IsotopicDistribution {
        IsotopicDistribution::get_distribution_with(
            formula,
            DEFAULT_FINE_RESOLUTION,
            DEFAULT_MIN_PROBABILITY,
            DEFAULT_MOLECULAR_WEIGHT_RESOLUTION,
        )
    }

    /// Compute the fine-grained isotopic distribution with explicit resolution parameters.
    ///
    /// Mirrors the four-argument C# `GetDistribution(formula, fineResolution, minProbability,
    /// molecularWeightResolution)`.
    pub fn get_distribution_with(
        formula: &ChemicalFormula,
        fine_resolution: f64,
        min_probability: f64,
        molecular_weight_resolution: f64,
    ) -> IsotopicDistribution {
        let (new_fine_resolution, merge_fine_resolution) =
            get_new_fine_and_merge_resolutions(fine_resolution);

        let pt = periodic_table();
        let mut elemental_composition: Vec<Vec<Composition>> = Vec::new();

        // For every element (with its count) build one Composition per isotope, ordered by
        // atomic mass (C# `elementAndCount.Key.Isotopes.OrderBy(iso => iso.AtomicMass)`).
        for (atomic_number, count) in formula.elements() {
            let element = pt
                .element_by_number(atomic_number)
                .expect("formula element references a loaded element");
            let mut isotopes: Vec<_> = element.isotopes.iter().collect();
            isotopes.sort_by(|a, b| {
                a.atomic_mass
                    .partial_cmp(&b.atomic_mass)
                    .expect("isotope masses are finite")
            });

            let mut isotope_composition = Vec::with_capacity(isotopes.len());
            for isotope in isotopes {
                isotope_composition.push(Composition {
                    atoms: count,
                    molecular_weight: isotope.atomic_mass,
                    power: isotope.atomic_mass,
                    probability: isotope.relative_abundance,
                    log_probability: 0.0,
                });
            }
            elemental_composition.push(isotope_composition);
        }

        // Normalise each element's isotope probabilities to sum to 1, then recompute
        // log-probability and snap the power onto the molecular-weight grid.
        for compositions in elemental_composition.iter_mut() {
            let sum_prob: f64 = compositions.iter().map(|c| c.probability).sum();
            for composition in compositions.iter_mut() {
                composition.probability /= sum_prob;
                composition.log_probability = composition.probability.ln();
                composition.power =
                    (composition.molecular_weight / molecular_weight_resolution + 0.5).floor();
            }
        }

        let mut dist = calculate_fine_grain(
            &elemental_composition,
            molecular_weight_resolution,
            merge_fine_resolution,
            new_fine_resolution,
            min_probability,
        );

        // Isotope-specified atoms (e.g. an explicit ¹³C) contribute a fixed mass offset.
        let mut additional_mass = 0.0;
        for ((atomic_number, mass_number), count) in formula.isotopes() {
            let isotope = pt
                .element_by_number(atomic_number)
                .and_then(|e| e.isotope_by_mass_number(mass_number))
                .expect("formula isotope references a loaded isotope");
            additional_mass += isotope.atomic_mass * count as f64;
        }
        for m in dist.masses.iter_mut() {
            *m += additional_mass;
        }

        dist
    }
}

/// Port of `CalculateFineGrain`: expand, merge, then read off `(mass, intensity)` pairs.
fn calculate_fine_grain(
    elemental_composition: &[Vec<Composition>],
    mw_resolution: f64,
    merge_fine_resolution: f64,
    fine_resolution: f64,
    fine_min_prob: f64,
) -> IsotopicDistribution {
    let mut f_polynomial =
        multiply_fine_polynomial(elemental_composition, fine_resolution, mw_resolution, fine_min_prob);
    let f_polynomial = merge_fine_polynomial(&mut f_polynomial, mw_resolution, merge_fine_resolution);

    // Convert polynomial terms to a spectrum. (C# also tracks totalProbability/basePeak but
    // never uses them for the returned arrays, so they are omitted here.)
    let count = f_polynomial.len();
    let mut masses = Vec::with_capacity(count);
    let mut intensities = Vec::with_capacity(count);
    for polynomial in &f_polynomial {
        masses.push(polynomial.power * mw_resolution);
        intensities.push(polynomial.probability);
    }
    IsotopicDistribution { masses, intensities }
}

/// Port of `MultiplyFinePolynomial`: build each element's polynomial, then convolve them.
fn multiply_fine_polynomial(
    elemental_composition: &[Vec<Composition>],
    fine_resolution: f64,
    mw_resolution: f64,
    fine_min_prob: f64,
) -> Vec<Polynomial> {
    const NC: i32 = 10;
    const NC_ADD_VALUE: i32 = 1;
    const N_ATOMS: i32 = 200;

    // Number of elements that actually contribute at least one isotope.
    let n = elemental_composition
        .iter()
        .filter(|c| !c.is_empty())
        .count();

    let mut f_polynomial: Vec<Vec<Polynomial>> = (0..n).map(|_| Vec::new()).collect();

    for k in 0..n {
        let composition = &elemental_composition[k];
        let size = composition.len();
        let atoms = composition[0].atoms;
        let nc_add = if atoms < N_ATOMS { 10 } else { NC_ADD_VALUE };

        if size == 1 {
            let probability = composition[0].probability;
            let n1 = (atoms as f64 * probability) as i32;
            let prob = factor_ln(atoms) - factor_ln(n1) + n1 as f64 * composition[0].log_probability;
            let prob = prob.exp();
            f_polynomial[k].push(Polynomial {
                power: n1 as f64 * composition[0].power,
                probability: prob,
            });
        } else {
            let mut means = vec![0i32; size];
            let mut stds = vec![0i32; size];
            // nPolynomialTerms is accumulated in C# but never used; omitted.
            for i in 0..size {
                let n1 = (composition[0].atoms as f64 * composition[i].probability) as i32;
                let s1 = (nc_add as f64
                    + NC as f64
                        * (composition[i].atoms as f64
                            * composition[i].probability
                            * (1.0 - composition[i].probability))
                            .sqrt())
                .ceil() as i32;
                means[i] = n1;
                stds[i] = s1;
            }

            // Recursion ranges over all isotopes except the last (whose count is the residual).
            let mut mins = vec![0i32; size - 1];
            let mut maxs = vec![0i32; size - 1];
            let mut indices = vec![0i32; size - 1];
            for i in 0..size - 1 {
                let max = (means[i] - stds[i]).max(0);
                indices[i] = max;
                mins[i] = max;
                maxs[i] = means[i] + stds[i];
            }
            let max_value = means[size - 1] + stds[size - 1];

            multiple_fine_polynomial_recursive_helper(
                &mins,
                &maxs,
                &mut indices,
                0,
                &mut f_polynomial[k],
                composition,
                atoms,
                fine_min_prob,
                max_value,
            );
        }
    }

    let mut t_polynomial = f_polynomial[0].clone();

    if n <= 1 {
        return t_polynomial;
    }

    let mut fgid_polynomial: Vec<Polynomial> = Vec::new();
    for k in 1..n {
        multiply_fine_final_polynomial(
            &mut t_polynomial,
            &f_polynomial[k],
            &mut fgid_polynomial,
            fine_resolution,
            mw_resolution,
            fine_min_prob,
        );
    }

    t_polynomial
}

/// Port of `MultiplyFineFinalPolynomial`: convolve `t_polynomial` with `f_polynomial`, binning
/// the products onto the shared `fgid_polynomial` grid, then compact back into `t_polynomial`.
fn multiply_fine_final_polynomial(
    t_polynomial: &mut Vec<Polynomial>,
    f_polynomial: &[Polynomial],
    fgid_polynomial: &mut Vec<Polynomial>,
    fine_resolution: f64,
    mw_resolution: f64,
    fine_min_prob: f64,
) {
    let i_count = t_polynomial.len();
    let j_count = f_polynomial.len();
    if i_count == 0 || j_count == 0 {
        return;
    }

    let delta_mass = fine_resolution / mw_resolution;
    let min_probability = fine_min_prob;

    let min_mass = t_polynomial[0].power + f_polynomial[0].power;
    let max_mass = t_polynomial[i_count - 1].power + f_polynomial[j_count - 1].power;

    let max_index = ((max_mass - min_mass).abs() / delta_mass + 0.5) as i32 as usize;
    if max_index >= fgid_polynomial.len() {
        let extra = max_index - fgid_polynomial.len();
        for _ in 0..=extra {
            fgid_polynomial.push(Polynomial::NAN);
        }
    }

    for t in 0..t_polynomial.len() {
        for f in 0..f_polynomial.len() {
            let prob = t_polynomial[t].probability * f_polynomial[f].probability;
            if prob <= min_probability {
                continue;
            }

            let power = t_polynomial[t].power + f_polynomial[f].power;
            let indext = ((power - min_mass).abs() / delta_mass + 0.5) as i32 as usize;

            let temp = fgid_polynomial[indext];
            if temp.power.is_nan() || prob.is_nan() {
                fgid_polynomial[indext] = Polynomial {
                    power: power * prob,
                    probability: prob,
                };
            } else {
                fgid_polynomial[indext] = Polynomial {
                    power: temp.power + power * prob,
                    probability: temp.probability + prob,
                };
            }
        }
    }

    // Compact the populated grid bins back into t_polynomial, reusing existing slots first.
    let index = t_polynomial.len();
    let mut j = 0usize;
    for i in 0..fgid_polynomial.len() {
        if !fgid_polynomial[i].probability.is_nan() {
            let term = Polynomial {
                power: fgid_polynomial[i].power / fgid_polynomial[i].probability,
                probability: fgid_polynomial[i].probability,
            };
            if j < index {
                t_polynomial[j] = term;
                j += 1;
            } else {
                t_polynomial.push(term);
            }
        }
        fgid_polynomial[i] = Polynomial::NAN;
    }

    if j < index {
        t_polynomial.truncate(j);
    }
}

/// Port of `MultipleFinePolynomialRecursiveHelper`: enumerate the multinomial index grid for one
/// element, accumulating each combination's `(power, probability)` into `f_polynomial`.
#[allow(clippy::too_many_arguments)]
fn multiple_fine_polynomial_recursive_helper(
    mins: &[i32],
    maxs: &[i32],
    indices: &mut [i32],
    index: usize,
    f_polynomial: &mut Vec<Polynomial>,
    elemental_composition: &[Composition],
    atoms: i32,
    min_prob: f64,
    max_value: i32,
) {
    indices[index] = mins[index];
    while indices[index] <= maxs[index] {
        if index < mins.len() - 1 {
            multiple_fine_polynomial_recursive_helper(
                mins,
                maxs,
                indices,
                index + 1,
                f_polynomial,
                elemental_composition,
                atoms,
                min_prob,
                max_value,
            );
        } else {
            let l = atoms - indices.iter().sum::<i32>();
            if !(l < 0 || l > max_value) {
                let last = elemental_composition.len() - 1;
                let mut prob = factor_ln(atoms) - factor_ln(l)
                    + l as f64 * elemental_composition[last].log_probability;
                let mut power = l as f64 * elemental_composition[last].power;
                for i in 0..=last - 1 {
                    let index_value = indices[i];
                    let t_composition = &elemental_composition[i];
                    prob -= factor_ln(index_value);
                    prob += index_value as f64 * t_composition.log_probability;
                    power += index_value as f64 * t_composition.power;
                }

                let prob = prob.exp();
                if prob >= min_prob {
                    f_polynomial.push(Polynomial { probability: prob, power });
                }
            }
        }
        indices[index] += 1;
    }
}

/// Port of `MergeFinePolynomial`: 9 passes that coalesce adjacent terms whose masses fall
/// within a `k`-scaled threshold, keeping each merged bin's abundance-weighted centroid.
fn merge_fine_polynomial(
    t_polynomial: &mut Vec<Polynomial>,
    mw_resolution: f64,
    merge_fine_resolution: f64,
) -> Vec<Polynomial> {
    // Sort by mass (power); no NaN terms exist yet, so the comparator is total here.
    t_polynomial.sort_by(|a, b| {
        a.power
            .partial_cmp(&b.power)
            .expect("polynomial powers are finite before merging")
    });

    let count = t_polynomial.len();

    for k in 1..=9 {
        for i in 0..count {
            let power = t_polynomial[i].power;
            if power.is_nan() {
                continue;
            }

            let probability = t_polynomial[i].probability;
            // tempPolynomial accumulates the centroid numerator (power*prob) and total prob.
            let mut temp_power = power * probability;
            let mut temp_prob = probability;

            for j in (i + 1)..count {
                // Reads t_polynomial[i].power *fresh* each iteration: it is the running merged
                // centroid, updated in the branch below — this ordering matters for parity.
                let value = (t_polynomial[i].power * mw_resolution
                    - t_polynomial[j].power * mw_resolution)
                    .abs();
                let threshold = if k <= 8 {
                    k as f64 * merge_fine_resolution / 8.0
                } else {
                    merge_fine_resolution + merge_fine_resolution / 100.0
                };

                if value <= threshold {
                    temp_power += t_polynomial[j].power * t_polynomial[j].probability;
                    temp_prob += t_polynomial[j].probability;
                    t_polynomial[i] = Polynomial {
                        power: temp_power / temp_prob,
                        probability: temp_prob,
                    };
                    t_polynomial[j] = Polynomial::NAN;
                } else {
                    break;
                }
            }

            t_polynomial[i] = Polynomial {
                power: temp_power / temp_prob,
                probability: temp_prob,
            };
        }
    }

    // Return only the non-empty (non-NaN) terms.
    t_polynomial
        .iter()
        .filter(|p| !p.power.is_nan())
        .copied()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sum of intensities, and the (mass, intensity) of the most abundant peak.
    fn summarise(dist: &IsotopicDistribution) -> (f64, f64, f64) {
        let total: f64 = dist.intensities.iter().sum();
        let mut base_mass = 0.0;
        let mut base_int = f64::NEG_INFINITY;
        for (&m, &it) in dist.masses.iter().zip(dist.intensities.iter()) {
            if it > base_int {
                base_int = it;
                base_mass = m;
            }
        }
        (total, base_mass, base_int)
    }

    #[test]
    fn water_distribution_is_nonempty_and_monotone_in_mass() {
        let water = ChemicalFormula::parse_formula("H2O").unwrap();
        let dist = IsotopicDistribution::get_distribution(&water);

        assert_eq!(dist.masses.len(), dist.intensities.len());
        assert!(!dist.masses.is_empty(), "water must produce an envelope");

        // Merge sorts terms by mass; the emitted arrays preserve that ascending order.
        for w in dist.masses.windows(2) {
            assert!(w[0] <= w[1], "masses must be ascending: {} > {}", w[0], w[1]);
        }

        // The monoisotopic (lightest) peak of water should sit near 18.0106 Da and be the base.
        let (total, base_mass, base_int) = summarise(&dist);
        assert!(total > 0.0, "total intensity must be positive");
        assert!(
            (dist.masses[0] - 18.0106).abs() < 0.02,
            "lightest peak near 18.0106 Da, got {}",
            dist.masses[0]
        );
        assert!(
            (base_mass - 18.0106).abs() < 0.02,
            "base peak should be the monoisotopic peak for water, got {base_mass}"
        );
        assert!(base_int > 0.0);
    }

    #[test]
    fn carbon_isotopes_give_two_resolved_peaks() {
        // C100 has a strong M and M+1 (from ~1.1% ¹³C) separated by ~1 neutron mass.
        let c100 = ChemicalFormula::parse_formula("C100").unwrap();
        let dist = IsotopicDistribution::get_distribution(&c100);
        assert!(dist.masses.len() >= 2, "expected at least M and M+1 peaks");

        // First two peaks ~1 Da apart (the ¹³C−¹²C mass difference is ~1.0034 Da).
        let spacing = dist.masses[1] - dist.masses[0];
        assert!(
            (spacing - 1.0034).abs() < 0.05,
            "M+1 spacing should be ~1.0034 Da, got {spacing}"
        );

        // For 100 carbons the M+1 peak (~100 * 1.1%) should exceed the monoisotopic peak.
        assert!(
            dist.intensities[1] > dist.intensities[0],
            "M+1 should dominate M for C100: {} !> {}",
            dist.intensities[1],
            dist.intensities[0]
        );
    }

    #[test]
    fn isotope_label_adds_fixed_mass_offset() {
        // Replacing one of glucose's six carbons with an explicit ¹³C goes through the
        // `additionalMass` path: every peak shifts up by (¹³C − ¹²C) ≈ 1.00336 Da and the
        // peak count is unchanged (the labelled atom is no longer part of the polynomial).
        let light = ChemicalFormula::parse_formula("C6H12O6").unwrap();
        let labelled = ChemicalFormula::parse_formula("C5C{13}1H12O6").unwrap();

        let dl = IsotopicDistribution::get_distribution(&light);
        let dh = IsotopicDistribution::get_distribution(&labelled);

        assert!(!dl.masses.is_empty() && !dh.masses.is_empty());
        // The lightest peak is the all-light-isotope term: swapping one ¹²C for ¹³C lifts it
        // by exactly (¹³C − ¹²C) ≈ 1.00336 Da via the `additionalMass` offset.
        let shift = 13.00335483507 - 12.0;
        assert!(
            (dh.masses[0] - dl.masses[0] - shift).abs() < 1e-6,
            "lightest labelled peak {} should be {shift} heavier than {}",
            dh.masses[0],
            dl.masses[0]
        );
    }

    #[test]
    fn peptide_envelope_base_peak_is_near_monoisotopic_mass() {
        // A real peptide formula (PEPTIDE = C34H53N7O15): the base peak should be within a
        // couple of Da of the monoisotopic mass and the envelope should have several peaks.
        let peptide = ChemicalFormula::parse_formula("C34H53N7O15").unwrap();
        let dist = IsotopicDistribution::get_distribution(&peptide);
        assert!(dist.masses.len() >= 3, "peptide envelope should be multi-peak");
        let (_total, base_mass, _) = summarise(&dist);
        let mono = peptide.monoisotopic_mass();
        assert!(
            (base_mass - mono).abs() < 2.5,
            "base peak {base_mass} should be within ~2.5 Da of mono mass {mono}"
        );
    }

    #[test]
    fn factor_ln_matches_direct_log_sum() {
        // FactorLn(n) == sum_{k=2}^{n} ln(k), bit-identical to the cumulative cache.
        for n in [0, 1, 2, 5, 10, 50] {
            let mut expected = 0.0;
            for k in 2..=n {
                expected += (k as f64).ln();
            }
            assert_eq!(factor_ln(n), expected, "factor_ln({n})");
        }
    }
}
