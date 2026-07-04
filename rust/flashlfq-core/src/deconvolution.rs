//! Classic isotopic-envelope deconvolution — faithful port of mzLib's
//! `ClassicDeconvolutionAlgorithm` (plus the pieces it depends on: the `Averagine` model,
//! the deconvolution-flavoured `IsotopicEnvelope`, `ClassicDeconvolutionParameters`, and the
//! `Chemistry/ClassExtensions` helpers `ToMz`/`ToMass`/`GetClosestIndex`).
//!
//! This is a **parity port**: the C# control flow, comparison operators, integer division,
//! rounding modes, summation order and sort order are transcribed verbatim so a future golden
//! (C#-generated) parity harness can match bit-for-bit. It is deliberately *not* idiomatic
//! where idiom would change the arithmetic.
//!
//! ## Naming
//! The C# type produced by deconvolution is `MassSpectrometry.IsotopicEnvelope`. That name is
//! already taken in this crate by the FlashLFQ tracing envelope
//! ([`crate::isotopic_envelope::IsotopicEnvelope`], a *different* class), so the deconvolution
//! envelope is named [`DeconEnvelope`] here.
//!
//! ## Model
//! [`Averagine`] precomputes 1500 theoretical isotopic envelopes (one per half-averagine step),
//! matching `Averagine`'s static constructor. Each entry is sorted by intensity **descending**
//! (the C# does `Array.Sort(intensities, masses)` ascending then `Array.Reverse` on both). The
//! most-intense mass and its gap to the monoisotopic mass are cached per entry.
//!
//! ## Faithful-port judgment calls (flagged for a future golden reviewer)
//! - **Averagine intensity-tie order**: C# `Array.Sort` is an *unstable* introsort, so the order
//!   of equal-intensity peaks is unspecified. We sort a `(intensity, mass)` pair ascending with a
//!   stable sort then reverse (mirroring "sort then reverse"); equal-intensity ties therefore come
//!   out in reversed original order. A golden could disagree only on the ordering of exactly-equal
//!   intensities, which does not affect `masses[0]`/`DiffToMonoisotopic`.
//! - **`HashSet<int>` charge-state iteration order**: C# enumerates an add-only `HashSet<int>` in
//!   insertion order (int hashcode == value, slots fill sequentially). We replicate that with an
//!   insertion-ordered dedup `Vec<i32>`.
//! - **Sample standard deviation**: `List<double>.StandardDeviation()` dispatches to MathNet
//!   `StreamingStatistics.StandardDeviation` (sample, `n-1`; `NaN` for `n < 2`). The exact
//!   streaming recurrence is reproduced in [`standard_deviation`].
//! - **Median**: `.Median()` sorts and averages the two middle order statistics for even counts.
//! - **`HashSet<double>` final dedup**: C# compares raw `double` values by exact equality; we use
//!   `f64::to_bits()` in a `HashSet<u64>` (the compared m/z are literally the same array elements,
//!   so bitwise equality matches).

use std::cmp::Ordering;
use std::collections::HashSet;
use std::sync::LazyLock;

use crate::chemical_formula::ChemicalFormula;
use crate::isotopic_distribution::{IsotopicDistribution, DEFAULT_MOLECULAR_WEIGHT_RESOLUTION};
use crate::isotopic_envelope::{C13_MINUS_C12, PROTON_MASS};
use crate::periodic_table::periodic_table;
use crate::theoretical_isotope_distribution::{FINE_RESOLUTION, MIN_PROBABILITY};

// ---------------------------------------------------------------------------
// Chemistry/ClassExtensions helpers (ToMz / ToMass / GetClosestIndex)
// ---------------------------------------------------------------------------

/// `mass.ToMz(charge)` — m/z of a neutral mass, assuming charge from (de)protonation.
/// Mirrors `ClassExtensions.ToMz(this double mass, int charge)`:
/// `mass / |charge| + sign(charge) * ProtonMass`.
#[inline]
fn to_mz(mass: f64, charge: i32) -> f64 {
    mass / (charge.abs() as f64) + (charge.signum() as f64) * PROTON_MASS
}

/// `mz.ToMass(charge)` — neutral mass from an m/z, assuming charge from a proton.
/// Mirrors `ClassExtensions.ToMass(this double mz, int charge)`:
/// `|charge| * mz - charge * ProtonMass`.
#[inline]
fn to_mass(mz: f64, charge: i32) -> f64 {
    (charge.abs() as f64) * mz - (charge as f64) * PROTON_MASS
}

/// Which neighbour `get_closest_index` returns when the target is between two array elements.
/// Mirrors `Chemistry.ArraySearchOption`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ArraySearchOption {
    Previous,
    Closest,
    Next,
}

/// `Array.BinarySearch` semantics over an ascending `f64` slice: returns the index of an exact
/// match, or the bitwise complement (`!insertionPoint`) of the first element greater than
/// `target`. The midpoint is `lo + ((hi - lo) >> 1)`, matching .NET.
fn binary_search(arr: &[f64], target: f64) -> isize {
    let mut lo: isize = 0;
    let mut hi: isize = arr.len() as isize - 1;
    while lo <= hi {
        let mid = lo + ((hi - lo) >> 1);
        match arr[mid as usize].partial_cmp(&target).unwrap() {
            Ordering::Equal => return mid,
            Ordering::Less => lo = mid + 1,
            Ordering::Greater => hi = mid - 1,
        }
    }
    !lo
}

/// Faithful port of `ClassExtensions.GetClosestIndex`. `sorted_array` is ascending. Returns the
/// index of the array element closest to (or previous/next of) `target`. Empty arrays are not
/// expected here (the algorithm guards on `Size == 0` before any call).
fn get_closest_index(sorted_array: &[f64], target: f64, search_type: ArraySearchOption) -> usize {
    let idx = binary_search(sorted_array, target);
    if idx >= 0 {
        // exact match
        return idx as usize;
    }
    let idx = (!idx) as usize; // insertion point: first element greater than target
    if idx == sorted_array.len() {
        // target greater than the largest element
        return sorted_array.len() - 1;
    }
    if idx == 0 {
        // target less than the smallest element
        return 0;
    }
    match search_type {
        ArraySearchOption::Previous => idx - 1,
        ArraySearchOption::Closest => {
            let backwards_diff = target - sorted_array[idx - 1];
            let forwards_diff = sorted_array[idx] - target;
            if backwards_diff < forwards_diff {
                idx - 1
            } else {
                idx
            }
        }
        ArraySearchOption::Next => idx,
    }
}

/// `spectrum.GetClosestPeakIndex(x)` == `GetClosestIndex(x, Closest)` over the ascending m/z array.
#[inline]
fn get_closest_peak_index(mz: &[f64], target: f64) -> usize {
    get_closest_index(mz, target, ArraySearchOption::Closest)
}

// ---------------------------------------------------------------------------
// Statistics helpers (MathNet parity)
// ---------------------------------------------------------------------------

/// Sample standard deviation, reproducing MathNet `StreamingStatistics.StandardDeviation`
/// (`Statistics.StandardDeviation` on a `List<double>` dispatches here, since a `List<double>`
/// is not a `double[]`). Divides by `n - 1`; returns `NaN` for fewer than two samples.
fn standard_deviation(samples: &[f64]) -> f64 {
    let mut variance = 0.0_f64;
    let mut t = 0.0_f64;
    let mut j: u64 = 0;
    let mut iter = samples.iter();
    if let Some(&first) = iter.next() {
        j += 1;
        t = first;
    }
    for &xi in iter {
        j += 1;
        t += xi;
        let diff = (j as f64) * xi - t;
        variance += (diff * diff) / ((j as f64) * (j as f64 - 1.0));
    }
    if j > 1 {
        (variance / (j as f64 - 1.0)).sqrt()
    } else {
        f64::NAN
    }
}

/// MathNet `.Median()`: sort ascending, take the middle element for odd counts, or the average
/// of the two middle elements for even counts. Empty input returns `NaN` (not reached in the
/// algorithm, which always seeds at least one prediction).
fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Ionization polarity. Mirrors `MassSpectrometry.Polarity` for the two cases the algorithm uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Polarity {
    Positive,
    Negative,
}

impl Default for Polarity {
    fn default() -> Self {
        Polarity::Positive
    }
}

/// Parameters for classic deconvolution — the target (non-decoy) subset of
/// `ClassicDeconvolutionParameters` / `DeconvolutionParameters`. Decoy averagine, the generic
/// scorer, and equality/hashing are intentionally out of scope for this port.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClassicDeconvolutionParameters {
    pub min_assumed_charge_state: i32,
    pub max_assumed_charge_state: i32,
    pub deconvolution_tolerance_ppm: f64,
    pub intensity_ratio_limit: f64,
    pub polarity: Polarity,
    /// The expected spacing between isotope peaks (Da). Defaults to `C13MinusC12`. Retained for
    /// parity with the C# constructor's default even though the classic algorithm derives spacing
    /// from the averagine model rather than reading this field directly.
    pub expected_isotope_spacing: f64,
}

impl ClassicDeconvolutionParameters {
    /// Construct classic deconvolution parameters. Mirrors the C# constructor
    /// `ClassicDeconvolutionParameters(minCharge, maxCharge, deconPpm, intensityRatio, polarity,
    /// ..., expectedIsotopeSpacing = Constants.C13MinusC12)`.
    pub fn new(
        min_charge: i32,
        max_charge: i32,
        decon_ppm: f64,
        intensity_ratio: f64,
        polarity: Polarity,
    ) -> Self {
        ClassicDeconvolutionParameters {
            min_assumed_charge_state: min_charge,
            max_assumed_charge_state: max_charge,
            deconvolution_tolerance_ppm: decon_ppm,
            intensity_ratio_limit: intensity_ratio,
            polarity,
            expected_isotope_spacing: C13_MINUS_C12,
        }
    }
}

impl Default for ClassicDeconvolutionParameters {
    /// A convenient default for bottom-up work (min charge 1, max charge 6, 4 ppm,
    /// intensity ratio 3, positive). Tryptic peptides rarely exceed charge 6, so the
    /// charge range is capped there rather than at the top-down-oriented 60.
    fn default() -> Self {
        ClassicDeconvolutionParameters::new(1, 6, 4.0, 3.0, Polarity::Positive)
    }
}

// ---------------------------------------------------------------------------
// Deconvolution envelope (MassSpectrometry.IsotopicEnvelope, renamed DeconEnvelope)
// ---------------------------------------------------------------------------

/// A deconvoluted isotopic envelope. Port of `MassSpectrometry.IsotopicEnvelope` (the
/// deconvolution class, *not* the FlashLFQ tracing envelope of the same name).
#[derive(Debug, Clone, PartialEq)]
pub struct DeconEnvelope {
    /// Observed `(mz, intensity)` peaks, in the order they were matched (most intense first, then
    /// heavier isotopes as found).
    pub peaks: Vec<(f64, f64)>,
    /// Monoisotopic mass. Refined to the median of per-peak predictions when the envelope is
    /// accepted (see `set_median_monoisotopic_mass`).
    pub monoisotopic_mass: f64,
    /// Charge state.
    pub charge: i32,
    /// Summed intensity of the observed peaks.
    pub total_intensity: f64,
    /// Algorithm-specific score (empirical heuristic; may be aggregated across adjacent charge
    /// states via `aggregate_charge_state_score`).
    pub score: f64,
    /// Most abundant observed isotopic peak as `mz * |charge|` (NOT proton-corrected to a neutral
    /// mass). Mirrors `MostAbundantObservedIsotopicMass`.
    pub most_abundant_observed_isotopic_mass: f64,
    /// The per-peak monoisotopic-mass predictions gathered across this and adjacent charge states
    /// (the raw values the median in [`monoisotopic_mass`](Self::monoisotopic_mass) collapses);
    /// empty until the envelope is accepted.
    ///
    /// This is an **additive, parity-safe** extension over the C# `IsotopicEnvelope`: the charge-
    /// state consensus step in [`crate::feature_refinement`] needs the raw prediction list (the
    /// true monoisotopic mass recurs across charges while spurious ±1 Da predictions do not), which
    /// the median discards. It is populated at the acceptance point in [`classic_deconvolute`]
    /// (`set_median_monoisotopic_mass`) and does not affect any C#-matched field's value.
    pub monoisotopic_mass_predictions: Vec<f64>,
}

impl DeconEnvelope {
    /// Mirrors the C# constructor `IsotopicEnvelope(bestListOfPeaks, bestMonoisotopicMass,
    /// bestChargeState, bestTotalIntensity, bestStDev)`: computes `MostAbundantObservedIsotopicMass`
    /// from the most intense peak and the score from the ratio standard deviation.
    fn new(
        peaks: Vec<(f64, f64)>,
        monoisotopic_mass: f64,
        charge: i32,
        total_intensity: f64,
        st_dev: f64,
    ) -> Self {
        // MaxBy(p => p.intensity): first peak with the maximum intensity.
        let most_intense_mz = peaks
            .iter()
            .copied()
            .fold(None::<(f64, f64)>, |acc, p| match acc {
                Some(cur) if cur.1 >= p.1 => Some(cur),
                _ => Some(p),
            })
            .map(|p| p.0)
            .unwrap_or(0.0);
        let most_abundant_observed_isotopic_mass = most_intense_mz * (charge.abs() as f64);
        let score = score_isotope_envelope(peaks.len(), total_intensity, charge, st_dev);
        DeconEnvelope {
            peaks,
            monoisotopic_mass,
            charge,
            total_intensity,
            score,
            most_abundant_observed_isotopic_mass,
            monoisotopic_mass_predictions: Vec::new(),
        }
    }

    /// `AggregateChargeStateScore`: additively fold another charge state's score into this one.
    fn aggregate_charge_state_score(&mut self, other: &DeconEnvelope) {
        self.score += other.score;
    }

    /// `SetMedianMonoisotopicMass`: fine-tune the monoisotopic mass to the median of all per-peak
    /// predictions gathered across this and adjacent charge states.
    ///
    /// Additionally retains the raw prediction list on
    /// [`monoisotopic_mass_predictions`](Self::monoisotopic_mass_predictions) for the charge-state
    /// consensus step (parity-safe: `monoisotopic_mass` is set exactly as the C# does).
    fn set_median_monoisotopic_mass(&mut self, predictions: &[f64]) {
        self.monoisotopic_mass = median(predictions);
        self.monoisotopic_mass_predictions = predictions.to_vec();
    }
}

/// `ScoreIsotopeEnvelope`: `count >= 2 ? total / stDev^0.13 * count^0.4 / |charge|^0.06 : 0`.
fn score_isotope_envelope(peak_count: usize, total_intensity: f64, charge: i32, st_dev: f64) -> f64 {
    if peak_count >= 2 {
        total_intensity / st_dev.powf(0.13) * (peak_count as f64).powf(0.4)
            / ((charge.abs() as f64).powf(0.06))
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Averagine model
// ---------------------------------------------------------------------------

/// Number of averagine envelopes to precompute (`AverageResidue.NumAveraginesToGenerate`).
const NUM_AVERAGINES_TO_GENERATE: usize = 1500;

/// The precomputed Averagine model. Port of `MassSpectrometry.Averagine`'s static tables.
struct Averagine {
    all_masses: Vec<Vec<f64>>,
    all_intensities: Vec<Vec<f64>>,
    most_intense_masses: Vec<f64>,
    diff_to_monoisotopic: Vec<f64>,
    /// Monoisotopic mass of each entry (`most_intense − diff_to_monoisotopic`), ascending. Not part
    /// of the C# model; added so an arbitrary mono mass can be looked up (for the off-by-one
    /// corrector), which the intensity-keyed `most_intense_masses` table cannot do directly.
    monoisotopic_masses: Vec<f64>,
}

impl Averagine {
    /// Build the 1500 envelopes, mirroring the C# static constructor exactly.
    fn build() -> Self {
        // Magic numbers from https://pmc.ncbi.nlm.nih.gov/articles/PMC6166224/
        const AVERAGE_C: f64 = 4.9384;
        const AVERAGE_H: f64 = 7.7583;
        const AVERAGE_O: f64 = 1.4773;
        const AVERAGE_N: f64 = 1.3577;
        const AVERAGE_S: f64 = 0.0417;

        let pt = periodic_table();
        let an = |symbol: &str| {
            pt.element_by_symbol(symbol)
                .unwrap_or_else(|| panic!("averagine element {symbol} missing from periodic table"))
                .atomic_number
        };
        let (c, h, o, n, s) = (an("C"), an("H"), an("O"), an("N"), an("S"));

        let mut all_masses = Vec::with_capacity(NUM_AVERAGINES_TO_GENERATE);
        let mut all_intensities = Vec::with_capacity(NUM_AVERAGINES_TO_GENERATE);
        let mut most_intense_masses = vec![0.0; NUM_AVERAGINES_TO_GENERATE];
        let mut diff_to_monoisotopic = vec![0.0; NUM_AVERAGINES_TO_GENERATE];

        for i in 0..NUM_AVERAGINES_TO_GENERATE {
            let averagine_multiplier = (i as f64 + 1.0) / 2.0;
            let mut formula = ChemicalFormula::new();
            // C# `Convert.ToInt32` is round-half-to-even (banker's rounding).
            formula.add_element(c, (AVERAGE_C * averagine_multiplier).round_ties_even() as i32);
            formula.add_element(h, (AVERAGE_H * averagine_multiplier).round_ties_even() as i32);
            formula.add_element(o, (AVERAGE_O * averagine_multiplier).round_ties_even() as i32);
            formula.add_element(n, (AVERAGE_N * averagine_multiplier).round_ties_even() as i32);
            formula.add_element(s, (AVERAGE_S * averagine_multiplier).round_ties_even() as i32);

            let dist = IsotopicDistribution::get_distribution_with(
                &formula,
                FINE_RESOLUTION,
                MIN_PROBABILITY,
                DEFAULT_MOLECULAR_WEIGHT_RESOLUTION,
            );

            // C#: Array.Sort(intensities, masses) [ascending by intensity], then Reverse both →
            // sorted by intensity DESCENDING. Reproduce "sort ascending then reverse" with a
            // stable sort over (intensity, mass) pairs.
            let mut pairs: Vec<(f64, f64)> = dist
                .intensities
                .iter()
                .copied()
                .zip(dist.masses.iter().copied())
                .collect();
            pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            pairs.reverse();

            let masses: Vec<f64> = pairs.iter().map(|p| p.1).collect();
            let intensities: Vec<f64> = pairs.iter().map(|p| p.0).collect();

            most_intense_masses[i] = masses[0];
            diff_to_monoisotopic[i] = masses[0] - formula.monoisotopic_mass();
            all_masses.push(masses);
            all_intensities.push(intensities);
        }

        let monoisotopic_masses = most_intense_masses
            .iter()
            .zip(diff_to_monoisotopic.iter())
            .map(|(mi, d)| mi - d)
            .collect();

        Averagine {
            all_masses,
            all_intensities,
            most_intense_masses,
            diff_to_monoisotopic,
            monoisotopic_masses,
        }
    }

    /// `GetMostIntenseMassIndex`: closest index in the (ascending) most-intense-mass table.
    fn get_most_intense_mass_index(&self, test_mass: f64) -> usize {
        get_closest_index(&self.most_intense_masses, test_mass, ArraySearchOption::Closest)
    }

    fn get_all_theoretical_masses(&self, index: usize) -> &[f64] {
        &self.all_masses[index]
    }

    fn get_all_theoretical_intensities(&self, index: usize) -> &[f64] {
        &self.all_intensities[index]
    }

    fn get_diff_to_monoisotopic(&self, index: usize) -> f64 {
        self.diff_to_monoisotopic[index]
    }
}

/// Process-wide, lazily built averagine model (mirrors the C# static tables).
static AVERAGINE: LazyLock<Averagine> = LazyLock::new(Averagine::build);

/// Averagine isotope envelope as per-isotope **comb weights**, for the untargeted trace kernel's
/// `CombWeightModel::Averagine`. Given the neutral mass of the **most intense** isotope peak (the
/// trace-kernel seed), returns intensities indexed from the monoisotope (`k = 0`, then `k = 1` at
/// `+ (C13 − C12)`, …), summed per isotope index and normalized so the maximum weight is `1.0`.
///
/// This is the table-driven analogue of the closed-form Poisson comb
/// ([`crate::trace_kernel::poisson_comb_weights`]). It is more accurate near the mass where the
/// envelope mode shifts off the monoisotope (~1.8 kDa) — exactly where a Poisson `i*` can misplace
/// the monoisotope by one ¹³C unit — because it uses the real averagine peak intensities rather than
/// a single-parameter approximation. The vector is truncated once a post-mode weight falls below
/// `min_weight` (relative to the max), or at `max_isotopes` entries.
pub fn averagine_comb_weights(most_intense_mass: f64, min_weight: f64, max_isotopes: usize) -> Vec<f64> {
    let model = &*AVERAGINE;
    let idx = model.get_most_intense_mass_index(most_intense_mass);
    let mono = model.most_intense_masses[idx] - model.get_diff_to_monoisotopic(idx);
    averagine_envelope_by_index(idx, mono, min_weight, max_isotopes)
}

/// Averagine isotope envelope for a peptide of monoisotopic mass `mono_mass`, returned as per-isotope
/// intensities indexed from the monoisotope (`k = 0`), normalized so the maximum weight is `1.0`.
///
/// The mono-keyed counterpart of [`averagine_comb_weights`] (which is keyed by the *most-intense*
/// mass). Used by the off-by-one corrector to compare the observed envelope against the expected
/// averagine envelope anchored at several candidate monoisotopes.
pub fn averagine_intensities_from_mono(mono_mass: f64, min_weight: f64, max_isotopes: usize) -> Vec<f64> {
    let model = &*AVERAGINE;
    let idx = get_closest_index(&model.monoisotopic_masses, mono_mass, ArraySearchOption::Closest);
    let mono = model.monoisotopic_masses[idx];
    averagine_envelope_by_index(idx, mono, min_weight, max_isotopes)
}

/// Shared core of the two averagine-envelope accessors: bins averagine entry `idx`'s
/// intensity-sorted peaks back into isotope indices relative to `mono`, normalizes to max 1.0, and
/// trims the descending high-mass tail below `min_weight` (keeping every tooth up to and including
/// the mode — the low-mass teeth are what place the monoisotope).
fn averagine_envelope_by_index(idx: usize, mono: f64, min_weight: f64, max_isotopes: usize) -> Vec<f64> {
    let model = &*AVERAGINE;
    let masses = model.get_all_theoretical_masses(idx);
    let intensities = model.get_all_theoretical_intensities(idx);

    let mut weights: Vec<f64> = Vec::new();
    for (&m, &inten) in masses.iter().zip(intensities.iter()) {
        let k = ((m - mono) / C13_MINUS_C12).round();
        if k < 0.0 {
            continue;
        }
        let k = k as usize;
        if k >= weights.len() {
            weights.resize(k + 1, 0.0);
        }
        weights[k] += inten;
    }
    if weights.is_empty() {
        return weights;
    }

    let max_w = weights.iter().copied().fold(0.0_f64, f64::max);
    if max_w > 0.0 {
        for w in weights.iter_mut() {
            *w /= max_w;
        }
    }

    weights.truncate(max_isotopes);
    let mode = weights
        .iter()
        .enumerate()
        .fold(0usize, |best, (i, &w)| if w > weights[best] { i } else { best });
    let mut end = weights.len();
    while end > mode + 1 && weights[end - 1] < min_weight {
        end -= 1;
    }
    weights.truncate(end);
    weights
}

// ---------------------------------------------------------------------------
// The algorithm
// ---------------------------------------------------------------------------

/// Faithful port of `ClassicDeconvolutionAlgorithm.Deconvolute` over an ascending `(mz, intensity)`
/// spectrum. `range_min`/`range_max` bound the m/z window to deconvolute (the C# `MzRange`).
///
/// Returns the accepted [`DeconEnvelope`]s in the same order the C# iterator yields them
/// (descending score, first-claim-wins peak deduplication).
pub fn classic_deconvolute(
    mz: &[f64],
    intensity: &[f64],
    range_min: f64,
    range_max: f64,
    params: &ClassicDeconvolutionParameters,
) -> Vec<DeconEnvelope> {
    assert_eq!(mz.len(), intensity.len(), "mz and intensity must be parallel");

    // if no peaks, stop
    if mz.is_empty() {
        return Vec::new();
    }

    let model = &*AVERAGINE;
    let mut isolated_masses_and_charges: Vec<DeconEnvelope> = Vec::new();

    let (start, end) = extract_indices(mz, range_min, range_max);

    // find the most intense peak in the range
    let mut max_intensity = 0.0_f64;
    let mut index = start;
    while index < end {
        if intensity[index as usize] > max_intensity {
            max_intensity = intensity[index as usize];
        }
        index += 1;
    }

    // go through each peak in the range, assuming it is the most intense peak of its envelope
    let mut candidate = start;
    while candidate < end {
        let candidate_idx = candidate as usize;
        candidate += 1;

        let candidate_intensity = intensity[candidate_idx];
        // ignore peaks over 100x less intense than the most intense peak in the range
        if candidate_intensity * 100.0 < max_intensity {
            continue;
        }

        let mut best_envelope: Option<DeconEnvelope> = None;
        let candidate_mz = mz[candidate_idx];

        // find candidate charge states from the spacing of nearby higher-m/z peaks
        let mut all_possible_charge_states: Vec<i32> = Vec::new();
        for i in (candidate_idx + 1)..mz.len() {
            let delta_mass = mz[i] - candidate_mz;
            if delta_mass < 1.1 {
                // lower-bound charge state
                let mut charge = if params.polarity == Polarity::Negative {
                    (-1.0 / delta_mass).floor() as i32
                } else {
                    (1.0 / delta_mass).floor() as i32
                };
                if charge >= params.min_assumed_charge_state
                    && charge <= params.max_assumed_charge_state
                    && !all_possible_charge_states.contains(&charge)
                {
                    all_possible_charge_states.push(charge);
                }
                // upper-bound charge state
                charge += 1;
                if charge >= params.min_assumed_charge_state
                    && charge <= params.max_assumed_charge_state
                    && !all_possible_charge_states.contains(&charge)
                {
                    all_possible_charge_states.push(charge);
                }
            } else {
                break;
            }
        }

        // investigate the putative charge states
        for &charge_state in &all_possible_charge_states {
            // mass of this peak assuming this charge
            let test_most_intense_mass = to_mass(candidate_mz, charge_state);
            // averagine envelope index closest in mass
            let mass_index = model.get_most_intense_mass_index(test_most_intense_mass);

            let mut mono_predictions: Vec<f64> = Vec::new();

            // look for other isotopes at the assumed charge
            let putative = find_isotopic_envelope(
                mz,
                intensity,
                model,
                mass_index,
                candidate_mz,
                candidate_intensity,
                test_most_intense_mass,
                charge_state,
                params.deconvolution_tolerance_ppm,
                params.intensity_ratio_limit,
                &mut mono_predictions,
            );

            if putative.peaks.len() >= 2 {
                // NOTE: this mutates `putative`'s score (via aggregate), so it must run before the
                // score comparison below (matching the C# ordering).
                let mut putative = putative;
                let num_other_charge_states_observed = observe_adjacent_charge_states(
                    mz,
                    intensity,
                    model,
                    &mut putative,
                    candidate_mz,
                    mass_index,
                    params.deconvolution_tolerance_ppm,
                    params.intensity_ratio_limit,
                    params.min_assumed_charge_state,
                    params.max_assumed_charge_state,
                    &mut mono_predictions,
                );

                let better = match &best_envelope {
                    None => true,
                    Some(best) => putative.score > best.score,
                };
                // INTEGER division: charge / 5 in i32
                if better && (putative.charge / 5 <= num_other_charge_states_observed) {
                    putative.set_median_monoisotopic_mass(&mono_predictions);
                    best_envelope = Some(putative);
                }
            }
        }

        if let Some(best) = best_envelope {
            isolated_masses_and_charges.push(best);
        }
    }

    // Final dedup: order by score descending (stable), first-claim-wins over peak m/z.
    isolated_masses_and_charges
        .sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal));

    let mut seen: HashSet<u64> = HashSet::new();
    let mut result: Vec<DeconEnvelope> = Vec::new();
    for ok in isolated_masses_and_charges {
        // seen.Overlaps(ok.Peaks.Select(mz))
        if ok.peaks.iter().any(|p| seen.contains(&p.0.to_bits())) {
            continue;
        }
        for p in &ok.peaks {
            seen.insert(p.0.to_bits());
        }
        result.push(ok);
    }
    result
}

/// `FindIsotopicEnvelope`: starting from the assumed-most-intense peak, walk the averagine
/// envelope's remaining theoretical peaks (most→least intense) and collect experimental matches
/// that pass the ppm tolerance, the intensity-ratio test, and the not-already-observed test.
#[allow(clippy::too_many_arguments)]
fn find_isotopic_envelope(
    mz: &[f64],
    intensity: &[f64],
    model: &Averagine,
    mass_index: usize,
    candidate_mz: f64,
    candidate_intensity: f64,
    test_most_intense_mass: f64,
    charge_state: i32,
    deconvolution_tolerance_ppm: f64,
    intensity_ratio_limit: f64,
    mono_predictions: &mut Vec<f64>,
) -> DeconEnvelope {
    let theoretical_masses = model.get_all_theoretical_masses(mass_index);
    let theoretical_intensities = model.get_all_theoretical_intensities(mass_index);

    // add "most intense peak"
    let mut observed_peaks: Vec<(f64, f64)> = vec![(candidate_mz, candidate_intensity)];
    let mut ratios: Vec<f64> = vec![theoretical_intensities[0] / candidate_intensity];

    // mass difference (actual - theoretical) for the tallest peak
    let difference_between_theor_and_actual_mass = test_most_intense_mass - theoretical_masses[0];
    let mut total_intensity = candidate_intensity;
    let monoisotopic_mass = test_most_intense_mass - model.get_diff_to_monoisotopic(mass_index);
    mono_predictions.push(monoisotopic_mass);

    for index_to_look_at in 1..theoretical_intensities.len() {
        let theor_mass_that_trying_to_find =
            theoretical_masses[index_to_look_at] + difference_between_theor_and_actual_mass;
        let closest_peak_index =
            get_closest_peak_index(mz, to_mz(theor_mass_that_trying_to_find, charge_state));
        let closest_peak_mz = mz[closest_peak_index];
        let closest_peak_intensity = intensity[closest_peak_index];
        let closest_peak_mass = to_mass(closest_peak_mz, charge_state);

        if (closest_peak_mass - theor_mass_that_trying_to_find).abs()
            / theor_mass_that_trying_to_find
            * 1e6
            <= deconvolution_tolerance_ppm
            && peak2_satisfies_ratio(
                theoretical_intensities[0],
                theoretical_intensities[index_to_look_at],
                candidate_intensity,
                closest_peak_intensity,
                intensity_ratio_limit,
            )
            && !observed_peaks.contains(&(closest_peak_mz, closest_peak_intensity))
        {
            observed_peaks.push((closest_peak_mz, closest_peak_intensity));
            total_intensity += closest_peak_intensity;
            ratios.push(theoretical_intensities[index_to_look_at] / closest_peak_intensity);
            let monoisotopic_mass_from_this_peak =
                monoisotopic_mass + closest_peak_mass - theor_mass_that_trying_to_find;
            mono_predictions.push(monoisotopic_mass_from_this_peak);
        } else {
            break;
        }
    }

    DeconEnvelope::new(
        observed_peaks,
        monoisotopic_mass,
        charge_state,
        total_intensity,
        standard_deviation(&ratios),
    )
}

/// `ObserveAdjacentChargeStates`: look for the proposed neutral mass at lower and higher charge
/// states (stopping at the first miss in each direction), aggregating their scores into
/// `original_envelope` and collecting their monoisotopic-mass predictions.
#[allow(clippy::too_many_arguments)]
fn observe_adjacent_charge_states(
    mz: &[f64],
    intensity: &[f64],
    model: &Averagine,
    original_envelope: &mut DeconEnvelope,
    most_intense_peak_mz: f64,
    mass_index: usize,
    deconvolution_tolerance_ppm: f64,
    intensity_ratio_limit: f64,
    min_charge_to_look_for: i32,
    max_charge_to_look_for: i32,
    mono_predictions: &mut Vec<f64>,
) -> i32 {
    let mut num_adjacent_charge_states_observed = 0;
    let original_z = original_envelope.charge;
    let most_abundant_neutral_isotope = to_mass(most_intense_peak_mz, original_z);

    // lower charge states
    let mut lower_z = original_z - 1;
    while lower_z >= min_charge_to_look_for {
        if find_charge_state_of_mass(
            mz,
            intensity,
            model,
            original_envelope,
            lower_z,
            most_abundant_neutral_isotope,
            mass_index,
            deconvolution_tolerance_ppm,
            intensity_ratio_limit,
            mono_predictions,
        ) {
            num_adjacent_charge_states_observed += 1;
        } else {
            break;
        }
        lower_z -= 1;
    }

    // higher charge states
    let mut higher_z = original_z + 1;
    while higher_z <= max_charge_to_look_for {
        if find_charge_state_of_mass(
            mz,
            intensity,
            model,
            original_envelope,
            higher_z,
            most_abundant_neutral_isotope,
            mass_index,
            deconvolution_tolerance_ppm,
            intensity_ratio_limit,
            mono_predictions,
        ) {
            num_adjacent_charge_states_observed += 1;
        } else {
            break;
        }
        higher_z += 1;
    }

    num_adjacent_charge_states_observed
}

/// `FindChargeStateOfMass`: test whether the given neutral mass appears at `z_to_investigate`. On
/// a successful envelope (score != 0), aggregate its score into `original_envelope` and keep its
/// predictions; on failure, roll back the single prediction the inner `find_isotopic_envelope`
/// seeded.
#[allow(clippy::too_many_arguments)]
fn find_charge_state_of_mass(
    mz: &[f64],
    intensity: &[f64],
    model: &Averagine,
    original_envelope: &mut DeconEnvelope,
    z_to_investigate: i32,
    most_abundant_neutral_isotope_to_investigate: f64,
    mass_index: usize,
    deconvolution_tolerance_ppm: f64,
    intensity_ratio_limit: f64,
    mono_predictions: &mut Vec<f64>,
) -> bool {
    let most_abundant_isotope_mz_theoretical =
        to_mz(most_abundant_neutral_isotope_to_investigate, z_to_investigate);
    let observed_peak_index = get_closest_peak_index(mz, most_abundant_isotope_mz_theoretical);

    let most_abundant_isotope_mz_observed = mz[observed_peak_index];
    let most_abundant_isotope_mass_observed =
        to_mass(most_abundant_isotope_mz_observed, z_to_investigate);

    if (most_abundant_isotope_mass_observed - most_abundant_neutral_isotope_to_investigate).abs()
        / most_abundant_neutral_isotope_to_investigate
        * 1e6
        <= deconvolution_tolerance_ppm
    {
        let test = find_isotopic_envelope(
            mz,
            intensity,
            model,
            mass_index,
            most_abundant_isotope_mz_observed,
            intensity[observed_peak_index],
            most_abundant_isotope_mass_observed,
            z_to_investigate,
            deconvolution_tolerance_ppm,
            intensity_ratio_limit,
            mono_predictions,
        );

        if test.score != 0.0 {
            original_envelope.aggregate_charge_state_score(&test);
            true
        } else {
            // remove the test mass from this failed isotopic envelope
            mono_predictions.pop();
            false
        }
    } else {
        false
    }
}

/// `ExtractIndices`: `(1, 0)` (an empty range) when the window misses the spectrum entirely,
/// otherwise `(GetClosestIndex(min, Next), GetClosestIndex(max, Previous))`. Returned as `i64` so
/// the empty case (`start > end`) yields a non-iterating range.
fn extract_indices(mz: &[f64], min_x: f64, max_x: f64) -> (i64, i64) {
    if *mz.last().unwrap() < min_x || *mz.first().unwrap() > max_x {
        return (1, 0);
    }
    (
        get_closest_index(mz, min_x, ArraySearchOption::Next) as i64,
        get_closest_index(mz, max_x, ArraySearchOption::Previous) as i64,
    )
}

/// `Peak2satisfiesRatio`: the second peak's observed intensity must be within `intensity_ratio` of
/// the intensity implied by scaling peak 1's observed/theoretical ratio onto peak 2's theoretical
/// intensity.
fn peak2_satisfies_ratio(
    peak1_theor_intensity: f64,
    peak2_theor_intensity: f64,
    peak1_intensity: f64,
    peak2_intensity: f64,
    intensity_ratio: f64,
) -> bool {
    let compared_should_be = peak1_intensity / peak1_theor_intensity * peak2_theor_intensity;
    !(peak2_intensity < compared_should_be / intensity_ratio
        || peak2_intensity > compared_should_be * intensity_ratio)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic spectrum from the averagine envelope at `mass_index`, placing each
    /// theoretical isotope peak at `theorMass.to_mz(charge)` with the theoretical intensity, and
    /// returning ascending-by-m/z parallel arrays plus the true monoisotopic mass.
    fn synthetic_spectrum(mass_index: usize, charge: i32) -> (Vec<f64>, Vec<f64>, f64) {
        let model = &*AVERAGINE;
        let masses = model.get_all_theoretical_masses(mass_index);
        let intensities = model.get_all_theoretical_intensities(mass_index);

        let mut peaks: Vec<(f64, f64)> = masses
            .iter()
            .zip(intensities.iter())
            .map(|(&m, &int)| (to_mz(m, charge), int))
            .collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let mz: Vec<f64> = peaks.iter().map(|p| p.0).collect();
        let inten: Vec<f64> = peaks.iter().map(|p| p.1).collect();
        // truth mono = most_intense_mass - diff_to_monoisotopic == formula.MonoisotopicMass
        let truth_mono =
            model.most_intense_masses[mass_index] - model.diff_to_monoisotopic[mass_index];
        (mz, inten, truth_mono)
    }

    fn ppm_diff(a: f64, b: f64) -> f64 {
        (a - b).abs() / b * 1e6
    }

    #[test]
    fn recovers_single_envelope_charge_and_mass() {
        let model = &*AVERAGINE;
        let charge = 2;
        // pick the averagine entry closest to ~2000 Da
        let mass_index = model.get_most_intense_mass_index(2000.0);
        let (mz, inten, truth_mono) = synthetic_spectrum(mass_index, charge);

        let params = ClassicDeconvolutionParameters::new(1, 10, 20.0, 3.0, Polarity::Positive);
        let range_min = *mz.first().unwrap();
        let range_max = *mz.last().unwrap();
        let envelopes = classic_deconvolute(&mz, &inten, range_min, range_max, &params);

        assert!(!envelopes.is_empty(), "should recover at least one envelope");
        let found = envelopes
            .iter()
            .find(|e| e.charge == charge)
            .expect("an envelope at the true charge state should be recovered");
        assert!(
            ppm_diff(found.monoisotopic_mass, truth_mono) <= 20.0,
            "mono mass {} not within 20 ppm of truth {} (ppm={})",
            found.monoisotopic_mass,
            truth_mono,
            ppm_diff(found.monoisotopic_mass, truth_mono)
        );
        assert!(found.peaks.len() >= 2, "envelope should have >= 2 peaks");
        assert!(found.score > 0.0, "envelope score should be positive");
    }

    #[test]
    fn recovers_charge_three_envelope() {
        let model = &*AVERAGINE;
        let charge = 3;
        let mass_index = model.get_most_intense_mass_index(3000.0);
        let (mz, inten, truth_mono) = synthetic_spectrum(mass_index, charge);

        let params = ClassicDeconvolutionParameters::new(1, 10, 20.0, 3.0, Polarity::Positive);
        let range_min = *mz.first().unwrap();
        let range_max = *mz.last().unwrap();
        let envelopes = classic_deconvolute(&mz, &inten, range_min, range_max, &params);

        let found = envelopes
            .iter()
            .find(|e| e.charge == charge)
            .expect("an envelope at charge 3 should be recovered");
        assert!(ppm_diff(found.monoisotopic_mass, truth_mono) <= 20.0);
    }

    #[test]
    fn two_charge_states_exercise_adjacent_observation() {
        // Same neutral mass present at charge 2 and charge 3 in one spectrum.
        let model = &*AVERAGINE;
        let mass_index = model.get_most_intense_mass_index(2500.0);
        let truth_mono =
            model.most_intense_masses[mass_index] - model.diff_to_monoisotopic[mass_index];

        let masses = model.get_all_theoretical_masses(mass_index);
        let intensities = model.get_all_theoretical_intensities(mass_index);

        let mut peaks: Vec<(f64, f64)> = Vec::new();
        for charge in [2, 3] {
            for (&m, &int) in masses.iter().zip(intensities.iter()) {
                peaks.push((to_mz(m, charge), int));
            }
        }
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mz: Vec<f64> = peaks.iter().map(|p| p.0).collect();
        let inten: Vec<f64> = peaks.iter().map(|p| p.1).collect();

        let params = ClassicDeconvolutionParameters::new(1, 10, 20.0, 3.0, Polarity::Positive);
        let range_min = *mz.first().unwrap();
        let range_max = *mz.last().unwrap();
        let envelopes = classic_deconvolute(&mz, &inten, range_min, range_max, &params);

        assert!(
            !envelopes.is_empty(),
            "should recover envelope(s) from the two-charge-state spectrum"
        );
        // At least one recovered envelope should match the true neutral mass at one of the charges.
        assert!(
            envelopes.iter().any(|e| {
                (e.charge == 2 || e.charge == 3) && ppm_diff(e.monoisotopic_mass, truth_mono) <= 20.0
            }),
            "expected an envelope matching the true mass at charge 2 or 3"
        );
    }

    #[test]
    fn accepted_envelope_surfaces_monoisotopic_predictions() {
        // The additive `monoisotopic_mass_predictions` list must be populated on acceptance and be
        // bit-equal to the median that `monoisotopic_mass` collapses to (that is exactly how it is
        // set), and must not perturb any C#-matched field. Reuse the single-envelope fixture.
        let model = &*AVERAGINE;
        let charge = 2;
        let mass_index = model.get_most_intense_mass_index(2000.0);
        let (mz, inten, _truth_mono) = synthetic_spectrum(mass_index, charge);

        let params = ClassicDeconvolutionParameters::new(1, 10, 20.0, 3.0, Polarity::Positive);
        let range_min = *mz.first().unwrap();
        let range_max = *mz.last().unwrap();
        let envelopes = classic_deconvolute(&mz, &inten, range_min, range_max, &params);

        let found = envelopes
            .iter()
            .find(|e| e.charge == charge)
            .expect("an envelope at the true charge state should be recovered");
        assert!(
            !found.monoisotopic_mass_predictions.is_empty(),
            "accepted envelope should surface its raw monoisotopic predictions"
        );
        // median(predictions) is exactly how monoisotopic_mass was set → bit-equal.
        assert_eq!(
            median(&found.monoisotopic_mass_predictions),
            found.monoisotopic_mass
        );
    }

    #[test]
    fn averagine_comb_weights_normalize_and_shift_mode() {
        // Light peptide: the monoisotope is the tallest tooth (mode at k=0); max normalized to 1.
        let light = averagine_comb_weights(1000.0, 1e-3, 20);
        assert!(!light.is_empty());
        let max_l = light.iter().copied().fold(0.0_f64, f64::max);
        assert!((max_l - 1.0).abs() < 1e-12, "weights normalized to max 1");
        let mode_l = light
            .iter()
            .enumerate()
            .fold(0usize, |b, (i, &w)| if w > light[b] { i } else { b });
        assert_eq!(mode_l, 0, "light mass: monoisotope is the tallest tooth");

        // Heavy peptide: the envelope mode shifts *off* the monoisotope (k > 0) — precisely the
        // off-by-one regime averagine models better than a single-parameter Poisson.
        let heavy = averagine_comb_weights(3500.0, 1e-3, 30);
        let mode_h = heavy
            .iter()
            .enumerate()
            .fold(0usize, |b, (i, &w)| if w > heavy[b] { i } else { b });
        assert!(mode_h >= 1, "heavy mass: mode above the monoisotope, got {mode_h}");
        assert!(heavy[0] > 0.0, "monoisotope tooth retained (needed to place the mono)");
    }

    #[test]
    fn empty_spectrum_returns_empty() {
        let params = ClassicDeconvolutionParameters::default();
        let out = classic_deconvolute(&[], &[], 0.0, 5000.0, &params);
        assert!(out.is_empty());
    }

    #[test]
    fn single_lone_peak_yields_nothing() {
        // One peak: no higher-m/z neighbour within 1.1 Th, so no charge states → no envelopes.
        let mz = vec![1001.007];
        let inten = vec![1.0e6];
        let params = ClassicDeconvolutionParameters::new(1, 10, 20.0, 3.0, Polarity::Positive);
        let out = classic_deconvolute(&mz, &inten, 0.0, 5000.0, &params);
        assert!(out.is_empty(), "a lone peak should not yield an envelope");
    }

    #[test]
    fn two_distant_peaks_yield_nothing() {
        // Two peaks more than 1.1 Th apart: the spacing loop breaks immediately, no charge states.
        let mz = vec![1000.0, 1005.0];
        let inten = vec![1.0e6, 9.0e5];
        let params = ClassicDeconvolutionParameters::new(1, 10, 20.0, 3.0, Polarity::Positive);
        let out = classic_deconvolute(&mz, &inten, 0.0, 5000.0, &params);
        assert!(out.is_empty(), "two distant peaks should not form an envelope");
    }

    #[test]
    fn get_closest_index_matches_binary_search_semantics() {
        let arr = [1.0, 3.0, 5.0, 7.0, 9.0];
        // exact match
        assert_eq!(get_closest_index(&arr, 5.0, ArraySearchOption::Closest), 2);
        // below the minimum → 0
        assert_eq!(get_closest_index(&arr, -10.0, ArraySearchOption::Previous), 0);
        // above the maximum → last
        assert_eq!(get_closest_index(&arr, 100.0, ArraySearchOption::Next), 4);
        // between 3 and 5, closer to 3
        assert_eq!(get_closest_index(&arr, 3.4, ArraySearchOption::Closest), 1);
        // between 3 and 5, closer to 5
        assert_eq!(get_closest_index(&arr, 4.6, ArraySearchOption::Closest), 2);
        // Previous / Next around 4.0 (equidistant): Closest tie → forward (idx), Previous → idx-1
        assert_eq!(get_closest_index(&arr, 4.0, ArraySearchOption::Previous), 1);
        assert_eq!(get_closest_index(&arr, 4.0, ArraySearchOption::Next), 2);
        assert_eq!(get_closest_index(&arr, 4.0, ArraySearchOption::Closest), 2);
    }

    #[test]
    fn standard_deviation_is_sample_and_nan_for_singletons() {
        assert!(standard_deviation(&[42.0]).is_nan());
        // sample stddev of [2,4,4,4,5,5,7,9] is 2.13809...
        let s = standard_deviation(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]);
        assert!((s - 2.138089935299395).abs() < 1e-12, "got {s}");
    }

    #[test]
    fn median_odd_and_even() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
    }

    #[test]
    fn to_mz_to_mass_round_trip() {
        let mass = 1998.6;
        let mz = to_mz(mass, 3);
        assert!((to_mass(mz, 3) - mass).abs() < 1e-9);
    }
}
