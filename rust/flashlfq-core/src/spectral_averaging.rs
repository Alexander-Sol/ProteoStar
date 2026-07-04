//! Spectral averaging — faithful port of mzLib's `SpectralAveraging` project, **default-config
//! subset**.
//!
//! Purpose in the untargeted feature-detection pipeline: once the trace kernel has located a
//! candidate feature and derived an RT window, we average that window's MS1 scans into a single
//! composite spectrum (SNR ~√N) before the final [`crate::deconvolution`] pass. The **m/z binning**
//! is the load-bearing step — it is what gives the composite its mass accuracy — so it is ported
//! verbatim (see `agent_info/Feature-Detection-Design.md`, "Window + averaging").
//!
//! ## Scope (default-config subset)
//! mzLib's full project supports several outlier-rejection algorithms, several weighting schemes,
//! and file-level scan windowing. This port covers the **default** configuration faithfully and
//! stubs the rest, per the design decision to "port binning + the default rejection/weighting
//! config faithfully (golden-gated), stub the alternate rejection algorithms until needed":
//!
//! | Knob | mzLib default | Ported here |
//! |------|---------------|-------------|
//! | `SpectralAveragingType` | `MzBinning` (only variant) | ✅ full |
//! | `NormalizationType` | `RelativeToTics` | ✅ all four variants (trivial) |
//! | `SpectraWeightingType` | `WeightEvenly` | ✅ `WeightEvenly` + `TicValue`; `MrsNoiseEstimation` panics |
//! | `OutlierRejectionType` | `NoRejection` | ✅ `NoRejection`; the six clipping variants panic |
//! | `BinSize` | `0.01` | ✅ |
//!
//! The file-level windowing (`SpectraFileAveraging`, `AverageEverynScansWithOverlap`, scan overlap,
//! output type, thread count) is **not** ported: feature detection supplies its own scan window, so
//! the relevant surface is just `AverageSpectra(double[][] xArrays, double[][] yArrays, params)`.
//!
//! ## Parity
//! Control flow, summation order, the `floor((x - minX) / binSize)` bin index, the zero-padding of
//! bins missing a spectrum, and the "divide the bin's summed intensity by the *spectrum* count (not
//! the present-peak count)" behaviour are transcribed verbatim so a future C# golden matches at the
//! standard tolerance (counts exact, floats rel-1e-6). Where mzLib mutates the caller's `yArrays`
//! in place during normalization, this port normalizes an internal clone instead — the arithmetic
//! is identical; only the (undesirable) side effect on the caller is dropped.

// ---------------------------------------------------------------------------
// Configuration enums (mirroring SpectralAveraging/DataStructures/Enums)
// ---------------------------------------------------------------------------

/// `SpectralAveraging.OutlierRejectionType`. Only [`OutlierRejectionType::NoRejection`] is ported;
/// the clipping variants are carried for config/parity completeness but panic if dispatched
/// (see module scope).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutlierRejectionType {
    NoRejection,
    MinMaxClipping,
    PercentileClipping,
    SigmaClipping,
    WinsorizedSigmaClipping,
    AveragedSigmaClipping,
    BelowThresholdRejection,
}

/// `SpectralAveraging.SpectraWeightingType`. [`SpectraWeightingType::WeightEvenly`] and
/// [`SpectraWeightingType::TicValue`] are ported; `MrsNoiseEstimation` panics (it needs the
/// MRS noise estimator + biweight midvariance, out of the default-config subset).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectraWeightingType {
    WeightEvenly,
    MrsNoiseEstimation,
    TicValue,
}

/// `SpectralAveraging.NormalizationType`. All four variants are ported (they are trivial).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalizationType {
    NoNormalization,
    RelativeToTics,
    AbsoluteToTic,
    RelativeIntensity,
}

/// `SpectralAveraging.SpectralAveragingType`. mzLib defines only one variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectralAveragingType {
    MzBinning,
}

// ---------------------------------------------------------------------------
// Parameters (mirroring SpectralAveragingParameters, averaging-only subset)
// ---------------------------------------------------------------------------

/// The subset of `SpectralAveragingParameters` that governs [`average_spectra`]. The file-level
/// windowing fields (`SpectraFileAveragingType`, `NumberOfScansToAverage`, `ScanOverlap`,
/// `OutputType`, `MaxThreadsToUsePerFile`) are intentionally omitted — feature detection supplies
/// its own scan window. `MinSigmaValue`/`MaxSigmaValue`/`Percentile` are carried so the struct can
/// still describe the (currently panicking) clipping configs without changing shape later.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralAveragingParameters {
    pub outlier_rejection_type: OutlierRejectionType,
    pub spectral_weighting_type: SpectraWeightingType,
    pub spectral_averaging_type: SpectralAveragingType,
    pub normalization_type: NormalizationType,
    pub bin_size: f64,
    pub percentile: f64,
    pub min_sigma_value: f64,
    pub max_sigma_value: f64,
}

impl Default for SpectralAveragingParameters {
    /// Mirrors `SpectralAveragingParameters.SetDefaultValues()` (averaging-relevant fields only):
    /// `NoRejection`, `WeightEvenly`, `MzBinning`, `RelativeToTics`, bin size `0.01`.
    fn default() -> Self {
        SpectralAveragingParameters {
            outlier_rejection_type: OutlierRejectionType::NoRejection,
            spectral_weighting_type: SpectraWeightingType::WeightEvenly,
            spectral_averaging_type: SpectralAveragingType::MzBinning,
            normalization_type: NormalizationType::RelativeToTics,
            bin_size: 0.01,
            percentile: 0.1,
            min_sigma_value: 1.5,
            max_sigma_value: 1.5,
        }
    }
}

// ---------------------------------------------------------------------------
// BinnedPeak (mirroring DataStructures/BinnedPeak.cs)
// ---------------------------------------------------------------------------

/// One peak assigned to an m/z bin, tagged with the spectrum it came from. Mirrors the internal
/// `BinnedPeak` record struct. Zero-intensity instances are synthesized to pad bins that a given
/// spectrum did not contribute a real peak to (see [`get_bins`]).
#[derive(Debug, Clone, Copy)]
struct BinnedPeak {
    mz: f64,
    intensity: f64,
    spectra_id: usize,
}

// ---------------------------------------------------------------------------
// Public entry point (mirroring SpectraAveraging.AverageSpectra)
// ---------------------------------------------------------------------------

/// Averages a group of spectra into one composite `(mz, intensity)` spectrum.
///
/// Faithful port of `SpectraAveraging.AverageSpectra` → `MzBinning`. `x_arrays[i]` / `y_arrays[i]`
/// are the m/z and intensity arrays of the *i*-th spectrum (each internally ascending in m/z, as
/// mzML centroided scans are). Returns `(mz, intensity)` for the averaged spectrum, ascending in
/// m/z, with zero-intensity bins dropped.
///
/// Panics if `x_arrays`/`y_arrays` differ in outer length, if any paired inner arrays differ in
/// length, or if the input is empty — matching the fact that the C# would throw / divide-by-zero
/// on those degenerate inputs rather than return a meaningful spectrum.
pub fn average_spectra(
    x_arrays: &[Vec<f64>],
    y_arrays: &[Vec<f64>],
    parameters: &SpectralAveragingParameters,
) -> (Vec<f64>, Vec<f64>) {
    match parameters.spectral_averaging_type {
        SpectralAveragingType::MzBinning => mz_binning(x_arrays, y_arrays, parameters),
    }
}

/// Faithful port of `SpectraAveraging.MzBinning`.
fn mz_binning(
    x_arrays: &[Vec<f64>],
    y_arrays: &[Vec<f64>],
    parameters: &SpectralAveragingParameters,
) -> (Vec<f64>, Vec<f64>) {
    assert_eq!(
        x_arrays.len(),
        y_arrays.len(),
        "x and y arrays must have the same number of spectra"
    );
    assert!(!x_arrays.is_empty(), "cannot average zero spectra");
    for (x, y) in x_arrays.iter().zip(y_arrays.iter()) {
        assert_eq!(x.len(), y.len(), "each spectrum's x and y arrays must match in length");
    }

    // get tics (from the *original* intensities, before normalization mutates them)
    let tics: Vec<f64> = y_arrays.iter().map(|y| sum(y)).collect();
    let average_tic = sum(&tics) / tics.len() as f64;

    // normalize spectra — mzLib mutates the caller's arrays here; we normalize an internal clone.
    let mut y_norm: Vec<Vec<f64>> = y_arrays.to_vec();
    normalize_spectra(&mut y_norm, parameters.normalization_type);

    // get bins
    let bins = get_bins(x_arrays, &y_norm, parameters.bin_size);

    // get weights
    let weights = calculate_spectra_weights(x_arrays, &y_norm, parameters.spectral_weighting_type);

    // reject outliers and average bins
    let mut averaged_peaks: Vec<(f64, f64)> = Vec::with_capacity(bins.len());
    for peaks_from_bin in &bins {
        let kept = reject_outliers(peaks_from_bin, parameters);
        if kept.is_empty() {
            continue;
        }
        averaged_peaks.push(average_bin(&kept, &weights));
    }

    // return averaged: drop zero-intensity bins, order by m/z; AbsoluteToTic re-scales by averageTic
    let mut ordered: Vec<(f64, f64)> =
        averaged_peaks.into_iter().filter(|p| p.1 != 0.0).collect();
    ordered.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mzs: Vec<f64> = ordered.iter().map(|p| p.0).collect();
    let intensities: Vec<f64> = ordered
        .iter()
        .map(|p| match parameters.normalization_type {
            NormalizationType::AbsoluteToTic => p.1 * average_tic,
            _ => p.1,
        })
        .collect();
    (mzs, intensities)
}

// ---------------------------------------------------------------------------
// Binning (mirroring SpectraAveraging.GetBins / AverageBin)
// ---------------------------------------------------------------------------

/// Faithful port of `SpectraAveraging.GetBins`. Sorts every peak of every spectrum into an m/z bin
/// (`floor((mz - minX) / binSize)`), then, for each bin, pads a **zero-intensity** peak for every
/// spectrum that did not contribute a real peak — so outlier rejection and averaging see one peak
/// per spectrum per bin.
///
/// Returned as a `Vec` of bins rather than a keyed map: the bin *index* never affects the output
/// (the final spectrum is re-sorted by m/z), so only the per-bin peak groups matter. Bins are
/// produced in ascending bin-index order.
fn get_bins(x_arrays: &[Vec<f64>], y_arrays: &[Vec<f64>], bin_size: f64) -> Vec<Vec<BinnedPeak>> {
    let num_spectra = x_arrays.len();
    let min_x_value = x_arrays
        .iter()
        .flat_map(|x| x.iter().copied())
        .fold(f64::INFINITY, f64::min);

    // Sort all peaks into bins keyed by bin index. BTreeMap gives deterministic ascending-index
    // iteration; the C# uses a Dictionary whose order is irrelevant post-sort.
    use std::collections::BTreeMap;
    let mut bins: BTreeMap<i64, Vec<BinnedPeak>> = BTreeMap::new();
    for i in 0..num_spectra {
        for j in 0..x_arrays[i].len() {
            let mz = x_arrays[i][j];
            let bin_index = ((mz - min_x_value) / bin_size).floor() as i64;
            bins.entry(bin_index).or_default().push(BinnedPeak {
                mz,
                intensity: y_arrays[i][j],
                spectra_id: i,
            });
        }
    }

    // Pad each bin with zero-intensity peaks for absent spectra. The padded m/z is the running
    // average of the bin's current peaks — evaluated after each insertion, exactly as mzLib does
    // (which, because a zero peak is added at the current mean, leaves the mean unchanged).
    for bin in bins.values_mut() {
        let spectra_in_bin: Vec<usize> = bin.iter().map(|p| p.spectra_id).collect();
        for i in 0..num_spectra {
            if !spectra_in_bin.contains(&i) {
                let mz = mean(bin.iter().map(|p| p.mz));
                bin.push(BinnedPeak { mz, intensity: 0.0, spectra_id: i });
            }
        }
    }

    bins.into_values().collect()
}

/// Faithful port of `SpectraAveraging.AverageBin`. Weighted mean intensity (numerator
/// `Σ intensity·weight`, denominator `Σ weight` over **all** peaks including the zero-intensity
/// padding), and the plain arithmetic mean of the peaks' m/z.
fn average_bin(peaks_in_bin: &[BinnedPeak], weights: &[f64]) -> (f64, f64) {
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for peak in peaks_in_bin {
        numerator += peak.intensity * weights[peak.spectra_id];
        denominator += weights[peak.spectra_id];
    }
    let mz = mean(peaks_in_bin.iter().map(|p| p.mz));
    let intensity = numerator / denominator;
    (mz, intensity)
}

// ---------------------------------------------------------------------------
// Normalization (mirroring SpectraNormalization)
// ---------------------------------------------------------------------------

/// Faithful port of `SpectraNormalization.NormalizeSpectra`. Mutates `y_arrays` in place.
fn normalize_spectra(y_arrays: &mut [Vec<f64>], normalization_type: NormalizationType) {
    match normalization_type {
        NormalizationType::NoNormalization => {}
        NormalizationType::AbsoluteToTic => normalize_absolute_to_tic(y_arrays),
        NormalizationType::RelativeToTics => normalize_relative_to_tics(y_arrays),
        NormalizationType::RelativeIntensity => to_relative_intensity(y_arrays),
    }
}

/// Divide each y by its own TIC (sum-to-one per spectrum). `NormalizeAbsoluteToTic`.
fn normalize_absolute_to_tic(y_arrays: &mut [Vec<f64>]) {
    for y in y_arrays.iter_mut() {
        let mut total_ion_current = sum(y);
        if total_ion_current == 0.0 {
            total_ion_current = 1.0;
        }
        for v in y.iter_mut() {
            *v /= total_ion_current;
        }
    }
}

/// Divide each y by its own TIC then multiply by the average TIC. `NormalizeRelativeToTics`
/// — the default. Puts every spectrum on the same total-intensity scale while preserving the
/// overall magnitude.
fn normalize_relative_to_tics(y_arrays: &mut [Vec<f64>]) {
    let tics: Vec<f64> = y_arrays.iter().map(|y| sum(y)).collect();
    let average_tic = sum(&tics) / tics.len() as f64;
    for i in 0..y_arrays.len() {
        let tic = if tics[i] == 0.0 { 1.0 } else { tics[i] };
        for v in y_arrays[i].iter_mut() {
            *v = *v / tic * average_tic;
        }
    }
}

/// Divide each y by the spectrum's maximum intensity. `ToRelativeIntensity`.
fn to_relative_intensity(y_arrays: &mut [Vec<f64>]) {
    for y in y_arrays.iter_mut() {
        let max_value = y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        for v in y.iter_mut() {
            *v /= max_value;
        }
    }
}

// ---------------------------------------------------------------------------
// Weighting (mirroring SpectralWeighting)
// ---------------------------------------------------------------------------

/// Faithful port of `SpectralWeighting.CalculateSpectraWeights`. Returns weights indexed by
/// spectrum id (mzLib returns a `Dictionary<int,double>`; the keys are a dense `0..count`, so a
/// `Vec` indexed by id is equivalent).
fn calculate_spectra_weights(
    x_arrays: &[Vec<f64>],
    y_arrays: &[Vec<f64>],
    spectra_weighting_type: SpectraWeightingType,
) -> Vec<f64> {
    match spectra_weighting_type {
        SpectraWeightingType::WeightEvenly => vec![1.0; x_arrays.len()],
        SpectraWeightingType::TicValue => weight_by_tic_value(y_arrays),
        SpectraWeightingType::MrsNoiseEstimation => unimplemented!(
            "MrsNoiseEstimation weighting is outside the default-config subset (needs the MRS \
             noise estimator + biweight midvariance); port it when a config requires it"
        ),
    }
}

/// Weight each spectrum by `tic_i / max_tic`. `WeightByTicValue`.
fn weight_by_tic_value(y_arrays: &[Vec<f64>]) -> Vec<f64> {
    let tics: Vec<f64> = y_arrays.iter().map(|y| sum(y)).collect();
    let max_tic = tics.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    tics.iter().map(|&t| t / max_tic).collect()
}

// ---------------------------------------------------------------------------
// Outlier rejection (mirroring OutlierRejection.RejectOutliers)
// ---------------------------------------------------------------------------

/// Faithful port of the `OutlierRejection.RejectOutliers(List<BinnedPeak>, ...)` overload for the
/// **default** `NoRejection` config. The six clipping variants are stubbed (they belong to the
/// not-yet-ported alternate-config work) and panic if dispatched.
fn reject_outliers(
    peaks: &[BinnedPeak],
    parameters: &SpectralAveragingParameters,
) -> Vec<BinnedPeak> {
    match parameters.outlier_rejection_type {
        OutlierRejectionType::NoRejection => peaks.to_vec(),
        other => unimplemented!(
            "outlier rejection {:?} is outside the default-config subset; only NoRejection is \
             ported (see module scope)",
            other
        ),
    }
}

// ---------------------------------------------------------------------------
// Small numeric helpers (kept explicit for summation-order parity)
// ---------------------------------------------------------------------------

/// Left-to-right sum, matching `IEnumerable<double>.Sum()` accumulation order.
#[inline]
fn sum(values: &[f64]) -> f64 {
    values.iter().copied().sum()
}

/// Arithmetic mean, matching `IEnumerable<double>.Average()` (sum then divide by count).
#[inline]
fn mean(values: impl Iterator<Item = f64>) -> f64 {
    let mut total = 0.0;
    let mut count = 0usize;
    for v in values {
        total += v;
        count += 1;
    }
    total / count as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() <= 1e-9 * b.abs().max(1.0), "expected {b}, got {a}");
    }

    /// Two spectra, one shared bin (0.01 wide), default config. Hand-computed:
    /// tics 10/20 → averageTic 15 → both normalize to 15 → bin mean intensity 15,
    /// bin m/z = mean(100.00, 100.005) = 100.0025.
    #[test]
    fn two_spectra_single_shared_bin() {
        let x = vec![vec![100.000], vec![100.005]];
        let y = vec![vec![10.0], vec![20.0]];
        let (mz, inten) = average_spectra(&x, &y, &SpectralAveragingParameters::default());
        assert_eq!(mz.len(), 1);
        approx(mz[0], 100.0025);
        approx(inten[0], 15.0);
    }

    /// Padding path: spectrum B lacks the 200 m/z bin, so it is padded with a zero peak, and the
    /// bin intensity divides by the spectrum count (2), not the present-peak count (1).
    /// A: x=[100,200] y=[10,30] tic 40; B: x=[100] y=[20] tic 20; averageTic 30.
    /// Normalized A=[7.5,22.5], B=[30]. Bin100: (7.5+30)/2=18.75. Bin200: (22.5+0)/2=11.25.
    #[test]
    fn padding_divides_by_spectrum_count() {
        let x = vec![vec![100.000, 200.000], vec![100.000]];
        let y = vec![vec![10.0, 30.0], vec![20.0]];
        let (mz, inten) = average_spectra(&x, &y, &SpectralAveragingParameters::default());
        assert_eq!(mz.len(), 2);
        approx(mz[0], 100.0);
        approx(inten[0], 18.75);
        approx(mz[1], 200.0);
        approx(inten[1], 11.25);
    }

    /// With `NoNormalization` + `WeightEvenly`, a single spectrum passes through as-is (bins of one
    /// real peak each, no padding), only re-sorted by m/z with zero-intensity peaks dropped.
    #[test]
    fn no_normalization_single_spectrum_passthrough() {
        let params = SpectralAveragingParameters {
            normalization_type: NormalizationType::NoNormalization,
            ..SpectralAveragingParameters::default()
        };
        let x = vec![vec![300.0, 100.0, 200.0]];
        let y = vec![vec![3.0, 1.0, 2.0]];
        let (mz, inten) = average_spectra(&x, &y, &params);
        assert_eq!(mz, vec![100.0, 200.0, 300.0]);
        assert_eq!(inten, vec![1.0, 2.0, 3.0]);
    }

    /// Zero-intensity bins are dropped from the output entirely.
    #[test]
    fn zero_intensity_bins_dropped() {
        let params = SpectralAveragingParameters {
            normalization_type: NormalizationType::NoNormalization,
            ..SpectralAveragingParameters::default()
        };
        let x = vec![vec![100.0, 200.0]];
        let y = vec![vec![0.0, 5.0]];
        let (mz, inten) = average_spectra(&x, &y, &params);
        assert_eq!(mz, vec![200.0]);
        assert_eq!(inten, vec![5.0]);
    }

    /// `TicValue` weighting: two spectra sharing a bin, weights tic_i / max_tic.
    /// A tic 10 → w 0.5; B tic 20 → w 1.0. NoNormalization keeps raw intensities.
    /// intensity = (10·0.5 + 20·1.0)/(0.5+1.0) = 25/1.5 = 16.666…
    #[test]
    fn tic_value_weighting() {
        let params = SpectralAveragingParameters {
            normalization_type: NormalizationType::NoNormalization,
            spectral_weighting_type: SpectraWeightingType::TicValue,
            ..SpectralAveragingParameters::default()
        };
        let x = vec![vec![100.000], vec![100.005]];
        let y = vec![vec![10.0], vec![20.0]];
        let (mz, inten) = average_spectra(&x, &y, &params);
        assert_eq!(mz.len(), 1);
        approx(inten[0], 25.0 / 1.5);
    }

    #[test]
    #[should_panic(expected = "outside the default-config subset")]
    fn sigma_clipping_panics_as_stub() {
        let params = SpectralAveragingParameters {
            outlier_rejection_type: OutlierRejectionType::SigmaClipping,
            ..SpectralAveragingParameters::default()
        };
        let x = vec![vec![100.0], vec![100.0]];
        let y = vec![vec![10.0], vec![20.0]];
        let _ = average_spectra(&x, &y, &params);
    }
}
