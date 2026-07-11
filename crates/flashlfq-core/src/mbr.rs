//! Match-between-runs: retention-time alignment spline (PLAN.md P3.1).
//!
//! Ports `FlashLfqEngine.GetRtCalSpline` (+ the `ChooseBestPeak` helper and the
//! `RetentionTimeCalibDataPoint` it builds) — the first piece of FlashLFQ's MBR algorithm.
//! Given the error-checked `ChromatographicPeak` lists of two runs (a *donor* and an
//! *acceptor*), it builds the **RT calibration spline**: the set of anchor peptides MS2-detected
//! in *both* runs, each contributing a `(donor apex RT, acceptor apex RT)` data point, sorted by
//! donor apex RT. Later MBR stages (P3.2, `PredictRetentionTime`) binary-search this spline to
//! predict where a donor peptide should elute in the acceptor run.
//!
//! ## Faithfulness to `FlashLfqEngine.cs`
//!
//! [`get_rt_cal_spline`] mirrors `GetRtCalSpline` (`FlashLfqEngine.cs:559`):
//! 1. From each run keep only peaks that are **unambiguous, MS2-detected, actually traced, and
//!    confidently identified** — `NumIdentificationsByFullSeq == 1`, `DetectionType == MSMS`,
//!    `IsotopicEnvelopes.Any()`, and `Identifications.Min(QValue) < DonorQValueThreshold`.
//! 2. Group those by modified sequence and pick one [`choose_best_peak`] per sequence
//!    (`ChooseBestPeak`, `FlashLfqEngine.cs:670`).
//! 3. For every acceptor best-peak whose sequence also has a donor best-peak, emit a
//!    [`RetentionTimeCalibDataPoint`]; collect `donorApexRT - acceptorApexRT` for the anchor
//!    peptides (apex RTs both `> 0`) — the distribution P3.3's scorer consumes.
//! 4. Order the donor best-peaks by peakfinding mass (the `GetRandomPeak` lookup list, P3.2) and
//!    sort the calibration curve by donor apex RT.
//!
//! The `scorer.AddRtPredErrorDistribution(...)` call in the C# method only *stores* the anchor
//! RT-diff distribution on the (P3.3) `MbrScorer`; here [`get_rt_cal_spline`] returns
//! `anchor_rt_diffs` so a later task feeds it to the scorer. The peaks themselves are genuine
//! engine output (the L5/L6-parity `run_msms` peaks), so the spline is built over real data.

use crate::chromatographic_peak::ChromatographicPeak;
use crate::detection_type::DetectionType;

/// Default donor q-value threshold (`FlashLfqParameters.DonorQValueThreshold`, set to `0.01` in
/// the parameters constructor, `FlashLfqParameters.cs:36`). An anchor peptide's best identification
/// must have a q-value strictly below this to seed the alignment.
pub const DONOR_Q_VALUE_THRESHOLD: f64 = 0.01;

/// Number of anchor peptides used on each side for local alignment when predicting acceptor RTs
/// (`FlashLfqEngine.NumberOfAnchorPeptidesForMbr`, `FlashLfqEngine.cs:35`). Unused by the spline
/// construction itself; carried here for the P3.2 `PredictRetentionTime` port.
pub const NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR: usize = 3;

/// Window (minutes) in which the `Neighbors` donor criterion counts neighbouring peaks
/// (`FlashLfqEngine.MbrAlignmentWindow`, `FlashLfqEngine.cs:36`).
pub const MBR_ALIGNMENT_WINDOW: f64 = 2.5;

/// Maximum half-window (minutes) the predicted-RT search window is ever allowed to span
/// (`FlashLfqParameters.MaxMbrRtWindow`, set to `1.0` in the parameters constructor,
/// `FlashLfqParameters.cs:32`). [`predict_retention_time`] clamps the `6·stddev` RT range to this.
pub const MAX_MBR_RT_WINDOW: f64 = 1.0;

/// Ppm tolerance used for MBR acceptor-peak matching (`FlashLfqParameters.MbrPpmTolerance`,
/// `FlashLfqParameters.cs:31`). Carried here for the P3.2 acceptor search; unused by RT prediction.
pub const MBR_PPM_TOLERANCE: f64 = 10.0;

/// How the single best peak is chosen among the candidates for one modified sequence.
/// Port of `FlashLFQ.DonorCriterion` (`FlashLfqParameters.cs:98`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DonorCriterion {
    /// Select the peak with the highest-scoring identification (the parameters-constructor
    /// default, `FlashLfqParameters.cs:35`). Falls back to [`DonorCriterion::Intensity`] when the
    /// chosen peak's first identification has a non-positive PSM score.
    Score,
    /// Select the peak with the most neighbouring peaks (distinct sequences) within
    /// [`MBR_ALIGNMENT_WINDOW`]. Not used by the default engine path; see
    /// [`choose_best_peak`] for the (currently unsupported) status.
    Neighbors,
    /// Select the most intense peak.
    Intensity,
}

/// One anchor point of the RT calibration spline: a peptide MS2-detected in both runs, pairing the
/// donor-run peak with the acceptor-run peak. Port of `FlashLFQ.RetentionTimeCalibDataPoint`.
///
/// The C# class holds object references to the two `ChromatographicPeak`s; here we keep owned
/// clones (the peaks are cheap value types in this port, and the spline outlives the per-file peak
/// vectors in the MBR loop). `rt_diff = acceptorApexRT - donorApexRT` matches the C# field, and is
/// `NaN` when either peak is absent (only the donor-present, acceptor-absent test-point shape).
#[derive(Debug, Clone)]
pub struct RetentionTimeCalibDataPoint {
    /// The peak in the donor run.
    pub donor_peak: ChromatographicPeak,
    /// The peak in the acceptor run, if any (the spline only stores paired points; a `None`
    /// acceptor matches the C# `new RetentionTimeCalibDataPoint(donorPeak, null)` test point used
    /// by the P3.2 binary search).
    pub acceptor_peak: Option<ChromatographicPeak>,
    /// `acceptorApexRT - donorApexRT` (C# `RtDiff`); `NaN` when `acceptor_peak` is `None`.
    pub rt_diff: f64,
}

impl RetentionTimeCalibDataPoint {
    /// Builds a data point, computing `rt_diff` exactly as the C# constructor does.
    pub fn new(
        donor_peak: ChromatographicPeak,
        acceptor_peak: Option<ChromatographicPeak>,
    ) -> RetentionTimeCalibDataPoint {
        let rt_diff = match &acceptor_peak {
            Some(acc) => acc.apex_retention_time() - donor_peak.apex_retention_time(),
            None => f64::NAN,
        };
        RetentionTimeCalibDataPoint {
            donor_peak,
            acceptor_peak,
            rt_diff,
        }
    }

    /// The donor peak's apex retention time (C# `DonorFilePeak.Apex.IndexedPeak.RetentionTime`).
    pub fn donor_apex_rt(&self) -> f64 {
        self.donor_peak.apex_retention_time()
    }
}

/// The products of [`get_rt_cal_spline`], mirroring `GetRtCalSpline`'s return value plus its two
/// `out`/side-channel outputs.
#[derive(Debug, Clone, Default)]
pub struct RtCalSpline {
    /// The RT calibration spline: paired donor/acceptor anchor points, sorted by donor apex RT
    /// (the C# `return rtCalibrationCurve.OrderBy(p => p.DonorFilePeak.Apex.IndexedPeak.RetentionTime)`).
    pub calibration_curve: Vec<RetentionTimeCalibDataPoint>,
    /// `donorApexRT - acceptorApexRT` for every anchor peptide with both apex RTs `> 0` — the
    /// distribution C# hands to `scorer.AddRtPredErrorDistribution`.
    pub anchor_rt_diffs: Vec<f64>,
    /// Donor best-peaks ordered by their first identification's peakfinding mass (the C#
    /// `donorFileBestMsmsPeaksOrderedByMass` out-parameter, used by `GetRandomPeak` in P3.2).
    pub donor_best_peaks_ordered_by_mass: Vec<ChromatographicPeak>,
}

/// Whether a peak qualifies as an MBR anchor candidate: a single unambiguous full sequence,
/// MS2-detected, actually traced (has envelopes), and confidently identified
/// (`min QValue < donor_q_value_threshold`). Mirrors the `.Where(...)` filter applied to both
/// `_results.Peaks[donor]` and `_results.Peaks[acceptor]` in `GetRtCalSpline`.
fn is_anchor_candidate(peak: &ChromatographicPeak, donor_q_value_threshold: f64) -> bool {
    if peak.num_identifications_by_full_seq != 1
        || peak.detection_type != DetectionType::MSMS
        || peak.isotopic_envelopes.is_empty()
    {
        return false;
    }
    // C# `Identifications.Min(id => id.QValue)`; every anchor candidate has ids here.
    let min_q = peak
        .identifications
        .iter()
        .map(|id| id.q_value)
        .fold(f64::INFINITY, f64::min);
    min_q < donor_q_value_threshold
}

/// Groups the anchor-candidate peaks of one run by their first identification's modified sequence,
/// preserving first-seen key order and within-group order (C# LINQ `GroupBy`, which is order
/// preserving), then picks one [`choose_best_peak`] per sequence.
///
/// Returns the `(modified_sequence, best_peak)` pairs in group-key order (the C# dictionary
/// insertion order), so the acceptor-side iteration that builds the calibration curve matches C#.
fn best_msms_peaks_by_sequence<'a>(
    peaks: &'a [ChromatographicPeak],
    donor_q_value_threshold: f64,
    criterion: DonorCriterion,
) -> Vec<(String, &'a ChromatographicPeak)> {
    // Ordered grouping: keys in first-seen order, each holding its peaks in input order.
    let mut keys: Vec<String> = Vec::new();
    let mut groups: Vec<Vec<&ChromatographicPeak>> = Vec::new();

    for peak in peaks {
        if !is_anchor_candidate(peak, donor_q_value_threshold) {
            continue;
        }
        let Some(seq) = peak.first_modified_sequence() else {
            continue;
        };
        match keys.iter().position(|k| k == seq) {
            Some(i) => groups[i].push(peak),
            None => {
                keys.push(seq.to_string());
                groups.push(vec![peak]);
            }
        }
    }

    let mut out: Vec<(String, &ChromatographicPeak)> = Vec::with_capacity(keys.len());
    for (key, group) in keys.into_iter().zip(groups.into_iter()) {
        // C# `if (!peaksForPeptide.Any()) continue;` — never empty here, kept for fidelity.
        if group.is_empty() {
            continue;
        }
        if let Some(best) = choose_best_peak(&group, criterion) {
            out.push((key, best));
        }
    }
    out
}

/// Picks the single best peak among candidates for one modified sequence. Port of
/// `FlashLfqEngine.ChooseBestPeak` (`FlashLfqEngine.cs:670`).
///
/// - [`DonorCriterion::Score`]: the peak whose highest-scoring identification is greatest (LINQ
///   `MaxBy` → first peak attaining the max). If that peak's *first* identification has a positive
///   PSM score the result stands; otherwise (every score `≤ 0`) it falls through to the intensity
///   rule — exactly the C# `goto default`.
/// - [`DonorCriterion::Intensity`] (and the fall-through): the most intense peak.
/// - [`DonorCriterion::Neighbors`]: requires the whole file's peak list to count neighbours and is
///   not used by the default engine path; it is intentionally not implemented here (P3.2 will add
///   it if a non-default criterion is ever exercised). It currently behaves like `Intensity`'s
///   fall-through is *not* taken — instead it panics to make an unported path obvious.
///
/// `peaks` is assumed non-empty (the caller skips empty groups), matching the C# precondition.
pub fn choose_best_peak<'a>(
    peaks: &[&'a ChromatographicPeak],
    criterion: DonorCriterion,
) -> Option<&'a ChromatographicPeak> {
    match criterion {
        DonorCriterion::Score => {
            let best = max_by_key_first(peaks, |p| max_psm_score(p))?;
            let first_score = best
                .identifications
                .first()
                .map(|id| id.score)
                .unwrap_or(0.0);
            if first_score > 0.0 {
                Some(best)
            } else {
                // C# `goto default;` — fall through to the intensity rule.
                max_by_key_first(peaks, |p| p.intensity)
            }
        }
        DonorCriterion::Intensity => max_by_key_first(peaks, |p| p.intensity),
        DonorCriterion::Neighbors => {
            unimplemented!(
                "DonorCriterion::Neighbors needs the full file peak list; the default engine \
                 path uses Score, so this is unported (PLAN.md P3.2)"
            )
        }
    }
}

/// The maximum PSM score among a peak's identifications (C# `peak.Identifications.Max(id => id.PsmScore)`).
/// `NEG_INFINITY` for a peak with no identifications (never the case for engine peaks).
fn max_psm_score(peak: &ChromatographicPeak) -> f64 {
    peak.identifications
        .iter()
        .map(|id| id.score)
        .fold(f64::NEG_INFINITY, f64::max)
}

/// LINQ `MaxBy` semantics: returns the **first** element attaining the maximum key, or `None` for
/// an empty slice. The key is `f64`; ties keep the earlier element (stable, like C# `MaxBy`).
fn max_by_key_first<'a>(
    peaks: &[&'a ChromatographicPeak],
    key: impl Fn(&ChromatographicPeak) -> f64,
) -> Option<&'a ChromatographicPeak> {
    let mut best: Option<(&ChromatographicPeak, f64)> = None;
    for &p in peaks {
        let k = key(p);
        match best {
            None => best = Some((p, k)),
            Some((_, bk)) if k > bk => best = Some((p, k)),
            _ => {}
        }
    }
    best.map(|(p, _)| p)
}

/// Builds the RT calibration spline between a donor and an acceptor run. Port of
/// `FlashLfqEngine.GetRtCalSpline` (`FlashLfqEngine.cs:559`).
///
/// `donor_peaks` / `acceptor_peaks` are the error-checked `ChromatographicPeak` lists for the two
/// runs (`_results.Peaks[donor]` / `_results.Peaks[acceptor]`). `donor_q_value_threshold` is
/// [`DONOR_Q_VALUE_THRESHOLD`] by default; `criterion` is [`DonorCriterion::Score`] by default.
///
/// Returns the sorted calibration curve, the anchor RT-diff distribution, and the donor best-peaks
/// ordered by peakfinding mass (see [`RtCalSpline`]).
pub fn get_rt_cal_spline(
    donor_peaks: &[ChromatographicPeak],
    acceptor_peaks: &[ChromatographicPeak],
    donor_q_value_threshold: f64,
    criterion: DonorCriterion,
) -> RtCalSpline {
    let donor_best = best_msms_peaks_by_sequence(donor_peaks, donor_q_value_threshold, criterion);
    let acceptor_best =
        best_msms_peaks_by_sequence(acceptor_peaks, donor_q_value_threshold, criterion);

    // create RT calibration curve: for each acceptor best-peak whose sequence also has a donor
    // best-peak, pair them. Iterate acceptor best-peaks in group order (C# dictionary order).
    let mut calibration_curve: Vec<RetentionTimeCalibDataPoint> = Vec::new();
    let mut anchor_rt_diffs: Vec<f64> = Vec::new();
    for (seq, acceptor_peak) in &acceptor_best {
        if let Some((_, donor_peak)) = donor_best.iter().find(|(k, _)| k == seq) {
            calibration_curve.push(RetentionTimeCalibDataPoint::new(
                (*donor_peak).clone(),
                Some((*acceptor_peak).clone()),
            ));
            let donor_rt = donor_peak.apex_retention_time();
            let acceptor_rt = acceptor_peak.apex_retention_time();
            if donor_rt > 0.0 && acceptor_rt > 0.0 {
                anchor_rt_diffs.push(donor_rt - acceptor_rt);
            }
        }
    }

    // donorFileBestMsmsPeaksOrderedByMass = donor best-peaks ordered by first id's peakfinding mass.
    let mut donor_best_peaks_ordered_by_mass: Vec<ChromatographicPeak> =
        donor_best.iter().map(|(_, p)| (*p).clone()).collect();
    donor_best_peaks_ordered_by_mass.sort_by(|a, b| {
        first_peakfinding_mass(a)
            .partial_cmp(&first_peakfinding_mass(b))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // sort the calibration curve by donor apex RT (stable, matching C# OrderBy).
    calibration_curve.sort_by(|a, b| {
        a.donor_apex_rt()
            .partial_cmp(&b.donor_apex_rt())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    RtCalSpline {
        calibration_curve,
        anchor_rt_diffs,
        donor_best_peaks_ordered_by_mass,
    }
}

/// The first identification's peakfinding mass for a peak (C# `Identifications.First().PeakfindingMass`),
/// read from the parallel `identification_peakfinding_masses` list. `NaN` if absent.
fn first_peakfinding_mass(peak: &ChromatographicPeak) -> f64 {
    peak.identification_peakfinding_masses
        .first()
        .copied()
        .unwrap_or(f64::NAN)
}

/// The predicted location + width of an acceptor peak, mirroring `FlashLFQ.RtInfo`
/// (`MBR/RtInfo.cs`). [`predict_retention_time`] produces it; the P3.2 acceptor search consumes
/// `rt_start_hypothesis`/`rt_end_hypothesis` to bound the MS1 scan region to look in. `width` is
/// mutable in C# (the search-widening loop grows it); kept public here for the same use.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RtInfo {
    /// The predicted acceptor-run retention time of the donor peptide (C# `PredictedRt`).
    pub predicted_rt: f64,
    /// The full width (minutes) of the predicted RT window (C# `Width`).
    pub width: f64,
}

impl RtInfo {
    /// Builds an `RtInfo` (C# `new RtInfo(predictedRt, width)`).
    pub fn new(predicted_rt: f64, width: f64) -> RtInfo {
        RtInfo {
            predicted_rt,
            width,
        }
    }

    /// Lower bound of the predicted RT window, `PredictedRt - Width/2` (C# `RtStartHypothesis`).
    pub fn rt_start_hypothesis(&self) -> f64 {
        self.predicted_rt - (self.width / 2.0)
    }

    /// Upper bound of the predicted RT window, `PredictedRt + Width/2` (C# `RtEndHypothesis`).
    pub fn rt_end_hypothesis(&self) -> f64 {
        self.predicted_rt + (self.width / 2.0)
    }
}

/// Median of the samples, matching MathNet `Statistics.Median` (the order statistic; average of the
/// two central values for an even count). `NaN` for an empty slice. No NaN inputs are expected (the
/// RT diffs are finite), so the sort uses a total order via `partial_cmp().unwrap()`.
pub(crate) fn median(samples: &[f64]) -> f64 {
    let n = samples.len();
    if n == 0 {
        return f64::NAN;
    }
    let mut s = samples.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).expect("RT diffs are finite"));
    if n % 2 == 1 {
        s[(n - 1) / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

/// Sample standard deviation (Bessel-corrected, `N-1` denominator), matching MathNet
/// `Statistics.StandardDeviation` for an `IEnumerable<double>` — i.e. the streaming variance of
/// `StreamingStatistics.Variance`, then `sqrt`. `NaN` for fewer than two samples. The running
/// formula (and thus the floating-point accumulation order) is reproduced verbatim so the result
/// matches the C# `List<double>.StandardDeviation()` used by `PredictRetentionTime`.
pub(crate) fn standard_deviation(samples: &[f64]) -> f64 {
    let n = samples.len();
    if n <= 1 {
        return f64::NAN;
    }
    // StreamingStatistics.Variance: j is the running 1-based count, t the running sum.
    let mut variance = 0.0_f64;
    let mut t = samples[0];
    for i in 1..n {
        let xi = samples[i];
        t += xi;
        let j = (i + 1) as f64;
        let diff = j * xi - t;
        variance += (diff * diff) / (j * (j - 1.0));
    }
    (variance / (n as f64 - 1.0)).sqrt()
}

/// Sample variance (Bessel-corrected, `N-1` denominator), matching MathNet
/// `Statistics.Variance` for a `List<double>` (i.e. `StreamingStatistics.Variance`). `NaN` for
/// fewer than two samples. Uses the identical running accumulation as [`standard_deviation`]
/// (which is its `sqrt`), so the floating-point order matches the C# used by the MBR scorer's
/// isotopic-correlation gamma fit. Consumed by [`crate::mbr_scorer`].
pub(crate) fn variance(samples: &[f64]) -> f64 {
    let n = samples.len();
    if n <= 1 {
        return f64::NAN;
    }
    let mut variance = 0.0_f64;
    let mut t = samples[0];
    for i in 1..n {
        let xi = samples[i];
        t += xi;
        let j = (i + 1) as f64;
        let diff = j * xi - t;
        variance += (diff * diff) / (j * (j - 1.0));
    }
    variance / (n as f64 - 1.0)
}

/// Incremental arithmetic mean, matching MathNet `Statistics.Mean` for a `List<double>`
/// (`StreamingStatistics.Mean`): `mean += (x - mean) / count`. `NaN` for an empty slice. Distinct
/// from [`average`] (LINQ `Enumerable.Average`, a plain sum/count) — the scorer uses **both**, so
/// the two are kept separate to preserve the C# floating-point results bit-for-bit.
pub(crate) fn mean(samples: &[f64]) -> f64 {
    if samples.is_empty() {
        return f64::NAN;
    }
    let mut m = 0.0_f64;
    let mut count = 0.0_f64;
    for &d in samples {
        count += 1.0;
        m += (d - m) / count;
    }
    m
}

/// LINQ `Enumerable.Average`: a plain running sum divided by the count, in iteration order. `NaN`
/// for an empty slice (the C# overload throws on empty; the scorer guards `Count > 0` upstream).
/// Distinct from [`mean`] (MathNet's incremental mean) — see that function's note.
pub(crate) fn average(samples: &[f64]) -> f64 {
    let mut sum = 0.0_f64;
    for &d in samples {
        sum += d;
    }
    sum / samples.len() as f64
}

/// The `tau`-quantile of the samples, matching MathNet `ArrayStatistics.QuantileInplace` — the
/// default (type-8, median-unbiased) definition `h = (n + 1/3)·tau + 1/3`, linearly interpolating
/// between the two bracketing order statistics. `NaN` for an empty slice or `tau` outside `[0, 1]`.
/// MathNet's partial `SelectInplace` is replaced by a full sort (same value at each order rank).
pub(crate) fn quantile(samples: &[f64], tau: f64) -> f64 {
    let n = samples.len();
    if !(0.0..=1.0).contains(&tau) || n == 0 {
        return f64::NAN;
    }
    let mut s = samples.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).expect("scorer samples are finite"));
    let h = (n as f64 + 1.0 / 3.0) * tau + 1.0 / 3.0;
    let hf = h.trunc() as i64; // C# (int)h, truncation toward zero (h > 0 here, so floor)
    if hf <= 0 || tau == 0.0 {
        return s[0];
    }
    if hf >= n as i64 || tau == 1.0 {
        return s[n - 1];
    }
    let a = s[(hf - 1) as usize];
    let b = s[hf as usize];
    a + (h - hf as f64) * (b - a)
}

/// Interquartile range, matching MathNet `Statistics.InterquartileRange`:
/// [`quantile`]`(0.75) -` [`quantile`]`(0.25)`. Used by the MBR scorer's distribution fits.
pub(crate) fn interquartile_range(samples: &[f64]) -> f64 {
    quantile(samples, 0.75) - quantile(samples, 0.25)
}

/// Three-way comparison matching C# `double.CompareTo` (NaN sorts below everything, NaN==NaN).
/// Returns `-1`/`0`/`1`. Used only for finite RTs in the binary search.
fn compare_f64(a: f64, b: f64) -> i32 {
    if a < b {
        -1
    } else if a > b {
        1
    } else if a == b {
        0
    } else {
        // NaN handling: a==NaN sorts as less than a non-NaN b; two NaNs are equal.
        match (a.is_nan(), b.is_nan()) {
            (true, true) => 0,
            (true, false) => -1,
            (false, true) => 1,
            (false, false) => 0, // unreachable
        }
    }
}

/// Faithful port of .NET `Array.BinarySearch` over the calibration curve, comparing each point's
/// donor apex RT against `key_rt` (C# `RetentionTimeCalibDataPoint.CompareTo`, which compares donor
/// apex RTs). Returns the index of a matching element, or the bitwise complement `~insertionPoint`
/// when not found — exactly the .NET `lo/hi` midpoint algorithm, so the chosen index on duplicate
/// keys matches C#.
fn binary_search_donor_rt(curve: &[RetentionTimeCalibDataPoint], key_rt: f64) -> i64 {
    let mut lo: i64 = 0;
    let mut hi: i64 = curve.len() as i64 - 1;
    while lo <= hi {
        let i = lo + ((hi - lo) >> 1);
        let order = compare_f64(curve[i as usize].donor_apex_rt(), key_rt);
        if order == 0 {
            return i;
        }
        if order < 0 {
            lo = i + 1;
        } else {
            hi = i - 1;
        }
    }
    !lo // ~lo == -lo - 1
}

/// Predicts where a donor peptide elutes in the acceptor run by **local RT alignment** against the
/// calibration spline. Port of `FlashLfqEngine.PredictRetentionTime` (`FlashLfqEngine.cs:716`).
///
/// Binary-searches `calibration_curve` (sorted by donor apex RT) for `donor_peak`, then gathers up
/// to [`NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR`] anchor points on each side — forward from the insertion
/// point, then backward — stopping early once a neighbour's donor-RT difference exceeds `0.5`
/// minutes (no longer "local"). From the donor−acceptor RT differences of those anchors it returns:
/// - no anchors: the donor peak's own RT with a `0.25`-minute width;
/// - one anchor: `donorRT − medianDiff` with a `0.25`-minute width;
/// - otherwise: `donorRT − medianDiff` with width `min(6·stddev, max_mbr_rt_window)`.
///
/// The anchors are appended forward-then-backward exactly as C# does, so the
/// [`standard_deviation`] accumulation order (and therefore the width) matches bit-for-bit.
///
/// **Fraction gating omitted on purpose:** the C# method's leading early-return (skip donor/acceptor
/// pairs more than one fraction apart) depends on `SpectraFileInfo.Fraction`, which the core
/// [`ChromatographicPeak`] does not model. That gate belongs to the P3.2 MBR orchestration (where
/// file metadata lives); this function ports the alignment math, which always yields an `RtInfo`.
pub fn predict_retention_time(
    calibration_curve: &[RetentionTimeCalibDataPoint],
    donor_peak: &ChromatographicPeak,
    max_mbr_rt_window: f64,
    number_of_anchor_peptides: usize,
) -> RtInfo {
    let donor_rt = donor_peak.apex_retention_time();

    // binary search for this donor peak in the calibration spline (by donor apex RT).
    let mut index = binary_search_donor_rt(calibration_curve, donor_rt);
    if index < 0 {
        index = !index; // ~index
    }
    let len = calibration_curve.len() as i64;
    if index >= len && index >= 1 {
        index = len - 1;
    }

    // Helper: does this curve point have a usable acceptor apex (non-null, apex RT > 0)?
    let acceptor_rt_if_usable = |p: &RetentionTimeCalibDataPoint| -> Option<f64> {
        p.acceptor_peak.as_ref().and_then(|acc| {
            let rt = acc.apex_retention_time();
            (rt > 0.0).then_some(rt)
        })
    };

    let mut nearby: Vec<RetentionTimeCalibDataPoint> = Vec::new();

    // forward anchors: r = index+1 .. len
    let mut forward_anchors = 0usize;
    let mut r = index + 1;
    while r < len {
        let point = &calibration_curve[r as usize];
        let rt_diff = point.donor_apex_rt() - donor_rt;
        if acceptor_rt_if_usable(point).is_some() {
            if rt_diff.abs() > 0.5 {
                break;
            }
            nearby.push(point.clone());
            forward_anchors += 1;
            if forward_anchors >= number_of_anchor_peptides {
                break;
            }
        }
        r += 1;
    }

    // backward anchors: r = index-1 .. 0
    let mut backward_anchors = 0usize;
    let mut r = index - 1;
    while r >= 0 {
        let point = &calibration_curve[r as usize];
        let rt_diff = point.donor_apex_rt() - donor_rt;
        if acceptor_rt_if_usable(point).is_some() {
            if rt_diff.abs() > 0.5 {
                break;
            }
            nearby.push(point.clone());
            backward_anchors += 1;
            if backward_anchors >= number_of_anchor_peptides {
                break;
            }
        }
        r -= 1;
    }

    if nearby.is_empty() {
        // No nearby calibration points: donor RT with a 15-second (0.25 min) width.
        return RtInfo::new(donor_rt, 0.25);
    }

    // donor − acceptor RT difference for each anchor (C# `DonorFilePeak.ApexRT - AcceptorFilePeak.ApexRT`).
    let rt_diffs: Vec<f64> = nearby
        .iter()
        .map(|p| {
            p.donor_apex_rt()
                - p.acceptor_peak
                    .as_ref()
                    .expect("nearby anchors have an acceptor peak")
                    .apex_retention_time()
        })
        .collect();

    let median_rt_diff = median(&rt_diffs);
    if rt_diffs.len() == 1 {
        return RtInfo::new(donor_rt - median_rt_diff, 0.25);
    }

    let rt_range = (standard_deviation(&rt_diffs) * 6.0).min(max_mbr_rt_window);
    RtInfo::new(donor_rt - median_rt_diff, rt_range)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isotopic_envelope::IsotopicEnvelope;
    use crate::peak_indexing::IndexedMassSpectralPeak;
    use crate::psm_tsv::Identification;

    fn id(modseq: &str, base: &str, charge: i32, score: f64, q_value: f64) -> Identification {
        Identification {
            file_name: "f".to_string(),
            base_sequence: base.to_string(),
            modified_sequence: modseq.to_string(),
            monoisotopic_mass: 1000.0,
            ms2_retention_time_in_minutes: 1.0,
            precursor_charge_state: charge,
            score,
            q_value,
            is_decoy: false,
        }
    }

    /// Envelope at the given scan/RT/intensity (charge 1).
    fn env(scan: i32, rt: f64, intensity: f64) -> IsotopicEnvelope {
        let peak = IndexedMassSpectralPeak::new(1000.0, intensity, scan, rt);
        IsotopicEnvelope::new(peak, 1, intensity, 0.95)
    }

    /// A traced MSMS peak for one identification, with a single envelope fixing its apex RT.
    fn peak(modseq: &str, score: f64, q_value: f64, apex_rt: f64, intensity: f64) -> ChromatographicPeak {
        let mut p =
            ChromatographicPeak::from_identification(id(modseq, modseq, 2, score, q_value), 1000.0);
        p.isotopic_envelopes = vec![env(0, apex_rt, intensity)];
        p.calculate_intensity_for_this_feature(false);
        p
    }

    #[test]
    fn anchor_candidate_filter_rejects_high_qvalue_and_empty() {
        let good = peak("SEQ", 10.0, 0.005, 30.0, 100.0);
        assert!(is_anchor_candidate(&good, DONOR_Q_VALUE_THRESHOLD));

        let high_q = peak("SEQ", 10.0, 0.02, 30.0, 100.0);
        assert!(!is_anchor_candidate(&high_q, DONOR_Q_VALUE_THRESHOLD));

        let mut no_env = peak("SEQ", 10.0, 0.005, 30.0, 100.0);
        no_env.isotopic_envelopes.clear();
        no_env.calculate_intensity_for_this_feature(false);
        assert!(!is_anchor_candidate(&no_env, DONOR_Q_VALUE_THRESHOLD));
    }

    #[test]
    fn choose_best_peak_score_then_intensity_fallback() {
        // Two peaks for one sequence: lower score but higher intensity vs higher score.
        let hi_score = peak("SEQ", 20.0, 0.001, 30.0, 10.0);
        let hi_int = peak("SEQ", 5.0, 0.001, 31.0, 999.0);
        let cands = vec![&hi_int, &hi_score];
        let best = choose_best_peak(&cands, DonorCriterion::Score).unwrap();
        assert_eq!(best.apex_retention_time(), 30.0); // the higher-score peak

        // All scores zero -> fall through to intensity.
        let z1 = peak("SEQ", 0.0, 0.001, 30.0, 10.0);
        let z2 = peak("SEQ", 0.0, 0.001, 31.0, 999.0);
        let cands = vec![&z1, &z2];
        let best = choose_best_peak(&cands, DonorCriterion::Score).unwrap();
        assert_eq!(best.apex_retention_time(), 31.0); // the more intense peak
    }

    #[test]
    fn spline_pairs_shared_sequences_and_sorts_by_donor_rt() {
        // donor has A (rt 30), B (rt 10), C (rt 20); acceptor has A (rt 32), B (rt 12), D (rt 99).
        let donor = vec![
            peak("A", 10.0, 0.001, 30.0, 100.0),
            peak("B", 10.0, 0.001, 10.0, 100.0),
            peak("C", 10.0, 0.001, 20.0, 100.0),
        ];
        let acceptor = vec![
            peak("A", 10.0, 0.001, 32.0, 100.0),
            peak("B", 10.0, 0.001, 12.0, 100.0),
            peak("D", 10.0, 0.001, 99.0, 100.0),
        ];
        let spline =
            get_rt_cal_spline(&donor, &acceptor, DONOR_Q_VALUE_THRESHOLD, DonorCriterion::Score);

        // Only A and B are shared -> two anchor points, sorted by donor apex RT (B=10, A=30).
        assert_eq!(spline.calibration_curve.len(), 2);
        assert_eq!(spline.calibration_curve[0].donor_apex_rt(), 10.0);
        assert_eq!(spline.calibration_curve[1].donor_apex_rt(), 30.0);
        // rt_diff = acceptor - donor.
        assert!((spline.calibration_curve[0].rt_diff - (12.0 - 10.0)).abs() < 1e-9);
        assert!((spline.calibration_curve[1].rt_diff - (32.0 - 30.0)).abs() < 1e-9);
        // anchor diffs = donor - acceptor for both points.
        assert_eq!(spline.anchor_rt_diffs.len(), 2);
        // donor best-peaks ordered by mass = all three donor sequences (same mass -> stable).
        assert_eq!(spline.donor_best_peaks_ordered_by_mass.len(), 3);
    }

    #[test]
    fn spline_excludes_high_qvalue_anchors() {
        let donor = vec![peak("A", 10.0, 0.02, 30.0, 100.0)]; // q-value too high
        let acceptor = vec![peak("A", 10.0, 0.001, 32.0, 100.0)];
        let spline =
            get_rt_cal_spline(&donor, &acceptor, DONOR_Q_VALUE_THRESHOLD, DonorCriterion::Score);
        assert!(spline.calibration_curve.is_empty());
    }

    #[test]
    fn rt_info_hypothesis_bounds() {
        let info = RtInfo::new(50.0, 0.5);
        assert_eq!(info.rt_start_hypothesis(), 49.75);
        assert_eq!(info.rt_end_hypothesis(), 50.25);
    }

    #[test]
    fn median_matches_mathnet() {
        assert!(median(&[]).is_nan());
        assert_eq!(median(&[5.0]), 5.0);
        // odd count -> middle order statistic.
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        // even count -> average of the two central order statistics.
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
    }

    #[test]
    fn standard_deviation_is_sample_sd() {
        assert!(standard_deviation(&[]).is_nan());
        assert!(standard_deviation(&[7.0]).is_nan());
        // Sample SD of [2,4,4,4,5,5,7,9] is 2.13808993...; population SD would be 2.0.
        let sd = standard_deviation(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]);
        assert!((sd - 2.138089935299395).abs() < 1e-12, "sd = {sd}");
        // Two identical samples -> zero spread.
        assert_eq!(standard_deviation(&[3.0, 3.0]), 0.0);
    }

    /// Build a calibration point pairing a donor peak at `donor_rt` with an acceptor peak at
    /// `acceptor_rt` (both traced, single id).
    fn cal_point(donor_rt: f64, acceptor_rt: f64) -> RetentionTimeCalibDataPoint {
        RetentionTimeCalibDataPoint::new(
            peak("D", 10.0, 0.001, donor_rt, 100.0),
            Some(peak("A", 10.0, 0.001, acceptor_rt, 100.0)),
        )
    }

    #[test]
    fn predict_rt_no_nearby_points_returns_donor_rt() {
        let curve: Vec<RetentionTimeCalibDataPoint> = Vec::new();
        let donor = peak("X", 10.0, 0.001, 42.0, 100.0);
        let info = predict_retention_time(&curve, &donor, MAX_MBR_RT_WINDOW, NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR);
        assert_eq!(info.predicted_rt, 42.0);
        assert_eq!(info.width, 0.25);
    }

    #[test]
    fn predict_rt_single_anchor_shifts_by_median_diff() {
        // Donor at 30.0 sits below both curve points, so the binary search insertion index is 0.
        // The element AT that index (30.2) is the pivot and is excluded from both directions
        // (faithful to .NET BinarySearch + the C# index+1/index-1 gathers); forward gathering
        // starts at index 1, so only the 30.3 point is collected -> a single anchor.
        let curve = vec![cal_point(30.2, 30.0), cal_point(30.3, 30.1)];
        let donor = peak("X", 10.0, 0.001, 30.0, 100.0);
        let info = predict_retention_time(&curve, &donor, MAX_MBR_RT_WINDOW, NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR);
        // single anchor: predicted = donorRt - (donorRt_anchor - acceptorRt_anchor) = 30.0 - 0.2.
        // RTs round-trip through f32 (envelope storage), so allow ~1e-4 slack.
        assert!((info.predicted_rt - 29.8).abs() < 1e-4, "predicted = {}", info.predicted_rt);
        assert_eq!(info.width, 0.25);
    }

    #[test]
    fn predict_rt_multi_anchor_uses_stddev_width() {
        // Donor at 30.0 below all points -> insertion index 0, pivot 30.2 excluded, forward
        // gathers 30.3 (diff 0.2) and 30.4 (diff 0.3).
        let curve = vec![
            cal_point(30.2, 30.0), // pivot (index 0) -> skipped
            cal_point(30.3, 30.1), // forward anchor, donor-acceptor diff 0.2
            cal_point(30.4, 30.1), // forward anchor, donor-acceptor diff 0.3
        ];
        let donor = peak("X", 10.0, 0.001, 30.0, 100.0);
        let info = predict_retention_time(&curve, &donor, MAX_MBR_RT_WINDOW, NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR);
        // median of {0.2, 0.3} = 0.25 -> predicted = 30.0 - 0.25 = 29.75. (RTs widen from f32.)
        assert!((info.predicted_rt - 29.75).abs() < 1e-4, "predicted = {}", info.predicted_rt);
        // width = min(6 * sampleSD, MAX); sampleSD of {0.2,0.3} = sqrt(0.005) -> 6*… ≈ 0.4242641.
        assert!((info.width - 0.424264068711929).abs() < 1e-4, "width = {}", info.width);
    }

    #[test]
    fn predict_rt_clamps_width_to_max_window() {
        // Wide spread of donor-acceptor diffs -> 6*stddev exceeds the 1.0-min cap.
        let curve = vec![
            cal_point(30.2, 30.0), // pivot (index 0) -> skipped
            cal_point(30.3, 29.9), // forward anchor, diff 0.4
            cal_point(30.4, 30.35), // forward anchor, diff 0.05
        ];
        let donor = peak("X", 10.0, 0.001, 30.0, 100.0);
        let info = predict_retention_time(&curve, &donor, MAX_MBR_RT_WINDOW, NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR);
        assert_eq!(info.width, MAX_MBR_RT_WINDOW);
    }
}
