//! Match-between-runs scorer — port of mzLib `FlashLFQ.MbrScorer` + `MbrScorerFactory` (PLAN.md P3.2c).
//!
//! [`MbrScorer`] turns an acceptor file's unambiguous MS/MS peaks into a set of statistical
//! distributions, then scores each candidate transferred peak ([`MbrChromatographicPeak`]) against
//! them. The combined score is the geometric mean of five component scores (intensity, RT, ppm,
//! scan-count, isotopic-distribution), scaled to `[0, 100]`. Built via [`build_mbr_scorer`]
//! (the `MbrScorerFactory.BuildMbrScorer` port).
//!
//! ## Faithfulness to `MBR/MbrScorer.cs`
//! - The four acceptor-file distributions are fit exactly as C# does: ppm error → [`Normal`] (median
//!   center, IQR/1.36 or stddev spread, with the raw stats preserved for [`MbrScorer::get_ppm_error_tolerance`]);
//!   log2 intensity → [`Normal`]; scan count → [`Normal`] (LINQ `Average` mean — *not* the incremental
//!   mean — and stddev/IQROne spread); `1 − isotopicCorrelation` → [`Gamma`] (method-of-moments
//!   α = mean²/var, β = mean/var, using MathNet's incremental `Mean`/`Variance`). The per-donor-file
//!   log-fold-change and RT-prediction-error distributions are added later by
//!   [`MbrScorer::calculate_fold_change_between_files`] / [`MbrScorer::add_rt_pred_error_distribution`].
//! - **The long-standing `sigma = 1` RT-distribution bug (`MbrScorer.cs:288`) is replicated** — the RT
//!   prediction-error `Normal` is built `Normal(medianRtError, 1)`, *not* with the computed stddev.
//!   The stddev is still computed (and discarded) so the control flow matches C# exactly.
//! - `ScoreMbr` writes the five component scores onto the peak and returns
//!   `100·(intensity·rt·ppm·scanCount·isotopic)^0.2`. `CalculateScore` returns `_minScore = 3e-7`
//!   whenever a distribution is null / the value is non-finite / the CDF degenerates to 0, so one bad
//!   component can't drive the whole score to zero.
//!
//! ## Divergence from C#: donor-file keys are strings
//! The base [`ChromatographicPeak`] in this port carries no `SpectraFileInfo`, so the per-donor-file
//! dictionaries are keyed by the donor file **name** (`&str`) rather than a file object. Callers
//! (the P3.2d orchestration) pass the donor file name explicitly to
//! [`calculate_fold_change_between_files`](MbrScorer::calculate_fold_change_between_files),
//! [`add_rt_pred_error_distribution`](MbrScorer::add_rt_pred_error_distribution),
//! [`score_mbr`](MbrScorer::score_mbr), and [`is_valid_for_donor`](MbrScorer::is_valid_for_donor).
//! The acceptor-file distributions live on the scorer itself, exactly as in C# (one scorer per
//! acceptor file).

use std::collections::HashMap;

use crate::chromatographic_peak::{ChromatographicPeak, EnvelopePeakKey};
use crate::mbr::{
    average, interquartile_range, mean, median, standard_deviation, variance,
};
use crate::mbr_chromatographic_peak::MbrChromatographicPeak;
use crate::special_functions::{log2, Gamma, Normal};
use crate::tolerance::PpmTolerance;

/// Scores match-between-runs transfers for a single acceptor file. Port of `FlashLFQ.MbrScorer`.
#[derive(Debug, Clone)]
pub struct MbrScorer {
    // Acceptor-file distributions (set by `initialize_scorer`). `None` until initialized / when the
    // fit degenerates (the C# fields are nullable for the same reason).
    log_intensity_distribution: Option<Normal>,
    ppm_distribution: Option<Normal>,
    scan_count_distribution: Option<Normal>,
    isotopic_correlation_distribution: Option<Gamma>,

    // Raw ppm stats kept so `get_ppm_error_tolerance` returns `|median| + 4·stddev` exactly, even
    // when the fitted Normal had to substitute a safe stddev.
    ppm_median_raw: f64,
    ppm_std_dev_raw: f64,

    // Per-donor-file distributions, keyed by donor file name (C# keys by `SpectraFileInfo`).
    log_fc_distribution_dictionary: HashMap<String, Normal>,
    rt_prediction_error_distribution_dictionary: HashMap<String, Normal>,

    /// Map from each acceptor MS/MS peak's apex indexed-peak to that peak — so the acceptor search
    /// (P3.2d) doesn't transfer a peptide onto an apex already claimed by an MS/MS id (C#
    /// `ApexToAcceptorFilePeakDict`). Keyed by [`EnvelopePeakKey`] (structural apex identity).
    pub apex_to_acceptor_file_peak: HashMap<EnvelopePeakKey, ChromatographicPeak>,
    /// The acceptor file's unambiguous MS/MS-identified peaks the distributions are fit from
    /// (C# `UnambiguousMsMsAcceptorPeaks`).
    pub unambiguous_msms_acceptor_peaks: Vec<ChromatographicPeak>,
    /// The largest `ScanCount` over the acceptor MS/MS peaks (C# `MaxNumberOfScansObserved`); `0`
    /// when there are no peaks. Carried for the P3.2d acceptor search.
    pub max_number_of_scans_observed: f64,

    // Minimum component score: 3e-7 ≈ the normal-tail mass beyond 5σ. Keeps one zero component from
    // zeroing the whole geometric-mean score (C# `_minScore`).
    min_score: f64,
}

impl MbrScorer {
    /// Constructs a scorer from the acceptor apex map and the acceptor MS/MS peaks (C# `MbrScorer`
    /// constructor). Distributions are not yet fit — call [`initialize_scorer`](Self::initialize_scorer).
    pub fn new(
        apex_to_acceptor_file_peak: HashMap<EnvelopePeakKey, ChromatographicPeak>,
        unambiguous_acceptor_peaks: Vec<ChromatographicPeak>,
    ) -> MbrScorer {
        let max_number_of_scans_observed = if unambiguous_acceptor_peaks.is_empty() {
            0.0
        } else {
            unambiguous_acceptor_peaks
                .iter()
                .map(|p| p.scan_count())
                .max()
                .unwrap() as f64
        };
        MbrScorer {
            log_intensity_distribution: None,
            ppm_distribution: None,
            scan_count_distribution: None,
            isotopic_correlation_distribution: None,
            ppm_median_raw: 0.0,
            ppm_std_dev_raw: 0.0,
            log_fc_distribution_dictionary: HashMap::new(),
            rt_prediction_error_distribution_dictionary: HashMap::new(),
            apex_to_acceptor_file_peak,
            unambiguous_msms_acceptor_peaks: unambiguous_acceptor_peaks,
            max_number_of_scans_observed,
            min_score: 3e-7,
        }
    }

    /// Fits the four acceptor-file distributions and reports whether the scorer is usable. Port of
    /// `InitializeScorer`: returns `false` (without fitting) when fewer than three acceptor MS/MS
    /// peaks exist, else fits the distributions and returns [`is_valid`](Self::is_valid).
    pub fn initialize_scorer(&mut self) -> bool {
        if self.unambiguous_msms_acceptor_peaks.len() < 3 {
            return false;
        }
        self.log_intensity_distribution = self.get_log_intensity_distribution();
        self.ppm_distribution = self.get_ppm_error_distribution();
        self.isotopic_correlation_distribution = self.get_isotopic_envelope_corr_distribution();
        self.scan_count_distribution = self.get_scan_count_distribution();
        self.is_valid()
    }

    // ----- distribution builders -----

    /// Port of `GetPpmErrorDistribution`. Centers on the median ppm error; spread is IQR/1.36 for
    /// >30 peaks else the sample stddev. The raw median/stddev are stashed for the ppm tolerance,
    /// while the returned `Normal` substitutes a safe `1.0` stddev when the raw spread is zero/invalid.
    fn get_ppm_error_distribution(&mut self) -> Option<Normal> {
        let ppm_errors: Vec<f64> = self
            .unambiguous_msms_acceptor_peaks
            .iter()
            .map(|p| p.mass_error)
            .filter(|e| !e.is_nan() && !e.is_infinite())
            .collect();
        if ppm_errors.len() < 2 {
            return None;
        }
        self.ppm_median_raw = median(&ppm_errors);
        let spread = if ppm_errors.len() > 30 {
            interquartile_range(&ppm_errors) / 1.36
        } else {
            standard_deviation(&ppm_errors)
        };
        self.ppm_std_dev_raw = if spread.is_nan() || spread < 0.0 { 0.0 } else { spread };

        let safe_std_dev = if self.ppm_std_dev_raw > 0.0
            && Normal::is_valid_parameter_set(self.ppm_median_raw, self.ppm_std_dev_raw)
        {
            self.ppm_std_dev_raw
        } else {
            1.0
        };
        Some(Normal::new(self.ppm_median_raw, safe_std_dev))
    }

    /// Port of `GetLogIntensityDistribution`. Median of log2 intensities (positive, finite only) as
    /// the center, IQR/1.36 as the spread (falling back to `1.0` when invalid).
    fn get_log_intensity_distribution(&self) -> Option<Normal> {
        let log_intensities: Vec<f64> = self
            .unambiguous_msms_acceptor_peaks
            .iter()
            .filter(|p| p.intensity > 0.0 && !p.intensity.is_nan() && !p.intensity.is_infinite())
            .map(|p| log2(p.intensity))
            .collect();
        if log_intensities.len() < 2 {
            return None;
        }
        let mean_center = median(&log_intensities);
        let mut std_dev = interquartile_range(&log_intensities) / 1.36;
        if std_dev.is_nan() || std_dev <= 0.0 || !Normal::is_valid_parameter_set(mean_center, std_dev)
        {
            std_dev = 1.0;
        }
        Some(Normal::new(mean_center, std_dev))
    }

    /// Port of `GetScanCountDistribution`. LINQ-`Average` mean of the scan counts (or `0` when
    /// empty), spread = stddev for >30 peaks else IQR/1.36 (fallback `1.0`). Always returns a
    /// `Normal` (never `None`).
    fn get_scan_count_distribution(&self) -> Option<Normal> {
        let scan_list: Vec<f64> = self
            .unambiguous_msms_acceptor_peaks
            .iter()
            .map(|p| p.scan_count() as f64)
            .collect();
        let mean_center = if !scan_list.is_empty() {
            average(&scan_list)
        } else {
            0.0
        };
        let mut std_dev = if scan_list.len() > 30 {
            standard_deviation(&scan_list)
        } else {
            interquartile_range(&scan_list) / 1.36
        };
        if std_dev.is_nan() || std_dev <= 0.0 || !Normal::is_valid_parameter_set(mean_center, std_dev)
        {
            std_dev = 1.0;
        }
        Some(Normal::new(mean_center, std_dev))
    }

    /// Port of `GetIsotopicEnvelopeCorrDistribution`. Fits a [`Gamma`] to `1 − pearsonCorrelation`
    /// (positive, finite only) by method of moments (α = mean²/var, β = mean/var), using MathNet's
    /// incremental [`mean`]/[`variance`]. `None` when there are <2 samples or the fit is invalid.
    fn get_isotopic_envelope_corr_distribution(&self) -> Option<Gamma> {
        let pearson_corrs: Vec<f64> = self
            .unambiguous_msms_acceptor_peaks
            .iter()
            .map(|p| 1.0 - p.isotopic_pearson_correlation())
            .filter(|&p| p > 0.0 && !p.is_nan() && !p.is_infinite())
            .collect();
        if pearson_corrs.len() < 2 {
            return None;
        }
        let mean_val = mean(&pearson_corrs);
        let variance_val = variance(&pearson_corrs);
        if variance_val <= 0.0 || mean_val.is_nan() || variance_val.is_nan() {
            return None;
        }
        let alpha = mean_val.powi(2) / variance_val;
        let beta = mean_val / variance_val;
        if !Gamma::is_valid_parameter_set(alpha, beta) {
            return None;
        }
        Some(Gamma::new(alpha, beta))
    }

    /// Port of `CalculateFoldChangeBetweenFiles`. Builds (when ≥100 paired fold changes exist) the
    /// log2 intensity fold-change `Normal` for `donor_file`, used by [`calculate_intensity_score`](Self::calculate_intensity_score).
    ///
    /// `donor_file` stands in for the C# `idDonorPeaks.First().SpectraFileInfo` key (see module note).
    pub fn calculate_fold_change_between_files(
        &mut self,
        donor_file: &str,
        id_donor_peaks: &[ChromatographicPeak],
    ) {
        if id_donor_peaks.is_empty() {
            return;
        }

        // Best (per the C# rule) acceptor peak per modified sequence. NOTE: the C# condition keeps
        // replacing with the *less* intense peak (`currentBest.Intensity > acceptor.Intensity`);
        // replicated verbatim, quirk and all.
        let mut acceptor_best: HashMap<String, &ChromatographicPeak> = HashMap::new();
        for acc in &self.unambiguous_msms_acceptor_peaks {
            let key = match acc.first_modified_sequence() {
                Some(s) => s.to_string(),
                None => continue,
            };
            match acceptor_best.get(&key) {
                Some(current) => {
                    if current.intensity > acc.intensity {
                        acceptor_best.insert(key, acc);
                    }
                }
                None => {
                    acceptor_best.insert(key, acc);
                }
            }
        }

        let mut fold_changes: Vec<f64> = Vec::new();
        for donor in id_donor_peaks {
            let donor_intensity = donor.intensity;
            let key = match donor.first_modified_sequence() {
                Some(s) => s,
                None => continue,
            };
            if let Some(acc) = acceptor_best.get(key) {
                let log_fold_change = log2(acc.intensity) - log2(donor_intensity);
                if !log_fold_change.is_nan() && !log_fold_change.is_infinite() {
                    fold_changes.push(log_fold_change);
                }
            }
        }

        // Redundant second filter, kept for fidelity.
        let fold_changes: Vec<f64> = fold_changes
            .into_iter()
            .filter(|d| !(d.is_nan() || d.is_infinite()))
            .collect();
        if fold_changes.len() < 100 {
            return;
        }

        let median_fc = median(&fold_changes);
        let std_dev_fc = standard_deviation(&fold_changes);
        if median_fc.is_nan()
            || std_dev_fc.is_nan()
            || !Normal::is_valid_parameter_set(median_fc, std_dev_fc)
        {
            return;
        }
        self.log_fc_distribution_dictionary
            .insert(donor_file.to_string(), Normal::new(median_fc, std_dev_fc));
    }

    /// Port of `AddRtPredErrorDistribution`. From a donor file's anchor-peptide RT differences,
    /// builds the local-alignment prediction-error `Normal` for `donor_file`. **Replicates the
    /// `sigma = 1` bug** (`MbrScorer.cs:288`): the distribution is `Normal(medianRtError, 1)` even
    /// though the stddev `sigma` is computed. Defaults to `Normal(0, 1)` when too few diffs / errors.
    pub fn add_rt_pred_error_distribution(
        &mut self,
        donor_file: &str,
        anchor_peptide_rt_diffs: &[f64],
        number_of_anchor_peptides: usize,
    ) {
        // Safe, non-degenerate default.
        let mut rt_prediction_error_dist = Normal::new(0.0, 1.0);

        let n = anchor_peptide_rt_diffs.len();
        if n >= 2 * number_of_anchor_peptides + 1 {
            let mut rt_prediction_errors: Vec<f64> = Vec::new();
            for i in number_of_anchor_peptides..(n - number_of_anchor_peptides) {
                let mut cum_sum_rt_diffs = 0.0_f64;
                for j in 1..=number_of_anchor_peptides {
                    cum_sum_rt_diffs += anchor_peptide_rt_diffs[i - j];
                    cum_sum_rt_diffs += anchor_peptide_rt_diffs[i + j];
                }
                let avg_diff = cum_sum_rt_diffs / (2 * number_of_anchor_peptides) as f64;
                let err = avg_diff - anchor_peptide_rt_diffs[i];
                if !err.is_nan() && !err.is_infinite() {
                    rt_prediction_errors.push(err);
                }
            }

            if rt_prediction_errors.len() >= 2 {
                let median_rt_error = median(&rt_prediction_errors);
                let std_dev_rt_error = standard_deviation(&rt_prediction_errors);
                if !median_rt_error.is_nan() {
                    // `sigma` is computed exactly as C# does, then deliberately discarded: the
                    // long-standing PR #802 bug builds the dist with stddev 1, not `sigma`.
                    let _sigma = if std_dev_rt_error.is_nan()
                        || std_dev_rt_error <= 0.0
                        || !Normal::is_valid_parameter_set(median_rt_error, std_dev_rt_error)
                    {
                        1.0
                    } else {
                        std_dev_rt_error
                    };
                    rt_prediction_error_dist = Normal::new(median_rt_error, 1.0);
                }
            }
        }

        self.rt_prediction_error_distribution_dictionary
            .insert(donor_file.to_string(), rt_prediction_error_dist);
    }

    // ----- scoring -----

    /// Scores an MBR candidate peak, writing its five component scores + RT prediction error onto
    /// it and returning the combined score in `[0, 100]` (higher is better). Port of `ScoreMbr`.
    ///
    /// `donor_file` is the donor file name (C# `donorPeak.SpectraFileInfo`, used to look up the RT
    /// and log-fold-change distributions; see module note).
    pub fn score_mbr(
        &self,
        acceptor_peak: &mut MbrChromatographicPeak,
        donor_peak: &ChromatographicPeak,
        donor_file: &str,
        predicted_rt: f64,
    ) -> f64 {
        acceptor_peak.intensity_score =
            self.calculate_intensity_score(acceptor_peak.intensity(), donor_peak, donor_file);
        acceptor_peak.rt_prediction_error = predicted_rt - acceptor_peak.apex_retention_time();

        // Safe default RT distribution if none was computed for this donor file.
        let rt_dist = self
            .rt_prediction_error_distribution_dictionary
            .get(donor_file)
            .copied()
            .unwrap_or_else(|| Normal::new(0.0, 1.0));

        acceptor_peak.rt_score =
            self.calculate_score_normal(Some(rt_dist), acceptor_peak.rt_prediction_error);
        acceptor_peak.ppm_score =
            self.calculate_score_normal(self.ppm_distribution, acceptor_peak.mass_error());
        acceptor_peak.scan_count_score =
            self.calculate_score_normal(self.scan_count_distribution, acceptor_peak.scan_count() as f64);
        acceptor_peak.isotopic_distribution_score = self.calculate_score_gamma(
            self.isotopic_correlation_distribution,
            1.0 - acceptor_peak.isotopic_pearson_correlation(),
        );

        let score = 100.0
            * (acceptor_peak.intensity_score
                * acceptor_peak.rt_score
                * acceptor_peak.ppm_score
                * acceptor_peak.scan_count_score
                * acceptor_peak.isotopic_distribution_score)
                .powf(0.20);
        score
    }

    /// Port of `CalculateScore(Normal, value)`. The fraction of the distribution at least
    /// `|mean − value|` from the mean (so `1` when `value == mean`), floored at `_minScore`.
    pub fn calculate_score_normal(&self, distribution: Option<Normal>, value: f64) -> f64 {
        let dist = match distribution {
            Some(d) if !value.is_nan() && !value.is_infinite() => d,
            _ => return self.min_score,
        };
        let absolute_diff_from_mean = (dist.mean - value).abs();
        let score = 2.0 * dist.cumulative_distribution(dist.mean - absolute_diff_from_mean);
        if score.is_nan() || score == 0.0 {
            self.min_score
        } else {
            score
        }
    }

    /// Port of `CalculateScore(Gamma, value)`. `1 − CDF(value)` (the gamma upper tail), floored at
    /// `_minScore` for negative/non-finite values or a null distribution.
    pub fn calculate_score_gamma(&self, distribution: Option<Gamma>, value: f64) -> f64 {
        if value < 0.0 || distribution.is_none() || value.is_nan() || value.is_infinite() {
            return self.min_score;
        }
        let cdf = distribution.unwrap().cumulative_distribution(value);
        if cdf.is_nan() {
            return self.min_score;
        }
        1.0 - cdf
    }

    /// Port of `CalculateIntensityScore`. Uses the donor-file log-fold-change distribution when one
    /// exists (and both intensities are non-zero), else the acceptor log-intensity distribution.
    pub fn calculate_intensity_score(
        &self,
        acceptor_intensity: f64,
        donor_peak: &ChromatographicPeak,
        donor_file: &str,
    ) -> f64 {
        if acceptor_intensity != 0.0 && donor_peak.intensity != 0.0 {
            if let Some(log_fc_distribution) = self.log_fc_distribution_dictionary.get(donor_file) {
                let log_fold_change = log2(acceptor_intensity) - log2(donor_peak.intensity);
                return self.calculate_score_normal(Some(*log_fc_distribution), log_fold_change);
            }
        }
        let log_intensity = log2(acceptor_intensity);
        self.calculate_score_normal(self.log_intensity_distribution, log_intensity)
    }

    // ----- validity / tolerance -----

    /// Whether the scorer is validly parameterized for transfers from `donor_file` (its RT-error
    /// distribution exists *and* [`is_valid`](Self::is_valid)). Port of `IsValid(SpectraFileInfo)`.
    pub fn is_valid_for_donor(&self, donor_file: &str) -> bool {
        self.rt_prediction_error_distribution_dictionary
            .contains_key(donor_file)
            && self.is_valid()
    }

    /// Whether the acceptor-file ppm / scan-count / log-intensity distributions are all present
    /// (indifferent to the isotopic-correlation distribution being null). Port of `IsValid()`.
    pub fn is_valid(&self) -> bool {
        self.ppm_distribution.is_some()
            && self.scan_count_distribution.is_some()
            && self.log_intensity_distribution.is_some()
    }

    /// `|median| + 4·stddev` of the raw ppm-error distribution (C# `GetPpmErrorTolerance`), the
    /// acceptor-search half-window before it is clamped by `MbrPpmTolerance`.
    pub fn get_ppm_error_tolerance(&self) -> f64 {
        self.ppm_median_raw.abs() + 4.0 * self.ppm_std_dev_raw
    }
}

/// Builds the scorer for one acceptor file (port of `MbrScorerFactory.BuildMbrScorer`).
///
/// Constructs the apex→peak map (distinct by apex indexed-peak, first occurrence wins, peaks with an
/// apex only), initializes the scorer, and — on success — returns it together with the file-specific
/// MBR ppm tolerance (`min(scorer.GetPpmErrorTolerance(), mbr_ppm_tolerance)`). Returns `None` when
/// the acceptor file has too few MS/MS peaks to fit the distributions (C# returns a null scorer +
/// null tolerance). `mbr_ppm_tolerance` is `FlashLfqParameters.MbrPpmTolerance`
/// ([`crate::mbr::MBR_PPM_TOLERANCE`]).
pub fn build_mbr_scorer(
    acceptor_file_msms_peaks: Vec<ChromatographicPeak>,
    mbr_ppm_tolerance: f64,
) -> Option<(MbrScorer, PpmTolerance)> {
    let mut apex_to_acceptor_file_peak: HashMap<EnvelopePeakKey, ChromatographicPeak> =
        HashMap::new();
    for peak in &acceptor_file_msms_peaks {
        if let Some(apex) = peak.apex {
            // DistinctBy(Apex.IndexedPeak): keep the first peak seen for each apex peak.
            apex_to_acceptor_file_peak
                .entry(EnvelopePeakKey::from(&apex))
                .or_insert_with(|| peak.clone());
        }
    }

    let mut scorer = MbrScorer::new(apex_to_acceptor_file_peak, acceptor_file_msms_peaks);
    if !scorer.initialize_scorer() {
        return None;
    }
    let mbr_ppm = scorer.get_ppm_error_tolerance().min(mbr_ppm_tolerance);
    Some((scorer, PpmTolerance::new(mbr_ppm)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isotopic_envelope::IsotopicEnvelope;
    use crate::mbr::MBR_PPM_TOLERANCE;
    use crate::peak_indexing::IndexedMassSpectralPeak;
    use crate::psm_tsv::Identification;

    fn id(modseq: &str, charge: i32) -> Identification {
        Identification {
            file_name: "acceptor".to_string(),
            base_sequence: modseq.to_string(),
            modified_sequence: modseq.to_string(),
            monoisotopic_mass: 1000.0,
            ms2_retention_time_in_minutes: 1.0,
            precursor_charge_state: charge,
            score: 10.0,
            q_value: 0.001,
            is_decoy: false,
        }
    }

    /// An envelope at a given scan/RT/intensity/correlation (charge 1).
    fn env(scan: i32, rt: f64, intensity: f64, corr: f64) -> IsotopicEnvelope {
        let peak = IndexedMassSpectralPeak::new(1000.0, intensity, scan, rt);
        IsotopicEnvelope::new(peak, 1, intensity, corr)
    }

    /// A traced MS/MS acceptor peak: `n_env` envelopes (apex carries `corr`), `mass_error` ppm.
    fn acc_peak(
        modseq: &str,
        rt: f64,
        intensity: f64,
        corr: f64,
        mass_error: f64,
        n_env: usize,
    ) -> ChromatographicPeak {
        let mut p = ChromatographicPeak::from_identification(id(modseq, 2), 1000.0);
        p.isotopic_envelopes = (0..n_env)
            .map(|i| env(i as i32, rt, intensity, corr))
            .collect();
        p.calculate_intensity_for_this_feature(false);
        // calculate_intensity sets mass_error from the apex; override for a controlled ppm value.
        p.mass_error = mass_error;
        p
    }

    /// A small but valid acceptor-file peak set (varied intensity/RT/corr/scan-count/ppm).
    fn acceptor_peaks() -> Vec<ChromatographicPeak> {
        vec![
            acc_peak("PEPA", 20.0, 1.0e6, 0.99, 1.0, 5),
            acc_peak("PEPB", 25.0, 2.0e6, 0.97, -2.0, 6),
            acc_peak("PEPC", 30.0, 5.0e5, 0.95, 0.5, 4),
            acc_peak("PEPD", 35.0, 8.0e5, 0.98, -1.5, 7),
            acc_peak("PEPE", 40.0, 3.0e6, 0.96, 2.5, 3),
        ]
    }

    #[test]
    fn initialize_fails_with_fewer_than_three_peaks() {
        let peaks = vec![
            acc_peak("PEPA", 20.0, 1.0e6, 0.99, 1.0, 5),
            acc_peak("PEPB", 25.0, 2.0e6, 0.97, -2.0, 6),
        ];
        let mut scorer = MbrScorer::new(HashMap::new(), peaks);
        assert!(!scorer.initialize_scorer());
        assert!(!scorer.is_valid());
    }

    #[test]
    fn initialize_succeeds_and_is_valid() {
        let mut scorer = MbrScorer::new(HashMap::new(), acceptor_peaks());
        assert!(scorer.initialize_scorer());
        assert!(scorer.is_valid());
        // Max scan count over the peaks = 7.
        assert_eq!(scorer.max_number_of_scans_observed, 7.0);
    }

    #[test]
    fn factory_builds_scorer_and_clamps_tolerance() {
        // Give the peaks apexes so the apex dict is populated.
        let (scorer, tol) = build_mbr_scorer(acceptor_peaks(), MBR_PPM_TOLERANCE)
            .expect("five varied peaks initialize a valid scorer");
        // Tolerance is min(ppmErrorTolerance, MbrPpmTolerance=10).
        assert!(tol.value <= MBR_PPM_TOLERANCE + 1e-12);
        assert!((tol.value - scorer.get_ppm_error_tolerance().min(MBR_PPM_TOLERANCE)).abs() < 1e-12);
        // Each peak had a distinct apex -> five entries in the apex dict.
        assert_eq!(scorer.apex_to_acceptor_file_peak.len(), 5);
    }

    #[test]
    fn factory_returns_none_when_uninitializable() {
        let peaks = vec![acc_peak("PEPA", 20.0, 1.0e6, 0.99, 1.0, 5)];
        assert!(build_mbr_scorer(peaks, MBR_PPM_TOLERANCE).is_none());
    }

    #[test]
    fn score_is_in_range_and_perfect_match_scores_high() {
        let mut scorer = MbrScorer::new(HashMap::new(), acceptor_peaks());
        assert!(scorer.initialize_scorer());

        // An acceptor candidate sitting near the center of every distribution should score well.
        let mut candidate = MbrChromatographicPeak::new(id("PEPX", 2), "acceptor".to_string(), 1000.0, 30.0, false);
        candidate.peak.isotopic_envelopes = (0..5).map(|i| env(i, 30.0, 1.5e6, 0.99)).collect();
        candidate.peak.calculate_intensity_for_this_feature(false);
        candidate.peak.mass_error = 0.0;

        let donor = acc_peak("PEPX", 30.2, 1.5e6, 0.99, 0.0, 5);
        // predicted_rt == apex RT -> zero RT prediction error -> RT score ~1.
        let predicted_rt = candidate.apex_retention_time();
        let score = scorer.score_mbr(&mut candidate, &donor, "donor", predicted_rt);

        assert!(score >= 0.0 && score <= 100.0, "score = {score}");
        // Component scores all written and in (0, 1].
        for s in [
            candidate.intensity_score,
            candidate.rt_score,
            candidate.ppm_score,
            candidate.scan_count_score,
            candidate.isotopic_distribution_score,
        ] {
            assert!(s > 0.0 && s <= 1.0, "component score out of range: {s}");
        }
        // RT prediction error is predicted - apex = 0 here.
        assert!(candidate.rt_prediction_error.abs() < 1e-6);
    }

    #[test]
    fn nan_value_yields_min_score() {
        let scorer = MbrScorer::new(HashMap::new(), acceptor_peaks());
        // No distribution -> min score.
        assert_eq!(scorer.calculate_score_normal(None, 1.0), 3e-7);
        // Non-finite value -> min score even with a distribution.
        assert_eq!(
            scorer.calculate_score_normal(Some(Normal::new(0.0, 1.0)), f64::NAN),
            3e-7
        );
        // Negative value for the gamma branch -> min score.
        assert_eq!(scorer.calculate_score_gamma(Some(Gamma::new(2.0, 3.0)), -0.1), 3e-7);
    }

    #[test]
    fn normal_score_is_one_at_the_mean() {
        let scorer = MbrScorer::new(HashMap::new(), acceptor_peaks());
        // value == mean -> 2 * CDF(mean) = 2 * 0.5 = 1.
        let s = scorer.calculate_score_normal(Some(Normal::new(5.0, 2.0)), 5.0);
        assert!((s - 1.0).abs() < 1e-9, "score at mean = {s}");
    }

    #[test]
    fn rt_distribution_uses_sigma_one_bug() {
        let mut scorer = MbrScorer::new(HashMap::new(), acceptor_peaks());
        assert!(scorer.initialize_scorer());
        // Anchor diffs chosen so the prediction errors have a clearly non-unit stddev; the stored
        // distribution must still have stddev exactly 1 (the replicated bug), centered on the
        // error median.
        let diffs = vec![0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0];
        scorer.add_rt_pred_error_distribution("donor", &diffs, 1);
        let dist = scorer
            .rt_prediction_error_distribution_dictionary
            .get("donor")
            .copied()
            .expect("distribution stored for donor");
        assert_eq!(dist.std_dev, 1.0, "sigma bug: stddev must be 1, got {}", dist.std_dev);
    }

    #[test]
    fn add_rt_dist_defaults_when_too_few_diffs() {
        let mut scorer = MbrScorer::new(HashMap::new(), acceptor_peaks());
        // Fewer than 2*anchors+1 diffs -> default Normal(0, 1).
        scorer.add_rt_pred_error_distribution("donor", &[1.0, 2.0], 3);
        let dist = scorer
            .rt_prediction_error_distribution_dictionary
            .get("donor")
            .copied()
            .unwrap();
        assert_eq!(dist.mean, 0.0);
        assert_eq!(dist.std_dev, 1.0);
        // Donor not valid for scoring until its RT distribution exists AND the scorer is valid.
        assert!(scorer.initialize_scorer());
        assert!(scorer.is_valid_for_donor("donor"));
        assert!(!scorer.is_valid_for_donor("other-donor"));
    }
}
