//! Match-between-runs peak — port of mzLib `FlashLFQ.MbrChromatographicPeak` (PLAN.md P3.2b).
//!
//! In C# `MbrChromatographicPeak : ChromatographicPeak` *extends* the base peak with the MBR
//! scoring fields. Rust has no inheritance, so [`MbrChromatographicPeak`] **composes** a base
//! [`ChromatographicPeak`] (its `detection_type` is always [`DetectionType::MBR`]) and adds the
//! MBR-exclusive fields: the predicted acceptor RT, the `random_rt` decoy flag, the MBR q-value /
//! PEP, the five component scores the [`crate::mbr`] scorer (P3.2c) writes, the combined
//! [`MbrChromatographicPeak::mbr_score`], and the [`MbrChromatographicPeak::charge_list`] the
//! acceptor search (P3.2d) sets.
//!
//! ## Faithfulness to `MBR/MbrChromatographicPeak.cs`
//! - The constructor mirrors `MbrChromatographicPeak(Identification, SpectraFileInfo, double, bool)`:
//!   it builds the base peak from the single identification with `DetectionType.MBR`, stores the
//!   predicted RT (C# `init`) and the `random_rt` flag (C# read-only `RandomRt`). All scores start
//!   at `0.0`, `mbr_q_value` at `0.0`, `mbr_pep` at `None`, `charge_list` empty — the same default
//!   state a freshly-constructed C# peak has before the scorer/search populate it.
//! - `RandomRt == true` marks an MBR **decoy** peak (the C# doc-comment: "implies that this peak is
//!   a decoy peak identified by the MBR algorithm"); exposed as [`MbrChromatographicPeak::is_decoy_peak`].
//! - The base peak carries no `SpectraFileInfo` in this port, so the acceptor file is tracked as
//!   [`MbrChromatographicPeak::spectra_file_name`] (the extension-stripped file name) — the key the
//!   P3.2d MBR loop groups acceptor peaks by, standing in for the C# `SpectraFileInfo` reference.
//! - The scorer reads `Intensity`, `ApexRetentionTime`, `MassError`, `ScanCount`, and
//!   `IsotopicPearsonCorrelation` off the (base) peak; those are reached through
//!   [`MbrChromatographicPeak::peak`] (or the delegating accessors here).

use crate::chromatographic_peak::ChromatographicPeak;
use crate::detection_type::DetectionType;
use crate::psm_tsv::Identification;

/// An MBR-transferred chromatographic peak: a base [`ChromatographicPeak`] (detection type
/// [`DetectionType::MBR`]) plus the match-between-runs scoring fields. Port of
/// `FlashLFQ.MbrChromatographicPeak`.
#[derive(Debug, Clone)]
pub struct MbrChromatographicPeak {
    /// The underlying chromatographic peak (C# base class). Holds the identifications, isotopic
    /// envelopes, intensity, apex, mass error, and distinct-sequence counts. Its `detection_type`
    /// is [`DetectionType::MBR`].
    pub peak: ChromatographicPeak,
    /// Extension-stripped name of the acceptor spectra file this peak belongs to. Stands in for the
    /// C# `SpectraFileInfo` reference (the base peak carries no file info in this port); used by the
    /// P3.2d loop to group acceptor peaks and by the scorer's per-file distributions.
    pub spectra_file_name: String,
    /// The retention time predicted for this peptide in the acceptor run by the RT-alignment spline
    /// (C# `PredictedRetentionTime`, an `init` field).
    pub predicted_retention_time: f64,
    /// Whether this peak's RT was randomized — i.e. it is an MBR **decoy** peak (C# read-only
    /// `RandomRt`).
    pub random_rt: bool,
    /// `predictedRt - apexRetentionTime`, set by the scorer (C# `RtPredictionError`).
    pub rt_prediction_error: f64,
    /// MBR q-value, assigned by the FDR step (C# internal `MbrQValue`; P3.4).
    pub mbr_q_value: f64,
    /// MBR posterior error probability, assigned by the Python ML model (C# `MbrPep`, nullable;
    /// P3.3). `None` until scored.
    pub mbr_pep: Option<f64>,
    /// The combined MBR score in `[0, 100]`; higher is better (C# `MbrScore`). Set by the scorer
    /// as `100·(intensity·rt·ppm·scanCount·isotopic)^0.2`.
    pub mbr_score: f64,
    /// Ppm component score (C# `PpmScore`).
    pub ppm_score: f64,
    /// Intensity component score (C# `IntensityScore`).
    pub intensity_score: f64,
    /// RT component score (C# `RtScore`).
    pub rt_score: f64,
    /// Scan-count component score (C# `ScanCountScore`).
    pub scan_count_score: f64,
    /// Isotopic-distribution component score (C# `IsotopicDistributionScore`).
    pub isotopic_distribution_score: f64,
    /// Charge states matched for this acceptor peak (C# base `ChargeList`, set by the MBR acceptor
    /// search in `FlashLfqEngine.cs:1232`). Empty until the P3.2d search populates it.
    pub charge_list: Vec<i32>,
}

impl MbrChromatographicPeak {
    /// Builds an MBR peak from one identification, mirroring the C# constructor
    /// `MbrChromatographicPeak(id, fileInfo, predictedRetentionTime, randomRt)` (which calls the
    /// base ctor with `detectionType: DetectionType.MBR`).
    ///
    /// `peakfinding_mass` is the identification's `PeakfindingMass` (which `Identification` does not
    /// itself store), forwarded to the base peak for its mass-error term — exactly as the Phase 1
    /// MSMS path passes it (see [`ChromatographicPeak::from_identification`]). `spectra_file_name`
    /// is the acceptor file the peak is searched in.
    pub fn new(
        id: Identification,
        spectra_file_name: String,
        peakfinding_mass: f64,
        predicted_retention_time: f64,
        random_rt: bool,
    ) -> MbrChromatographicPeak {
        let peak = ChromatographicPeak::from_identifications(
            vec![id],
            vec![peakfinding_mass],
            DetectionType::MBR,
        );
        MbrChromatographicPeak {
            peak,
            spectra_file_name,
            predicted_retention_time,
            random_rt,
            rt_prediction_error: 0.0,
            mbr_q_value: 0.0,
            mbr_pep: None,
            mbr_score: 0.0,
            ppm_score: 0.0,
            intensity_score: 0.0,
            rt_score: 0.0,
            scan_count_score: 0.0,
            isotopic_distribution_score: 0.0,
            charge_list: Vec::new(),
        }
    }

    /// Whether this is an MBR decoy peak — i.e. its RT was randomized (C# `RandomRt`; the
    /// doc-comment states a randomized RT "implies that this peak is a decoy peak identified by the
    /// MBR algorithm"). Distinct from [`ChromatographicPeak::decoy_peptide`], which reflects the
    /// identification being a decoy *peptide*; the FDR step (P3.4) keys on the `(decoy_peptide,
    /// random_rt)` pair.
    pub fn is_decoy_peak(&self) -> bool {
        self.random_rt
    }

    /// Feature intensity of the underlying peak (C# `Intensity`), read by the scorer.
    pub fn intensity(&self) -> f64 {
        self.peak.intensity
    }

    /// Apex retention time of the underlying peak (C# `ApexRetentionTime`), read by the scorer.
    pub fn apex_retention_time(&self) -> f64 {
        self.peak.apex_retention_time()
    }

    /// Apex mass error in ppm of the underlying peak (C# `MassError`), read by the scorer.
    pub fn mass_error(&self) -> f64 {
        self.peak.mass_error
    }

    /// Number of isotopic envelopes traced (C# `ScanCount`), read by the scorer.
    pub fn scan_count(&self) -> usize {
        self.peak.scan_count()
    }

    /// Apex isotopic Pearson correlation (C# `IsotopicPearsonCorrelation`), read by the scorer.
    pub fn isotopic_pearson_correlation(&self) -> f64 {
        self.peak.isotopic_pearson_correlation()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a minimal identification for constructor tests.
    fn id(modified: &str, decoy: bool) -> Identification {
        Identification {
            file_name: "acceptor".to_string(),
            base_sequence: "PEPTIDE".to_string(),
            modified_sequence: modified.to_string(),
            monoisotopic_mass: 799.36,
            ms2_retention_time_in_minutes: -1.0,
            precursor_charge_state: 2,
            score: 10.0,
            q_value: 0.001,
            is_decoy: decoy,
        }
    }

    #[test]
    fn new_sets_mbr_detection_type_and_defaults() {
        let peak = MbrChromatographicPeak::new(
            id("PEPTIDE", false),
            "acceptor".to_string(),
            799.36,
            42.5,
            false,
        );

        // Base peak built with the MBR detection type and the single id resolved.
        assert_eq!(peak.peak.detection_type, DetectionType::MBR);
        assert_eq!(peak.peak.identifications.len(), 1);
        assert_eq!(peak.peak.num_identifications_by_full_seq, 1);
        assert_eq!(peak.peak.num_identifications_by_base_seq, 1);

        // MBR fields take their fresh-construction defaults.
        assert_eq!(peak.predicted_retention_time, 42.5);
        assert!(!peak.random_rt);
        assert_eq!(peak.rt_prediction_error, 0.0);
        assert_eq!(peak.mbr_q_value, 0.0);
        assert!(peak.mbr_pep.is_none());
        assert_eq!(peak.mbr_score, 0.0);
        assert_eq!(peak.ppm_score, 0.0);
        assert_eq!(peak.intensity_score, 0.0);
        assert_eq!(peak.rt_score, 0.0);
        assert_eq!(peak.scan_count_score, 0.0);
        assert_eq!(peak.isotopic_distribution_score, 0.0);
        assert!(peak.charge_list.is_empty());
        assert_eq!(peak.spectra_file_name, "acceptor");
    }

    #[test]
    fn random_rt_marks_decoy_peak() {
        let target = MbrChromatographicPeak::new(
            id("PEPTIDE", false),
            "acceptor".to_string(),
            799.36,
            1.0,
            false,
        );
        assert!(!target.is_decoy_peak());

        let decoy = MbrChromatographicPeak::new(
            id("PEPTIDE", false),
            "acceptor".to_string(),
            799.36,
            1.0,
            true,
        );
        // A randomized RT marks an MBR decoy peak, independent of whether the peptide id is a decoy.
        assert!(decoy.is_decoy_peak());
        assert!(!decoy.peak.decoy_peptide());
    }

    #[test]
    fn accessors_reflect_empty_base_peak() {
        // With no isotopic envelopes the base peak has no apex: scan count 0, apex RT/correlation -1.
        let peak = MbrChromatographicPeak::new(
            id("PEPTIDE", false),
            "acceptor".to_string(),
            799.36,
            1.0,
            false,
        );
        assert_eq!(peak.scan_count(), 0);
        assert_eq!(peak.apex_retention_time(), -1.0);
        assert_eq!(peak.isotopic_pearson_correlation(), -1.0);
        assert_eq!(peak.intensity(), 0.0);
        assert!(peak.mass_error().is_nan());
    }
}
