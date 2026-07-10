//! Peak detection type — verbatim port of mzLib `FlashLFQ.DetectionType`.
//!
//! Records how a [`crate::chromatographic_peak::ChromatographicPeak`] was found. Phase 1 only
//! produces `MSMS` (and, after charge aggregation, `MSMSAmbiguousPeakfinding`) and reports
//! `MSMSIdentifiedButNotQuantified` / `NotDetected` at the peptide level; the MBR / IsoTracker
//! variants are carried for fidelity and used in later phases.

/// How a chromatographic peak was detected. Verbatim from `FlashLFQ/DetectionType.cs`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DetectionType {
    /// The peak is detected from an MS2 identification.
    MSMS,
    /// The peak is detected from match-between-runs.
    MBR,
    /// The peak is detected from more than one MS2 identification.
    MSMSAmbiguousPeakfinding,
    /// The peak is detected from an MS2 identification by IsoTrack (priority IsoTrack_MSMS > MSMS).
    IsoTrackMSMS,
    /// The peak is detected from MBR by IsoTrack.
    IsoTrackMBR,
    /// Multiple peptides are mapped to this peak by IsoTracker.
    IsoTrackAmbiguous,
    /// We have an MS2 identification but no peak for quantification.
    MSMSIdentifiedButNotQuantified,
    /// Peptide was not detected by MS2 or MBR (only used in peptide-level quantification).
    NotDetected,
}

impl DetectionType {
    /// The detection type's name as written in FlashLFQ output (the verbatim C# enum member
    /// names, including the underscores on the IsoTrack variants). Used by the Parquet/TSV
    /// results writers (P1.16).
    pub fn as_str(self) -> &'static str {
        match self {
            DetectionType::MSMS => "MSMS",
            DetectionType::MBR => "MBR",
            DetectionType::MSMSAmbiguousPeakfinding => "MSMSAmbiguousPeakfinding",
            DetectionType::IsoTrackMSMS => "IsoTrack_MSMS",
            DetectionType::IsoTrackMBR => "IsoTrack_MBR",
            DetectionType::IsoTrackAmbiguous => "IsoTrack_Ambiguous",
            DetectionType::MSMSIdentifiedButNotQuantified => "MSMSIdentifiedButNotQuantified",
            DetectionType::NotDetected => "NotDetected",
        }
    }
}
