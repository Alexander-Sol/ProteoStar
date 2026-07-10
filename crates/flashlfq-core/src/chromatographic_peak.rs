//! Chromatographic peak + `CutPeak` — partial port of mzLib `FlashLFQ.ChromatographicPeak`
//! (the `IsotopicEnvelope` container and intensity recomputation) and the
//! `FlashLfqEngine.CutPeak` boundary trimming.
//!
//! A [`ChromatographicPeak`] holds the [`IsotopicEnvelope`]s found for one identification across
//! charge states (the P1.13 output). [`ChromatographicPeak::calculate_intensity_for_this_feature`]
//! ports `ChromatographicPeak.CalculateIntensityForThisFeature` (apex = most intense envelope,
//! intensity = sum or apex, mass error, distinct charge-state count) and
//! [`ChromatographicPeak::cut_peak`] ports `FlashLfqEngine.CutPeak` (recursive valley splitting).
//!
//! ## Fidelity notes
//! - `DiscriminationFactorToCutPeak = 0.6` (mzLib `FlashLfqEngine` const).
//! - `CutPeak` assumes the envelopes are already sorted by MS1 scan number (the caller's
//!   responsibility, exactly as in C#).
//! - The apex is `IsotopicEnvelopes.MaxBy(Intensity)` — LINQ `MaxBy` returns the **first**
//!   element attaining the maximum, so [`max_by_intensity`] keeps the first on ties.
//! - C#'s `timePointsForApexZ.IndexOf(apex)`/`IndexOf(valleyEnvelope)` use reference equality;
//!   here the valley's index is the loop index `i` directly (faithful: the valley *is* that
//!   element), and the apex's index is found by value equality (the global-max envelope is unique
//!   by scan index, so no earlier value-equal element can shadow it).
//! - `RemoveAll` filters **all** envelopes (every charge state) by retention time relative to the
//!   valley, even though the valley was chosen from the apex-charge subset — matching C#.
//! - Mass error mirrors `CalculateIntensityForThisFeature`: for each identification peakfinding
//!   mass, `((ToMass(Apex.M, Apex.ChargeState) - pfm) / pfm) * 1e6`, keeping the smallest-magnitude
//!   value (seeded from `NaN`). `ToMass` is the `f32` overload (see [`mz_to_mass_f32`]).

use std::collections::HashSet;

use crate::detection_type::DetectionType;
use crate::isotopic_envelope::{mz_to_mass_f32, IsotopicEnvelope};
use crate::psm_tsv::Identification;

/// Determines splitting behaviour for bimodal peaks. mzLib `FlashLfqEngine`
/// `DiscriminationFactorToCutPeak`.
pub const DISCRIMINATION_FACTOR_TO_CUT_PEAK: f64 = 0.6;

/// A chromatographic peak: the isotopic envelopes found for one identification across charge
/// states, plus the derived feature summary. Partial port of `FlashLFQ.ChromatographicPeak`
/// limited to what the P1.14 tracing path (`CutPeak` + intensity recomputation) needs.
///
/// The full peak (identifications, protein groups, detection type, output formatting) is built up
/// in later tasks; here we keep the envelopes, the precomputed peakfinding masses of the
/// identifications (for the mass-error term), and the fields `CutPeak` reads/writes.
#[derive(Debug, Clone)]
pub struct ChromatographicPeak {
    /// The isotopic envelopes, expected to be sorted by MS1 scan number (the caller's invariant,
    /// matching the C# `CutPeak` precondition).
    pub isotopic_envelopes: Vec<IsotopicEnvelope>,
    /// The peakfinding masses of the identifications mapped to this peak (C#
    /// `Identifications.Select(id => id.PeakfindingMass)`), used only for the mass-error term.
    /// Kept parallel to [`identifications`](Self::identifications); these are the per-id
    /// `PeakfindingMass` values, which `Identification` does not itself store.
    pub identification_peakfinding_masses: Vec<f64>,
    /// The identifications mapped to this peak (C# `Identifications`). Empty for the bare
    /// CutPeak fixtures built via [`ChromatographicPeak::new`].
    pub identifications: Vec<Identification>,
    /// How this peak was detected (C# `DetectionType`). Phase 1 produces `MSMS`, promoted to
    /// `MSMSAmbiguousPeakfinding` when merged across ambiguous identifications.
    pub detection_type: DetectionType,
    /// Summed (or apex) intensity of the feature (C# `Intensity`).
    pub intensity: f64,
    /// The most intense envelope (C# `Apex`); `None` when there are no envelopes.
    pub apex: Option<IsotopicEnvelope>,
    /// Retention time of the split valley, if the peak was cut (C# `SplitRT`).
    pub split_rt: f64,
    /// Smallest-magnitude apex mass error in ppm across the identifications (C# `MassError`);
    /// `NaN` when there are no envelopes.
    pub mass_error: f64,
    /// Number of distinct charge states observed across the envelopes (C#
    /// `NumChargeStatesObserved`).
    pub num_charge_states_observed: usize,
    /// Number of distinct base sequences among the identifications (C#
    /// `NumIdentificationsByBaseSeq`).
    pub num_identifications_by_base_seq: usize,
    /// Number of distinct modified (full) sequences among the identifications (C#
    /// `NumIdentificationsByFullSeq`).
    pub num_identifications_by_full_seq: usize,
}

impl ChromatographicPeak {
    /// Builds an empty peak from peakfinding masses only (no identifications), matching the C#
    /// constructor state (`SplitRT = 0`, `MassError = NaN`, no envelopes). Used by the CutPeak
    /// fixtures, which only need the mass-error term; the identification-bearing constructors
    /// ([`from_identification`](Self::from_identification) /
    /// [`from_identifications`](Self::from_identifications)) are used by the aggregation path.
    pub fn new(identification_peakfinding_masses: Vec<f64>) -> ChromatographicPeak {
        ChromatographicPeak {
            isotopic_envelopes: Vec::new(),
            identification_peakfinding_masses,
            identifications: Vec::new(),
            detection_type: DetectionType::MSMS,
            intensity: 0.0,
            apex: None,
            split_rt: 0.0,
            mass_error: f64::NAN,
            num_charge_states_observed: 0,
            num_identifications_by_base_seq: 0,
            num_identifications_by_full_seq: 0,
        }
    }

    /// Builds a peak for a single identification, mirroring the C#
    /// `ChromatographicPeak(Identification, SpectraFileInfo, DetectionType)` constructor:
    /// `SplitRT = 0`, `MassError = NaN`, no envelopes, then `ResolveIdentifications`.
    ///
    /// `peakfinding_mass` is the identification's `PeakfindingMass`
    /// (`MonoisotopicMass + mostAbundantIsotopeShift`), which `Identification` itself does not
    /// store — it is computed from the theoretical distribution (see
    /// [`crate::theoretical_isotope_distribution::peakfinding_mass`]).
    pub fn from_identification(id: Identification, peakfinding_mass: f64) -> ChromatographicPeak {
        Self::from_identifications(vec![id], vec![peakfinding_mass], DetectionType::MSMS)
    }

    /// Builds a peak for a set of identifications (the C# multi-id / isobaric-ambiguity
    /// constructor), running [`resolve_identifications`](Self::resolve_identifications) so the
    /// distinct-sequence counts are populated from the start.
    pub fn from_identifications(
        identifications: Vec<Identification>,
        identification_peakfinding_masses: Vec<f64>,
        detection_type: DetectionType,
    ) -> ChromatographicPeak {
        let mut peak = ChromatographicPeak {
            isotopic_envelopes: Vec::new(),
            identification_peakfinding_masses,
            identifications,
            detection_type,
            intensity: 0.0,
            apex: None,
            split_rt: 0.0,
            mass_error: f64::NAN,
            num_charge_states_observed: 0,
            num_identifications_by_base_seq: 0,
            num_identifications_by_full_seq: 0,
        };
        peak.resolve_identifications();
        peak
    }

    /// Apex retention time, or `-1.0` when there is no apex (C# `ApexRetentionTime`).
    pub fn apex_retention_time(&self) -> f64 {
        match self.apex {
            Some(apex) => apex.indexed_peak.retention_time as f64,
            None => -1.0,
        }
    }

    /// Number of isotopic envelopes traced for this peak (C# `ScanCount =>
    /// IsotopicEnvelopes.Count`). Consumed by the MBR scorer (P3.2c) as the scan-count feature.
    pub fn scan_count(&self) -> usize {
        self.isotopic_envelopes.len()
    }

    /// Pearson correlation of the apex envelope to the theoretical isotope abundances, or `-1.0`
    /// when there is no apex (C# `IsotopicPearsonCorrelation => Apex?.PearsonCorrelation ?? -1`).
    /// Consumed by the MBR scorer (P3.2c) as the isotopic-distribution feature.
    pub fn isotopic_pearson_correlation(&self) -> f64 {
        match self.apex {
            Some(apex) => apex.pearson_correlation,
            None => -1.0,
        }
    }

    /// Whether the peak's first identification is a decoy (C# `DecoyPeptide`). Returns `false`
    /// when there are no identifications (the bare CutPeak fixtures).
    pub fn decoy_peptide(&self) -> bool {
        self.identifications.first().map(|id| id.is_decoy).unwrap_or(false)
    }

    /// The modified (full) sequence of the first identification, if any (C#
    /// `Identifications.First().ModifiedSequence`).
    pub fn first_modified_sequence(&self) -> Option<&str> {
        self.identifications
            .first()
            .map(|id| id.modified_sequence.as_str())
    }

    /// Recomputes `NumIdentificationsByBaseSeq` / `NumIdentificationsByFullSeq` as the number of
    /// distinct base / modified sequences among the identifications. Port of
    /// `ChromatographicPeak.ResolveIdentifications`.
    pub fn resolve_identifications(&mut self) {
        let base: HashSet<&str> = self
            .identifications
            .iter()
            .map(|id| id.base_sequence.as_str())
            .collect();
        let full: HashSet<&str> = self
            .identifications
            .iter()
            .map(|id| id.modified_sequence.as_str())
            .collect();
        self.num_identifications_by_base_seq = base.len();
        self.num_identifications_by_full_seq = full.len();
    }

    /// Merges `other` into this peak by unioning identifications and isotopic envelopes, then
    /// recomputing the feature intensity. Port of `ChromatographicPeak.MergeFeatureWith`.
    ///
    /// `other` is not modified. Identifications already present (by value equality, mirroring the
    /// C# reference-equality `Union().Distinct()`) are not duplicated, and envelopes whose indexed
    /// peak is already present are skipped. The peakfinding-mass list is unioned in lockstep with
    /// the identifications so the mass-error term still sees every id.
    pub fn merge_feature_with(&mut self, other: &ChromatographicPeak, integrate: bool) {
        // C# guards `otherFeature != this`; an explicit self-merge is a no-op here too.
        if std::ptr::eq(self, other) {
            return;
        }

        let existing_peaks: HashSet<EnvelopePeakKey> = self
            .isotopic_envelopes
            .iter()
            .map(|e| EnvelopePeakKey::from(e))
            .collect();

        for (id, &pfm) in other
            .identifications
            .iter()
            .zip(other.identification_peakfinding_masses.iter())
        {
            if !self.identifications.iter().any(|existing| existing == id) {
                self.identifications.push(id.clone());
                self.identification_peakfinding_masses.push(pfm);
            }
        }
        self.resolve_identifications();

        for &env in &other.isotopic_envelopes {
            if !existing_peaks.contains(&EnvelopePeakKey::from(&env)) {
                self.isotopic_envelopes.push(env);
            }
        }

        self.calculate_intensity_for_this_feature(integrate);

        // A merge onto an IsoTracker peak promotes its detection type (C# `MergeFeatureWith`).
        if matches!(
            self.detection_type,
            DetectionType::IsoTrackMBR | DetectionType::IsoTrackMSMS
        ) {
            self.detection_type = DetectionType::MSMSAmbiguousPeakfinding;
        }
    }

    /// Recomputes apex, intensity, mass error, and the distinct-charge-state count from the current
    /// envelopes. Port of `ChromatographicPeak.CalculateIntensityForThisFeature`.
    ///
    /// When `integrate` is true the intensity is the sum of all envelope intensities; otherwise it
    /// is the apex envelope's intensity.
    pub fn calculate_intensity_for_this_feature(&mut self, integrate: bool) {
        if let Some(apex) = max_by_intensity(&self.isotopic_envelopes) {
            self.apex = Some(apex);

            self.intensity = if integrate {
                self.isotopic_envelopes.iter().map(|p| p.intensity).sum()
            } else {
                apex.intensity
            };

            self.mass_error = f64::NAN;
            let apex_mass = mz_to_mass_f32(apex.indexed_peak.m(), apex.charge_state);
            for &pfm in &self.identification_peakfinding_masses {
                let mass_error_for_id = ((apex_mass - pfm) / pfm) * 1e6;
                if self.mass_error.is_nan() || mass_error_for_id.abs() < self.mass_error.abs() {
                    self.mass_error = mass_error_for_id;
                }
            }

            let charge_states: HashSet<i32> = self
                .isotopic_envelopes
                .iter()
                .map(|p| p.charge_state)
                .collect();
            self.num_charge_states_observed = charge_states.len();
        } else {
            self.intensity = 0.0;
            self.mass_error = f64::NAN;
            self.num_charge_states_observed = 0;
            self.apex = None;
        }
    }

    /// Recursively trims envelopes occurring before/after "valleys" surrounding the
    /// identification's MS2 retention time, recomputing intensity after each cut. Port of
    /// `FlashLfqEngine.CutPeak`.
    ///
    /// Assumes the envelopes are already sorted by MS1 scan number. `identification_time` is the
    /// MS2 scan retention time; `integrate` is forwarded to the intensity recomputation.
    pub fn cut_peak(&mut self, identification_time: f64, integrate: bool) {
        // this method assumes the isotope envelopes are already sorted by MS1 scan number
        let mut cut_this_peak = false;

        if self.isotopic_envelopes.len() < 5 {
            return;
        }

        let apex = match self.apex {
            Some(a) => a,
            None => return,
        };

        // Ordered list of all time points where the apex charge state had a valid envelope.
        let time_points: Vec<IsotopicEnvelope> = self
            .isotopic_envelopes
            .iter()
            .filter(|p| p.charge_state == apex.charge_state)
            .copied()
            .collect();
        let scan_numbers: HashSet<i32> = time_points
            .iter()
            .map(|p| p.indexed_peak.zero_based_scan_index)
            .collect();
        let apex_index = time_points
            .iter()
            .position(|p| *p == apex)
            .expect("apex is one of the apex-charge envelopes") as i32;
        let mut valley_envelope: Option<IsotopicEnvelope> = None;

        // +1 checks the right side, -1 checks the left side.
        for &direction in &[1i32, -1i32] {
            valley_envelope = None;
            let mut index_of_valley: i32 = 0;

            let mut i = apex_index + direction;
            while i < time_points.len() as i32 && i >= 0 {
                let timepoint = time_points[i as usize];

                // Valley envelope is the lowest-intensity point encountered thus far.
                if valley_envelope.is_none()
                    || timepoint.intensity < valley_envelope.unwrap().intensity
                {
                    valley_envelope = Some(timepoint);
                    index_of_valley = i;
                }

                let valley = valley_envelope.unwrap();
                let mut discrimination_factor =
                    (timepoint.intensity - valley.intensity) / timepoint.intensity;

                // If the time point is at least discriminationFactor times more intense than the
                // valley, check whether it is also more intense than the point next to the valley.
                if discrimination_factor > DISCRIMINATION_FACTOR_TO_CUT_PEAK
                    && index_of_valley + direction < time_points.len() as i32
                    && index_of_valley + direction >= 0
                {
                    let second_valley_timepoint =
                        time_points[(index_of_valley + direction) as usize];

                    discrimination_factor = (timepoint.intensity - second_valley_timepoint.intensity)
                        / timepoint.intensity;

                    // Cut if the timepoint is more intense than the second valley, or if no
                    // envelope is observed in the scan immediately after the valley.
                    if discrimination_factor > DISCRIMINATION_FACTOR_TO_CUT_PEAK
                        || !scan_numbers
                            .contains(&(valley.indexed_peak.zero_based_scan_index + direction))
                    {
                        cut_this_peak = true;
                        break;
                    }
                }

                i += direction;
            }

            if cut_this_peak {
                break;
            }
        }

        // cut
        if cut_this_peak {
            let valley = valley_envelope.expect("a valley was found when cutting");
            let valley_rt = valley.indexed_peak.retention_time as f64;

            if identification_time > valley_rt {
                // MS2 id is right of the valley; remove all envelopes left of (and at) the valley.
                self.isotopic_envelopes
                    .retain(|p| !(p.indexed_peak.retention_time as f64 <= valley_rt));
            } else {
                // MS2 id is left of the valley; remove all envelopes right of (and at) the valley.
                self.isotopic_envelopes
                    .retain(|p| !(p.indexed_peak.retention_time as f64 >= valley_rt));
            }

            // recalculate intensity for the peak
            self.calculate_intensity_for_this_feature(integrate);
            self.split_rt = valley_rt;

            // recursively cut
            self.cut_peak(identification_time, integrate);
        }
    }
}

/// A structural identity key for an envelope's indexed peak, standing in for the C# reference
/// equality of `IIndexedPeak` (C# threads the same peak object; our peaks are `Copy` values, so
/// two envelopes referring to the same physical peak share scan index + m/z + intensity bits).
/// Used to dedup envelopes on merge and to group peaks by apex in `RunErrorChecking`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EnvelopePeakKey {
    /// Zero-based MS1 scan index of the indexed peak.
    pub zero_based_scan_index: i32,
    /// Bit pattern of the peak m/z (`f32`), for exact structural equality.
    pub mz_bits: u32,
    /// Bit pattern of the peak intensity (`f32`), for exact structural equality.
    pub intensity_bits: u32,
}

impl EnvelopePeakKey {
    /// Builds the key from an envelope's indexed peak.
    pub fn from(envelope: &IsotopicEnvelope) -> EnvelopePeakKey {
        let peak = envelope.indexed_peak;
        EnvelopePeakKey {
            zero_based_scan_index: peak.zero_based_scan_index,
            mz_bits: peak.mz.to_bits(),
            intensity_bits: peak.intensity.to_bits(),
        }
    }
}

/// Returns the first envelope attaining the maximum intensity (LINQ `MaxBy` semantics), or `None`
/// for an empty list.
fn max_by_intensity(envelopes: &[IsotopicEnvelope]) -> Option<IsotopicEnvelope> {
    let mut best: Option<IsotopicEnvelope> = None;
    for &env in envelopes {
        match best {
            None => best = Some(env),
            Some(b) if env.intensity > b.intensity => best = Some(env),
            _ => {}
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::peak_indexing::IndexedMassSpectralPeak;

    /// Builds an envelope at the given scan/RT with the given intensity (charge 1, the
    /// `IsotopicEnvelope::new` divisor of 1 leaves `intensity` unchanged).
    fn env(scan: i32, rt: f64, intensity: f64) -> IsotopicEnvelope {
        env_z(scan, rt, intensity, 1)
    }

    fn env_z(scan: i32, rt: f64, intensity: f64, charge: i32) -> IsotopicEnvelope {
        let peak = IndexedMassSpectralPeak::new(1000.0, intensity, scan, rt);
        // pass intensity * charge so that after the /charge in `new` the stored value is `intensity`
        IsotopicEnvelope::new(peak, charge, intensity * charge as f64, 0.95)
    }

    fn peak_from(envelopes: Vec<IsotopicEnvelope>, integrate: bool) -> ChromatographicPeak {
        let mut peak = ChromatographicPeak::new(vec![1000.0]);
        peak.isotopic_envelopes = envelopes;
        peak.calculate_intensity_for_this_feature(integrate);
        peak
    }

    #[test]
    fn discrimination_factor_constant() {
        assert_eq!(DISCRIMINATION_FACTOR_TO_CUT_PEAK, 0.6);
    }

    #[test]
    fn calculate_intensity_apex_and_sum() {
        let envs = vec![env(0, 1.0, 10.0), env(1, 1.1, 50.0), env(2, 1.2, 20.0)];
        let mut peak = peak_from(envs.clone(), false);
        assert_eq!(peak.apex.unwrap().indexed_peak.zero_based_scan_index, 1);
        assert!((peak.intensity - 50.0).abs() < 1e-9);
        peak.calculate_intensity_for_this_feature(true);
        assert!((peak.intensity - 80.0).abs() < 1e-9);
        assert_eq!(peak.num_charge_states_observed, 1);
    }

    #[test]
    fn calculate_intensity_empty_resets() {
        let mut peak = ChromatographicPeak::new(vec![1000.0]);
        peak.calculate_intensity_for_this_feature(true);
        assert!(peak.apex.is_none());
        assert_eq!(peak.intensity, 0.0);
        assert!(peak.mass_error.is_nan());
        assert_eq!(peak.num_charge_states_observed, 0);
    }

    #[test]
    fn cut_peak_no_cut_when_fewer_than_five_envelopes() {
        let envs = vec![env(0, 1.0, 1.0), env(1, 1.1, 5.0), env(2, 1.2, 1.0)];
        let mut peak = peak_from(envs, true);
        let before = peak.isotopic_envelopes.len();
        peak.cut_peak(1.1, true);
        assert_eq!(peak.isotopic_envelopes.len(), before);
        assert_eq!(peak.split_rt, 0.0);
    }

    #[test]
    fn cut_peak_no_cut_on_clean_unimodal_peak() {
        // A single clean elution profile (no deep valley) must not be split.
        let envs = vec![
            env(0, 1.0, 1.0),
            env(1, 1.1, 3.0),
            env(2, 1.2, 6.0),
            env(3, 1.3, 3.0),
            env(4, 1.4, 1.0),
        ];
        let mut peak = peak_from(envs, true);
        peak.cut_peak(1.2, true);
        assert_eq!(peak.isotopic_envelopes.len(), 5);
    }

    #[test]
    fn cut_peak_splits_bimodal_peak_keeping_id_side() {
        // Two elution humps separated by a deep valley at scan 4 (RT 1.4). The apex is the
        // left hump (intensity 100 at scan 2). The MS2 id RT (1.2) is left of the valley, so the
        // right hump (scans 5..8) must be removed.
        let envs = vec![
            env(0, 1.0, 10.0),
            env(1, 1.1, 60.0),
            env(2, 1.2, 100.0), // apex (left hump)
            env(3, 1.3, 60.0),
            env(4, 1.4, 5.0),  // valley
            env(5, 1.5, 10.0), // shallow rising edge past the valley
            env(6, 1.6, 40.0),
            env(7, 1.7, 90.0), // right hump
            env(8, 1.8, 8.0),
        ];
        let mut peak = peak_from(envs, true);
        peak.cut_peak(1.2, true);

        // All surviving envelopes are at or left of the valley RT.
        assert!(peak
            .isotopic_envelopes
            .iter()
            .all(|p| (p.indexed_peak.retention_time as f64) < 1.4));
        assert_eq!(peak.split_rt, 1.4_f32 as f64);
        // Apex preserved on the kept side.
        assert_eq!(peak.apex.unwrap().indexed_peak.zero_based_scan_index, 2);
    }

    #[test]
    fn cut_peak_splits_keeping_right_side_when_id_is_right_of_valley() {
        // Mirror image: apex on the right hump, id RT (1.6) right of the valley -> left hump cut.
        let envs = vec![
            env(0, 1.0, 8.0),
            env(1, 1.1, 40.0),
            env(2, 1.2, 90.0), // left hump
            env(3, 1.3, 10.0), // shallow rising edge toward the valley
            env(4, 1.4, 5.0),  // valley
            env(5, 1.5, 60.0),
            env(6, 1.6, 100.0), // apex (right hump)
            env(7, 1.7, 60.0),
            env(8, 1.8, 10.0),
        ];
        let mut peak = peak_from(envs, true);
        peak.cut_peak(1.6, true);

        assert!(peak
            .isotopic_envelopes
            .iter()
            .all(|p| (p.indexed_peak.retention_time as f64) > 1.4));
        assert_eq!(peak.apex.unwrap().indexed_peak.zero_based_scan_index, 6);
    }
}
