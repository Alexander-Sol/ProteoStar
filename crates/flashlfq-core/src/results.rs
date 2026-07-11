//! Charge aggregation + peptide-level results (PLAN.md P1.15).
//!
//! Ports the post-peakfinding aggregation FlashLFQ performs once every MS2 identification has a
//! [`ChromatographicPeak`] (built in P1.13/P1.14):
//!
//! 1. [`run_error_checking`] — port of the MSMS path of `FlashLfqEngine.RunErrorChecking`. After
//!    peak finding, several identifications can resolve to peaks sharing the same apex MS1 peak
//!    (typically the same peptide found at multiple charge states, or near-isobaric peptides).
//!    Those duplicates are merged so one apex yields one feature; conflicting peaks are kept or
//!    overwritten according to whether their modified sequence is in the quantify set.
//! 2. [`calculate_peptide_results`] — port of the MSMS path of
//!    `FlashLfqResults.CalculatePeptideResults`. For each file, peaks are grouped by modified
//!    sequence and the most intense unambiguous peak fixes the peptide's intensity / apex RT /
//!    detection type. Ambiguous peaks (>1 full sequence) zero out a sequence whose intensity is
//!    mostly explained by the ambiguous feature.
//!
//! ## Phase-1 scope (no MBR, no IsoTracker, no fractions)
//!
//! The full C# routines also resolve MBR-vs-MSMS conflicts, IsoTracker ambiguity, and
//! fraction-level ambiguity. Phase 1 quantifies mzML with MS2 ids only, so every peak is
//! `DetectionType::MSMS`; the MBR / IsoTracker branches are unreachable and are intentionally not
//! ported here (they land with Phase 3). `HandleAmbiguityInFractions` is likewise deferred — it
//! only changes results for fractionated samples, which Phase 1 does not target. These omissions
//! are noted so a later phase knows exactly what is missing.
//!
//! ## FDR (P2.2 — decoy / quantify-set hardening)
//!
//! The MS2 path computes no per-peak q-value: a peak's q-value is the identification's search
//! q-value, and FDR control is exercised entirely through **two** modified-sequence sets that the
//! C# engine path threads separately. Phase 1 conflated them (legal only because the corpus is all
//! targets); P2.2 splits them faithfully so the decoy path is correct:
//!
//! - The **quantify set** ([`engine_quantify_set`]) — the membership gate used by
//!   [`run_error_checking`] (merge / overwrite decisions) and by the unambiguous / ambiguous
//!   filters in [`calculate_peptide_results`]. The C# `FlashLfqEngine` constructor defaults this to
//!   `allIdentifications.Select(id => id.ModifiedSequence).ToHashSet()` — **including decoy modified
//!   sequences** (`FlashLfqEngine.cs:94`). When MetaMorpheus supplies its own peptide-level
//!   q<0.01 set, that supplied set replaces the default; that injection is the MS2 FDR knob.
//! - The **output peptide set** ([`output_peptide_sequences`]) — the modified sequences that
//!   actually get a row in the peptide×file table. C# builds it from
//!   `identifications.Where(id => !id.IsDecoy && quantifySet.Contains(id.ModifiedSequence))`
//!   (`FlashLfqResults.cs:42`): **non-decoy** ids whose sequence is in the quantify set. Decoy-only
//!   sequences therefore never appear in the results even though they sit in the quantify set and
//!   can still influence error-checking merges. Decoy-bearing peaks are additionally filtered out
//!   of [`calculate_peptide_results`] (`!IsDecoy`), matching C#.
//!
//! [`default_quantify_set`] is retained: it mirrors the *`FlashLfqResults` constructor fallback*
//! (`?? identifications.Where(!IsDecoy)...`, `FlashLfqResults.cs:33`), reached only when a caller
//! constructs results with a null set — never in the default engine path, which always passes the
//! engine set. The MBR q-value FDR (`CalculateFdrForMbrPeaks` / `CorrectQValues`) is Phase 3.

use std::collections::{HashMap, HashSet};

use crate::chromatographic_peak::{ChromatographicPeak, EnvelopePeakKey};
use crate::detection_type::DetectionType;
use crate::psm_tsv::Identification;

/// The engine-default **quantify set**: every distinct modified sequence among the
/// identifications, **including decoys**. Mirrors the C# `FlashLfqEngine` constructor default
/// `allIdentifications.Select(id => id.ModifiedSequence).ToHashSet()` (`FlashLfqEngine.cs:94`),
/// which is reached whenever no external `peptideSequencesToQuantify` is supplied.
///
/// This is the membership gate used by [`run_error_checking`] and [`calculate_peptide_results`].
/// Decoy sequences belong here (so they participate in error-checking merges exactly as C# does);
/// they are excluded from the *output* by [`output_peptide_sequences`], not by this set.
pub fn engine_quantify_set(identifications: &[Identification]) -> HashSet<String> {
    identifications
        .iter()
        .map(|id| id.modified_sequence.clone())
        .collect()
}

/// The set of modified sequences that get a row in the peptide×file results: non-decoy
/// identifications whose modified sequence is in `quantify_set`. Mirrors the C#
/// `identifications.Where(id => !id.IsDecoy && _peptideModifiedSequencesToQuantify.Contains(
/// id.ModifiedSequence))` loop that seeds `FlashLfqResults.PeptideModifiedSequences`
/// (`FlashLfqResults.cs:42`). Decoy-only sequences are absent even when they sit in `quantify_set`.
pub fn output_peptide_sequences(
    identifications: &[Identification],
    quantify_set: &HashSet<String>,
) -> HashSet<String> {
    identifications
        .iter()
        .filter(|id| !id.is_decoy && quantify_set.contains(&id.modified_sequence))
        .map(|id| id.modified_sequence.clone())
        .collect()
}

/// The `FlashLfqResults` constructor *fallback* quantify set: every non-decoy modified sequence.
/// Mirrors `peptideModifiedSequencesToQuantify ?? identifications.Where(id => !id.IsDecoy)
/// .Select(id => id.ModifiedSequence).ToHashSet()` (`FlashLfqResults.cs:33`). This is reached only
/// when results are constructed with a null set; the default engine path passes
/// [`engine_quantify_set`] instead, so this is kept for fidelity to the `FlashLfqResults` API
/// rather than used by [`crate::engine::run_msms`].
pub fn default_quantify_set(identifications: &[Identification]) -> HashSet<String> {
    identifications
        .iter()
        .filter(|id| !id.is_decoy)
        .map(|id| id.modified_sequence.clone())
        .collect()
}

/// Merges duplicate peaks that share an apex MS1 peak and recomputes feature intensities. Port of
/// the MSMS path of `FlashLfqEngine.RunErrorChecking` for one spectra file.
///
/// Steps, faithful to C#:
/// - drop `null`/empty MBR peaks (Phase 1: just keep every peak, all non-null);
/// - for each peak (stable order — C# orders by `DetectionType == MBR`, all `false` here),
///   recompute its intensity and resolve its identification counts;
/// - a peak with no apex (no envelopes) survives directly (the MBR-skip branch is unreachable);
/// - peaks are grouped by apex peak (structural [`EnvelopePeakKey`], standing in for C#
///   reference equality). When two MSMS peaks land on the same apex:
///   - if both modified sequences are in `quantify_set`, the stored peak absorbs the new one via
///     [`ChromatographicPeak::merge_feature_with`];
///   - if only the new peak's sequence is in the set, it overwrites the stored peak;
///   - otherwise the new peak is dropped (its sequence isn't quantified, the stored one is).
///
/// Returns the error-checked peaks in first-seen apex order (matching C# `Dictionary` iteration
/// order, on which mzLib's output ordering relies), with apex-less peaks appended in input order
/// ahead of the grouped peaks — exactly as C# builds `errorCheckedPeaks` then appends the grouped
/// values.
pub fn run_error_checking(
    peaks: Vec<ChromatographicPeak>,
    quantify_set: &HashSet<String>,
    integrate: bool,
) -> Vec<ChromatographicPeak> {
    // C# first adds apex-less (non-MBR) peaks straight to errorCheckedPeaks, then appends the
    // grouped values. We accumulate the two lists separately and concatenate at the end.
    let mut apexless: Vec<ChromatographicPeak> = Vec::new();
    // Insertion-ordered map: parallel key list + peak store, mirroring C# Dictionary order.
    let mut grouped_keys: Vec<EnvelopePeakKey> = Vec::new();
    let mut grouped: HashMap<EnvelopePeakKey, ChromatographicPeak> = HashMap::new();

    for mut try_peak in peaks {
        try_peak.calculate_intensity_for_this_feature(integrate);
        try_peak.resolve_identifications();

        let apex = match try_peak.apex {
            Some(apex) => apex,
            None => {
                // DetectionType::MBR would `continue` here; Phase 1 peaks are MSMS, so keep it.
                apexless.push(try_peak);
                continue;
            }
        };

        let apex_key = EnvelopePeakKey::from(&apex);

        match grouped.get_mut(&apex_key) {
            Some(stored) => {
                // Both peaks are MSMS in Phase 1 (the MBR branches are unreachable).
                let try_in_set = try_peak
                    .first_modified_sequence()
                    .is_some_and(|s| quantify_set.contains(s));
                if try_in_set {
                    let stored_in_set = stored
                        .first_modified_sequence()
                        .is_some_and(|s| quantify_set.contains(s));
                    if stored_in_set {
                        stored.merge_feature_with(&try_peak, integrate);
                    } else {
                        // Stored peak's sequence isn't quantified; overwrite it in place,
                        // preserving its position in the iteration order (C# `dict[key] = peak`).
                        *stored = try_peak;
                    }
                }
                // else: neither replaces — the new peak is dropped (its sequence isn't quantified).
            }
            None => {
                grouped_keys.push(apex_key);
                grouped.insert(apex_key, try_peak);
            }
        }
    }

    let mut error_checked = apexless;
    for key in grouped_keys {
        if let Some(peak) = grouped.remove(&key) {
            error_checked.push(peak);
        }
    }
    error_checked
}

/// A per-(file, peptide) quantification cell: intensity, apex retention time, and detection type.
/// The peptide-level analogue of one cell of the C# `Peptide` intensity/RT/detection tables.
#[derive(Debug, Clone, PartialEq)]
pub struct PeptideQuant {
    /// Reported peptide intensity for this file (the most intense unambiguous peak's intensity,
    /// or `0.0` when only ambiguous / not-quantified).
    pub intensity: f64,
    /// Apex retention time of the chosen peak, or `0.0` / the ambiguous peak's apex RT.
    pub retention_time: f64,
    /// How the peptide was quantified in this file.
    pub detection_type: DetectionType,
}

impl PeptideQuant {
    /// The C# initial state set for every (peptide, file) before peaks are considered:
    /// `NotDetected`, intensity 0, RT 0.
    fn not_detected() -> PeptideQuant {
        PeptideQuant {
            intensity: 0.0,
            retention_time: 0.0,
            detection_type: DetectionType::NotDetected,
        }
    }
}

/// Peptide-level results: for each modified sequence, a per-file [`PeptideQuant`]. The Phase-1
/// (peptide × file → intensity) table that the L6 parity gate will check.
#[derive(Debug, Clone, Default)]
pub struct PeptideResults {
    /// `modified_sequence -> (file_name -> quant cell)`.
    pub quant: HashMap<String, HashMap<String, PeptideQuant>>,
}

impl PeptideResults {
    /// The quant cell for a (modified sequence, file), or `None` if the sequence isn't quantified.
    pub fn get(&self, modified_sequence: &str, file_name: &str) -> Option<&PeptideQuant> {
        self.quant.get(modified_sequence)?.get(file_name)
    }

    /// The reported intensity for a (modified sequence, file), `0.0` if absent.
    pub fn intensity(&self, modified_sequence: &str, file_name: &str) -> f64 {
        self.get(modified_sequence, file_name)
            .map(|q| q.intensity)
            .unwrap_or(0.0)
    }
}

/// Aggregates error-checked peaks into peptide-level results per file. Port of the MSMS path of
/// `FlashLfqResults.CalculatePeptideResults` (without fraction-ambiguity handling).
///
/// `peaks_by_file` is `file_name -> error-checked peaks` (the output of [`run_error_checking`] per
/// file). `quantify_set` is the membership gate (the C# `_peptideModifiedSequencesToQuantify`,
/// decoys included — see [`engine_quantify_set`]); `peptide_sequences` is the set of sequences that
/// get an output row (non-decoy, see [`output_peptide_sequences`], = the C#
/// `PeptideModifiedSequences` keys). Every sequence in `peptide_sequences` gets a `NotDetected`
/// cell for every file first, then:
/// - **unambiguous** peaks (`NumIdentificationsByFullSeq == 1`, non-decoy) are grouped by modified
///   sequence; the most intense one fixes the cell. `MSMS` if intensity > 0, else
///   `MSMSIdentifiedButNotQuantified`. (C# `intensity = group.Max(Intensity)`, best peak = first
///   with that intensity.)
/// - **ambiguous** peaks (`NumIdentificationsByFullSeq > 1`, non-decoy): for each of their
///   identifications whose sequence is in the quantify set, if the ambiguous feature explains more
///   than 30% of the already-recorded intensity, the cell is zeroed and marked
///   `MSMSAmbiguousPeakfinding` (the `quantifyAmbiguousPeptides = false` branch).
pub fn calculate_peptide_results(
    peaks_by_file: &HashMap<String, Vec<ChromatographicPeak>>,
    quantify_set: &HashSet<String>,
    peptide_sequences: &HashSet<String>,
) -> PeptideResults {
    let mut results = PeptideResults::default();

    // Initialize every output sequence (non-decoy, in the quantify set) to NotDetected for every
    // file. Decoy-only sequences in `quantify_set` get no row, matching C# `PeptideModifiedSequences`.
    let files: Vec<&String> = peaks_by_file.keys().collect();
    for sequence in peptide_sequences {
        let per_file = results.quant.entry(sequence.clone()).or_default();
        for file in &files {
            per_file.insert((*file).clone(), PeptideQuant::not_detected());
        }
    }

    for (file_name, file_peaks) in peaks_by_file {
        // --- unambiguous peaks ---
        // group non-decoy single-full-seq peaks by modified sequence, keeping only quantified ones
        let mut grouped: HashMap<&str, Vec<&ChromatographicPeak>> = HashMap::new();
        for peak in file_peaks {
            if peak.num_identifications_by_full_seq != 1 || peak.decoy_peptide() {
                continue;
            }
            let Some(seq) = peak.first_modified_sequence() else {
                continue;
            };
            if !quantify_set.contains(seq) {
                continue;
            }
            grouped.entry(seq).or_default().push(peak);
        }

        for (seq, group) in grouped {
            // intensity = max; best peak = first peak attaining it (C# `First(p => p.Intensity == intensity)`)
            let intensity = group
                .iter()
                .map(|p| p.intensity)
                .fold(f64::NEG_INFINITY, f64::max);
            let best_peak = group
                .iter()
                .find(|p| p.intensity == intensity)
                .expect("max came from the group");

            let detection_type = if intensity > 0.0 {
                DetectionType::MSMS
            } else {
                DetectionType::MSMSIdentifiedButNotQuantified
            };

            if let Some(cell) = results
                .quant
                .get_mut(seq)
                .and_then(|m| m.get_mut(file_name))
            {
                cell.intensity = intensity;
                cell.retention_time = best_peak.apex_retention_time();
                cell.detection_type = detection_type;
            }
        }

        // --- ambiguous peaks (>1 full sequence) ---
        for peak in file_peaks {
            if peak.num_identifications_by_full_seq <= 1 || peak.decoy_peptide() {
                continue;
            }
            for id in peak.identifications.iter().filter(|id| !id.is_decoy) {
                let seq = id.modified_sequence.as_str();
                if !quantify_set.contains(seq) {
                    continue;
                }
                let already = results.intensity(seq, file_name);
                let fraction_ambiguous = peak.intensity / (already + peak.intensity);
                if fraction_ambiguous > 0.3 {
                    if let Some(cell) =
                        results.quant.get_mut(seq).and_then(|m| m.get_mut(file_name))
                    {
                        cell.detection_type = DetectionType::MSMSAmbiguousPeakfinding;
                        cell.intensity = 0.0;
                        cell.retention_time = peak.apex_retention_time();
                    }
                }
            }
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isotopic_envelope::IsotopicEnvelope;
    use crate::peak_indexing::IndexedMassSpectralPeak;

    fn id(modseq: &str, base: &str, charge: i32, decoy: bool) -> Identification {
        Identification {
            file_name: "f1".to_string(),
            base_sequence: base.to_string(),
            modified_sequence: modseq.to_string(),
            monoisotopic_mass: 1000.0,
            ms2_retention_time_in_minutes: 1.2,
            precursor_charge_state: charge,
            score: 10.0,
            q_value: 0.001,
            is_decoy: decoy,
        }
    }

    /// Envelope at the given scan/RT/intensity/charge (intensity*charge so the /charge in
    /// `IsotopicEnvelope::new` leaves the stored value at `intensity`).
    fn env(scan: i32, rt: f64, intensity: f64, charge: i32, mz: f32) -> IsotopicEnvelope {
        let peak = IndexedMassSpectralPeak::new(mz as f64, intensity, scan, rt);
        IsotopicEnvelope::new(peak, charge, intensity * charge as f64, 0.95)
    }

    fn peak_with(id: Identification, envelopes: Vec<IsotopicEnvelope>) -> ChromatographicPeak {
        let mut peak = ChromatographicPeak::from_identification(id, 1000.0);
        peak.isotopic_envelopes = envelopes;
        peak.calculate_intensity_for_this_feature(false);
        peak
    }

    #[test]
    fn resolve_identifications_counts_distinct_sequences() {
        let mut peak = ChromatographicPeak::from_identifications(
            vec![
                id("PEP[ox]TIDE", "PEPTIDE", 2, false),
                id("PEPTIDE", "PEPTIDE", 3, false), // same base, different full seq
                id("PEPTIDE", "PEPTIDE", 2, false), // duplicate full seq
            ],
            vec![1000.0, 1000.0, 1000.0],
            DetectionType::MSMS,
        );
        peak.resolve_identifications();
        assert_eq!(peak.num_identifications_by_base_seq, 1);
        assert_eq!(peak.num_identifications_by_full_seq, 2);
    }

    #[test]
    fn decoy_peptide_reads_first_identification() {
        let target = ChromatographicPeak::from_identification(id("PEPTIDE", "PEPTIDE", 2, false), 1000.0);
        let decoy = ChromatographicPeak::from_identification(id("PEPTIDE", "PEPTIDE", 2, true), 1000.0);
        assert!(!target.decoy_peptide());
        assert!(decoy.decoy_peptide());
    }

    #[test]
    fn merge_feature_unions_ids_and_envelopes() {
        let mut a = peak_with(
            id("PEPTIDE", "PEPTIDE", 2, false),
            vec![env(0, 1.0, 10.0, 2, 500.0), env(1, 1.1, 30.0, 2, 500.0)],
        );
        // b: another charge state of the same peptide, overlapping scan 1, new scan 2.
        let b = peak_with(
            id("PEPTIDE", "PEPTIDE", 3, false),
            vec![env(1, 1.1, 30.0, 2, 500.0), env(2, 1.2, 5.0, 3, 333.0)],
        );
        a.merge_feature_with(&b, true);

        // ids unioned (two distinct identifications, same full seq -> full-seq count 1)
        assert_eq!(a.identifications.len(), 2);
        assert_eq!(a.num_identifications_by_full_seq, 1);
        // envelopes: 2 original + 1 new (the scan-1 duplicate is skipped by structural key)
        assert_eq!(a.isotopic_envelopes.len(), 3);
        // integrated intensity = 10 + 30 + 5 = 45
        assert!((a.intensity - 45.0).abs() < 1e-9);
    }

    #[test]
    fn run_error_checking_merges_peaks_sharing_an_apex() {
        // Two peaks (different charge states of one peptide) whose apex MS1 peak is identical.
        let shared_apex = env(5, 1.5, 100.0, 2, 500.0);
        let mut p1 = ChromatographicPeak::from_identification(id("SEQ", "SEQ", 2, false), 1000.0);
        p1.isotopic_envelopes = vec![env(4, 1.4, 20.0, 2, 500.0), shared_apex];
        let mut p2 = ChromatographicPeak::from_identification(id("SEQ", "SEQ", 3, false), 1000.0);
        p2.isotopic_envelopes = vec![shared_apex, env(6, 1.6, 10.0, 3, 333.0)];

        let quantify: HashSet<String> = ["SEQ".to_string()].into_iter().collect();
        let out = run_error_checking(vec![p1, p2], &quantify, true);

        assert_eq!(out.len(), 1, "peaks sharing an apex should merge to one");
        // merged envelopes: scan 4, 5 (shared, deduped), 6 = 3 envelopes
        assert_eq!(out[0].isotopic_envelopes.len(), 3);
        assert_eq!(out[0].identifications.len(), 2);
    }

    #[test]
    fn run_error_checking_keeps_distinct_apex_peaks() {
        let mut p1 = ChromatographicPeak::from_identification(id("A", "A", 2, false), 1000.0);
        p1.isotopic_envelopes = vec![env(5, 1.5, 100.0, 2, 500.0)];
        let mut p2 = ChromatographicPeak::from_identification(id("B", "B", 2, false), 1000.0);
        p2.isotopic_envelopes = vec![env(9, 2.5, 80.0, 2, 600.0)];

        let quantify: HashSet<String> = ["A".to_string(), "B".to_string()].into_iter().collect();
        let out = run_error_checking(vec![p1, p2], &quantify, true);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn run_error_checking_overwrites_unquantified_stored_peak() {
        // Stored peak's sequence is NOT in the quantify set; the conflicting peak's sequence IS,
        // so it overwrites the stored peak (C# `errorCheckedPeaksGroupedByApex[apex] = tryPeak`).
        let shared_apex = env(5, 1.5, 100.0, 2, 500.0);
        let mut stored = ChromatographicPeak::from_identification(id("UNWANTED", "UNWANTED", 2, false), 1000.0);
        stored.isotopic_envelopes = vec![shared_apex];
        let mut wanted = ChromatographicPeak::from_identification(id("WANTED", "WANTED", 3, false), 1000.0);
        wanted.isotopic_envelopes = vec![shared_apex];

        let quantify: HashSet<String> = ["WANTED".to_string()].into_iter().collect();
        let out = run_error_checking(vec![stored, wanted], &quantify, true);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].first_modified_sequence(), Some("WANTED"));
    }

    #[test]
    fn calculate_peptide_results_picks_most_intense_peak_per_file() {
        // Two MSMS peaks for the same peptide in one file (e.g. an MS2 id at two RTs); the more
        // intense one fixes the intensity. A decoy peak and an unquantified sequence are ignored.
        let mut hi = peak_with(id("SEQ", "SEQ", 2, false), vec![env(5, 1.5, 100.0, 2, 500.0)]);
        hi.detection_type = DetectionType::MSMS;
        let mut lo = peak_with(id("SEQ", "SEQ", 2, false), vec![env(8, 2.0, 40.0, 2, 500.0)]);
        lo.detection_type = DetectionType::MSMS;
        let decoy = peak_with(id("SEQ", "SEQ", 2, true), vec![env(9, 2.1, 999.0, 2, 500.0)]);

        let mut by_file: HashMap<String, Vec<ChromatographicPeak>> = HashMap::new();
        by_file.insert("fileA".to_string(), vec![hi, lo, decoy]);

        let quantify: HashSet<String> = ["SEQ".to_string()].into_iter().collect();
        let results = calculate_peptide_results(&by_file, &quantify, &quantify);

        let cell = results.get("SEQ", "fileA").expect("SEQ quantified");
        assert!((cell.intensity - 100.0).abs() < 1e-9);
        assert_eq!(cell.detection_type, DetectionType::MSMS);
        assert!((cell.retention_time - 1.5).abs() < 1e-6);
    }

    #[test]
    fn calculate_peptide_results_marks_undetected_sequences() {
        // A sequence in the quantify set with no peaks in the file stays NotDetected at intensity 0.
        let by_file: HashMap<String, Vec<ChromatographicPeak>> = {
            let mut m = HashMap::new();
            m.insert("fileA".to_string(), Vec::new());
            m
        };
        let quantify: HashSet<String> = ["GHOST".to_string()].into_iter().collect();
        let results = calculate_peptide_results(&by_file, &quantify, &quantify);

        let cell = results.get("GHOST", "fileA").expect("GHOST has a cell");
        assert_eq!(cell.intensity, 0.0);
        assert_eq!(cell.detection_type, DetectionType::NotDetected);
    }

    #[test]
    fn calculate_peptide_results_zeroes_ambiguous_sequence() {
        // An ambiguous peak (two full sequences) explaining all of a sequence's intensity zeroes
        // that sequence and marks it MSMSAmbiguousPeakfinding (quantifyAmbiguousPeptides = false).
        let ambiguous = {
            let mut p = ChromatographicPeak::from_identifications(
                vec![
                    id("AMB1", "AMB", 2, false),
                    id("AMB2", "AMB", 2, false),
                ],
                vec![1000.0, 1000.0],
                DetectionType::MSMS,
            );
            p.isotopic_envelopes = vec![env(5, 1.5, 50.0, 2, 500.0)];
            p.calculate_intensity_for_this_feature(false);
            p
        };

        let mut by_file: HashMap<String, Vec<ChromatographicPeak>> = HashMap::new();
        by_file.insert("fileA".to_string(), vec![ambiguous]);

        let quantify: HashSet<String> =
            ["AMB1".to_string(), "AMB2".to_string()].into_iter().collect();
        let results = calculate_peptide_results(&by_file, &quantify, &quantify);

        // No unambiguous intensity was recorded, so fraction = 50/(0+50) = 1 > 0.3 -> zeroed.
        let cell = results.get("AMB1", "fileA").expect("AMB1 has a cell");
        assert_eq!(cell.intensity, 0.0);
        assert_eq!(cell.detection_type, DetectionType::MSMSAmbiguousPeakfinding);
    }

    #[test]
    fn default_quantify_set_excludes_decoys() {
        let ids = vec![
            id("T1", "T1", 2, false),
            id("T2", "T2", 2, false),
            id("D1", "D1", 2, true),
        ];
        let set = default_quantify_set(&ids);
        assert!(set.contains("T1"));
        assert!(set.contains("T2"));
        assert!(!set.contains("D1"));
        assert_eq!(set.len(), 2);
    }

    // --- P2.2: decoy / quantify-set hardening ---

    #[test]
    fn engine_quantify_set_includes_decoys() {
        // The C# FlashLfqEngine default includes every modified sequence, decoys included
        // (FlashLfqEngine.cs:94) — this is the membership gate, not the output set.
        let ids = vec![
            id("T1", "T1", 2, false),
            id("D1", "D1", 2, true),
            id("D1", "D1", 3, true), // same decoy modseq, second charge -> one entry
        ];
        let set = engine_quantify_set(&ids);
        assert!(set.contains("T1"));
        assert!(set.contains("D1"), "decoy modseq belongs in the quantify set");
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn output_peptide_sequences_excludes_decoys_and_unquantified() {
        // Output rows = non-decoy ids whose modseq is in the quantify set (FlashLfqResults.cs:42).
        let ids = vec![
            id("T1", "T1", 2, false),
            id("T2", "T2", 2, false), // non-decoy but absent from the quantify set below
            id("D1", "D1", 2, true),  // decoy: in the quantify set, but never an output row
        ];
        let quantify: HashSet<String> = ["T1".to_string(), "D1".to_string()].into_iter().collect();
        let out = output_peptide_sequences(&ids, &quantify);
        assert!(out.contains("T1"));
        assert!(!out.contains("T2"), "T2 isn't in the quantify set");
        assert!(!out.contains("D1"), "decoy-only modseq is never an output row");
        assert_eq!(out.len(), 1);
    }

    #[test]
    fn calculate_peptide_results_omits_decoy_only_sequence_row() {
        // A decoy modseq sits in the quantify set (engine default) but is NOT an output sequence,
        // and its decoy peak is filtered out — so it produces no row at all, even NotDetected.
        let decoy_peak = peak_with(id("DEC", "DEC", 2, true), vec![env(5, 1.5, 80.0, 2, 500.0)]);
        let target_peak = peak_with(id("TGT", "TGT", 2, false), vec![env(5, 1.5, 60.0, 2, 510.0)]);

        let mut by_file: HashMap<String, Vec<ChromatographicPeak>> = HashMap::new();
        by_file.insert("fileA".to_string(), vec![decoy_peak, target_peak]);

        // quantify set includes the decoy (engine default); output set excludes it.
        let quantify: HashSet<String> =
            ["DEC".to_string(), "TGT".to_string()].into_iter().collect();
        let output: HashSet<String> = ["TGT".to_string()].into_iter().collect();
        let results = calculate_peptide_results(&by_file, &quantify, &output);

        assert!(results.get("DEC", "fileA").is_none(), "no row for a decoy-only sequence");
        let tgt = results.get("TGT", "fileA").expect("target quantified");
        assert!((tgt.intensity - 60.0).abs() < 1e-9);
        assert_eq!(tgt.detection_type, DetectionType::MSMS);
    }

    #[test]
    fn run_error_checking_decoy_in_quantify_set_overwrites_unquantified_target() {
        // A decoy peak shares an apex with a target whose sequence is NOT in the quantify set.
        // C# checks only quantify-set membership in the MSMS branch (FlashLfqEngine.cs:1333) — not
        // DecoyPeptide — so a decoy whose modseq is in the set overwrites the unquantified target.
        let shared_apex = env(5, 1.5, 100.0, 2, 500.0);
        let mut unquantified = ChromatographicPeak::from_identification(
            id("UNWANTED", "UNWANTED", 2, false),
            1000.0,
        );
        unquantified.isotopic_envelopes = vec![shared_apex];
        let mut decoy = ChromatographicPeak::from_identification(id("DEC", "DEC", 3, true), 1000.0);
        decoy.isotopic_envelopes = vec![shared_apex];

        // Decoy modseq is in the quantify set (engine default); the target sequence is not.
        let quantify: HashSet<String> = ["DEC".to_string()].into_iter().collect();
        let out = run_error_checking(vec![unquantified, decoy], &quantify, true);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].first_modified_sequence(), Some("DEC"));
        assert!(out[0].decoy_peptide());
    }
}
