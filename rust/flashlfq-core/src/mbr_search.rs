//! Match-between-runs acceptor search + orchestration → feature table (PLAN.md P3.2d).
//!
//! Ports the remaining MBR pieces of `FlashLfqEngine`:
//! - [`find_peptide_donor_files`] = `FindPeptideDonorFiles` (`FlashLfqEngine.cs:636`): pick one donor
//!   peak per quantifiable peptide and bucket it under the file it came from.
//! - [`get_random_peak`] = `GetRandomPeak` (`:822`): draw a *decoy* donor peak — a different peptide
//!   whose mass is a few hydrogens away and whose RT is far enough from the real donor.
//! - [`find_individual_acceptor_peak`] = `FindIndividualAcceptorPeak` (`:1250`): seed an acceptor
//!   peak from one isotopic envelope, trace + cut it, and score it.
//! - [`find_all_acceptor_peaks`] = `FindAllAcceptorPeaks` (`:1164`): scan the predicted-RT window of
//!   the acceptor file across every candidate charge state, returning the best-scoring acceptor peak.
//! - [`run_mbr`] = the `QuantifyMatchBetweenRunsPeaks` driver (`:873`): for every acceptor file build
//!   the [`MbrScorer`], map each donor file's peptides onto it (target + decoy hypotheses, with the
//!   window-widening retry loop), dedup, merge charge states, and emit **one feature row per
//!   transferred peak** — the table the P3.3 Python model trains on.
//!
//! ## Scope: the default, single-condition MS2 path
//!
//! This mirrors a default-`FlashLfqParameters` run (the C# golden generator's setup): one condition,
//! one biological replicate, **unfractionated** files. Consequently the fraction gate
//! (`PredictRetentionTime`'s early return) and the cross-condition fold-change branch
//! (`CalculateFoldChangeBetweenFiles`, guarded by `conditions.Distinct().Count() > 1`) never fire, and
//! `RequireMsmsIdInCondition` is off — exactly the behaviour of the two-K562-file corpus the rest of
//! Phase 3 is gated against. Those branches depend on `SpectraFileInfo` metadata the core
//! [`ChromatographicPeak`] does not model; they are intentionally omitted (documented divergences).
//!
//! ## Faithfulness notes & documented divergences
//! - The `Parallel.ForEach` over donor peaks is run **sequentially**, matching the golden generator's
//!   `MaxThreads = 1` (so the pseudo-random decoy draw and peak ordering are reproducible).
//! - The C# bug at `FlashLfqEngine.cs:986`/`:1004` — the decoy `FindAllAcceptorPeaks` call passes the
//!   *target* `rtInfo` (its `Width`) together with the decoy's `randomRt` centre — is replicated.
//! - The C# `sigma = 1` RT-distribution bug is already replicated in [`MbrScorer`].
//! - `AddPeakToConcurrentDict` keys by the peak's apex; peaks with **no apex** (no envelopes) cannot
//!   be keyed by a C# `Dictionary` either, so they are dropped here rather than throwing.
//! - The C# driver adds the winning peak to `_results.Peaks` **twice** (lines 1092 + 1110) when charge
//!   states are merged; `RunErrorChecking` dedups it by apex. We add it once.
//! - The final per-acceptor `RunErrorChecking` MBR/MSMS conflict pass is **not** re-run here: the
//!   `ApexToAcceptorFilePeakDict` guard (in [`find_individual_acceptor_peak`]) and the `msmsImsPeaks`
//!   guard (in the best-result loop) already prevent transferring onto an MS/MS-claimed apex, which
//!   are the conflicts that pass removes. The emitted feature rows are the candidate transferred peaks.

use std::collections::{HashMap, HashSet};

use crate::chromatographic_peak::{ChromatographicPeak, EnvelopePeakKey};
use crate::engine::{ISOTOPE_PPM_TOLERANCE, MISSED_SCANS_ALLOWED, NUM_ISOTOPES_REQUIRED};
use crate::isotopic_envelope::{
    get_isotopic_envelopes, mass_to_mz_f64, IsotopicEnvelope,
};
use crate::mbr::{
    get_rt_cal_spline, predict_retention_time, DonorCriterion, RtInfo, DONOR_Q_VALUE_THRESHOLD,
    MAX_MBR_RT_WINDOW, MBR_PPM_TOLERANCE, NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR,
};
use crate::mbr_chromatographic_peak::MbrChromatographicPeak;
use crate::mbr_scorer::{build_mbr_scorer, MbrScorer};
use crate::peak_indexing::{IndexedMassSpectralPeak, PeakIndexingEngine};
use crate::periodic_table::periodic_table;
use crate::psm_tsv::Identification;
use crate::theoretical_isotope_distribution::{
    expected_isotope_peaks, most_abundant_isotope_shift, ExpectedIsotopePeak,
};
use crate::tolerance::PpmTolerance;

/// Whether to integrate envelope intensities (`FlashLfqParameters.Integrate`, default `false`).
const INTEGRATE: bool = false;
/// `maxPeakHalfWidth` default for `GetXic` (`int.MaxValue` in the C# default argument).
const MAX_PEAK_HALF_WIDTH: f64 = i32::MAX as f64;

/// One row of the MBR feature table — a single candidate transferred peak. Port of the per-peak
/// features the C# `RunPEPAnalysis` / FastTree model reads off each `MbrChromatographicPeak`. The
/// P3.3 Python model trains on these rows (targets vs. `random_rt` decoys).
#[derive(Debug, Clone, PartialEq)]
pub struct FeatureRow {
    /// Modified (full) sequence of the donor peptide this peak was transferred from.
    pub donor_modified_sequence: String,
    /// Base sequence of the donor peptide.
    pub donor_base_sequence: String,
    /// Acceptor file the peak was found in (extension-stripped name).
    pub acceptor_file: String,
    /// Retention time predicted for the peptide in the acceptor run.
    pub predicted_retention_time: f64,
    /// Apex retention time of the transferred peak (`-1` when it has no apex).
    pub apex_retention_time: f64,
    /// Feature intensity of the transferred peak.
    pub intensity: f64,
    /// Ppm component score.
    pub ppm_score: f64,
    /// Intensity component score.
    pub intensity_score: f64,
    /// RT component score.
    pub rt_score: f64,
    /// Scan-count component score.
    pub scan_count_score: f64,
    /// Isotopic-distribution component score.
    pub isotopic_distribution_score: f64,
    /// Combined MBR score in `[0, 100]`.
    pub mbr_score: f64,
    /// Apex mass error in ppm.
    pub mass_error: f64,
    /// Number of isotopic envelopes traced.
    pub scan_count: usize,
    /// Apex isotopic Pearson correlation (`-1` when there is no apex).
    pub isotopic_pearson_correlation: f64,
    /// `predictedRt - apexRt`, the scorer's RT prediction error.
    pub rt_prediction_error: f64,
    /// Whether this peak's RT was randomized — i.e. it is an MBR **decoy** peak.
    pub random_rt: bool,
    /// Whether the donor identification is a decoy *peptide* (distinct from `random_rt`).
    pub decoy_peptide: bool,
}

impl FeatureRow {
    /// Builds a feature row from a scored MBR peak.
    fn from_peak(peak: &MbrChromatographicPeak, acceptor_file: &str) -> FeatureRow {
        let id = peak.peak.identifications.first();
        FeatureRow {
            donor_modified_sequence: id
                .map(|i| i.modified_sequence.clone())
                .unwrap_or_default(),
            donor_base_sequence: id.map(|i| i.base_sequence.clone()).unwrap_or_default(),
            acceptor_file: acceptor_file.to_string(),
            predicted_retention_time: peak.predicted_retention_time,
            apex_retention_time: peak.apex_retention_time(),
            intensity: peak.intensity(),
            ppm_score: peak.ppm_score,
            intensity_score: peak.intensity_score,
            rt_score: peak.rt_score,
            scan_count_score: peak.scan_count_score,
            isotopic_distribution_score: peak.isotopic_distribution_score,
            mbr_score: peak.mbr_score,
            mass_error: peak.mass_error(),
            scan_count: peak.scan_count(),
            isotopic_pearson_correlation: peak.isotopic_pearson_correlation(),
            rt_prediction_error: peak.rt_prediction_error,
            random_rt: peak.random_rt,
            decoy_peptide: peak.peak.decoy_peptide(),
        }
    }
}

/// The products of [`run_mbr`]: the flat feature table plus the surviving MBR peaks per acceptor file.
#[derive(Debug, Clone, Default)]
pub struct MbrResult {
    /// One row per candidate transferred peak (targets + `random_rt` decoys), across all acceptors.
    pub feature_rows: Vec<FeatureRow>,
    /// `acceptor_file -> transferred MBR peaks` (the peaks the rows were derived from).
    pub mbr_peaks_by_file: HashMap<String, Vec<MbrChromatographicPeak>>,
}

/// The maximum PSM score among a peak's identifications (C# `Identifications.Max(id => id.PsmScore)`).
fn max_psm_score(peak: &ChromatographicPeak) -> f64 {
    peak.identifications
        .iter()
        .map(|id| id.score)
        .fold(f64::NEG_INFINITY, f64::max)
}

/// Picks the index of the best peak among `candidates` by the default [`DonorCriterion::Score`] rule
/// (LINQ `MaxBy` over the per-peak max PSM score, first on ties); falls through to the most intense
/// peak when the chosen peak's first identification has a non-positive score. Mirrors
/// `FlashLfqEngine.ChooseBestPeak`, but returns an index so the caller can recover the source file.
fn choose_best_index(candidates: &[(String, &ChromatographicPeak)]) -> usize {
    // Score: first peak attaining the maximum per-peak max-PSM-score.
    let mut best = 0usize;
    let mut best_key = max_psm_score(candidates[0].1);
    for (i, (_, p)) in candidates.iter().enumerate().skip(1) {
        let k = max_psm_score(p);
        if k > best_key {
            best_key = k;
            best = i;
        }
    }
    let first_score = candidates[best]
        .1
        .identifications
        .first()
        .map(|id| id.score)
        .unwrap_or(0.0);
    if first_score > 0.0 {
        return best;
    }
    // Fall through to Intensity (C# `goto default`).
    let mut bi = 0usize;
    let mut bk = candidates[0].1.intensity;
    for (i, (_, p)) in candidates.iter().enumerate().skip(1) {
        if p.intensity > bk {
            bk = p.intensity;
            bi = i;
        }
    }
    bi
}

/// Selects one donor peak per quantifiable peptide and buckets it under the file it was found in.
/// Port of `FlashLfqEngine.FindPeptideDonorFiles` (`FlashLfqEngine.cs:636`).
///
/// A peak qualifies as a donor candidate when it is unambiguous (`NumIdentificationsByFullSeq == 1`),
/// actually traced (`IsotopicEnvelopes.Any()`), confidently identified (`min QValue <`
/// [`DONOR_Q_VALUE_THRESHOLD`]), and its sequence is in `peptide_sequences_to_quantify`. Candidates are
/// grouped by modified sequence (first-seen order, like LINQ `GroupBy`), one [`choose_best_index`] per
/// group, and bucketed under the source file. Files are visited in sorted name order for determinism.
pub fn find_peptide_donor_files(
    peaks_by_file: &HashMap<String, Vec<ChromatographicPeak>>,
    peptide_sequences_to_quantify: &HashSet<String>,
) -> HashMap<String, Vec<ChromatographicPeak>> {
    let mut files: Vec<&String> = peaks_by_file.keys().collect();
    files.sort();

    let mut seq_order: Vec<String> = Vec::new();
    let mut groups: HashMap<String, Vec<(String, &ChromatographicPeak)>> = HashMap::new();

    for file in &files {
        for peak in &peaks_by_file[*file] {
            if peak.num_identifications_by_full_seq != 1 || peak.isotopic_envelopes.is_empty() {
                continue;
            }
            let min_q = peak
                .identifications
                .iter()
                .map(|id| id.q_value)
                .fold(f64::INFINITY, f64::min);
            if !(min_q < DONOR_Q_VALUE_THRESHOLD) {
                continue;
            }
            let seq = match peak.first_modified_sequence() {
                Some(s) => s.to_string(),
                None => continue,
            };
            if !peptide_sequences_to_quantify.contains(&seq) {
                continue;
            }
            if !groups.contains_key(&seq) {
                seq_order.push(seq.clone());
            }
            groups.entry(seq).or_default().push(((*file).clone(), peak));
        }
    }

    let mut donor_file_to_peaks: HashMap<String, Vec<ChromatographicPeak>> = HashMap::new();
    for seq in &seq_order {
        let cands = &groups[seq];
        if cands.is_empty() {
            continue;
        }
        let best = choose_best_index(cands);
        let (file, peak) = &cands[best];
        donor_file_to_peaks
            .entry(file.clone())
            .or_default()
            .push((*peak).clone());
    }
    donor_file_to_peaks
}

/// Draws a *decoy* donor peak: a different peptide whose peakfinding mass is between 5 and 11
/// hydrogens away from the donor's (widening by ×10 up to `1e5` if needed) and whose apex RT differs
/// from the donor's by more than `retention_time_min_diff`. Port of `FlashLfqEngine.GetRandomPeak`
/// (`FlashLfqEngine.cs:822`). Returns `None` when no candidate exists even after widening.
///
/// `peaks_ordered_by_mass` is the donor file's best peaks ordered by peakfinding mass (the spline's
/// `donor_best_peaks_ordered_by_mass`). The selection index is the verbatim pseudo-random draw
/// `(int)(1e5·(pfm mod 1)·(ms2Rt mod 1)) mod count`.
pub fn get_random_peak<'a>(
    peaks_ordered_by_mass: &'a [ChromatographicPeak],
    donor_peak_retention_time: f64,
    retention_time_min_diff: f64,
    donor_identification: &Identification,
    donor_peakfinding_mass: f64,
) -> Option<&'a ChromatographicPeak> {
    let h_mass = periodic_table()
        .element_by_symbol("H")
        .expect("hydrogen is in the periodic table")
        .principal_isotope()
        .atomic_mass;
    let min_diff = 5.0 * h_mass;
    let mut max_diff = 11.0 * h_mass;

    let candidates = |max_diff: f64| -> Vec<&ChromatographicPeak> {
        peaks_ordered_by_mass
            .iter()
            .filter(|p| {
                let apex_rt = p.apex_retention_time();
                let first = match p.identifications.first() {
                    Some(id) => id,
                    None => return false,
                };
                let pfm = match p.identification_peakfinding_masses.first() {
                    Some(&m) => m,
                    None => return false,
                };
                apex_rt > 0.0
                    && (apex_rt - donor_peak_retention_time).abs() > retention_time_min_diff
                    && first.base_sequence != donor_identification.base_sequence
                    && (pfm - donor_peakfinding_mass).abs() > min_diff
                    && (pfm - donor_peakfinding_mass).abs() < max_diff
            })
            .collect()
    };

    let mut random_peak_candidates = candidates(max_diff);
    while random_peak_candidates.is_empty() && max_diff < 1e5 {
        max_diff *= 10.0;
        random_peak_candidates = candidates(max_diff);
    }
    if random_peak_candidates.is_empty() {
        return None;
    }

    // Pseudo-random index from the donor peakfinding mass + RT (C# (int) truncation).
    let pseudo = (1e5
        * (donor_peakfinding_mass % 1.0)
        * (donor_identification.ms2_retention_time_in_minutes % 1.0)) as i64;
    let idx = (pseudo.rem_euclid(random_peak_candidates.len() as i64)) as usize;
    Some(random_peak_candidates[idx])
}

/// Seeds an acceptor peak from the first (least intense) envelope in `charge_envelopes`, traces +
/// cuts it, removes the claimed envelopes from `charge_envelopes`, and scores it. Port of
/// `FlashLfqEngine.FindIndividualAcceptorPeak` (`FlashLfqEngine.cs:1250`). Returns `None` when the
/// seed envelope's apex is already claimed by an MS/MS identification (the `ApexToAcceptorFilePeakDict`
/// guard). `charge_envelopes` is always shrunk (the claimed seed is removed) so the caller's loop
/// terminates even on a `None`.
#[allow(clippy::too_many_arguments)]
fn find_individual_acceptor_peak(
    engine: &PeakIndexingEngine,
    scorer: &MbrScorer,
    donor_peak: &ChromatographicPeak,
    donor_file: &str,
    acceptor_file: &str,
    mbr_tol: &PpmTolerance,
    rt_info: &RtInfo,
    z: i32,
    charge_envelopes: &mut Vec<IsotopicEnvelope>,
    random_rt: Option<f64>,
    dist: &[ExpectedIsotopePeak],
    monoisotopic_mass: f64,
    peakfinding_mass: f64,
) -> Option<MbrChromatographicPeak> {
    // donorId = Identifications.OrderBy(QValue).First()
    let donor_id = donor_peak
        .identifications
        .iter()
        .min_by(|a, b| {
            a.q_value
                .partial_cmp(&b.q_value)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .expect("donor peak has at least one identification");

    let predicted = random_rt.unwrap_or(rt_info.predicted_rt);
    let mut acceptor_peak = MbrChromatographicPeak::new(
        donor_id.clone(),
        acceptor_file.to_string(),
        peakfinding_mass,
        predicted,
        random_rt.is_some(),
    );

    let seed_env = charge_envelopes[0];
    let xic = engine.get_xic(
        mass_to_mz_f64(peakfinding_mass, z),
        seed_env.indexed_peak.retention_time as f64,
        mbr_tol,
        MISSED_SCANS_ALLOWED,
        MAX_PEAK_HALF_WIDTH,
        None,
    );
    let best_charge_envelopes = get_isotopic_envelopes(
        engine,
        &xic,
        dist,
        monoisotopic_mass,
        peakfinding_mass,
        z,
        NUM_ISOTOPES_REQUIRED,
        ISOTOPE_PPM_TOLERANCE,
    );
    acceptor_peak.peak.isotopic_envelopes.extend(best_charge_envelopes);
    acceptor_peak
        .peak
        .calculate_intensity_for_this_feature(INTEGRATE);
    acceptor_peak
        .peak
        .cut_peak(seed_env.indexed_peak.retention_time as f64, INTEGRATE);

    // Claimed peaks = this peak's envelopes + the seed (the seed prevents infinite loops).
    let mut claimed: HashSet<EnvelopePeakKey> = acceptor_peak
        .peak
        .isotopic_envelopes
        .iter()
        .map(EnvelopePeakKey::from)
        .collect();
    claimed.insert(EnvelopePeakKey::from(&seed_env));
    charge_envelopes.retain(|p| !claimed.contains(&EnvelopePeakKey::from(p)));

    // Peak already identified by MS/MS - skip it.
    if scorer
        .apex_to_acceptor_file_peak
        .contains_key(&EnvelopePeakKey::from(&seed_env))
    {
        return None;
    }

    acceptor_peak.mbr_score = scorer.score_mbr(&mut acceptor_peak, donor_peak, donor_file, predicted);
    Some(acceptor_peak)
}

/// Searches the acceptor file's predicted-RT window across every candidate charge state and returns
/// the best-scoring acceptor peak (or `None`). Port of `FlashLfqEngine.FindAllAcceptorPeaks`
/// (`FlashLfqEngine.cs:1164`).
///
/// `random_rt`, when `Some`, centres the RT window on the decoy's predicted RT (using `rt_info.width`
/// — replicating the C# call that passes the target `rtInfo`). `dist`/`monoisotopic_mass`/
/// `peakfinding_mass` are the donor peptide's theoretical envelope + masses.
#[allow(clippy::too_many_arguments)]
fn find_all_acceptor_peaks(
    engine: &PeakIndexingEngine,
    scorer: &MbrScorer,
    rt_info: &RtInfo,
    mbr_tol: &PpmTolerance,
    donor_peak: &ChromatographicPeak,
    donor_file: &str,
    acceptor_file: &str,
    random_rt: Option<f64>,
    dist: &[ExpectedIsotopePeak],
    monoisotopic_mass: f64,
    peakfinding_mass: f64,
) -> Option<MbrChromatographicPeak> {
    let scan_infos = engine.scan_info();
    if scan_infos.is_empty() {
        return None;
    }

    let rt_start = match random_rt {
        None => rt_info.rt_start_hypothesis(),
        Some(r) => r - rt_info.width / 2.0,
    };
    let rt_end = match random_rt {
        None => rt_info.rt_end_hypothesis(),
        Some(r) => r + rt_info.width / 2.0,
    };

    // Snip the MS1 scans to the region the analyte should appear in.
    let mut start = scan_infos[0];
    let mut end = scan_infos[scan_infos.len() - 1];
    for scan in scan_infos {
        if scan.retention_time <= rt_start {
            start = *scan;
        }
        if scan.retention_time >= rt_end {
            end = *scan;
            break;
        }
    }

    // Charges to match = distinct donor precursor charges + the donor apex charge.
    let mut charges_to_match: Vec<i32> = Vec::new();
    for id in &donor_peak.identifications {
        if !charges_to_match.contains(&id.precursor_charge_state) {
            charges_to_match.push(id.precursor_charge_state);
        }
    }
    if let Some(apex) = donor_peak.apex {
        if !charges_to_match.contains(&apex.charge_state) {
            charges_to_match.push(apex.charge_state);
        }
    }

    let mut best_acceptor: Option<MbrChromatographicPeak> = None;

    for &z in &charges_to_match {
        let mut charge_xic: Vec<IndexedMassSpectralPeak> = Vec::new();
        for j in start.zero_based_scan_index..=end.zero_based_scan_index {
            if let Some(peak) =
                engine.get_indexed_peak(mass_to_mz_f64(peakfinding_mass, z), j, mbr_tol)
            {
                charge_xic.push(*peak);
            }
        }
        if charge_xic.is_empty() {
            continue;
        }

        let mut charge_envelopes = get_isotopic_envelopes(
            engine,
            &charge_xic,
            dist,
            monoisotopic_mass,
            peakfinding_mass,
            z,
            NUM_ISOTOPES_REQUIRED,
            ISOTOPE_PPM_TOLERANCE,
        );
        // OrderBy(env => env.Intensity): ascending; First() is the least intense (C# quirk).
        charge_envelopes.sort_by(|a, b| {
            a.intensity
                .partial_cmp(&b.intensity)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        while !charge_envelopes.is_empty() {
            let acceptor_peak = find_individual_acceptor_peak(
                engine,
                scorer,
                donor_peak,
                donor_file,
                acceptor_file,
                mbr_tol,
                rt_info,
                z,
                &mut charge_envelopes,
                random_rt,
                dist,
                monoisotopic_mass,
                peakfinding_mass,
            );
            let acceptor_peak = match acceptor_peak {
                Some(p) => p,
                None => continue,
            };
            let better = match &best_acceptor {
                None => true,
                Some(b) => b.mbr_score < acceptor_peak.mbr_score,
            };
            if better {
                let mut p = acceptor_peak;
                p.charge_list = charges_to_match.clone();
                best_acceptor = Some(p);
            }
        }
    }

    best_acceptor
}

/// `modified_sequence -> (apex key -> candidate peaks)` accumulator (the C#
/// `matchBetweenRunsIdentifiedPeaks` concurrent dictionary).
type MbrPeakDict = HashMap<String, HashMap<EnvelopePeakKey, Vec<MbrChromatographicPeak>>>;

/// Adds a found peak to the accumulator, keyed by donor modified sequence and the peak's apex (C#
/// `AddPeakToConcurrentDict`). Peaks with no apex cannot be keyed and are dropped.
fn add_peak_to_dict(
    dict: &mut MbrPeakDict,
    modified_sequence: &str,
    peak: &Option<MbrChromatographicPeak>,
) {
    if let Some(peak) = peak {
        if let Some(apex) = peak.peak.apex {
            let key = EnvelopePeakKey::from(&apex);
            dict.entry(modified_sequence.to_string())
                .or_default()
                .entry(key)
                .or_default()
                .push(peak.clone());
        }
    }
}

/// Runs the full MBR transfer for one acceptor file. Port of
/// `FlashLfqEngine.QuantifyMatchBetweenRunsPeaks` (`FlashLfqEngine.cs:873`), single-condition path.
fn quantify_mbr_for_acceptor(
    acceptor_file: &str,
    peaks_by_file: &HashMap<String, Vec<ChromatographicPeak>>,
    engines_by_file: &HashMap<String, PeakIndexingEngine>,
    donor_file_to_peaks: &HashMap<String, Vec<ChromatographicPeak>>,
    peptide_sequences_to_quantify: &HashSet<String>,
) -> Vec<MbrChromatographicPeak> {
    let acceptor_peaks = match peaks_by_file.get(acceptor_file) {
        Some(p) => p,
        None => return Vec::new(),
    };
    let acceptor_engine = match engines_by_file.get(acceptor_file) {
        Some(e) => e,
        None => return Vec::new(),
    };

    // Acceptor peaks whose ids contain a quantifiable peptide (scorer is fit from these).
    let acceptor_identified_peaks: Vec<ChromatographicPeak> = acceptor_peaks
        .iter()
        .filter(|p| {
            p.identifications
                .iter()
                .any(|id| peptide_sequences_to_quantify.contains(&id.modified_sequence))
        })
        .cloned()
        .collect();

    // Sequences already confidently identified in the acceptor (no need to transfer them).
    let mut acceptor_identified_sequences: HashSet<String> = HashSet::new();
    for peak in &acceptor_identified_peaks {
        if peak.isotopic_envelopes.is_empty() {
            continue;
        }
        let min_q = peak
            .identifications
            .iter()
            .map(|id| id.q_value)
            .fold(f64::INFINITY, f64::min);
        if min_q < 0.01 {
            for id in &peak.identifications {
                acceptor_identified_sequences.insert(id.modified_sequence.clone());
            }
        }
    }

    let mut scorer = match build_mbr_scorer(acceptor_identified_peaks, MBR_PPM_TOLERANCE) {
        Some((s, _file_tol)) => s,
        None => return Vec::new(),
    };
    // C# overwrites the file-specific tolerance with the flat MbrPpmTolerance (FlashLfqEngine.cs:895).
    let mbr_tol = PpmTolerance::new(MBR_PPM_TOLERANCE);

    let mut dict: MbrPeakDict = HashMap::new();

    // Map each donor file onto this acceptor (sorted for determinism).
    let mut donor_files: Vec<&String> = donor_file_to_peaks.keys().collect();
    donor_files.sort();

    for donor_file in donor_files {
        if donor_file == acceptor_file {
            continue;
        }
        let donor_best_peaks = &donor_file_to_peaks[donor_file];

        // Donor peaks not already identified in the acceptor and still quantifiable.
        let id_donor_peaks: Vec<&ChromatographicPeak> = donor_best_peaks
            .iter()
            .filter(|p| {
                let seq = match p.first_modified_sequence() {
                    Some(s) => s,
                    None => return false,
                };
                !acceptor_identified_sequences.contains(seq)
                    && peptide_sequences_to_quantify.contains(seq)
            })
            .collect();
        if id_donor_peaks.is_empty() {
            continue;
        }

        // RT calibration spline over the full donor/acceptor file peaks (not just the best peaks).
        let donor_all = match peaks_by_file.get(donor_file) {
            Some(p) => p,
            None => continue,
        };
        let spline = get_rt_cal_spline(
            donor_all,
            acceptor_peaks,
            DONOR_Q_VALUE_THRESHOLD,
            DonorCriterion::Score,
        );
        scorer.add_rt_pred_error_distribution(
            donor_file,
            &spline.anchor_rt_diffs,
            NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR,
        );
        if !scorer.is_valid_for_donor(donor_file) {
            continue;
        }
        let donor_peaks_mass_ordered = &spline.donor_best_peaks_ordered_by_mass;

        for donor_peak in &id_donor_peaks {
            let first_id = donor_peak
                .identifications
                .first()
                .expect("donor peak has an identification");
            let modseq = first_id.modified_sequence.clone();
            let dist = expected_isotope_peaks(
                None,
                &first_id.base_sequence,
                first_id.monoisotopic_mass,
                NUM_ISOTOPES_REQUIRED,
            );
            let mono = first_id.monoisotopic_mass;
            let pfm = mono
                + most_abundant_isotope_shift(&dist)
                    .expect("theoretical distribution normalizes a peak to abundance 1.0");

            let mut rt_info = predict_retention_time(
                &spline.calibration_curve,
                donor_peak,
                MAX_MBR_RT_WINDOW,
                NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR,
            );

            // Target: search the predicted-RT window.
            let mut best_acceptor = find_all_acceptor_peaks(
                acceptor_engine,
                &scorer,
                &rt_info,
                &mbr_tol,
                donor_peak,
                donor_file,
                acceptor_file,
                None,
                &dist,
                mono,
                pfm,
            );
            add_peak_to_dict(&mut dict, &modseq, &best_acceptor);

            // Decoy: draw a random donor far enough away, predict its RT, search at that centre.
            let minimum_rt_difference = rt_info.width * 2.0;
            let random_donor = get_random_peak(
                donor_peaks_mass_ordered,
                donor_peak.apex_retention_time(),
                minimum_rt_difference,
                first_id,
                pfm,
            );
            let mut best_decoy: Option<MbrChromatographicPeak> = None;
            let mut decoy_rt_info: Option<RtInfo> = None;
            if let Some(random_donor) = random_donor {
                let drt = predict_retention_time(
                    &spline.calibration_curve,
                    random_donor,
                    MAX_MBR_RT_WINDOW,
                    NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR,
                );
                decoy_rt_info = Some(drt);
                // C# bug: passes the *target* rtInfo (its width) with the decoy's randomRt centre.
                best_decoy = find_all_acceptor_peaks(
                    acceptor_engine,
                    &scorer,
                    &rt_info,
                    &mbr_tol,
                    donor_peak,
                    donor_file,
                    acceptor_file,
                    Some(drt.predicted_rt),
                    &dist,
                    mono,
                    pfm,
                );
                add_peak_to_dict(&mut dict, &modseq, &best_decoy);
            }

            // Widen the window and retry while nothing was found.
            let mut window_width = 0.5_f64.max(rt_info.width);
            while best_acceptor.is_none() && best_decoy.is_none() {
                window_width = window_width.min(MAX_MBR_RT_WINDOW);
                rt_info.width = window_width;
                best_acceptor = find_all_acceptor_peaks(
                    acceptor_engine,
                    &scorer,
                    &rt_info,
                    &mbr_tol,
                    donor_peak,
                    donor_file,
                    acceptor_file,
                    None,
                    &dist,
                    mono,
                    pfm,
                );
                add_peak_to_dict(&mut dict, &modseq, &best_acceptor);

                if let Some(mut drt) = decoy_rt_info {
                    drt.width = window_width;
                    decoy_rt_info = Some(drt);
                    best_decoy = find_all_acceptor_peaks(
                        acceptor_engine,
                        &scorer,
                        &rt_info,
                        &mbr_tol,
                        donor_peak,
                        donor_file,
                        acceptor_file,
                        Some(drt.predicted_rt),
                        &dist,
                        mono,
                        pfm,
                    );
                    add_peak_to_dict(&mut dict, &modseq, &best_decoy);
                }

                if window_width >= MAX_MBR_RT_WINDOW {
                    break;
                } else {
                    window_width *= 2.0;
                }
            }
        }
    }

    finalize_acceptor_peaks(&mut dict, acceptor_peaks, &acceptor_identified_sequences, peptide_sequences_to_quantify)
}

/// Dedups, resolves MS/MS conflicts, merges charge states, and selects the surviving MBR peaks from
/// the accumulator. Port of the post-processing in `QuantifyMatchBetweenRunsPeaks`
/// (`FlashLfqEngine.cs:1022`–`:1112`).
fn finalize_acceptor_peaks(
    dict: &mut MbrPeakDict,
    acceptor_peaks: &[ChromatographicPeak],
    acceptor_identified_sequences: &HashSet<String>,
    peptide_sequences_to_quantify: &HashSet<String>,
) -> Vec<MbrChromatographicPeak> {
    // Dedup: within each apex, keep the single highest-scoring peak (peaks share one modseq here).
    for apex_map in dict.values_mut() {
        for peak_list in apex_map.values_mut() {
            if let Some(best_idx) = max_score_index(peak_list) {
                let best = peak_list[best_idx].clone();
                peak_list.clear();
                peak_list.push(best);
            }
        }
    }

    // msmsImsPeaks: scan index -> set of MS/MS apex peaks already claimed in the acceptor.
    let mut msms_ims: HashMap<i32, HashSet<EnvelopePeakKey>> = HashMap::new();
    for peak in acceptor_peaks {
        if peak.decoy_peptide() {
            continue;
        }
        let apex = match peak.apex {
            Some(a) => a,
            None => continue,
        };
        let quantifiable = peak
            .first_modified_sequence()
            .map(|s| peptide_sequences_to_quantify.contains(s))
            .unwrap_or(false);
        if !quantifiable {
            continue;
        }
        msms_ims
            .entry(apex.indexed_peak.zero_based_scan_index)
            .or_default()
            .insert(EnvelopePeakKey::from(&apex));
    }

    let is_msms_claimed = |peak: &MbrChromatographicPeak| -> bool {
        match peak.peak.apex {
            Some(apex) => msms_ims
                .get(&apex.indexed_peak.zero_based_scan_index)
                .map(|set| set.contains(&EnvelopePeakKey::from(&apex)))
                .unwrap_or(false),
            None => false,
        }
    };

    let mut result: Vec<MbrChromatographicPeak> = Vec::new();

    // Take the best result (per RandomRt group) for each transferred peptide.
    let mut modseqs: Vec<&String> = dict.keys().collect();
    modseqs.sort();
    for modseq in modseqs {
        if acceptor_identified_sequences.contains(modseq) {
            continue;
        }
        let apex_map = &dict[modseq];

        // Flatten all apexes' (deduped) best peaks, sort by score descending.
        let mut all_best: Vec<MbrChromatographicPeak> =
            apex_map.values().filter_map(|v| v.first().cloned()).collect();
        all_best.sort_by(|a, b| {
            b.mbr_score
                .partial_cmp(&a.mbr_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // GroupBy(RandomRt), preserving first-appearance order.
        let mut group_keys: Vec<bool> = Vec::new();
        let mut groups: HashMap<bool, Vec<MbrChromatographicPeak>> = HashMap::new();
        for p in all_best {
            if !groups.contains_key(&p.random_rt) {
                group_keys.push(p.random_rt);
            }
            groups.entry(p.random_rt).or_default().push(p);
        }

        for key in group_keys {
            let mut peak_hypotheses = groups.remove(&key).unwrap();
            // best = First; remove from hypotheses.
            let mut best: Option<MbrChromatographicPeak> = if peak_hypotheses.is_empty() {
                None
            } else {
                Some(peak_hypotheses.remove(0))
            };

            // Discard a best already claimed by an MS/MS apex, stepping to the next hypothesis.
            loop {
                match &best {
                    Some(b) if is_msms_claimed(b) => {
                        if peak_hypotheses.is_empty() {
                            best = None;
                            break;
                        }
                        best = Some(peak_hypotheses.remove(0));
                    }
                    _ => break,
                }
            }
            let mut best = match best {
                Some(b) => b,
                None => continue,
            };

            // Merge in compatible (different-charge, in-RT-range, unclaimed) charge states.
            if !peak_hypotheses.is_empty() {
                let start = best
                    .peak
                    .isotopic_envelopes
                    .iter()
                    .map(|e| e.indexed_peak.retention_time as f64)
                    .fold(f64::INFINITY, f64::min);
                let end = best
                    .peak
                    .isotopic_envelopes
                    .iter()
                    .map(|e| e.indexed_peak.retention_time as f64)
                    .fold(f64::NEG_INFINITY, f64::max);
                let best_charge = best.peak.apex.map(|a| a.charge_state);
                for peak in &peak_hypotheses {
                    let peak_apex = match peak.peak.apex {
                        Some(a) => a,
                        None => continue,
                    };
                    if Some(peak_apex.charge_state) == best_charge {
                        continue;
                    }
                    let rt = peak_apex.indexed_peak.retention_time as f64;
                    if rt >= start && rt <= end {
                        if is_msms_claimed(peak) {
                            continue;
                        }
                        best.peak.merge_feature_with(&peak.peak, INTEGRATE);
                    }
                }
            }

            result.push(best);
        }
    }

    result
}

/// Index of the highest-`mbr_score` peak (LINQ `MaxBy`, first on ties), or `None` for an empty slice.
fn max_score_index(peaks: &[MbrChromatographicPeak]) -> Option<usize> {
    if peaks.is_empty() {
        return None;
    }
    let mut best = 0usize;
    let mut best_key = peaks[0].mbr_score;
    for (i, p) in peaks.iter().enumerate().skip(1) {
        if p.mbr_score > best_key {
            best_key = p.mbr_score;
            best = i;
        }
    }
    Some(best)
}

/// Runs match-between-runs across every file in `peaks_by_file`, treating each in turn as the acceptor
/// and the others as donors, and returns the feature table + surviving MBR peaks. Driver entry for
/// PLAN.md P3.2d (the per-acceptor [`quantify_mbr_for_acceptor`] over all files).
///
/// `peptide_sequences_to_quantify` is the engine's `PeptideModifiedSequencesToQuantify`
/// (non-decoy modified sequences in the quantify set). Acceptor files are processed in sorted name
/// order for determinism.
pub fn run_mbr(
    peaks_by_file: &HashMap<String, Vec<ChromatographicPeak>>,
    engines_by_file: &HashMap<String, PeakIndexingEngine>,
    peptide_sequences_to_quantify: &HashSet<String>,
) -> MbrResult {
    let donor_file_to_peaks = find_peptide_donor_files(peaks_by_file, peptide_sequences_to_quantify);

    let mut feature_rows: Vec<FeatureRow> = Vec::new();
    let mut mbr_peaks_by_file: HashMap<String, Vec<MbrChromatographicPeak>> = HashMap::new();

    let mut files: Vec<&String> = peaks_by_file.keys().collect();
    files.sort();
    for acceptor in files {
        let peaks = quantify_mbr_for_acceptor(
            acceptor,
            peaks_by_file,
            engines_by_file,
            &donor_file_to_peaks,
            peptide_sequences_to_quantify,
        );
        for p in &peaks {
            feature_rows.push(FeatureRow::from_peak(p, acceptor));
        }
        mbr_peaks_by_file.insert(acceptor.clone(), peaks);
    }

    MbrResult {
        feature_rows,
        mbr_peaks_by_file,
    }
}

/// Error from [`apply_mbr_pep`] when the supplied pep slice does not line up with the feature table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MbrPepLengthMismatch {
    /// Number of peps supplied.
    pub peps: usize,
    /// Number of feature rows (== number of MBR peaks) in the result.
    pub rows: usize,
}

impl std::fmt::Display for MbrPepLengthMismatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "expected one pep per feature row: got {} peps for {} rows",
            self.peps, self.rows
        )
    }
}

impl std::error::Error for MbrPepLengthMismatch {}

/// Writes the Python PEP model's scores back onto the MBR peaks (PLAN.md P3.3 — "scores returned
/// to Rust"). `peps` must be in **canonical feature-table order** — i.e. parallel to the rows of
/// [`crate::parquet_output::feature_table_record_batch`], which is the exact table the binding
/// hands to the Python model — so each pep maps unambiguously back onto the peak its row came from.
/// Sets `MbrChromatographicPeak::mbr_pep` on every transferred peak (the value the C#
/// `Compute_PEP_For_All_Peaks` stores: `1 - probability`). The FDR step (P3.4) then reads it.
pub fn apply_mbr_pep(result: &mut MbrResult, peps: &[f64]) -> Result<(), MbrPepLengthMismatch> {
    if peps.len() != result.feature_rows.len() {
        return Err(MbrPepLengthMismatch {
            peps: peps.len(),
            rows: result.feature_rows.len(),
        });
    }

    // `peps[k]` is the score for the row at canonical-sorted position `k`; `order[k]` is the
    // feature-row index (== the flattened peak index) that landed there. Invert to a per-row pep.
    let order = crate::parquet_output::feature_table_sort_order(&result.feature_rows);
    let mut pep_for_row = vec![0.0f64; result.feature_rows.len()];
    for (k, &row_index) in order.iter().enumerate() {
        pep_for_row[row_index] = peps[k];
    }

    // The flattened peak order matches the feature-row order: `run_mbr` pushes one row per peak,
    // iterating acceptor files in sorted name order then each file's peak vec in order.
    let mut files: Vec<String> = result.mbr_peaks_by_file.keys().cloned().collect();
    files.sort();
    let mut row_index = 0usize;
    for file in &files {
        let peaks = result
            .mbr_peaks_by_file
            .get_mut(file)
            .expect("file key enumerated from the map");
        for peak in peaks.iter_mut() {
            peak.mbr_pep = Some(pep_for_row[row_index]);
            row_index += 1;
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------------------------
// MBR FDR (PLAN.md P3.4) — port of FlashLfqEngine.CalculateFdrForMbrPeaks / EstimateFdr /
// CorrectQValues (FlashLfqEngine.cs:1426–1528).
// ---------------------------------------------------------------------------------------------

/// `Math.Round(x, 6)` — round to 6 decimal places, half-to-even (the C# default `MidpointRounding`).
fn round6(x: f64) -> f64 {
    (x * 1e6).round_ties_even() / 1e6
}

/// Decoy-peptide error count `max(0, decoyPeptides - doubleDecoys)` (C# `EstimateDecoyPeptideErrors`,
/// `:1494`). `usize::saturating_sub` is exactly the `Math.Max(0, ...)`.
fn estimate_decoy_peptide_errors(decoy_peptides: usize, double_decoys: usize) -> usize {
    decoy_peptides.saturating_sub(double_decoys)
}

/// Running MBR FDR estimate `(1 + decoyPeaks + max(0, decoyPeptides - doubleDecoys)) / totalPeaks`
/// (C# `EstimateFdr`, `:1499`). Targets contribute nothing; `random_rt` decoy peaks count directly,
/// decoy-peptide peaks count via the double-decoy correction.
fn estimate_fdr(
    double_decoys: usize,
    decoy_peptides: usize,
    decoy_peaks: usize,
    total_peaks: usize,
) -> f64 {
    (1 + decoy_peaks + estimate_decoy_peptide_errors(decoy_peptides, double_decoys)) as f64
        / total_peaks as f64
}

/// Monotone q-value correction (C# `CorrectQValues`, `:1510`): walking from the bottom of the list up,
/// each corrected q-value is the minimum of its raw value and every q-value below it (so q-values only
/// increase or stay flat as score improves). An empty input yields an empty result.
fn correct_q_values(temp_qs: &[f64]) -> Vec<f64> {
    let n = temp_qs.len();
    let mut corrected = vec![0.0f64; n];
    if n == 0 {
        return corrected;
    }
    corrected[n - 1] = temp_qs[n - 1];
    for i in (0..n - 1).rev() {
        corrected[i] = if temp_qs[i] > corrected[i + 1] {
            corrected[i + 1]
        } else {
            temp_qs[i]
        };
    }
    corrected
}

/// Orders two MBR peaks by `MbrPep` ascending, then `MbrScore` descending — the C#
/// `OrderBy(MbrPep).ThenByDescending(MbrScore)` used in the `usePep` branch. A `None` pep (peak not
/// scored by the PEP model) sorts first, matching C#'s null-is-smallest `OrderBy`.
fn pep_then_score_cmp(a: &MbrChromatographicPeak, b: &MbrChromatographicPeak) -> std::cmp::Ordering {
    let pa = a.mbr_pep.unwrap_or(f64::NEG_INFINITY);
    let pb = b.mbr_pep.unwrap_or(f64::NEG_INFINITY);
    pa.partial_cmp(&pb)
        .unwrap_or(std::cmp::Ordering::Equal)
        .then_with(|| {
            b.mbr_score
                .partial_cmp(&a.mbr_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
}

/// Assigns MBR q-values to one acceptor file's transferred peaks. Port of
/// `FlashLfqEngine.CalculateFdrForMbrPeaks` (`FlashLfqEngine.cs:1426`).
///
/// - When `use_pep` is `true` (the PEP model ran successfully), the peaks are first **deduplicated to
///   the single best acceptor per donor** — grouped by the donor identification (the donor modified
///   sequence here, the available proxy for the C# `Identifications.First()` reference; in this
///   single-donor-per-peptide port that *is* the donor identification), keeping the
///   `OrderBy(MbrPep).ThenByDescending(MbrScore)` first — then the surviving peaks are sorted the same
///   way. The dropped peaks are removed from `peaks` (C#: `_results.Peaks[acceptorFile]` keeps only the
///   filtered MBR peaks plus the non-MBR peaks). Target and decoy (`random_rt`) hypotheses of the same
///   donor compete in this grouping, exactly as the C# comment "acceptor can be target or decoy!" notes.
/// - When `use_pep` is `false`, the peaks are merely sorted by `MbrScore` descending and none are
///   dropped (C#: "better to err on the safe side and not remove the decoys").
///
/// Then it walks the ordered list, accumulating `(decoy_peptide, random_rt)` counts, and assigns each
/// peak the corrected (monotone) q-value from [`estimate_fdr`] + [`correct_q_values`]. On return `peaks`
/// is the ordered (and, under `use_pep`, filtered) list with every `mbr_q_value` set.
pub fn calculate_fdr_for_mbr_peaks(peaks: &mut Vec<MbrChromatographicPeak>, use_pep: bool) {
    if peaks.is_empty() {
        return;
    }

    let mut ordered: Vec<MbrChromatographicPeak> = if use_pep {
        let taken = std::mem::take(peaks);
        // GroupBy(Identifications.First()) — keyed by the donor modified sequence (one donor per
        // peptide in this port). First-appearance group order is irrelevant: the surviving peaks are
        // re-sorted below, so a HashMap is fine.
        let mut groups: HashMap<String, Vec<MbrChromatographicPeak>> = HashMap::new();
        let mut group_order: Vec<String> = Vec::new();
        for p in taken {
            let key = p
                .peak
                .first_modified_sequence()
                .unwrap_or("")
                .to_string();
            if !groups.contains_key(&key) {
                group_order.push(key.clone());
            }
            groups.entry(key).or_default().push(p);
        }
        // Per group: OrderBy(MbrPep).ThenByDescending(MbrScore).First(). `min_by` returns the first of
        // equal minima, and the group preserves insertion order, so this matches the C# stable sort.
        let mut selected: Vec<MbrChromatographicPeak> = Vec::with_capacity(group_order.len());
        for key in &group_order {
            let group = groups.remove(key).expect("group key enumerated above");
            let best = group
                .into_iter()
                .min_by(pep_then_score_cmp)
                .expect("non-empty group");
            selected.push(best);
        }
        // Final OrderBy(MbrPep).ThenByDescending(MbrScore) over the survivors.
        selected.sort_by(pep_then_score_cmp);
        selected
    } else {
        let mut taken = std::mem::take(peaks);
        // OrderByDescending(MbrScore) — stable, like C#.
        taken.sort_by(|a, b| {
            b.mbr_score
                .partial_cmp(&a.mbr_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        taken
    };

    // Walk the ordered list accumulating target/decoy-peptide/decoy-peak/double-decoy counts and the
    // running, 6-decimal-rounded FDR estimate (C#'s `tempQs`).
    let mut temp_qs: Vec<f64> = Vec::with_capacity(ordered.len());
    let mut total_peaks = 0usize;
    let mut decoy_peptides = 0usize;
    let mut decoy_peaks = 0usize;
    let mut double_decoys = 0usize;
    for p in &ordered {
        total_peaks += 1;
        match (p.peak.decoy_peptide(), p.random_rt) {
            (false, false) => {}
            (true, false) => decoy_peptides += 1,
            (false, true) => decoy_peaks += 1,
            (true, true) => double_decoys += 1,
        }
        temp_qs.push(round6(estimate_fdr(
            double_decoys,
            decoy_peptides,
            decoy_peaks,
            total_peaks,
        )));
    }

    let corrected = correct_q_values(&temp_qs);
    for (p, q) in ordered.iter_mut().zip(corrected.iter()) {
        p.mbr_q_value = *q;
    }

    *peaks = ordered;
}

/// Whether the MBR PEP model would run (C# `RunPEPAnalysis`, `:1393`): more than 100 transferred peaks
/// **and** more than 20 decoy (`random_rt`) peaks across all acceptor files. The result decides the
/// `use_pep` flag passed to [`apply_mbr_fdr`] (C#: `pepSuccesful` from `RunPEPAnalysis`).
pub fn mbr_pep_analysis_succeeded(result: &MbrResult) -> bool {
    let mut total = 0usize;
    let mut decoys = 0usize;
    for peaks in result.mbr_peaks_by_file.values() {
        for p in peaks {
            total += 1;
            if p.random_rt {
                decoys += 1;
            }
        }
    }
    total > 100 && decoys > 20
}

/// Assigns MBR q-values across every acceptor file of `result` (C# loops
/// `CalculateFdrForMbrPeaks(spectraFile, pepSuccesful)` over each acceptor — `FlashLfqEngine.cs:291`).
/// Files are processed in sorted name order for determinism. After this call every surviving MBR peak
/// in [`MbrResult::mbr_peaks_by_file`] carries an [`MbrChromatographicPeak::mbr_q_value`].
///
/// Note: under `use_pep` the per-file peak lists may **shrink** (only the best acceptor per donor
/// survives). [`MbrResult::feature_rows`] is the *pre-FDR* table (already consumed by the PEP model) and
/// is left unchanged; q-values live on the peaks. `use_pep` should come from [`mbr_pep_analysis_succeeded`]
/// (or whether a PEP model was actually applied via [`apply_mbr_pep`]).
pub fn apply_mbr_fdr(result: &mut MbrResult, use_pep: bool) {
    let mut files: Vec<String> = result.mbr_peaks_by_file.keys().cloned().collect();
    files.sort();
    for file in files {
        if let Some(peaks) = result.mbr_peaks_by_file.get_mut(&file) {
            calculate_fdr_for_mbr_peaks(peaks, use_pep);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detection_type::DetectionType;

    fn id(modseq: &str, base: &str, charge: i32, score: f64, q_value: f64, rt: f64) -> Identification {
        Identification {
            file_name: "f".to_string(),
            base_sequence: base.to_string(),
            modified_sequence: modseq.to_string(),
            monoisotopic_mass: 1000.0,
            ms2_retention_time_in_minutes: rt,
            precursor_charge_state: charge,
            score,
            q_value,
            is_decoy: false,
        }
    }

    fn env(scan: i32, rt: f64, intensity: f64) -> IsotopicEnvelope {
        let peak = IndexedMassSpectralPeak::new(1000.0, intensity, scan, rt);
        IsotopicEnvelope::new(peak, 2, intensity, 0.95)
    }

    /// A traced MSMS peak for one identification (single envelope fixing its apex), with the given
    /// peakfinding mass.
    fn peak(
        modseq: &str,
        base: &str,
        score: f64,
        q_value: f64,
        apex_rt: f64,
        intensity: f64,
        pfm: f64,
    ) -> ChromatographicPeak {
        let mut p = ChromatographicPeak::from_identification(
            id(modseq, base, 2, score, q_value, apex_rt),
            pfm,
        );
        p.isotopic_envelopes = vec![env(0, apex_rt, intensity)];
        p.calculate_intensity_for_this_feature(false);
        p
    }

    #[test]
    fn find_peptide_donor_files_keys_by_source_file_and_picks_best_score() {
        let mut peaks_by_file: HashMap<String, Vec<ChromatographicPeak>> = HashMap::new();
        // Sequence A appears in two files; the higher-PSM-score peak (file_b) should win.
        peaks_by_file.insert(
            "file_a".to_string(),
            vec![peak("A", "A", 5.0, 0.001, 30.0, 100.0, 1000.0)],
        );
        peaks_by_file.insert(
            "file_b".to_string(),
            vec![
                peak("A", "A", 20.0, 0.001, 31.0, 50.0, 1000.0),
                peak("B", "B", 10.0, 0.001, 40.0, 200.0, 1100.0),
            ],
        );
        let quantify: HashSet<String> = ["A".to_string(), "B".to_string()].into_iter().collect();

        let donor = find_peptide_donor_files(&peaks_by_file, &quantify);
        // A's best peak is in file_b (score 20 > 5); B is only in file_b.
        let file_b = donor.get("file_b").expect("file_b has donor peaks");
        assert_eq!(file_b.len(), 2);
        // file_a contributed nothing (its A peak lost to file_b's).
        assert!(donor.get("file_a").is_none());
    }

    #[test]
    fn find_peptide_donor_files_excludes_high_qvalue_and_unquantified() {
        let mut peaks_by_file: HashMap<String, Vec<ChromatographicPeak>> = HashMap::new();
        peaks_by_file.insert(
            "f".to_string(),
            vec![
                peak("A", "A", 5.0, 0.5, 30.0, 100.0, 1000.0), // q too high
                peak("B", "B", 5.0, 0.001, 40.0, 100.0, 1100.0), // not in quantify set
            ],
        );
        let quantify: HashSet<String> = ["A".to_string()].into_iter().collect();
        let donor = find_peptide_donor_files(&peaks_by_file, &quantify);
        assert!(donor.is_empty());
    }

    #[test]
    fn get_random_peak_picks_distant_mass_and_sequence() {
        // donor pfm ~1000; H mass ~1.0078, so 5..11 H is ~5.04..11.09 away.
        let h = periodic_table()
            .element_by_symbol("H")
            .unwrap()
            .principal_isotope()
            .atomic_mass;
        let donor_pfm = 1000.0;
        // Candidate at +8 H (within window), different base sequence, far RT.
        let cand = peak("Z", "ZZ", 10.0, 0.001, 80.0, 100.0, donor_pfm + 8.0 * h);
        // Too-close-in-mass peak (within min_diff), should be excluded.
        let near = peak("Y", "YY", 10.0, 0.001, 80.0, 100.0, donor_pfm + 1.0 * h);
        let ordered = vec![near, cand];
        let donor_id = id("D", "DD", 2, 10.0, 0.001, 30.0);

        let chosen = get_random_peak(&ordered, 30.0, 1.0, &donor_id, donor_pfm)
            .expect("a distant-mass decoy candidate exists");
        assert_eq!(chosen.first_modified_sequence(), Some("Z"));
    }

    #[test]
    fn get_random_peak_returns_none_when_no_candidate() {
        let donor_pfm = 1000.0;
        // Only candidate shares the donor base sequence -> excluded.
        let same_seq = peak("D", "DD", 10.0, 0.001, 80.0, 100.0, donor_pfm + 8.0 * 1.00782503223);
        let ordered = vec![same_seq];
        let donor_id = id("D", "DD", 2, 10.0, 0.001, 30.0);
        assert!(get_random_peak(&ordered, 30.0, 1.0, &donor_id, donor_pfm).is_none());
    }

    #[test]
    fn add_peak_to_dict_skips_apexless_peaks() {
        let mut dict: MbrPeakDict = HashMap::new();
        // A peak with no envelopes has no apex; it must be dropped, not panic.
        let apexless = MbrChromatographicPeak::new(
            id("A", "A", 2, 10.0, 0.001, 30.0),
            "acc".to_string(),
            1000.0,
            30.0,
            false,
        );
        add_peak_to_dict(&mut dict, "A", &Some(apexless));
        assert!(dict.is_empty());

        // A peak with an apex is keyed.
        let mut withapex = MbrChromatographicPeak::new(
            id("A", "A", 2, 10.0, 0.001, 30.0),
            "acc".to_string(),
            1000.0,
            30.0,
            false,
        );
        withapex.peak.isotopic_envelopes = vec![env(0, 30.0, 100.0)];
        withapex.peak.calculate_intensity_for_this_feature(false);
        add_peak_to_dict(&mut dict, "A", &Some(withapex));
        assert_eq!(dict.len(), 1);
        assert_eq!(dict["A"].len(), 1);
    }

    #[test]
    fn detection_type_default_is_mbr_for_constructed_peaks() {
        let p = MbrChromatographicPeak::new(
            id("A", "A", 2, 10.0, 0.001, 30.0),
            "acc".to_string(),
            1000.0,
            30.0,
            false,
        );
        assert_eq!(p.peak.detection_type, DetectionType::MBR);
    }

    fn mbr_peak(modseq: &str, random_rt: bool, predicted_rt: f64) -> MbrChromatographicPeak {
        MbrChromatographicPeak::new(
            id(modseq, modseq, 2, 10.0, 0.001, 30.0),
            "acc".to_string(),
            1000.0,
            predicted_rt,
            random_rt,
        )
    }

    fn frow(modseq: &str, random_rt: bool, predicted_rt: f64) -> FeatureRow {
        FeatureRow {
            donor_modified_sequence: modseq.to_string(),
            donor_base_sequence: modseq.to_string(),
            acceptor_file: "acc".to_string(),
            predicted_retention_time: predicted_rt,
            apex_retention_time: -1.0,
            intensity: 0.0,
            ppm_score: 0.5,
            intensity_score: 0.5,
            rt_score: 0.5,
            scan_count_score: 0.5,
            isotopic_distribution_score: 0.5,
            mbr_score: 50.0,
            mass_error: 0.0,
            scan_count: 0,
            isotopic_pearson_correlation: -1.0,
            rt_prediction_error: 0.0,
            random_rt,
            decoy_peptide: false,
        }
    }

    #[test]
    fn apply_mbr_pep_maps_table_order_scores_back_onto_peaks() {
        // Peaks (and their feature rows) in natural emission order: B/target, A/target, A/decoy.
        let peaks = vec![
            mbr_peak("B", false, 10.0),
            mbr_peak("A", false, 20.0),
            mbr_peak("A", true, 20.0),
        ];
        let feature_rows = vec![
            frow("B", false, 10.0),
            frow("A", false, 20.0),
            frow("A", true, 20.0),
        ];
        let mut by_file: HashMap<String, Vec<MbrChromatographicPeak>> = HashMap::new();
        by_file.insert("acc".to_string(), peaks);
        let mut result = MbrResult {
            feature_rows,
            mbr_peaks_by_file: by_file,
        };

        // Canonical table order is (A/target, A/decoy, B/target) = natural indices [1, 2, 0].
        // Hand back a distinct pep per sorted row; it must land on the right natural peak.
        apply_mbr_pep(&mut result, &[0.1, 0.9, 0.5]).expect("lengths match");
        let peaks = &result.mbr_peaks_by_file["acc"];
        assert_eq!(peaks[0].mbr_pep, Some(0.5)); // B/target, sorted position 2
        assert_eq!(peaks[1].mbr_pep, Some(0.1)); // A/target, sorted position 0
        assert_eq!(peaks[2].mbr_pep, Some(0.9)); // A/decoy, sorted position 1
    }

    #[test]
    fn apply_mbr_pep_rejects_length_mismatch() {
        let mut by_file: HashMap<String, Vec<MbrChromatographicPeak>> = HashMap::new();
        by_file.insert("acc".to_string(), vec![mbr_peak("A", false, 20.0)]);
        let mut result = MbrResult {
            feature_rows: vec![frow("A", false, 20.0)],
            mbr_peaks_by_file: by_file,
        };
        let err = apply_mbr_pep(&mut result, &[0.1, 0.2]).unwrap_err();
        assert_eq!(err.peps, 2);
        assert_eq!(err.rows, 1);
    }

    // ----- P3.4: MBR FDR / q-values -----------------------------------------------------------

    /// Builds a scored MBR peak with controllable combined score, pep, decoy-peak (`random_rt`) and
    /// decoy-peptide (`is_decoy` on the id) flags — the four inputs `CalculateFdrForMbrPeaks` keys on.
    fn mbr_scored(
        modseq: &str,
        mbr_score: f64,
        mbr_pep: Option<f64>,
        random_rt: bool,
        decoy_peptide: bool,
    ) -> MbrChromatographicPeak {
        let mut identification = id(modseq, modseq, 2, 10.0, 0.001, 30.0);
        identification.is_decoy = decoy_peptide;
        let mut p = MbrChromatographicPeak::new(identification, "acc".to_string(), 1000.0, 30.0, random_rt);
        p.mbr_score = mbr_score;
        p.mbr_pep = mbr_pep;
        p
    }

    #[test]
    fn estimate_fdr_matches_csharp_formula() {
        // (1 + decoyPeaks + max(0, decoyPeptides - doubleDecoys)) / totalPeaks
        assert_eq!(estimate_fdr(0, 0, 0, 1), 1.0);
        assert_eq!(estimate_fdr(0, 0, 0, 2), 0.5);
        assert_eq!(estimate_fdr(0, 0, 1, 4), 0.5); // one decoy peak
        assert_eq!(estimate_fdr(0, 3, 0, 5), (1.0 + 3.0) / 5.0); // 3 decoy peptides
        // double-decoy correction: max(0, 2 - 5) = 0, so only the +1 numerator term remains.
        assert_eq!(estimate_fdr(5, 2, 0, 10), 1.0 / 10.0);
    }

    #[test]
    fn correct_q_values_is_monotone_nondecreasing() {
        // Raw q's dip and rise; the corrected list is the running min-from-below, hence non-decreasing.
        let corrected = correct_q_values(&[1.0, 0.5, 0.333333, 0.5, 0.6]);
        assert_eq!(corrected, vec![0.333333, 0.333333, 0.333333, 0.5, 0.6]);
        for w in corrected.windows(2) {
            assert!(w[0] <= w[1], "q-values must be non-decreasing down the list");
        }
        assert!(correct_q_values(&[]).is_empty());
    }

    #[test]
    fn calculate_fdr_no_pep_orders_by_score_and_assigns_qvalues() {
        // Out-of-order on input; no peaks are dropped without PEP.
        let mut peaks = vec![
            mbr_scored("P4", 60.0, None, true, false),  // decoy peak
            mbr_scored("P1", 90.0, None, false, false), // target
            mbr_scored("P5", 50.0, None, true, false),  // decoy peak
            mbr_scored("P3", 70.0, None, false, false), // target
            mbr_scored("P2", 80.0, None, false, false), // target
        ];
        calculate_fdr_for_mbr_peaks(&mut peaks, false);

        // Sorted by MbrScore descending, nothing dropped.
        let scores: Vec<f64> = peaks.iter().map(|p| p.mbr_score).collect();
        assert_eq!(scores, vec![90.0, 80.0, 70.0, 60.0, 50.0]);

        // q-values match the hand-computed corrected FDR walk and are non-decreasing.
        let qs: Vec<f64> = peaks.iter().map(|p| p.mbr_q_value).collect();
        let approx = |a: f64, b: f64| (a - b).abs() < 1e-9;
        // q-values are stored rounded to 6 decimals (C# Math.Round(_, 6)): 1/3 -> 0.333333.
        assert!(approx(qs[0], 0.333333));
        assert!(approx(qs[1], 0.333333));
        assert!(approx(qs[2], 0.333333));
        assert!(approx(qs[3], 0.5));
        assert!(approx(qs[4], 0.6));
        for w in qs.windows(2) {
            assert!(w[0] <= w[1] + 1e-12);
        }
    }

    #[test]
    fn calculate_fdr_use_pep_dedups_best_acceptor_per_donor() {
        // Donor A has a target (pep 0.2) and a decoy (pep 0.7) hypothesis; donor B a single target.
        let mut peaks = vec![
            mbr_scored("A", 80.0, Some(0.2), false, false),
            mbr_scored("A", 90.0, Some(0.7), true, false), // higher score but worse pep -> loses
            mbr_scored("B", 50.0, Some(0.1), false, false),
        ];
        calculate_fdr_for_mbr_peaks(&mut peaks, true);

        // One survivor per donor; the A decoy is dropped. Survivors sorted by MbrPep ascending.
        assert_eq!(peaks.len(), 2);
        assert_eq!(peaks[0].peak.first_modified_sequence(), Some("B")); // pep 0.1
        assert_eq!(peaks[1].peak.first_modified_sequence(), Some("A")); // pep 0.2 (target, not decoy)
        assert!(!peaks[1].random_rt, "the surviving A peak is the target, not the decoy");

        // Both survivors carry q-values (both targets here -> 0.5 after correction).
        for p in &peaks {
            assert!(p.mbr_q_value > 0.0);
        }
    }

    #[test]
    fn calculate_fdr_handles_empty_and_double_decoys() {
        let mut empty: Vec<MbrChromatographicPeak> = Vec::new();
        calculate_fdr_for_mbr_peaks(&mut empty, false); // must not panic
        assert!(empty.is_empty());

        // A double decoy (decoy peptide AND random RT) is counted in the double-decoy bucket, which the
        // EstimateDecoyPeptideErrors correction subtracts back out — so a lone double decoy after a
        // target does not inflate the FDR the way a plain decoy peak would.
        let mut peaks = vec![
            mbr_scored("T", 90.0, None, false, false), // target
            mbr_scored("D", 80.0, None, true, true),   // double decoy
        ];
        calculate_fdr_for_mbr_peaks(&mut peaks, false);
        // total=2: decoyPeaks=0, decoyPeptides=0, doubleDecoys=1 -> (1+0+0)/2 = 0.5 at the second peak.
        assert!((peaks[1].mbr_q_value - 0.5).abs() < 1e-9);
    }

    #[test]
    fn mbr_pep_analysis_succeeded_applies_count_thresholds() {
        let make = |total: usize, decoys: usize| {
            let mut v = Vec::new();
            for i in 0..total {
                v.push(mbr_scored("X", 50.0, None, i < decoys, false));
            }
            let mut by_file = HashMap::new();
            by_file.insert("acc".to_string(), v);
            MbrResult {
                feature_rows: Vec::new(),
                mbr_peaks_by_file: by_file,
            }
        };
        // Both thresholds are strict ">": exactly 100 peaks / 20 decoys is not enough.
        assert!(!mbr_pep_analysis_succeeded(&make(100, 21)));
        assert!(!mbr_pep_analysis_succeeded(&make(101, 20)));
        assert!(mbr_pep_analysis_succeeded(&make(101, 21)));
    }
}
