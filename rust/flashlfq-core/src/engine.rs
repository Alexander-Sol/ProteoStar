//! Engine driver — port of the MS2 quantification path of `FlashLfqEngine.Run`
//! (PLAN.md P1.17b). Wires together the already-ported pieces (theoretical isotope
//! distributions, the binned m/z index + XIC extraction, envelope search, peak cutting,
//! charge aggregation, error checking, and peptide-result aggregation) into the end-to-end
//! flow that produces the L5 (per-`ChromatographicPeak`) and L6 (peptide × file) tables.
//!
//! ## Scope: the default MS2-only path
//!
//! This mirrors `FlashLfqEngine.Run` run with `FlashLfqParameters` defaults — the MS2 path
//! only. Concretely: `Normalize=false`, `PpmTolerance=10`, `IsotopePpmTolerance=5`,
//! `Integrate=false`, `NumIsotopesRequired=2`, `IdSpecificChargeState=false`,
//! `QuantifyAmbiguousPeptides=false`, `MatchBetweenRuns=false`, `IsoTracker=false`. The
//! MBR / IsoTracker / normalization / Bayesian-protein branches of `Run` are intentionally
//! absent (Phase 2+).
//!
//! ## Flow (faithful to `FlashLfqEngine.cs`)
//!
//! 1. **`CalculateTheoreticalIsotopeDistributions`** (`FlashLfqEngine.cs:351`): one theoretical
//!    envelope per distinct modified sequence (`expected_isotope_peaks(None, base, mono, 2)` —
//!    `OptionalChemicalFormula` is never set by `MakeIdentifications`, so it is always `None`),
//!    the engine-global charge-state range `[min..=max]` over **all** ids, and a per-id
//!    peakfinding mass `mono + mostAbundantIsotopeShift`.
//! 2. Per file, **`QuantifyMs2IdentifiedPeptides`** (`FlashLfqEngine.cs:472`): build a
//!    `ChromatographicPeak` per id, trace each charge's XIC, mass-filter it, search isotopic
//!    envelopes, cut the peak, then trim to the precursor-charge scan range.
//! 3. Per file, **`RunErrorChecking`** = [`run_error_checking`].
//! 4. **`CalculatePeptideResults`** = [`calculate_peptide_results`] over all files.

use std::collections::HashMap;
use std::fmt;
use std::path::{Path, PathBuf};

use crate::chromatographic_peak::ChromatographicPeak;
use crate::isotopic_envelope::{
    get_isotopic_envelopes, mass_to_mz_f64, mz_to_mass_f32, IsotopicEnvelope,
};
use crate::parquet_output::write_peptide_results_parquet;
use crate::peak_indexing::PeakIndexingEngine;
use crate::psm_tsv::{file_name_without_extension, read_identifications, Identification, PsmTsvError};
use crate::results::{
    calculate_peptide_results, engine_quantify_set, output_peptide_sequences, run_error_checking,
    PeptideResults,
};
use crate::theoretical_isotope_distribution::{
    expected_isotope_peaks, most_abundant_isotope_shift, ExpectedIsotopePeak,
};
use crate::tolerance::PpmTolerance;

/// Peak-finding ppm tolerance (`FlashLfqEngine.PeakfindingPpmTolerance`).
pub const PEAKFINDING_PPM_TOLERANCE: f64 = 20.0;
/// The narrower mass tolerance (`FlashLfqParameters.PpmTolerance` default).
pub const PPM_TOLERANCE: f64 = 10.0;
/// Isotope-matching ppm tolerance (`FlashLfqParameters.IsotopePpmTolerance` default).
pub const ISOTOPE_PPM_TOLERANCE: f64 = 5.0;
/// Minimum number of accurate isotope peaks required (`FlashLfqParameters.NumIsotopesRequired`).
pub const NUM_ISOTOPES_REQUIRED: usize = 2;
/// Consecutive missed scans allowed during XIC tracing (`FlashLfqEngine.MissedScansAllowed`).
pub const MISSED_SCANS_ALLOWED: i32 = 1;
/// Whether to integrate (sum) envelope intensities; default `false` → apex intensity.
pub const INTEGRATE: bool = false;
/// Whether to restrict XIC tracing to the id's own charge state
/// (`FlashLfqParameters.IdSpecificChargeState`, default `false`).
pub const ID_SPECIFIC_CHARGE_STATE: bool = false;

/// `maxPeakHalfWidth` default for `GetXic` (`int.MaxValue` in the C# default argument).
const MAX_PEAK_HALF_WIDTH: f64 = i32::MAX as f64;

/// The output of the MS2 quantification path: the error-checked peaks per file (the L5 table)
/// and the aggregated peptide × file results (the L6 table).
pub struct Ms2QuantResult {
    /// `file_name -> error-checked chromatographic peaks` (post-`RunErrorChecking`).
    pub peaks_by_file: HashMap<String, Vec<ChromatographicPeak>>,
    /// Peptide-level results (the peptide × file intensity/RT/detection table).
    pub peptide_results: PeptideResults,
    /// `file_name -> built peak index` (C# `IndexingEngineDictionary`). Kept so the MBR acceptor
    /// search (P3.2d) can look up indexed peaks / trace XICs in each acceptor file without
    /// re-reading the spectra. One entry per file in `peaks_by_file`.
    pub engines_by_file: HashMap<String, PeakIndexingEngine>,
    /// Non-decoy modified sequences in the quantify set (C#
    /// `PeptideModifiedSequencesToQuantify`) — the donor-selection / acceptor-search gate used by
    /// the MBR driver.
    pub peptide_sequences_to_quantify: std::collections::HashSet<String>,
}

/// Theoretical-distribution state computed once for the whole id set, mirroring the products of
/// `CalculateTheoreticalIsotopeDistributions`.
struct TheoreticalState {
    /// `modified_sequence -> theoretical isotope peaks`.
    mod_seq_to_dist: HashMap<String, Vec<ExpectedIsotopePeak>>,
    /// `modified_sequence -> most-abundant isotope mass shift` (the peakfinding-mass offset).
    mod_seq_to_shift: HashMap<String, f64>,
    /// Engine-global charge-state range `[min..=max]` over all ids.
    charge_states: Vec<i32>,
}

/// Ports `CalculateTheoreticalIsotopeDistributions`: builds the per-modified-sequence theoretical
/// envelope cache, the per-modified-sequence peakfinding shift, and the engine-global charge range.
fn calculate_theoretical_isotope_distributions(ids: &[Identification]) -> TheoreticalState {
    let mut mod_seq_to_dist: HashMap<String, Vec<ExpectedIsotopePeak>> = HashMap::new();

    // First occurrence of each modified sequence builds the distribution (C# `ContainsKey` skip).
    for id in ids {
        if mod_seq_to_dist.contains_key(&id.modified_sequence) {
            continue;
        }
        let dist = expected_isotope_peaks(
            None,
            &id.base_sequence,
            id.monoisotopic_mass,
            NUM_ISOTOPES_REQUIRED,
        );
        mod_seq_to_dist.insert(id.modified_sequence.clone(), dist);
    }

    let mut mod_seq_to_shift: HashMap<String, f64> = HashMap::new();
    for (seq, dist) in &mod_seq_to_dist {
        let shift = most_abundant_isotope_shift(dist)
            .expect("expected_isotope_peaks always normalizes one peak to abundance 1.0");
        mod_seq_to_shift.insert(seq.clone(), shift);
    }

    let min_charge = ids
        .iter()
        .map(|id| id.precursor_charge_state)
        .min()
        .expect("identification set is non-empty");
    let max_charge = ids
        .iter()
        .map(|id| id.precursor_charge_state)
        .max()
        .expect("identification set is non-empty");
    let charge_states: Vec<i32> = (min_charge..=max_charge).collect();

    TheoreticalState {
        mod_seq_to_dist,
        mod_seq_to_shift,
        charge_states,
    }
}

/// Ports `QuantifyMs2IdentifiedPeptides` for one file: produces one `ChromatographicPeak` per id
/// (in id order), tracing every charge state's XIC, mass-filtering it, searching isotopic
/// envelopes, cutting the peak, and trimming to the precursor-charge scan range. Empty peaks are
/// retained (matching C#: the array slot is filled before the `continue` guards).
fn quantify_ms2_identified_peptides(
    engine: &PeakIndexingEngine,
    ids_for_file: &[&Identification],
    state: &TheoreticalState,
) -> Vec<ChromatographicPeak> {
    let peakfinding_tol = PpmTolerance::new(PEAKFINDING_PPM_TOLERANCE);
    let ppm_tolerance = PpmTolerance::new(PPM_TOLERANCE);

    let mut peaks: Vec<ChromatographicPeak> = Vec::with_capacity(ids_for_file.len());

    for id in ids_for_file {
        let shift = state.mod_seq_to_shift[&id.modified_sequence];
        let peakfinding_mass = id.monoisotopic_mass + shift;
        let dist = &state.mod_seq_to_dist[&id.modified_sequence];

        let mut peak = ChromatographicPeak::from_identification((*id).clone(), peakfinding_mass);

        for &charge in &state.charge_states {
            if ID_SPECIFIC_CHARGE_STATE && charge != id.precursor_charge_state {
                continue;
            }

            // get XIC (peak finding); get_xic already returns RT-ascending order.
            let mut xic = engine.get_xic(
                mass_to_mz_f64(peakfinding_mass, charge),
                id.ms2_retention_time_in_minutes,
                &peakfinding_tol,
                MISSED_SCANS_ALLOWED,
                MAX_PEAK_HALF_WIDTH,
                None,
            );

            // filter by the smaller mass tolerance
            xic.retain(|p| ppm_tolerance.within(mz_to_mass_f32(p.m(), charge), peakfinding_mass));

            // filter by isotopic distribution
            let envelopes = get_isotopic_envelopes(
                engine,
                &xic,
                dist,
                id.monoisotopic_mass,
                peakfinding_mass,
                charge,
                NUM_ISOTOPES_REQUIRED,
                ISOTOPE_PPM_TOLERANCE,
            );

            peak.isotopic_envelopes.extend(envelopes);
        }

        peak.calculate_intensity_for_this_feature(INTEGRATE);
        peak.cut_peak(id.ms2_retention_time_in_minutes, INTEGRATE);

        if !peak.isotopic_envelopes.is_empty() {
            let precursor_charge = id.precursor_charge_state;
            let precursor_xic: Vec<IsotopicEnvelope> = peak
                .isotopic_envelopes
                .iter()
                .filter(|e| e.charge_state == precursor_charge)
                .copied()
                .collect();

            if precursor_xic.is_empty() {
                peak.isotopic_envelopes.clear();
            } else {
                let min = precursor_xic
                    .iter()
                    .map(|e| e.indexed_peak.zero_based_scan_index)
                    .min()
                    .expect("precursor_xic is non-empty");
                let max = precursor_xic
                    .iter()
                    .map(|e| e.indexed_peak.zero_based_scan_index)
                    .max()
                    .expect("precursor_xic is non-empty");
                peak.isotopic_envelopes.retain(|e| {
                    let s = e.indexed_peak.zero_based_scan_index;
                    s >= min && s <= max
                });
                peak.calculate_intensity_for_this_feature(INTEGRATE);
            }
        }

        peaks.push(peak);
    }

    peaks
}

/// Runs the MS2 quantification path end-to-end (PLAN.md P1.17b).
///
/// `ids` are the parsed identifications (psmtsv order preserved); `file_to_mzml` maps each id's
/// `file_name` to the mzML path to index. Returns the error-checked peaks per file and the
/// aggregated peptide × file results.
///
/// Files are processed in first-appearance order of `ids` (the order is immaterial to the L5/L6
/// outputs, which are sorted before comparison, but keeps per-file peak order deterministic).
pub fn run_msms(
    ids: Vec<Identification>,
    file_to_mzml: &HashMap<String, PathBuf>,
) -> std::io::Result<Ms2QuantResult> {
    let state = calculate_theoretical_isotope_distributions(&ids);
    // Two faithful sets (P2.2): the quantify set (membership gate, decoys included — the C#
    // `FlashLfqEngine` default) and the output peptide set (non-decoy rows — the C#
    // `FlashLfqResults.PeptideModifiedSequences` keys).
    let quantify_set = engine_quantify_set(&ids);
    let peptide_sequences = output_peptide_sequences(&ids, &quantify_set);

    // Unique file names in first-appearance order.
    let mut file_order: Vec<String> = Vec::new();
    for id in &ids {
        if !file_order.iter().any(|f| f == &id.file_name) {
            file_order.push(id.file_name.clone());
        }
    }

    let mut peaks_by_file: HashMap<String, Vec<ChromatographicPeak>> = HashMap::new();
    let mut engines_by_file: HashMap<String, PeakIndexingEngine> = HashMap::new();

    for file_name in &file_order {
        let mzml_path = file_to_mzml
            .get(file_name)
            .unwrap_or_else(|| panic!("no mzML path provided for file '{file_name}'"));

        let engine = PeakIndexingEngine::from_spectra_file(mzml_path.clone())?
            .unwrap_or_else(|| panic!("file '{file_name}' contained no MS1 peaks"));

        let ids_for_file: Vec<&Identification> =
            ids.iter().filter(|id| &id.file_name == file_name).collect();

        let file_peaks = quantify_ms2_identified_peptides(&engine, &ids_for_file, &state);
        let error_checked = run_error_checking(file_peaks, &quantify_set, INTEGRATE);
        peaks_by_file.insert(file_name.clone(), error_checked);
        // Keep the index alive for the MBR acceptor search (P3.2d).
        engines_by_file.insert(file_name.clone(), engine);
    }

    let peptide_results =
        calculate_peptide_results(&peaks_by_file, &quantify_set, &peptide_sequences);

    Ok(Ms2QuantResult {
        peaks_by_file,
        peptide_results,
        engines_by_file,
        peptide_sequences_to_quantify: peptide_sequences,
    })
}

/// Errors raised by the high-level [`quant`] / [`quant_to_parquet`] entry points.
#[derive(Debug)]
pub enum QuantError {
    /// The psmtsv could not be parsed.
    PsmTsv(PsmTsvError),
    /// An mzML file could not be opened/read.
    Io(std::io::Error),
    /// An identification referenced a spectra file with no matching path in `raw_paths`.
    /// Carries the bare file name (psmtsv `File Name`, extension stripped) and the bare names
    /// that *were* supplied, so the caller can see the mismatch.
    MissingSpectraFile { file_name: String, supplied: Vec<String> },
    /// Writing the Parquet output failed.
    Parquet(String),
}

impl fmt::Display for QuantError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QuantError::PsmTsv(e) => write!(f, "failed to read identifications: {e}"),
            QuantError::Io(e) => write!(f, "mzML read error: {e}"),
            QuantError::MissingSpectraFile { file_name, supplied } => write!(
                f,
                "no spectra file supplied for identification file {file_name:?}; \
                 supplied bare names: {supplied:?}"
            ),
            QuantError::Parquet(e) => write!(f, "failed to write Parquet output: {e}"),
        }
    }
}

impl std::error::Error for QuantError {}

impl From<PsmTsvError> for QuantError {
    fn from(e: PsmTsvError) -> Self {
        QuantError::PsmTsv(e)
    }
}

impl From<std::io::Error> for QuantError {
    fn from(e: std::io::Error) -> Self {
        QuantError::Io(e)
    }
}

/// Maps each id's spectra-file name to a supplied path, keyed by the *bare* file name (the
/// psmtsv `File Name` transform — extension + directory stripped — applied to each raw path so a
/// `--raw .../foo.mzML` matches ids whose `file_name` is `foo`). Errors if any id references a
/// file no path was supplied for.
fn resolve_file_to_mzml(
    ids: &[Identification],
    raw_paths: &[PathBuf],
) -> Result<HashMap<String, PathBuf>, QuantError> {
    // bare name -> path, for the supplied raw files.
    let mut by_bare: HashMap<String, PathBuf> = HashMap::new();
    for path in raw_paths {
        let stem = path
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default();
        by_bare.insert(file_name_without_extension(&stem), path.clone());
    }

    let mut file_to_mzml: HashMap<String, PathBuf> = HashMap::new();
    for id in ids {
        if file_to_mzml.contains_key(&id.file_name) {
            continue;
        }
        match by_bare.get(&id.file_name) {
            Some(path) => {
                file_to_mzml.insert(id.file_name.clone(), path.clone());
            }
            None => {
                let mut supplied: Vec<String> = by_bare.keys().cloned().collect();
                supplied.sort();
                return Err(QuantError::MissingSpectraFile {
                    file_name: id.file_name.clone(),
                    supplied,
                });
            }
        }
    }
    Ok(file_to_mzml)
}

/// High-level entry point: parse a psmtsv, index the supplied mzML files, run the MS2
/// quantification path, and return the in-memory [`Ms2QuantResult`].
///
/// `results_path` is a MetaMorpheus `*.psmtsv`; `raw_paths` are the spectra files (mzML in
/// Phase 1). Each identification is matched to a spectra file by *bare* file name — the same
/// extension-and-directory-stripped name the psmtsv `File Name` column resolves to — so the
/// order of `raw_paths` is immaterial and extra paths are ignored. Returns [`QuantError`] when
/// parsing fails, a referenced spectra file is missing, or a file holds no MS1 peaks.
///
/// `params` is reserved for the parameter surface; Phase 1 only supports the default MS2 path,
/// so callers pass no tunables yet (see the module docs for the fixed defaults).
pub fn quant<P: AsRef<Path>>(
    results_path: P,
    raw_paths: &[PathBuf],
) -> Result<Ms2QuantResult, QuantError> {
    let ids = read_identifications(results_path)?;
    let file_to_mzml = resolve_file_to_mzml(&ids, raw_paths)?;
    // `file_to_mzml` covers every id file name, so `run_msms` will not hit its missing-path path.
    run_msms(ids, &file_to_mzml).map_err(QuantError::Io)
}

/// Like [`quant`], but writes the peptide × file results table to a Parquet file at
/// `output_path` (tidy long form, via [`write_peptide_results_parquet`]) and returns the number
/// of rows written. This is the end-to-end FlashLFQ entry the Python `quant()` binding wraps.
pub fn quant_to_parquet<P: AsRef<Path>, Q: AsRef<Path>>(
    results_path: P,
    raw_paths: &[PathBuf],
    output_path: Q,
) -> Result<usize, QuantError> {
    let result = quant(results_path, raw_paths)?;
    write_peptide_results_parquet(&result.peptide_results, output_path)
        .map_err(|e| QuantError::Parquet(e.to_string()))
}

/// Runs the MS2 quant ([`quant`]) and then the match-between-runs transfer ([`run_mbr`]),
/// returning the MBR feature table (PLAN.md P3.2e). This is the end-to-end FlashLFQ-with-MBR
/// entry the Python `quant(match_between_runs=True)` binding wraps.
///
/// `run_mbr` reuses the per-file `PeakIndexingEngine`s kept alive on the [`Ms2QuantResult`], so no
/// spectra are re-read.
///
/// [`run_mbr`]: crate::mbr_search::run_mbr
pub fn quant_mbr<P: AsRef<Path>>(
    results_path: P,
    raw_paths: &[PathBuf],
) -> Result<crate::mbr_search::MbrResult, QuantError> {
    let result = quant(results_path, raw_paths)?;
    Ok(crate::mbr_search::run_mbr(
        &result.peaks_by_file,
        &result.engines_by_file,
        &result.peptide_sequences_to_quantify,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn charge_states_span_min_to_max() {
        let make = |charge: i32| Identification {
            file_name: "f".to_string(),
            base_sequence: "PEPTIDE".to_string(),
            modified_sequence: "PEPTIDE".to_string(),
            monoisotopic_mass: 799.359964,
            ms2_retention_time_in_minutes: 1.0,
            precursor_charge_state: charge,
            score: 0.0,
            q_value: 0.0,
            is_decoy: false,
        };
        let ids = vec![make(2), make(4), make(3)];
        let state = calculate_theoretical_isotope_distributions(&ids);
        assert_eq!(state.charge_states, vec![2, 3, 4]);
        // One distribution cached for the single distinct modified sequence.
        assert_eq!(state.mod_seq_to_dist.len(), 1);
        assert!(state.mod_seq_to_shift.contains_key("PEPTIDE"));
    }

    fn id_for_file(file_name: &str) -> Identification {
        Identification {
            file_name: file_name.to_string(),
            base_sequence: "PEPTIDE".to_string(),
            modified_sequence: "PEPTIDE".to_string(),
            monoisotopic_mass: 799.359964,
            ms2_retention_time_in_minutes: 1.0,
            precursor_charge_state: 2,
            score: 0.0,
            q_value: 0.0,
            is_decoy: false,
        }
    }

    #[test]
    fn resolve_file_to_mzml_matches_ids_by_bare_name() {
        // ids reference the bare names; raw paths carry directories + extensions in any order.
        let ids = vec![
            id_for_file("K562_3"),
            id_for_file("K562_4"),
            id_for_file("K562_3"), // duplicate file name resolves once
        ];
        let raw_paths = vec![
            PathBuf::from("/data/runs/K562_4.mzML"),
            PathBuf::from(r"C:\spectra\K562_3.raw"),
            PathBuf::from("/data/runs/UNUSED_extra.mzML"), // extra path is ignored
        ];

        let map = resolve_file_to_mzml(&ids, &raw_paths).expect("all ids resolve");
        assert_eq!(map.len(), 2);
        assert_eq!(map["K562_3"], PathBuf::from(r"C:\spectra\K562_3.raw"));
        assert_eq!(map["K562_4"], PathBuf::from("/data/runs/K562_4.mzML"));
    }

    #[test]
    fn resolve_file_to_mzml_errors_on_missing_spectra_file() {
        let ids = vec![id_for_file("K562_3"), id_for_file("MISSING")];
        let raw_paths = vec![PathBuf::from("/data/K562_3.mzML")];

        match resolve_file_to_mzml(&ids, &raw_paths) {
            Err(QuantError::MissingSpectraFile { file_name, supplied }) => {
                assert_eq!(file_name, "MISSING");
                assert_eq!(supplied, vec!["K562_3".to_string()]);
            }
            other => panic!("expected MissingSpectraFile, got {other:?}"),
        }
    }
}
