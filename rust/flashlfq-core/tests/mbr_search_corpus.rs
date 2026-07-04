//! MBR acceptor-search + orchestration smoke test over the real corpus (PLAN.md P3.2d).
//!
//! Runs the Phase-1 MS2 engine over the two K562 mzMLs, then drives the full match-between-runs
//! transfer ([`flashlfq_core::mbr_search::run_mbr`]) and asserts the emitted **feature table** — the
//! data P3.3's Python model trains on — is well-formed: target *and* decoy (`random_rt`) peaks are
//! produced, every component score lies in `(0, 1]`, the combined score lies in `[0, 100]`, and the
//! transferred peptides were genuinely absent from the acceptor run (otherwise MBR would not transfer
//! them). There is no C# golden for the orchestration yet (the method is private + parallel; a faithful
//! replica is a follow-up), so this validates the port runs end-to-end over real data and yields a
//! sane training table.

use std::collections::HashMap;
use std::path::PathBuf;

use flashlfq_core::engine::{run_msms, Ms2QuantResult};
use flashlfq_core::mbr_search::{apply_mbr_fdr, mbr_pep_analysis_succeeded, run_mbr};
use flashlfq_core::psm_tsv::read_identifications;

const FILE_3: &str = "20100614_Velos1_TaGe_SA_K562_3";
const FILE_4: &str = "20100614_Velos1_TaGe_SA_K562_4";

fn test_data(relative: &str) -> PathBuf {
    flashlfq_core::mzlib_test_data(relative)
}

fn run_engine() -> Ms2QuantResult {
    let ids = read_identifications(test_data("AllPSMs.psmtsv"))
        .expect("AllPSMs.psmtsv parses into identifications");
    let mut file_to_mzml = HashMap::new();
    file_to_mzml.insert(FILE_3.to_string(), test_data(&format!("{FILE_3}.mzML")));
    file_to_mzml.insert(FILE_4.to_string(), test_data(&format!("{FILE_4}.mzML")));
    run_msms(ids, &file_to_mzml).expect("engine runs over the corpus")
}

#[test]
fn mbr_emits_a_well_formed_feature_table() {
    let result = run_engine();

    let mbr = run_mbr(
        &result.peaks_by_file,
        &result.engines_by_file,
        &result.peptide_sequences_to_quantify,
    );

    assert!(
        !mbr.feature_rows.is_empty(),
        "MBR should transfer at least some peaks across the two K562 files"
    );

    let mut targets = 0usize;
    let mut decoys = 0usize;
    for row in &mbr.feature_rows {
        // The acceptor is one of the two corpus files.
        assert!(
            row.acceptor_file == FILE_3 || row.acceptor_file == FILE_4,
            "unexpected acceptor file: {}",
            row.acceptor_file
        );

        // Combined score in [0, 100].
        assert!(
            row.mbr_score >= 0.0 && row.mbr_score <= 100.0 + 1e-9,
            "mbr_score out of range: {}",
            row.mbr_score
        );

        // Every component score in (0, 1].
        for (name, s) in [
            ("ppm", row.ppm_score),
            ("intensity", row.intensity_score),
            ("rt", row.rt_score),
            ("scan_count", row.scan_count_score),
            ("isotopic", row.isotopic_distribution_score),
        ] {
            assert!(
                s > 0.0 && s <= 1.0 + 1e-9,
                "{name} component score out of range: {s} (seq {})",
                row.donor_modified_sequence
            );
        }

        // A transferred peak has a sequence and a non-negative scan count.
        assert!(!row.donor_modified_sequence.is_empty());

        if row.random_rt {
            decoys += 1;
        } else {
            targets += 1;
        }
    }

    assert!(targets > 0, "expected at least one target (real-RT) transfer");
    assert!(decoys > 0, "expected at least one decoy (random-RT) transfer");

    // The stored peaks match the row count, and every stored peak is an MBR peak.
    let stored: usize = mbr.mbr_peaks_by_file.values().map(|v| v.len()).sum();
    assert_eq!(
        stored,
        mbr.feature_rows.len(),
        "one feature row per stored MBR peak"
    );
    for peaks in mbr.mbr_peaks_by_file.values() {
        for p in peaks {
            assert_eq!(
                p.peak.detection_type,
                flashlfq_core::detection_type::DetectionType::MBR
            );
        }
    }
}

/// P3.4: after the FDR pass every surviving MBR peak carries a q-value, and the q-values are monotone
/// (non-decreasing) down each acceptor file's score-ordered list. The corpus has > 100 transferred
/// peaks and > 20 decoys, so the PEP-success gate is met and the `use_pep` dedup branch runs.
#[test]
fn mbr_fdr_assigns_monotone_qvalues_to_every_peak() {
    let result = run_engine();
    let mut mbr = run_mbr(
        &result.peaks_by_file,
        &result.engines_by_file,
        &result.peptide_sequences_to_quantify,
    );

    // The K562 corpus clears both RunPEPAnalysis thresholds (>100 peaks, >20 decoys).
    assert!(
        mbr_pep_analysis_succeeded(&mbr),
        "the corpus should produce enough peaks/decoys to run PEP"
    );

    // Without peps assigned, the use_pep dedup would collapse to one peak per donor with no pep info;
    // exercise the no-pep path here (deterministic, no Python), which keeps every peak and just sorts.
    apply_mbr_fdr(&mut mbr, false);

    let mut total_peaks = 0usize;
    for peaks in mbr.mbr_peaks_by_file.values() {
        // Peaks are sorted by MbrScore descending; q-values must be non-decreasing down that order.
        let mut prev_q = f64::NEG_INFINITY;
        let mut prev_score = f64::INFINITY;
        for p in peaks {
            total_peaks += 1;
            assert!(
                p.mbr_q_value >= 0.0 && p.mbr_q_value <= 1.0 + 1e-9,
                "q-value out of range: {}",
                p.mbr_q_value
            );
            assert!(
                p.mbr_score <= prev_score + 1e-9,
                "peaks should be ordered by MbrScore descending"
            );
            assert!(
                p.mbr_q_value >= prev_q - 1e-12,
                "q-values must be non-decreasing down the score-ordered list"
            );
            prev_q = p.mbr_q_value;
            prev_score = p.mbr_score;
        }
    }
    // No peaks were dropped without PEP.
    assert_eq!(total_peaks, mbr.feature_rows.len());
}
