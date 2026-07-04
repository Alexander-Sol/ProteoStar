//! MBR RT-calibration-spline parity gate (PLAN.md P3.1).
//!
//! Drives the Rust MS2 engine ([`flashlfq_core::engine::run_msms`]) over the Phase-1 corpus
//! (`AllPSMs.psmtsv` + the two K562 mzMLs), then builds the **donor=K562_3 / acceptor=K562_4**
//! retention-time calibration spline with [`flashlfq_core::mbr::get_rt_cal_spline`] and diffs it
//! against the C# golden `parity/golden/MBR_rt_cal_spline.tsv`.
//!
//! The golden is produced by `parity/csharp_golden/` from the **real** `FlashLfqEngine` peaks
//! (`results.Peaks`, genuine C# ground truth); the spline-assembly logic on that side is a faithful
//! replica of the private `FlashLfqEngine.GetRtCalSpline` + `ChooseBestPeak` (the method is private
//! and depends on an `MbrScorer`, a P3.3 type, so it is replicated rather than invoked). Because the
//! peaks are real and the Rust port reproduces the identical filter/group/choose/pair/sort logic,
//! this validates the spline end-to-end over real data.
//!
//! Both sides sort the anchor points by `(donor_apex_rt, donor_modified_sequence)` before a
//! positional diff (a tie-robust order: two peptides can apex at the same MS1 scan RT, and C#'s
//! `OrderBy` is only stable in donor-RT, leaving such ties in dictionary-insertion order). The
//! modified sequence is compared exactly; all RT / mass / rt-diff cells relative `1e-6`.

use std::collections::HashMap;
use std::path::PathBuf;

use flashlfq_core::engine::{run_msms, Ms2QuantResult};
use flashlfq_core::mbr::{
    get_rt_cal_spline, predict_retention_time, DonorCriterion, DONOR_Q_VALUE_THRESHOLD,
    MAX_MBR_RT_WINDOW, NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR,
};
use flashlfq_core::psm_tsv::read_identifications;

const FILE_3: &str = "20100614_Velos1_TaGe_SA_K562_3";
const FILE_4: &str = "20100614_Velos1_TaGe_SA_K562_4";

fn test_data(relative: &str) -> PathBuf {
    flashlfq_core::mzlib_test_data(relative)
}

fn golden(relative: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("parity")
        .join("golden")
        .join(relative)
}

fn file_to_mzml() -> HashMap<String, PathBuf> {
    let mut m = HashMap::new();
    m.insert(FILE_3.to_string(), test_data(&format!("{FILE_3}.mzML")));
    m.insert(FILE_4.to_string(), test_data(&format!("{FILE_4}.mzML")));
    m
}

fn run_engine() -> Ms2QuantResult {
    let ids = read_identifications(test_data("AllPSMs.psmtsv"))
        .expect("AllPSMs.psmtsv parses into identifications");
    run_msms(ids, &file_to_mzml()).expect("engine runs over the corpus")
}

/// `|a-b| / max(|a|,|b|) < 1e-6`, with exact-equal (incl. both `0`) accepted.
fn floats_match(a: f64, b: f64) -> bool {
    if a == b {
        return true;
    }
    if !a.is_finite() || !b.is_finite() {
        return false;
    }
    let denom = a.abs().max(b.abs());
    (a - b).abs() / denom < 1e-6
}

/// One spline anchor row as compared on both sides.
#[derive(Debug, Clone)]
struct SplineRow {
    donor_modseq: String,
    donor_apex_rt: f64,
    acceptor_apex_rt: f64,
    rt_diff: f64,
    donor_peakfinding_mass: f64,
}

/// Sort by donor apex RT then donor modseq (tie-robust; see module docs).
fn sort_rows(rows: &mut [SplineRow]) {
    rows.sort_by(|a, b| {
        a.donor_apex_rt
            .partial_cmp(&b.donor_apex_rt)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.donor_modseq.cmp(&b.donor_modseq))
    });
}

fn read_golden_rows() -> Vec<SplineRow> {
    let path = golden("MBR_rt_cal_spline.tsv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("could not read golden {}: {e}", path.display()));
    let mut lines = text.lines();
    let header = lines.next().expect("golden has a header");
    assert_eq!(
        header,
        "donor_modified_sequence\tdonor_apex_rt\tacceptor_apex_rt\trt_diff\tdonor_peakfinding_mass"
    );
    let mut rows: Vec<SplineRow> = lines
        .filter(|l| !l.is_empty())
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            SplineRow {
                donor_modseq: c[0].to_string(),
                donor_apex_rt: c[1].parse().expect("donor_apex_rt"),
                acceptor_apex_rt: c[2].parse().expect("acceptor_apex_rt"),
                rt_diff: c[3].parse().expect("rt_diff"),
                donor_peakfinding_mass: c[4].parse().expect("donor_peakfinding_mass"),
            }
        })
        .collect();
    sort_rows(&mut rows);
    rows
}

#[test]
fn rt_cal_spline_matches_csharp_golden() {
    let result = run_engine();
    let donor_peaks = result
        .peaks_by_file
        .get(FILE_3)
        .expect("K562_3 peaks present");
    let acceptor_peaks = result
        .peaks_by_file
        .get(FILE_4)
        .expect("K562_4 peaks present");

    let spline = get_rt_cal_spline(
        donor_peaks,
        acceptor_peaks,
        DONOR_Q_VALUE_THRESHOLD,
        DonorCriterion::Score,
    );

    let mut rust_rows: Vec<SplineRow> = spline
        .calibration_curve
        .iter()
        .map(|p| SplineRow {
            donor_modseq: p
                .donor_peak
                .first_modified_sequence()
                .expect("donor peak has a modified sequence")
                .to_string(),
            donor_apex_rt: p.donor_peak.apex_retention_time(),
            acceptor_apex_rt: p
                .acceptor_peak
                .as_ref()
                .expect("paired acceptor peak")
                .apex_retention_time(),
            rt_diff: p.rt_diff,
            donor_peakfinding_mass: p
                .donor_peak
                .identification_peakfinding_masses
                .first()
                .copied()
                .expect("donor peak has a peakfinding mass"),
        })
        .collect();
    sort_rows(&mut rust_rows);

    let golden_rows = read_golden_rows();

    assert_eq!(
        rust_rows.len(),
        golden_rows.len(),
        "spline anchor-point count: rust {} vs golden {}",
        rust_rows.len(),
        golden_rows.len()
    );
    assert!(!golden_rows.is_empty(), "golden spline must be non-empty");

    // anchor_rt_diffs has one entry per curve point (both apex RTs are > 0 for traced peaks).
    assert_eq!(
        spline.anchor_rt_diffs.len(),
        spline.calibration_curve.len(),
        "every paired anchor contributes an RT diff"
    );

    for (i, (r, g)) in rust_rows.iter().zip(golden_rows.iter()).enumerate() {
        assert_eq!(
            r.donor_modseq, g.donor_modseq,
            "anchor {i}: donor modified sequence mismatch"
        );
        assert!(
            floats_match(r.donor_apex_rt, g.donor_apex_rt),
            "anchor {i} ({}): donor apex RT rust {} vs golden {}",
            r.donor_modseq,
            r.donor_apex_rt,
            g.donor_apex_rt
        );
        assert!(
            floats_match(r.acceptor_apex_rt, g.acceptor_apex_rt),
            "anchor {i} ({}): acceptor apex RT rust {} vs golden {}",
            r.donor_modseq,
            r.acceptor_apex_rt,
            g.acceptor_apex_rt
        );
        assert!(
            floats_match(r.rt_diff, g.rt_diff),
            "anchor {i} ({}): rt_diff rust {} vs golden {}",
            r.donor_modseq,
            r.rt_diff,
            g.rt_diff
        );
        assert!(
            floats_match(r.donor_peakfinding_mass, g.donor_peakfinding_mass),
            "anchor {i} ({}): donor peakfinding mass rust {} vs golden {}",
            r.donor_modseq,
            r.donor_peakfinding_mass,
            g.donor_peakfinding_mass
        );
    }
}

/// One predicted-RT row, keyed by donor modified sequence (unique per donor best-peak).
#[derive(Debug, Clone)]
struct PredictedRtRow {
    donor_modseq: String,
    predicted_rt: f64,
    width: f64,
}

fn read_predicted_rt_golden() -> Vec<PredictedRtRow> {
    let path = golden("MBR_predicted_rt.tsv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("could not read golden {}: {e}", path.display()));
    let mut lines = text.lines();
    let header = lines.next().expect("golden has a header");
    assert_eq!(header, "donor_modified_sequence\tpredicted_rt\twidth");
    let mut rows: Vec<PredictedRtRow> = lines
        .filter(|l| !l.is_empty())
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            PredictedRtRow {
                donor_modseq: c[0].to_string(),
                predicted_rt: c[1].parse().expect("predicted_rt"),
                width: c[2].parse().expect("width"),
            }
        })
        .collect();
    rows.sort_by(|a, b| a.donor_modseq.cmp(&b.donor_modseq));
    rows
}

/// PLAN.md P3.2a — `PredictRetentionTime` parity.
///
/// Builds the donor=K562_3 / acceptor=K562_4 spline, then predicts the acceptor-run RT of every
/// donor best-peak (ordered by peakfinding mass, the C# `donorPeaksMassOrdered` set) via
/// [`predict_retention_time`], and diffs against the C# golden `MBR_predicted_rt.tsv`. The golden's
/// `PredictRtReplica` is a faithful replica of the *internal* `FlashLfqEngine.PredictRetentionTime`
/// over real `RetentionTimeCalibDataPoint[]` objects built from the same real engine peaks. Both
/// sides key rows by donor modified sequence (unique per donor best-peak), so a positional diff
/// after sorting by sequence needs no tie handling. Predicted RT + width compared relative `1e-6`.
#[test]
fn predicted_rt_matches_csharp_golden() {
    let result = run_engine();
    let donor_peaks = result
        .peaks_by_file
        .get(FILE_3)
        .expect("K562_3 peaks present");
    let acceptor_peaks = result
        .peaks_by_file
        .get(FILE_4)
        .expect("K562_4 peaks present");

    let spline = get_rt_cal_spline(
        donor_peaks,
        acceptor_peaks,
        DONOR_Q_VALUE_THRESHOLD,
        DonorCriterion::Score,
    );

    let mut rust_rows: Vec<PredictedRtRow> = spline
        .donor_best_peaks_ordered_by_mass
        .iter()
        .map(|donor_peak| {
            let info = predict_retention_time(
                &spline.calibration_curve,
                donor_peak,
                MAX_MBR_RT_WINDOW,
                NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR,
            );
            PredictedRtRow {
                donor_modseq: donor_peak
                    .first_modified_sequence()
                    .expect("donor peak has a modified sequence")
                    .to_string(),
                predicted_rt: info.predicted_rt,
                width: info.width,
            }
        })
        .collect();
    rust_rows.sort_by(|a, b| a.donor_modseq.cmp(&b.donor_modseq));

    let golden_rows = read_predicted_rt_golden();

    assert_eq!(
        rust_rows.len(),
        golden_rows.len(),
        "predicted-RT row count: rust {} vs golden {}",
        rust_rows.len(),
        golden_rows.len()
    );
    assert!(!golden_rows.is_empty(), "golden predicted-RT table must be non-empty");

    for (i, (r, g)) in rust_rows.iter().zip(golden_rows.iter()).enumerate() {
        assert_eq!(
            r.donor_modseq, g.donor_modseq,
            "row {i}: donor modified sequence mismatch"
        );
        assert!(
            floats_match(r.predicted_rt, g.predicted_rt),
            "row {i} ({}): predicted RT rust {} vs golden {}",
            r.donor_modseq,
            r.predicted_rt,
            g.predicted_rt
        );
        assert!(
            floats_match(r.width, g.width),
            "row {i} ({}): width rust {} vs golden {}",
            r.donor_modseq,
            r.width,
            g.width
        );
    }
}
