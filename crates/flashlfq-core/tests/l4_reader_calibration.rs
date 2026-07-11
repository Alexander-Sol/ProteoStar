//! L4 — reader calibration (report-only, **not** a pass/fail parity gate).
//!
//! P1.12 quantifies how far the Rust `mzdata` reader's per-MS1-scan `(m/z, intensity)` peak
//! lists drift from the file's actual decoded contents for `sliced-mzml.mzML`, so later
//! (L5/L6) tracing parity tolerances are set with eyes open. This is reader fidelity
//! (base64 + zlib decode → float reinterpretation), distinct from the L0–L3 chemistry gates.
//!
//! Ground truth is `parity/golden/L4_reader_peaks.tsv`, produced by
//! `parity/dump_l4_reader_peaks.py` — an independent raw decode of the mzML binary data
//! arrays. mzLib's mzML reader performs the identical base64+zlib+float decode but then **drops
//! peaks with intensity `< 0.01`** ("Remove Zero Intensity Peaks", `Readers/MzML/Mzml.cs`). The
//! Rust reader ([`read_ms1_scans`]) now applies the same drop so its peak lists match mzLib's
//! processed spectra (this is what L5/L6 tracing actually consumes; keeping the raw zero-filled
//! profile made the XIC walk cross gaps mzLib stops at — see P1.17b). To compare apples-to-apples,
//! the raw-decode golden is filtered with the identical `< 0.01` rule in [`load_golden`] before
//! diffing, so the assertion is "Rust reader == mzLib reader", not "Rust reader == raw file".
//!
//! The test prints a per-scan report (visible with `--nocapture`) and records the worst
//! offenders. Assertions are deliberately structural + loose: scan count and per-scan peak
//! counts must match exactly (a real reader invariant), and the float deltas must stay under a
//! generous sanity bound so a catastrophic reader mismatch still fails — the *magnitude* of
//! the drift is the deliverable, recorded in PLAN.md, not a tight tolerance.

use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

use flashlfq_core::peak_indexing::read_ms1_scans;

/// A scan's decoded peaks from the golden file.
struct GoldenScan {
    one_based_scan_number: i32,
    rt_minutes: f64,
    mz: Vec<f64>,
    intensity: Vec<f64>,
}

fn test_data(relative: &str) -> PathBuf {
    flashlfq_core::mzlib_test_data(relative)
}

fn golden_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("parity")
        .join("golden")
        .join("L4_reader_peaks.tsv")
}

/// Parses the golden TSV into per-scan decoded peak lists.
fn load_golden() -> Vec<GoldenScan> {
    let text = fs::read_to_string(golden_path()).expect("golden L4 file must exist (run dump_l4_reader_peaks.py)");
    let mut scans: BTreeMap<i32, GoldenScan> = BTreeMap::new();
    let mut order: Vec<i32> = Vec::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "S" => {
                let idx: i32 = cols[1].parse().unwrap();
                let one_based: i32 = cols[2].parse().unwrap();
                let rt: f64 = cols[3].parse().unwrap();
                let count: usize = cols[4].parse().unwrap();
                order.push(idx);
                scans.insert(
                    idx,
                    GoldenScan {
                        one_based_scan_number: one_based,
                        rt_minutes: rt,
                        mz: Vec::with_capacity(count),
                        intensity: Vec::with_capacity(count),
                    },
                );
            }
            "P" => {
                let idx: i32 = cols[1].parse().unwrap();
                let mz: f64 = cols[3].parse().unwrap();
                let inten: f64 = cols[4].parse().unwrap();
                let s = scans.get_mut(&idx).expect("P record before its S record");
                s.mz.push(mz);
                s.intensity.push(inten);
            }
            other => panic!("unexpected golden record type {other:?}"),
        }
    }

    order
        .into_iter()
        .map(|i| {
            let mut s = scans.remove(&i).unwrap();
            apply_zero_filter(&mut s);
            s
        })
        .collect()
}

/// Drops peaks with intensity `< 0.01` from a golden scan, replicating mzLib's reader (`Mzml.cs`,
/// "Remove Zero Intensity Peaks") — including its guard that an all-"zero" scan is kept verbatim.
fn apply_zero_filter(scan: &mut GoldenScan) {
    const ZERO_EQUIVALENT_INTENSITY: f64 = 0.01;
    let total = scan.intensity.len();
    let zero_count = scan
        .intensity
        .iter()
        .filter(|&&v| v < ZERO_EQUIVALENT_INTENSITY)
        .count();
    if zero_count == 0 || zero_count == total {
        return;
    }
    let mut mz = Vec::with_capacity(total - zero_count);
    let mut intensity = Vec::with_capacity(total - zero_count);
    for (m, i) in scan.mz.iter().zip(scan.intensity.iter()) {
        if *i >= ZERO_EQUIVALENT_INTENSITY {
            mz.push(*m);
            intensity.push(*i);
        }
    }
    scan.mz = mz;
    scan.intensity = intensity;
}

#[test]
fn l4_reader_calibration_report() {
    let scans = read_ms1_scans(test_data("sliced-mzml.mzML")).expect("sliced-mzml.mzML must be readable");
    let golden = load_golden();

    eprintln!("=== L4 reader calibration: mzdata vs raw mzML decode (sliced-mzml.mzML) ===");
    eprintln!(
        "mzdata MS1 scans: {}   golden MS1 scans: {}",
        scans.len(),
        golden.len()
    );
    assert_eq!(
        scans.len(),
        golden.len(),
        "MS1 scan count must match between mzdata and the raw decode"
    );

    // Aggregate worst offenders across all scans.
    let mut worst_dmz = 0.0_f64;
    let mut worst_dmz_where = (0usize, 0usize);
    let mut worst_dint_rel = 0.0_f64;
    let mut worst_dint_rel_where = (0usize, 0usize);
    let mut worst_dint_abs = 0.0_f64;
    let mut worst_drt = 0.0_f64;
    let mut total_peaks = 0usize;
    let mut count_mismatches = 0usize;

    for (si, (rust, gold)) in scans.iter().zip(golden.iter()).enumerate() {
        let n_rust = rust.mz.len();
        let n_gold = gold.mz.len();

        // Per-scan peak count is a structural invariant of a faithful decode.
        if n_rust != n_gold {
            count_mismatches += 1;
            eprintln!("  scan {si}: PEAK COUNT MISMATCH  mzdata={n_rust}  golden={n_gold}");
            continue;
        }

        let drt = (rust.retention_time - gold.rt_minutes).abs();
        if drt > worst_drt {
            worst_drt = drt;
        }

        let mut scan_worst_dmz = 0.0_f64;
        let mut scan_worst_dint_rel = 0.0_f64;
        for j in 0..n_rust {
            let dmz = (rust.mz[j] - gold.mz[j]).abs();
            if dmz > scan_worst_dmz {
                scan_worst_dmz = dmz;
            }
            if dmz > worst_dmz {
                worst_dmz = dmz;
                worst_dmz_where = (si, j);
            }

            let a = rust.intensity[j];
            let b = gold.intensity[j];
            let dint = (a - b).abs();
            if dint > worst_dint_abs {
                worst_dint_abs = dint;
            }
            let denom = a.abs().max(b.abs());
            let rel = if denom > 0.0 { dint / denom } else { 0.0 };
            if rel > scan_worst_dint_rel {
                scan_worst_dint_rel = rel;
            }
            if rel > worst_dint_rel {
                worst_dint_rel = rel;
                worst_dint_rel_where = (si, j);
            }
        }
        total_peaks += n_rust;

        eprintln!(
            "  scan {si:>2} (one-based {:>3}): peaks={n_rust:>5}  worst|Δm/z|={scan_worst_dmz:.3e}  worst rel|Δint|={scan_worst_dint_rel:.3e}  |Δrt|={drt:.3e}",
            gold.one_based_scan_number
        );
    }

    eprintln!("--- summary over {total_peaks} MS1 peaks ---");
    eprintln!(
        "worst |Δm/z|         = {worst_dmz:.6e} Da   at (scan {}, peak {})",
        worst_dmz_where.0, worst_dmz_where.1
    );
    eprintln!(
        "worst rel |Δintensity| = {worst_dint_rel:.6e}    at (scan {}, peak {})",
        worst_dint_rel_where.0, worst_dint_rel_where.1
    );
    eprintln!("worst abs |Δintensity| = {worst_dint_abs:.6e}");
    eprintln!("worst |Δretention time| = {worst_drt:.6e} min");

    // Structural invariant: a faithful decode must reproduce the peak counts exactly.
    assert_eq!(
        count_mismatches, 0,
        "every MS1 scan's peak count must match the raw decode"
    );

    // Loose sanity bounds (report-only): catch a catastrophic reader mismatch, not tight parity.
    assert!(
        worst_dmz < 1e-2,
        "m/z drift unexpectedly large: {worst_dmz} Da"
    );
    assert!(
        worst_dint_rel < 1e-3,
        "intensity drift unexpectedly large: rel {worst_dint_rel}"
    );
}
