//! End-to-end untargeted feature-detection runner.
//!
//! Reads a spectra file (mzML or Thermo `.raw`), runs the full untargeted pipeline
//! (index → `detect_features` → `refine_feature` → `resolve_charge_state_consensus`), writes the
//! resolved peptide-level features to a human-readable TSV, and — if given a base-FlashLFQ
//! `AllQuantifiedPeaks.tsv` — reports how many of those PSM-based peaks the untargeted run
//! independently rediscovered (mass + RT match).
//!
//! Usage:
//!   cargo run --release --example detect_features_tsv -- <spectra_file> <out.tsv> [reference.tsv]

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use flashlfq_core::deconvolution::{ClassicDeconvolutionParameters, Polarity};
use flashlfq_core::feature_refinement::{
    refine_feature, resolve_charge_state_consensus, RefinedFeature, ResolvedFeature,
};
use flashlfq_core::isotopic_envelope::mass_to_mz_f64;
use flashlfq_core::peak_indexing::{read_ms1_scans, PeakIndexingEngine};
use flashlfq_core::trace_kernel::{
    detect_features, median_ms1_scan_spacing_minutes, DetectedFeature, TraceKernelParameters,
};

/// Derives a sibling output path from the final path: `out.tsv` + tag `detected` -> `out.detected.tsv`.
fn sibling(out: &str, tag: &str) -> String {
    match out.strip_suffix(".tsv") {
        Some(stem) => format!("{stem}.{tag}.tsv"),
        None => format!("{out}.{tag}.tsv"),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: detect_features_tsv <spectra_file> <out.tsv> [reference.tsv]");
        std::process::exit(2);
    }
    let spectra_path = &args[1];
    let out_path = &args[2];
    let reference_path = args.get(3);

    let detected_path = sibling(out_path, "detected");
    let refined_path = sibling(out_path, "refined");
    eprintln!("output files:");
    eprintln!("  detected (pre-refinement): {detected_path}");
    eprintln!("  refined (post-decon):      {refined_path}");
    eprintln!("  resolved (final):          {out_path}");

    // --- read + index --------------------------------------------------------------------------
    let t0 = Instant::now();
    eprintln!("reading MS1 scans from {spectra_path} ...");
    let scans = read_ms1_scans(spectra_path).expect("failed to read spectra file");
    let engine = PeakIndexingEngine::index_peaks(&scans).expect("no indexable MS1 peaks");
    let n_peaks: usize = scans.iter().map(|s| s.mz.len()).sum();
    let total_intensity: f64 = scans.iter().flat_map(|s| s.intensity.iter()).sum();
    eprintln!(
        "  {} MS1 scans, {} peaks, ΣTIC {:.3e}, median scan spacing {:.4} min  ({:?})",
        scans.len(),
        n_peaks,
        total_intensity,
        median_ms1_scan_spacing_minutes(engine.scan_info()),
        t0.elapsed()
    );

    // --- detect --------------------------------------------------------------------------------
    let params = TraceKernelParameters {
        ppm_tolerance: 10.0,
        min_seed_intensity: 1000.0,
        coverage_target: 0.75,
        ..TraceKernelParameters::default()
    }
    .with_rt_from_scans(engine.scan_info(), 36.0);
    eprintln!(
        "detecting (charge {}..={}, {} ppm, σ_rt {:.4} min, ±{} scans, seed floor {:.0}, coverage {:.0}%) ...",
        params.min_charge,
        params.max_charge,
        params.ppm_tolerance,
        params.rt_sigma_minutes,
        params.half_window_scans,
        params.min_seed_intensity,
        params.coverage_target * 100.0
    );
    let t1 = Instant::now();
    let detected = detect_features(&engine, &params);
    let detected_intensity: f64 = detected.iter().map(|f| f.summed_intensity).sum();
    eprintln!(
        "  {} features detected, explained {:.1}% of ΣTIC  ({:?})",
        detected.len(),
        100.0 * detected_intensity / total_intensity,
        t1.elapsed()
    );
    write_detected_tsv(&detected_path, &detected);
    eprintln!("  wrote {} detected features -> {detected_path}", detected.len());

    // --- refine --------------------------------------------------------------------------------
    let avg = flashlfq_core::spectral_averaging::SpectralAveragingParameters::default();
    let decon = ClassicDeconvolutionParameters::new(
        params.min_charge,
        params.max_charge,
        10.0,
        3.0,
        Polarity::Positive,
    );
    let t2 = Instant::now();
    let mut refined: Vec<RefinedFeature> = Vec::with_capacity(detected.len());
    let progress_every = 1000usize;
    for (i, f) in detected.iter().enumerate() {
        if let Some(r) = refine_feature(f, &scans, &avg, &decon) {
            refined.push(r);
        }
        if (i + 1) % progress_every == 0 || i + 1 == detected.len() {
            eprintln!(
                "    refined {}/{} ({} kept)  [{:?}]",
                i + 1,
                detected.len(),
                refined.len(),
                t2.elapsed()
            );
        }
    }
    eprintln!(
        "  {} / {} features refined against averaged composites  ({:?})",
        refined.len(),
        detected.len(),
        t2.elapsed()
    );
    write_refined_tsv(&refined_path, &refined);
    eprintln!("  wrote {} refined features -> {refined_path}", refined.len());

    // --- resolve charge-state consensus --------------------------------------------------------
    let t3 = Instant::now();
    let resolved = resolve_charge_state_consensus(&refined, 10.0, 0.1);
    eprintln!(
        "  {} peptide-level features after charge-state consensus  ({:?})",
        resolved.len(),
        t3.elapsed()
    );

    write_tsv(out_path, &resolved);
    eprintln!("wrote {} -> {}", resolved.len(), out_path);

    if let Some(ref_path) = reference_path {
        compare_to_reference(ref_path, &resolved);
    }
}

/// Writes the raw detected features (pre-refinement, straight from the trace kernel) to a TSV,
/// sorted by summed intensity descending. Lets you inspect what the detector alone produced.
fn write_detected_tsv(path: &str, detected: &[DetectedFeature]) {
    let mut rows: Vec<&DetectedFeature> = detected.iter().collect();
    rows.sort_by(|a, b| b.summed_intensity.total_cmp(&a.summed_intensity));
    let f = File::create(path).expect("cannot create detected tsv");
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "Monoisotopic Mass\tCharge\tMono m/z\tApex RT\tRT Start\tRT End\tSummed Intensity\t\
         Detector Score\tNum Isotopes\tNum Peaks"
    )
    .unwrap();
    for d in rows {
        writeln!(
            w,
            "{:.5}\t{}\t{:.5}\t{:.4}\t{:.4}\t{:.4}\t{:.4e}\t{:.4e}\t{}\t{}",
            d.monoisotopic_mass,
            d.charge,
            d.mono_mz,
            d.apex_rt,
            d.start_rt,
            d.end_rt,
            d.summed_intensity,
            d.score,
            d.num_isotopes_observed,
            d.peaks.len()
        )
        .unwrap();
    }
    w.flush().unwrap();
}

/// Writes the refined features (post composite-deconvolution, pre charge-consensus) to a TSV,
/// sorted by summed intensity descending.
fn write_refined_tsv(path: &str, refined: &[RefinedFeature]) {
    let mut rows: Vec<&RefinedFeature> = refined.iter().collect();
    rows.sort_by(|a, b| {
        b.detected
            .summed_intensity
            .total_cmp(&a.detected.summed_intensity)
    });
    let f = File::create(path).expect("cannot create refined tsv");
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "Refined Monoisotopic Mass\tCharge\tApex RT\tSummed Intensity\tDecon Score\t\
         Num Candidate Masses\tDetector Mono Mass"
    )
    .unwrap();
    for r in rows {
        writeln!(
            w,
            "{:.5}\t{}\t{:.4}\t{:.4e}\t{:.4e}\t{}\t{:.5}",
            r.refined_monoisotopic_mass,
            r.refined_charge,
            r.detected.apex_rt,
            r.detected.summed_intensity,
            r.decon_score,
            r.candidate_masses.len(),
            r.detected.monoisotopic_mass
        )
        .unwrap();
    }
    w.flush().unwrap();
}

/// Writes the resolved features to a human-readable TSV, sorted by summed intensity descending.
fn write_tsv(path: &str, resolved: &[ResolvedFeature]) {
    let mut rows: Vec<&ResolvedFeature> = resolved.iter().collect();
    rows.sort_by(|a, b| b.summed_intensity.total_cmp(&a.summed_intensity));

    let f = File::create(path).expect("cannot create output tsv");
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "Monoisotopic Mass\tCharge States\tNum Charge States\tPrimary Charge\tMono m/z (primary)\t\
         Apex RT\tRT Start\tRT End\tSummed Intensity\tCross-Charge Support\tNum Members"
    )
    .unwrap();
    for r in rows {
        // Primary charge = the tallest member's charge.
        let primary_charge = r
            .members
            .iter()
            .max_by(|a, b| {
                a.detected
                    .summed_intensity
                    .total_cmp(&b.detected.summed_intensity)
            })
            .map(|m| m.refined_charge)
            .unwrap_or(0);
        let mono_mz = if primary_charge != 0 {
            mass_to_mz_f64(r.monoisotopic_mass, primary_charge)
        } else {
            0.0
        };
        let charges: Vec<String> = r.charge_states.iter().map(|c| c.to_string()).collect();
        writeln!(
            w,
            "{:.5}\t{}\t{}\t{}\t{:.5}\t{:.4}\t{:.4}\t{:.4}\t{:.4e}\t{}\t{}",
            r.monoisotopic_mass,
            charges.join(";"),
            r.charge_states.len(),
            primary_charge,
            mono_mz,
            r.apex_rt,
            r.start_rt,
            r.end_rt,
            r.summed_intensity,
            r.cross_charge_support,
            r.members.len()
        )
        .unwrap();
    }
    w.flush().unwrap();
}

/// Loads the base-FlashLFQ `AllQuantifiedPeaks.tsv` and reports how many of its peaks a resolved
/// feature independently matches by monoisotopic mass (±20 ppm) and apex RT (±0.3 min). The
/// reference is the calibrated file, so a small mass/RT offset vs. the raw is expected — tolerances
/// are generous accordingly.
fn compare_to_reference(path: &str, resolved: &[ResolvedFeature]) {
    let text = match std::fs::read_to_string(path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("could not read reference {path}: {e}");
            return;
        }
    };
    let mut lines = text.lines();
    let header = match lines.next() {
        Some(h) => h,
        None => return,
    };
    let cols: HashMap<&str, usize> = header.split('\t').enumerate().map(|(i, c)| (c, i)).collect();
    let mass_i = cols["Peptide Monoisotopic Mass"];
    let rt_i = cols["Peak RT Apex"];
    let charge_i = cols["Peak Charge"];
    let seq_i = cols["Full Sequence"];

    // Reference rows: (mass, apex_rt, charge, sequence). Skip malformed / empty-mass rows.
    struct RefPeak {
        mass: f64,
        rt: f64,
        charge: i32,
        seq: String,
    }
    let mut refs: Vec<RefPeak> = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let (m, rt) = match (
            fields.get(mass_i).and_then(|s| s.parse::<f64>().ok()),
            fields.get(rt_i).and_then(|s| s.parse::<f64>().ok()),
        ) {
            (Some(m), Some(rt)) => (m, rt),
            _ => continue,
        };
        let charge = fields
            .get(charge_i)
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or(0);
        let seq = fields.get(seq_i).map(|s| s.to_string()).unwrap_or_default();
        refs.push(RefPeak { mass: m, rt, charge, seq });
    }

    const MASS_PPM: f64 = 20.0;
    const RT_MIN: f64 = 0.3;

    let mut matched = 0usize;
    let mut matched_with_charge = 0usize;
    let mut unmatched_examples: Vec<String> = Vec::new();
    for rp in &refs {
        let hit = resolved.iter().find(|f| {
            (f.monoisotopic_mass - rp.mass).abs() / rp.mass * 1e6 <= MASS_PPM
                && (f.apex_rt - rp.rt).abs() <= RT_MIN
        });
        match hit {
            Some(f) => {
                matched += 1;
                if f.charge_states.contains(&rp.charge) {
                    matched_with_charge += 1;
                }
            }
            None => {
                if unmatched_examples.len() < 15 {
                    unmatched_examples.push(format!(
                        "{:.4} Da  RT {:.3}  z{}  {}",
                        rp.mass, rp.rt, rp.charge, rp.seq
                    ));
                }
            }
        }
    }

    eprintln!("\n=== comparison to base FlashLFQ ({}) ===", path);
    eprintln!("  reference PSM-based peaks: {}", refs.len());
    eprintln!("  resolved untargeted features: {}", resolved.len());
    eprintln!(
        "  reference peaks rediscovered (±{:.0} ppm mass, ±{:.1} min RT): {} / {}  ({:.1}%)",
        MASS_PPM,
        RT_MIN,
        matched,
        refs.len(),
        100.0 * matched as f64 / refs.len().max(1) as f64
    );
    eprintln!(
        "    ...of which the charge state also matched: {} ({:.1}%)",
        matched_with_charge,
        100.0 * matched_with_charge as f64 / refs.len().max(1) as f64
    );
    if !unmatched_examples.is_empty() {
        eprintln!("  examples of reference peaks NOT rediscovered:");
        for e in &unmatched_examples {
            eprintln!("    - {e}");
        }
    }
}
