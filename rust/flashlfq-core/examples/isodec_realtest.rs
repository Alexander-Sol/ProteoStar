//! Real-data check of the native IsoDec charge predictor: detect features on a top-down raw, then feed
//! each detected feature's claimed peaks (its observed isotope envelope) to IsoDec and compare IsoDec's
//! charge to the detector's. High agreement = the port works on real chimeric envelopes.
//!
//! Usage: isodec_realtest <raw_or_mzml>   (defaults to golden.raw)

use std::collections::HashMap;

use flashlfq_core::isodec::default_model;
use flashlfq_core::peak_indexing::{read_ms1_scans, PeakIndexingEngine};
use flashlfq_core::trace_kernel::{detect_features, TraceKernelParameters};

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| r"D:\CD_FDR_MSV000082367\raw\golden.raw".to_string());
    eprintln!("reading {path} ...");
    let scans = read_ms1_scans(&path).expect("read scans");
    let engine = PeakIndexingEngine::index_peaks(&scans).expect("index");

    // Top-down-ish detector config (charge 1..60, long comb).
    let params = TraceKernelParameters {
        min_charge: 1,
        max_charge: 60,
        max_isotopes: 60,
        min_isotopes_observed: 3,
        min_seed_intensity: 5000.0, // higher floor → fewer, stronger features → fast real check
        ..TraceKernelParameters::default()
    }
    .with_rt_from_index(&engine, 36.0);

    eprintln!("detecting ...");
    let features = detect_features(&engine, &params);
    eprintln!("{} features; predicting IsoDec charge from each envelope ...", features.len());

    let model = default_model();
    let mut agree = 0usize;
    let mut total = 0usize;
    let mut agree_strong = 0usize;
    let mut total_strong = 0usize;
    let mut total_dbg = 0usize;
    let mut off: HashMap<i32, usize> = HashMap::new();
    let mut zero_call = 0usize;
    for f in &features {
        // IsoDec expects a single-scan m/z cluster of ONE envelope. The detected feature's `peaks` span
        // the whole RT window (each m/z repeats per scan), so restrict to the apex scan's claimed teeth.
        let mut cluster: Vec<(f64, f32)> = f
            .peaks
            .iter()
            .filter(|p| p.zero_based_scan_index == f.apex_scan_index)
            .map(|p| (p.mz as f64, p.intensity))
            .collect();
        if cluster.len() < 3 {
            continue;
        }
        cluster.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mz: Vec<f64> = cluster.iter().map(|(m, _)| *m).collect();
        let inten: Vec<f32> = cluster.iter().map(|(_, i)| *i).collect();
        let z_iso = model.predict_charge(&mz, &inten);
        // Infer the true charge from the actual median adjacent-peak spacing: z ≈ 1.0033 / gap.
        if z_iso != f.charge && f.charge > 40 && total_dbg < 20 {
            total_dbg += 1;
            let mut gaps: Vec<f64> = mz.windows(2).map(|w| w[1] - w[0]).filter(|&g| g > 1e-6).collect();
            gaps.sort_by(f64::total_cmp);
            let med = gaps.get(gaps.len() / 2).copied().unwrap_or(0.0);
            let z_spacing = if med > 0.0 { (1.0033548 / med).round() as i32 } else { 0 };
            eprintln!(
                "  feat z_det={:>2} z_iso={:>2} z_from_spacing={:>2}  npk={} medgap={:.4}",
                f.charge, z_iso, z_spacing, mz.len(), med
            );
        }
        total += 1;
        if z_iso == 0 {
            zero_call += 1;
        }
        if z_iso == f.charge {
            agree += 1;
        }
        // Stratify by envelope quality: a well-observed envelope has many isotope teeth.
        if cluster.len() >= 8 {
            total_strong += 1;
            if z_iso == f.charge {
                agree_strong += 1;
            }
        }
        *off.entry(z_iso - f.charge).or_insert(0) += 1;
    }
    println!("\nfeatures scored: {total}");
    println!("IsoDec == detector charge (all):        {agree} ({:.1}%)", 100.0 * agree as f64 / total.max(1) as f64);
    println!("IsoDec == detector charge (>=8 teeth):  {agree_strong}/{total_strong} ({:.1}%)", 100.0 * agree_strong as f64 / total_strong.max(1) as f64);
    println!("IsoDec no-call (class 0):  {zero_call}");
    let mut offs: Vec<(&i32, &usize)> = off.iter().collect();
    offs.sort_by_key(|(k, _)| **k);
    println!("(isodec - detector) charge offset histogram (top entries):");
    let mut shown: Vec<(i32, usize)> = offs.iter().map(|(k, v)| (**k, **v)).collect();
    shown.sort_by(|a, b| b.1.cmp(&a.1));
    for (k, v) in shown.into_iter().take(9) {
        println!("   off {k:+}: {v}");
    }
}
