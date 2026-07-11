//! Dump MS1 peak intensities for noise-floor / seed-floor analysis.
//!
//! Reads a spectra file (mzML or Thermo `.raw`) via the same `read_ms1_scans` the detector
//! uses, then writes:
//!   <out_prefix>_intensities.f32   — every MS1 peak intensity, little-endian f32 (flat)
//!   <out_prefix>_per_scan.tsv      — per-scan: rt_min, n_peaks, tic, max, p50, p05, p01
//!
//! Usage:
//!   cargo run --release --example dump_ms1_intensities -- <spectra_file> <out_prefix>

use std::fs::File;
use std::io::{BufWriter, Write};

use flashlfq_core::peak_indexing::read_ms1_scans;

fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = (((sorted.len() - 1) as f64) * p).round() as usize;
    sorted[idx]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: dump_ms1_intensities <spectra_file> <out_prefix>");
        std::process::exit(1);
    }
    let path = &args[1];
    let prefix = &args[2];

    eprintln!("reading MS1 scans from {path} ...");
    let scans = read_ms1_scans(path).expect("failed to read MS1 scans");
    eprintln!("{} MS1 scans", scans.len());

    let bin_path = format!("{prefix}_intensities.f32");
    let scan_path = format!("{prefix}_per_scan.tsv");
    let mut bin = BufWriter::new(File::create(&bin_path).expect("create f32"));
    let mut per = BufWriter::new(File::create(&scan_path).expect("create tsv"));
    writeln!(per, "scan_index\trt_min\tn_peaks\ttic\tmax\tp50\tp05\tp01").unwrap();

    let mut total_peaks: u64 = 0;
    for (si, scan) in scans.iter().enumerate() {
        // flat intensity dump
        for &i in &scan.intensity {
            bin.write_all(&(i as f32).to_le_bytes()).unwrap();
        }
        total_peaks += scan.intensity.len() as u64;

        // per-scan stats
        let mut v: Vec<f64> = scan.intensity.iter().copied().filter(|&x| x > 0.0).collect();
        v.sort_by(|a, b| a.total_cmp(b));
        let tic: f64 = v.iter().sum();
        let maxv = v.last().copied().unwrap_or(0.0);
        writeln!(
            per,
            "{}\t{:.4}\t{}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}",
            si,
            scan.retention_time,
            v.len(),
            tic,
            maxv,
            percentile(&v, 0.50),
            percentile(&v, 0.05),
            percentile(&v, 0.01),
        )
        .unwrap();
    }
    bin.flush().unwrap();
    per.flush().unwrap();
    eprintln!("wrote {total_peaks} peak intensities -> {bin_path}");
    eprintln!("wrote per-scan stats -> {scan_path}");
}
