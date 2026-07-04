//! Quick diagnostic: what fraction of the total ion current is *isotopically structured*?
//!
//! For every MS1 scan we keep a peak's intensity only if that peak belongs to a run of at least
//! three peaks (itself plus two others) spaced at `C13_MINUS_C12 / z` Th for some charge state
//! `z` in `[MIN_Z, MAX_Z]`. Summing the intensity of the kept peaks and dividing by ΣTIC yields
//! the fraction of signal that even *has* isotope-envelope structure — an upper bound on what any
//! deconvolution-based feature detector could ever attribute to real peptide isotope patterns.
//!
//! This is deliberately per-scan and cheap (no RT tracing, no deconvolution). It is a ceiling, not
//! a detection: a peak counts as isotopic if it lines up under *any* allowed charge, regardless of
//! envelope-shape plausibility, so the true explainable fraction is somewhat below this number.
//!
//! Usage:
//!   cargo run --release --example isotope_signal_fraction -- <spectra_file>
//!
//! Env knobs:
//!   PPM=10            m/z match tolerance for isotope spacing (default 10 ppm)
//!   MIN_Z=2 MAX_Z=6   charge-state range to test (default 2..=6)
//!   MIN_INTENSITY=0   ignore peaks below this intensity when forming/counting chains

use flashlfq_core::isotopic_envelope::C13_MINUS_C12;
use flashlfq_core::peak_indexing::{read_ms1_scans, Scan};

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}
fn env_i32(key: &str, default: i32) -> i32 {
    std::env::var(key).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}

/// Index of the peak whose m/z is closest to `target`, provided it lies within `tol_th` Th.
/// `mz` is assumed ascending (as centroided spectra are), so a binary search + two-neighbour check
/// suffices.
fn nearest_within(mz: &[f64], target: f64, tol_th: f64) -> Option<usize> {
    if mz.is_empty() {
        return None;
    }
    let idx = mz.partition_point(|&m| m < target);
    let mut best: Option<usize> = None;
    let mut best_d = tol_th;
    for c in [idx.wrapping_sub(1), idx] {
        if c < mz.len() {
            let d = (mz[c] - target).abs();
            if d <= best_d {
                best_d = d;
                best = Some(c);
            }
        }
    }
    best
}

/// Marks, for one scan, every peak that belongs to a >=3-peak isotope run under some allowed charge.
/// Returns a parallel `Vec<bool>` over `scan`'s peaks.
fn mark_isotopic(scan: &Scan, min_z: i32, max_z: i32, ppm: f64, min_intensity: f64) -> Vec<bool> {
    let n = scan.mz.len();
    let mut isotopic = vec![false; n];
    if n < 3 {
        return isotopic;
    }
    for z in min_z..=max_z {
        let delta = C13_MINUS_C12 / z as f64;
        // next[i] = index of the peak one isotope step up from peak i (usize::MAX = none).
        let mut next = vec![usize::MAX; n];
        let mut has_prev = vec![false; n];
        for i in 0..n {
            if scan.intensity[i] < min_intensity {
                continue;
            }
            let target = scan.mz[i] + delta;
            let tol = target * ppm / 1e6;
            if let Some(j) = nearest_within(&scan.mz, target, tol) {
                // `next` links strictly upward in m/z (delta > tol), so chains never cycle.
                if j != i && j != usize::MAX && scan.intensity[j] >= min_intensity {
                    next[i] = j;
                    has_prev[j] = true;
                }
            }
        }
        // Walk each chain from its head (a peak that nothing points to); mark all if it reaches 3.
        for i in 0..n {
            if has_prev[i] {
                continue;
            }
            let mut chain = Vec::new();
            let mut cur = i;
            while cur != usize::MAX {
                chain.push(cur);
                cur = next[cur];
            }
            if chain.len() >= 3 {
                for c in chain {
                    isotopic[c] = true;
                }
            }
        }
    }
    isotopic
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: isotope_signal_fraction <spectra_file>");
        std::process::exit(2);
    }
    let spectra_path = &args[1];

    let ppm = env_f64("PPM", 10.0);
    let min_z = env_i32("MIN_Z", 2);
    let max_z = env_i32("MAX_Z", 6);
    let min_intensity = env_f64("MIN_INTENSITY", 0.0);

    eprintln!("reading MS1 scans from {spectra_path} ...");
    let scans = read_ms1_scans(spectra_path).expect("failed to read spectra file");

    let mut total_tic = 0.0f64;
    let mut isotopic_tic = 0.0f64;
    let mut total_peaks = 0usize;
    let mut isotopic_peaks = 0usize;

    for scan in &scans {
        let flags = mark_isotopic(scan, min_z, max_z, ppm, min_intensity);
        for (i, &inten) in scan.intensity.iter().enumerate() {
            total_tic += inten;
            total_peaks += 1;
            if flags[i] {
                isotopic_tic += inten;
                isotopic_peaks += 1;
            }
        }
    }

    let tic_frac = if total_tic > 0.0 { 100.0 * isotopic_tic / total_tic } else { 0.0 };
    let peak_frac = if total_peaks > 0 {
        100.0 * isotopic_peaks as f64 / total_peaks as f64
    } else {
        0.0
    };

    println!("spectra file:        {spectra_path}");
    println!("MS1 scans:           {}", scans.len());
    println!("charge range tested: {min_z}..={max_z}   tol: {ppm} ppm   min intensity: {min_intensity}");
    println!("ΣTIC (all peaks):    {total_tic:.4e}");
    println!("ΣTIC (isotopic):     {isotopic_tic:.4e}");
    println!(
        "isotope-structured signal: {tic_frac:.1}% of ΣTIC   ({isotopic_peaks} / {total_peaks} peaks = {peak_frac:.1}%)"
    );
}
