//! Measures the TRUE chromatographic FWHM distribution from a spectra file's real XIC intensity
//! profiles (half-max crossings), to ground the trace kernel's assumed-FWHM parameter in data.
//!
//! This is deliberately independent of the detector: it traces XICs with the production
//! `get_all_xics`, then for each clean (rise-then-fall, internal-apex) XIC interpolates the RT where
//! intensity crosses half the apex on each side and reports FWHM = right − left. The detector's
//! claimed-peak "width" is the ±2σ *window*, NOT the peak — this measures the peak itself.
//!
//! Usage: cargo run --release --example fwhm_probe -- <spectra_file>

use flashlfq_core::peak_indexing::{read_ms1_scans, ExtractedIonChromatogram, PeakIndexingEngine};
use flashlfq_core::tolerance::PpmTolerance;

/// FWHM (minutes) of one XIC via linear-interpolated half-max crossings. `None` if the apex sits at
/// an edge or the profile does not fall below half-max on both sides (a truncated / monotonic trace).
fn fwhm_minutes(xic: &ExtractedIonChromatogram) -> Option<f64> {
    let pts: Vec<(f64, f64)> = xic
        .peaks
        .iter()
        .map(|p| (p.retention_time as f64, p.intensity as f64))
        .collect();
    let n = pts.len();
    if n < 3 {
        return None;
    }
    // Apex = max-intensity sample; require it to be internal (a real rise-then-fall peak).
    let mut ai = 0usize;
    for i in 1..n {
        if pts[i].1 > pts[ai].1 {
            ai = i;
        }
    }
    if ai == 0 || ai == n - 1 {
        return None;
    }
    let half = pts[ai].1 / 2.0;
    if half <= 0.0 {
        return None;
    }
    // Left crossing: nearest sample left of apex at or below half, interpolate to `half`.
    let mut left = None;
    for i in (0..ai).rev() {
        if pts[i].1 <= half {
            let (t0, y0) = pts[i];
            let (t1, y1) = pts[i + 1];
            left = Some(if y1 != y0 {
                t0 + (half - y0) * (t1 - t0) / (y1 - y0)
            } else {
                t0
            });
            break;
        }
    }
    // Right crossing.
    let mut right = None;
    for i in (ai + 1)..n {
        if pts[i].1 <= half {
            let (t0, y0) = pts[i - 1];
            let (t1, y1) = pts[i];
            right = Some(if y1 != y0 {
                t0 + (half - y0) * (t1 - t0) / (y1 - y0)
            } else {
                t1
            });
            break;
        }
    }
    match (left, right) {
        (Some(l), Some(r)) if r > l => Some(r - l),
        _ => None,
    }
}

fn pct(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let i = ((p / 100.0) * (sorted.len() as f64 - 1.0)).round() as usize;
    sorted[i.min(sorted.len() - 1)]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: fwhm_probe <spectra_file>");
        std::process::exit(2);
    }
    let path = &args[1];
    eprintln!("reading + indexing {path} ...");
    let scans = read_ms1_scans(path).expect("read spectra");
    let engine = PeakIndexingEngine::index_peaks(&scans).expect("index");
    eprintln!("  {} MS1 scans", engine.scan_info().len());

    // Trace every XIC (>=5 peaks) with the production tracer; 1 missed scan, generous RT cap so a
    // real peak is never clipped before its half-max shoulders.
    let ppm = PpmTolerance::new(10.0);
    eprintln!("tracing XICs ...");
    let xics = engine.get_all_xics(&ppm, 1, 2.0, 5);
    eprintln!("  {} XICs (>=5 peaks)", xics.len());

    // Measure FWHM on clean peaks; keep apex intensity for intensity-weighted views.
    let mut all: Vec<f64> = Vec::new();
    let mut strong: Vec<f64> = Vec::new(); // apex in the top intensity decile
    let mut apex_int: Vec<f64> = Vec::new();
    let mut measured: Vec<(f64, f64)> = Vec::new(); // (fwhm_min, apex_intensity)
    for x in &xics {
        if let Some(w) = fwhm_minutes(x) {
            let ai = x.apex_peak.intensity as f64;
            all.push(w);
            apex_int.push(ai);
            measured.push((w, ai));
        }
    }
    all.sort_by(|a, b| a.partial_cmp(b).unwrap());
    apex_int.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let int_cut = pct(&apex_int, 90.0);
    for (w, ai) in &measured {
        if *ai >= int_cut {
            strong.push(*w);
        }
    }
    strong.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let show = |label: &str, v: &[f64]| {
        if v.is_empty() {
            println!("{label}: n=0");
            return;
        }
        let s2 = |m: f64| m * 60.0; // minutes -> seconds
        println!(
            "{label}: n={}  FWHM sec  p10={:.1} p25={:.1} p50={:.1} p75={:.1} p90={:.1}  (median {:.2} s)",
            v.len(),
            s2(pct(v, 10.0)),
            s2(pct(v, 25.0)),
            s2(pct(v, 50.0)),
            s2(pct(v, 75.0)),
            s2(pct(v, 90.0)),
            s2(pct(v, 50.0)),
        );
    };
    println!("\n=== measured chromatographic FWHM (real XIC half-max) ===");
    println!("XICs measured (clean rise-then-fall): {} of {}", all.len(), xics.len());
    show("ALL measured    ", &all);
    show("STRONG (top 10% apex)", &strong);
}
