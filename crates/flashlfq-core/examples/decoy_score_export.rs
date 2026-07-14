//! Score every refined feature against a given isotope-envelope model with BOTH a 2-D (apex-scan,
//! m/z-only) cosine and a 3-D (m/z × RT window) cosine, and emit a tidy scores table for the decoy
//! analysis. Used to compare how well each score separates real (target) features from decoys.
//!
//! The 2-D score is the union-grid envelope-fit cosine at the single apex scan (the same shape the
//! detector's `Decon Score` uses). The 3-D score extends the template across the feature's RT window:
//! `T[k,s] = w_k · g(s)` with `g` a Gaussian elution weight, scored against the observed intensity at
//! each isotope tooth `k` and scan `s`, plus an unexplained-peak penalty per scan (the union grid).
//! Both use the SAME envelope model + tooth lattice, so a decoy feature is scored against the very
//! decoy envelope it was detected with.
//!
//! Usage:
//!   decoy_score_export <raw> <refined.tsv> <out_scores.tsv> <model> <spacing_scale>
//! where <model> ∈ averagine|chlorinated|rotated|cbp|custom|shuffled|hybrid and <spacing_scale> is the
//! tooth-spacing multiplier (1.0 for the physical ¹³C lattice, e.g. 0.9368 for the 0.94-Da shifted
//! comb). Emits `mono<TAB>z<TAB>rt<TAB>score2d<TAB>score3d`, one row per refined feature.

use std::io::Write;

use flashlfq_core::deconvolution::EnvelopeModel;
use flashlfq_core::isotopic_envelope::{mass_to_mz_f64, C13_MINUS_C12};
use flashlfq_core::peak_indexing::{read_ms1_scans, Scan};
use flashlfq_core::trace_kernel::{tooth_offsets, LatticeMode};

const TEMPLATE_MIN_WEIGHT: f64 = 1e-3;
/// Default template length. Bottom-up envelopes fit inside ~24 teeth; heavy top-down proteoforms span
/// 40–60 significant peaks with the mode well above tooth 24, so `SCORE_MAX_ISOTOPES` raises this for
/// top-down (the ≥`MIN_REL` band is trimmed per feature anyway, so an oversized cap is harmless).
const TEMPLATE_MAX_ISOTOPES_DEFAULT: usize = 24;
const TOL_PPM: f64 = 20.0;
const MIN_REL: f64 = 0.2;
/// Default RT Gaussian σ (minutes) for the 3-D elution template — ~0.15 min ≈ 21 s FWHM, in the CA/Lumos
/// peak-width range. Top-down species elute broader (trace half-width ~60 s), so `SCORE_RT_SIGMA` widens
/// it. The same weighting is laid for target and decoy, so it does not bias the contrast.
const RT_SIGMA_MIN_DEFAULT: f64 = 0.15;

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}
fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}

fn model_from_name(name: &str) -> EnvelopeModel {
    match name {
        "averagine" => EnvelopeModel::Averagine,
        "chlorinated" | "decoy" => EnvelopeModel::ChlorinatedDecoy,
        "rotated" => EnvelopeModel::RotatedDecoy,
        "cbp" | "chloroboro" => EnvelopeModel::ChloroBoroPhosphateDecoy,
        "custom" | "weird" => EnvelopeModel::CustomDecoy,
        "shuffled" => EnvelopeModel::ShuffledDecoy,
        "hybrid" => EnvelopeModel::HybridDecoy,
        other => panic!("unknown model '{other}'"),
    }
}

struct Feature {
    mono: f64,
    z: i32,
    rt: f64,
}

/// Read (Refined Monoisotopic Mass, Charge, Apex RT) from a runner refined TSV.
fn read_refined(path: &str) -> Vec<Feature> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|_| panic!("read {path}"));
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name);
    let (cm, cz, crt) = (
        col("Refined Monoisotopic Mass"),
        col("Charge"),
        col("Apex RT"),
    );
    let mut out = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |c: Option<usize>| c.and_then(|i| f.get(i)).and_then(|s| s.parse::<f64>().ok());
        if let (Some(mono), Some(z), Some(rt)) = (get(cm), get(cz), get(crt)) {
            out.push(Feature { mono, z: z.round() as i32, rt });
        }
    }
    out
}

/// Index of the scan whose retention time is nearest `rt`.
fn nearest_scan(scans: &[Scan], rt: f64) -> usize {
    scans
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| (a.retention_time - rt).abs().total_cmp(&(b.retention_time - rt).abs()))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Nearest peak intensity within `tol_ppm` of `target_mz` in one scan's ascending m/z array, and the
/// index of the matched peak (so the caller can exclude it from the unexplained-peak set). `None` if
/// no peak is within tolerance.
fn nearest_peak(mz: &[f64], intensity: &[f64], target_mz: f64, tol_ppm: f64) -> Option<(usize, f64)> {
    if mz.is_empty() {
        return None;
    }
    let mut lo = mz.partition_point(|&m| m < target_mz);
    // Check the two candidates straddling target_mz.
    let mut best: Option<(usize, f64)> = None;
    for &j in &[lo.wrapping_sub(1), lo] {
        if j < mz.len() {
            let d = (mz[j] - target_mz).abs();
            if best.map_or(true, |(_, bd)| d < bd) {
                best = Some((j, d));
            }
        }
    }
    let _ = &mut lo;
    match best {
        Some((j, d)) if d / target_mz * 1e6 <= tol_ppm => Some((j, intensity[j])),
        _ => None,
    }
}

fn gaussian(delta: f64, sigma: f64) -> f64 {
    (-0.5 * (delta / sigma).powi(2)).exp()
}

/// Accumulate the template/observed contributions of ONE scan onto `(dot, nt2, no2)` for the union-grid
/// cosine: predicted teeth (matched or absent) plus unexplained in-window peaks. `g` is the scan's RT
/// weight (1.0 for the 2-D apex-only score).
#[allow(clippy::too_many_arguments)]
fn accumulate_scan(
    scan: &Scan,
    mono_mz: f64,
    offsets: &[f64],
    sig: &[(usize, f64)],
    kmax: usize,
    step: f64,
    g: f64,
    tol_ppm: f64,
    acc: &mut (f64, f64, f64),
) {
    let win_lo = mono_mz - 0.5 * step;
    let win_hi = mono_mz + offsets[kmax] + step;
    let lo = scan.mz.partition_point(|&m| m < win_lo);
    let hi = scan.mz.partition_point(|&m| m <= win_hi);
    let mut matched = vec![false; hi.saturating_sub(lo)];

    for &(k, w) in sig {
        let target = mono_mz + offsets[k];
        let o = match nearest_peak(&scan.mz[lo..hi], &scan.intensity[lo..hi], target, tol_ppm) {
            Some((j, inten)) => {
                matched[j] = true;
                inten
            }
            None => 0.0,
        };
        let t = w * g;
        acc.0 += t * o;
        acc.1 += t * t;
        acc.2 += o * o;
    }
    // Unexplained observed peaks in the window (template 0) → penalise fraction-explained.
    for j in lo..hi {
        if !matched[j - lo] && scan.intensity[j] > 0.0 {
            acc.2 += scan.intensity[j] * scan.intensity[j];
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 6 {
        eprintln!("usage: decoy_score_export <raw> <refined.tsv> <out.tsv> <model> <spacing_scale>");
        std::process::exit(2);
    }
    let raw = &args[1];
    let refined_path = &args[2];
    let out_path = &args[3];
    let model = model_from_name(&args[4]);
    let spacing_scale: f64 = args[5].parse().expect("spacing_scale must be a float");
    let lattice = if (spacing_scale - 1.0).abs() < 1e-9 {
        LatticeMode::Uniform
    } else {
        LatticeMode::Scaled(spacing_scale)
    };

    let template_max_isotopes = env_usize("SCORE_MAX_ISOTOPES", TEMPLATE_MAX_ISOTOPES_DEFAULT);
    let rt_sigma_min = env_f64("SCORE_RT_SIGMA", RT_SIGMA_MIN_DEFAULT);
    let rt_half_window_min = 2.0 * rt_sigma_min;

    let scans = read_ms1_scans(raw.clone()).expect("read raw");
    let feats = read_refined(refined_path);
    eprintln!(
        "scoring {} features from {refined_path} against {:?} (scale {spacing_scale}, max_iso {template_max_isotopes}, rt_sigma {rt_sigma_min})",
        feats.len(),
        model
    );

    let mut out = std::fs::File::create(out_path).expect("create out");
    writeln!(out, "mono\tz\trt\tscore2d\tscore3d").unwrap();

    for f in &feats {
        if f.z == 0 {
            continue;
        }
        let mono_mz = mass_to_mz_f64(f.mono, f.z);
        let template = model.intensities_from_mono(f.mono, TEMPLATE_MIN_WEIGHT, template_max_isotopes);
        if template.is_empty() {
            continue;
        }
        let mode_w = template.iter().cloned().fold(0.0_f64, f64::max).max(1e-12);
        let sig: Vec<(usize, f64)> = template
            .iter()
            .enumerate()
            .filter(|(_, &w)| w / mode_w >= MIN_REL)
            .map(|(k, &w)| (k, w))
            .collect();
        if sig.is_empty() {
            continue;
        }
        let kmax = sig.iter().map(|(k, _)| *k).max().unwrap();
        let offsets = tooth_offsets(f.z, template.len(), lattice);
        let step = C13_MINUS_C12 * spacing_scale / f.z.abs() as f64;

        // 2-D: single apex scan, g = 1.
        let apex = nearest_scan(&scans, f.rt);
        let mut acc2 = (0.0, 0.0, 0.0);
        accumulate_scan(&scans[apex], mono_mz, &offsets, &sig, kmax, step, 1.0, TOL_PPM, &mut acc2);
        let score2d = if acc2.1 > 0.0 && acc2.2 > 0.0 { acc2.0 / (acc2.1.sqrt() * acc2.2.sqrt()) } else { 0.0 };

        // 3-D: every scan in the RT window, weighted by the elution Gaussian.
        let mut acc3 = (0.0, 0.0, 0.0);
        let s_lo = scans.partition_point(|s| s.retention_time < f.rt - rt_half_window_min);
        let s_hi = scans.partition_point(|s| s.retention_time <= f.rt + rt_half_window_min);
        for scan in &scans[s_lo..s_hi] {
            let g = gaussian(scan.retention_time - f.rt, rt_sigma_min);
            if g < 1e-4 {
                continue;
            }
            accumulate_scan(scan, mono_mz, &offsets, &sig, kmax, step, g, TOL_PPM, &mut acc3);
        }
        let score3d = if acc3.1 > 0.0 && acc3.2 > 0.0 { acc3.0 / (acc3.1.sqrt() * acc3.2.sqrt()) } else { 0.0 };

        writeln!(out, "{:.5}\t{}\t{:.4}\t{:.5}\t{:.5}", f.mono, f.z, f.rt, score2d, score3d).unwrap();
    }
    eprintln!("wrote scores -> {out_path}");
}
