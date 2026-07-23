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

use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use flashlfq_core::deconvolution::{ClassicDeconvolutionParameters, EnvelopeModel, Polarity};
use flashlfq_core::feature_refinement::{
    four_way_decon, four_way_decon_detector_anchor, four_way_decon_gated, refine_feature,
    refine_feature_censored, refine_feature_multi, refine_feature_shift,
    refine_feature_shift_neighbor_with, refine_feature_shift_with, resolve_charge_state_consensus,
    resolve_consensus_by_apex, DeconView, FourWayDecon, NeighborIndex, RefinedFeature, ResolvedFeature,
};
use flashlfq_core::isotope_shift_decon::{envelope_fit_cosine, Deconvoluter};
use flashlfq_core::isotopic_envelope::{mass_to_mz_f64, C13_MINUS_C12};
use flashlfq_core::joint_fit::{joint_fit_target_shift, Component};
use flashlfq_core::peak_indexing::{
    read_ms1_scans, IndexedMassSpectralPeak, PeakIndexingEngine, PeakKey, Scan,
};
use flashlfq_core::trace_kernel::{
    detect_features, estimate_noise_floor, median_ms1_scan_spacing_minutes, CombWeightModel,
    DetectedFeature, LatticeMode, RtProfile, ScoreModel, TraceKernelParameters, FWHM_TO_SIGMA,
};

/// Neighbour isotope-m/z positions in feature `f`'s window that are **not** on its own grid — the peaks
/// to mask from its shift fit (`NEIGHBOR_REFINE`). Own grid = `mono_mz + k·spacing`, `k ∈ [-1, kmax+2]`.
/// Only neighbours at least `min_ratio ×` this feature's intensity contribute (defer to stronger
/// species only). `idx` may be built from detected OR refined features; querying by `f.apex_rt` avoids
/// any self-index alignment, and `min_ratio > 1` guarantees a feature never masks its own peaks.
fn neighbor_mask_for(idx: &NeighborIndex, f: &DetectedFeature, min_ratio: f64) -> Vec<f64> {
    const GRID_PPM: f64 = 15.0;
    let spacing = C13_MINUS_C12 / f.charge.max(1) as f64;
    let win_min = (f.mono_mz - 1.5).max(0.0);
    let win_max = f.mono_mz + (f.num_isotopes_observed as f64 + 3.0) * spacing + 1.0;
    let kmax = f.num_isotopes_observed as i32 + 2;
    let min_intensity = f.summed_intensity * min_ratio;
    idx.forbidden_positions_query(f.apex_rt, win_min, win_max, min_intensity)
        .into_iter()
        .filter(|&p| {
            !(-1..=kmax).any(|k| {
                let own = f.mono_mz + k as f64 * spacing;
                (p - own).abs() / p * 1e6 <= GRID_PPM
            })
        })
        .collect()
}

/// A growing, RT-bucketed set of *already-refined* feature grids, used by the iterative
/// (`NEIGHBOR_REFINE=iterative`) pass: features are refined in descending score order, and each locks
/// its corrected isotope grid here so subsequent (lower-scoring) features can mask its peaks.
struct LockedGrids {
    /// RT bin (`floor(apex_rt / rt_tol)`) → `(apex_rt, mono_mz, spacing, kmax)`.
    buckets: std::collections::HashMap<i64, Vec<(f64, f64, f64, i32)>>,
    rt_tol: f64,
}

impl LockedGrids {
    fn new(rt_tol: f64) -> Self {
        LockedGrids { buckets: std::collections::HashMap::new(), rt_tol }
    }
    fn bin(&self, rt: f64) -> i64 {
        (rt / self.rt_tol).floor() as i64
    }
    fn add(&mut self, apex_rt: f64, mono_mz: f64, charge: i32, kmax: i32) {
        let spacing = C13_MINUS_C12 / charge.max(1) as f64;
        let b = self.bin(apex_rt);
        self.buckets.entry(b).or_default().push((apex_rt, mono_mz, spacing, kmax));
    }
    /// Ascending isotope m/z of locked features co-eluting with `apex_rt` and in `[win_min, win_max]`.
    fn positions(&self, apex_rt: f64, win_min: f64, win_max: f64) -> Vec<f64> {
        let b = self.bin(apex_rt);
        let mut out = Vec::new();
        for bb in (b - 1)..=(b + 1) {
            let Some(v) = self.buckets.get(&bb) else { continue };
            for &(rt, mono_mz, spacing, kmax) in v {
                if (rt - apex_rt).abs() > self.rt_tol {
                    continue;
                }
                for k in 0..=kmax {
                    let m = mono_mz + k as f64 * spacing;
                    if m < win_min {
                        continue;
                    }
                    if m > win_max {
                        break;
                    }
                    out.push(m);
                }
            }
        }
        out.sort_by(f64::total_cmp);
        out
    }
}

/// Off-own-grid filter shared by the neighbour-mask paths: keep only positions not within `GRID_PPM`
/// of `f`'s own isotope grid (`mono_mz + k·spacing`, `k ∈ [-1, kmax+2]`).
fn filter_off_own_grid(raw: Vec<f64>, f: &DetectedFeature) -> Vec<f64> {
    const GRID_PPM: f64 = 15.0;
    let spacing = C13_MINUS_C12 / f.charge.max(1) as f64;
    let kmax = f.num_isotopes_observed as i32 + 2;
    raw.into_iter()
        .filter(|&p| {
            !(-1..=kmax).any(|k| {
                let own = f.mono_mz + k as f64 * spacing;
                (p - own).abs() / p * 1e6 <= GRID_PPM
            })
        })
        .collect()
}

/// Post-refine **joint linear-model** pass (`JOINT_FIT`): for each refined feature still scoring below
/// `max_score`, model its window as its own averagine envelope PLUS the overlapping co-eluting
/// neighbours and search the target's monoisotope over ¹³C shifts to maximise the joint (NNLS) fit —
/// correcting a residual mis-placement the single-envelope refine could not resolve because a
/// neighbour's peaks were confusing the score. Mutates `refined` in place; returns how many were moved.
fn apply_joint_fit_pass(
    refined: &mut [RefinedFeature],
    scans: &[Scan],
    max_score: f64,
    min_gain: f64,
) -> usize {
    let rt_tol = 0.05;
    let bin_of = |rt: f64| -> i64 { (rt / rt_tol).floor() as i64 };
    // Per-feature: (mono_mz, spacing, charge, apex_rt, intensity, kmax).
    let info: Vec<(f64, f64, i32, f64, f64, i32)> = refined
        .iter()
        .map(|r| {
            let z = r.refined_charge;
            let mono_mz = mass_to_mz_f64(r.refined_monoisotopic_mass, z);
            let spacing = C13_MINUS_C12 / z.max(1) as f64;
            (
                mono_mz,
                spacing,
                z,
                r.detected.apex_rt,
                r.detected.summed_intensity,
                r.detected.num_isotopes_observed as i32 + 2,
            )
        })
        .collect();
    let mut buckets: HashMap<i64, Vec<usize>> = HashMap::new();
    for (i, inf) in info.iter().enumerate() {
        buckets.entry(bin_of(inf.3)).or_default().push(i);
    }

    // Collect updates first (immutable borrow of refined via info), then apply.
    let mut updates: Vec<(usize, f64, f64)> = Vec::new(); // (idx, new_mono_mass, new_score)
    for i in 0..refined.len() {
        if refined[i].decon_score >= max_score {
            continue;
        }
        let (mono_mz, spacing, z, rt, _int, kmax) = info[i];
        let win_lo = mono_mz - 0.6 * spacing;
        let win_hi = mono_mz + (kmax as f64 + 1.0) * spacing;
        // Overlapping co-eluting neighbours (any tooth in the window), strongest first, up to 3.
        let b = bin_of(rt);
        let mut neigh: Vec<usize> = Vec::new();
        for bb in (b - 1)..=(b + 1) {
            let Some(v) = buckets.get(&bb) else { continue };
            for &j in v {
                if j == i || (info[j].3 - rt).abs() > rt_tol {
                    continue;
                }
                let (jm, js, _jz, _, _jint, jk) = info[j];
                let overlaps = (0..=jk).any(|k| {
                    let m = jm + k as f64 * js;
                    m >= win_lo && m <= win_hi
                });
                if overlaps {
                    neigh.push(j);
                }
            }
        }
        if neigh.is_empty() {
            continue;
        }
        neigh.sort_by(|&a, &c| info[c].4.total_cmp(&info[a].4));
        neigh.truncate(3);

        let ai = refined[i].detected.apex_scan_index.max(0) as usize;
        if ai >= scans.len() {
            continue;
        }
        let s = &scans[ai];
        let lo = s.mz.partition_point(|&m| m < win_lo - 0.1);
        let hi = s.mz.partition_point(|&m| m <= win_hi + 0.1);
        if hi <= lo + 1 {
            continue;
        }
        let (mzw, inw) = (&s.mz[lo..hi], &s.intensity[lo..hi]);

        let mut comps = vec![Component { mono_mz, charge: z }];
        for &j in &neigh {
            comps.push(Component { mono_mz: info[j].0, charge: info[j].2 });
        }
        let (shift, _jr) =
            joint_fit_target_shift(mzw, inw, &comps, &[-2, -1, 0, 1, 2], 20.0, 0.2, 0.0);
        if shift != 0 {
            let new_mass = refined[i].refined_monoisotopic_mass + shift as f64 * C13_MINUS_C12;
            let new_mono_mz = mass_to_mz_f64(new_mass, z);
            // Re-score the target alone at the corrected placement (single-envelope, for consensus).
            let new_score = envelope_fit_cosine(mzw, inw, new_mono_mz, z, 20.0, 0.2, 0.0);
            // Only accept a move that improves the target's own fit by a clear margin (guards against
            // over-correcting a feature onto a neighbour's peak for a marginal gain).
            if new_score > refined[i].decon_score + min_gain {
                updates.push((i, new_mass, new_score));
            }
        }
    }

    let n = updates.len();
    for (i, mass, score) in updates {
        refined[i].refined_monoisotopic_mass = mass;
        refined[i].candidate_masses = vec![mass];
        refined[i].decon_score = score;
    }
    n
}

/// Derives a sibling output path from the final path: `out.tsv` + tag `detected` -> `out.detected.tsv`.
fn sibling(out: &str, tag: &str) -> String {
    match out.strip_suffix(".tsv") {
        Some(stem) => format!("{stem}.{tag}.tsv"),
        None => format!("{out}.{tag}.tsv"),
    }
}

/// Opens an output writer, tolerating a locked target (e.g. the file is open in Excel for manual
/// validation): on failure it falls back to `<path>.new` and warns, rather than panicking and
/// discarding the whole run. Returns `None` only if even the fallback cannot be created.
fn open_out(path: &str) -> Option<BufWriter<File>> {
    match File::create(path) {
        Ok(f) => Some(BufWriter::new(f)),
        Err(e) => {
            let alt = format!("{path}.new");
            eprintln!("  WARN: could not write {path} ({e}) — is it open? writing {alt} instead");
            match File::create(&alt) {
                Ok(f) => Some(BufWriter::new(f)),
                Err(e2) => {
                    eprintln!("  WARN: fallback {alt} also failed ({e2}); skipping this file");
                    None
                }
            }
        }
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
    let log_path = match out_path.strip_suffix(".tsv") {
        Some(stem) => format!("{stem}.log"),
        None => format!("{out_path}.log"),
    };
    eprintln!("output files:");
    eprintln!("  detected (pre-refinement): {detected_path}");
    eprintln!("  refined (post-decon):      {refined_path}");
    eprintln!("  resolved (final):          {out_path}");
    eprintln!("  timing log:                {log_path}");

    // Per-step wall-clock timings, reported to the chat and written to the log file at the end.
    let run_start = Instant::now();
    let mut timings: Vec<(String, f64)> = Vec::new();

    // --- read + index --------------------------------------------------------------------------
    let t0 = Instant::now();
    eprintln!("reading MS1 scans from {spectra_path} ...");
    let scans = read_ms1_scans(spectra_path).expect("failed to read spectra file");
    let engine = PeakIndexingEngine::index_peaks(&scans).expect("no indexable MS1 peaks");
    let n_peaks: usize = scans.iter().map(|s| s.mz.len()).sum();
    let total_intensity: f64 = scans.iter().flat_map(|s| s.intensity.iter()).sum();
    let read_dur = t0.elapsed();
    timings.push(("read + index".into(), read_dur.as_secs_f64()));
    eprintln!(
        "  {} MS1 scans, {} peaks, ΣTIC {:.3e}, median scan spacing {:.4} min  ({:.1?})",
        scans.len(),
        n_peaks,
        total_intensity,
        median_ms1_scan_spacing_minutes(engine.scan_info()),
        read_dur
    );

    // --- detect --------------------------------------------------------------------------------
    // Comb-weight model selectable via COMB_MODEL=averagine|poisson (default averagine) for benchmarking.
    // The decoy models (target-decoy FDR) select a non-physical envelope whose detections are noise
    // coincidences — run a separate decoy pass per model and compare score histograms.
    let weight_model = match std::env::var("COMB_MODEL").as_deref() {
        Ok("poisson") => CombWeightModel::Poisson,
        // Chlorinated-averagine decoy (A+2-dominated, mode shifted off the mono).
        Ok("decoy") => CombWeightModel::Decoy,
        // Rotated-averagine decoy — real envelope, weights rotated by half (reversed-peptide analogue).
        Ok("rotated") => CombWeightModel::RotatedAveragine,
        // Chloro-boro-phosphate decoy — averagine + per-unit Cl/B/P (CBP_CL/CBP_B/CBP_P tune the load).
        Ok("cbp") | Ok("chloroboro") => CombWeightModel::ChloroBoroPhosphate,
        // Fully custom decoy — composition REPLACES the backbone, set via CUSTOM_AVERAGINE env
        // (default "P2 C1 N1 Cl1 Fe1" — the weird Fe+Cl averagine).
        Ok("custom") => CombWeightModel::Custom,
        // Shuffled-averagine decoy — real averagine weights randomly permuted (SHUFFLE_SEED varies it).
        Ok("shuffled") => CombWeightModel::ShuffledAveragine,
        // Hybrid — custom (CUSTOM_AVERAGINE) comb below HYBRID_MASS Da, shuffled averagine above.
        Ok("hybrid") => CombWeightModel::Hybrid,
        _ => CombWeightModel::Averagine,
    };
    // End-to-end envelope model for *refinement*, mirroring COMB_MODEL so a decoy is applied through
    // detection AND refinement (shift/recharge/walk-back/score all fit the decoy template, not the real
    // averagine — the fix for score "laundering" where refinement re-anchored decoys onto real peaks).
    let refine_model = match weight_model {
        CombWeightModel::Decoy => EnvelopeModel::ChlorinatedDecoy,
        CombWeightModel::RotatedAveragine => EnvelopeModel::RotatedDecoy,
        CombWeightModel::ChloroBoroPhosphate => EnvelopeModel::ChloroBoroPhosphateDecoy,
        CombWeightModel::Custom => EnvelopeModel::CustomDecoy,
        CombWeightModel::ShuffledAveragine => EnvelopeModel::ShuffledDecoy,
        CombWeightModel::Hybrid => EnvelopeModel::HybridDecoy,
        _ => EnvelopeModel::Averagine,
    };
    // Decoy comb *lattice* (tooth positions), independent of the weight model. DECOY_LATTICE=scaled with
    // DECOY_SPACING_SCALE=0.9368 lays a 0.94-Da off-lattice comb (teeth off the real isotope positions);
    // DECOY_LATTICE=mixed uses per-step mixed-charge spacings (globally impossible comb).
    let lattice_mode = match std::env::var("DECOY_LATTICE").as_deref() {
        Ok("mixed") => LatticeMode::MixedCharge,
        Ok("scaled") => {
            let s = std::env::var("DECOY_SPACING_SCALE")
                .ok()
                .and_then(|s| s.parse::<f64>().ok())
                .unwrap_or(0.5);
            LatticeMode::Scaled(s)
        }
        _ => LatticeMode::Uniform,
    };
    // Refine-stage deconvoluter: the decoy envelope model + (for a scaled lattice) the same spacing
    // scale, so detection and refinement share both the weight model AND the tooth positions. A mixed
    // lattice has no single scale, so refinement stays on the physical lattice there (weights still decoy).
    let refine_spacing_scale = match lattice_mode {
        LatticeMode::Scaled(s) => s,
        _ => 1.0,
    };
    let env_decon = Deconvoluter::new_with_spacing(refine_model, refine_spacing_scale);
    // Decoy elution profile (RT weighting). RT_PROFILE=uniform (flat) or inverted (1−gaussian, U-shaped)
    // for a target-decoy on the chromatographic axis; default gaussian (the real elution template).
    let rt_profile = match std::env::var("RT_PROFILE").as_deref() {
        Ok("uniform") => RtProfile::Uniform,
        Ok("inverted") => RtProfile::InvertedGaussian,
        _ => RtProfile::Gaussian,
    };
    // Coverage target (fraction of ΣTIC to explain) selectable via COVERAGE_TARGET=0.80|0.90|0.99…
    // Default 1.0 (uncapped): the detector runs to the seed-intensity floor. This is also what keeps the
    // default 2-D tiling path engaged — any coverage target < 1.0 needs a global running ΣTIC and forces
    // the serial fallback (see detect_features). Set COVERAGE_TARGET < 1 to cap coverage on the serial path.
    let coverage_target = std::env::var("COVERAGE_TARGET")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.0);
    // Assumed chromatographic FWHM (seconds) that sets the RT Gaussian σ and matched-filter window,
    // selectable via ASSUMED_FWHM_SEC (default 36). Real CA/Lumos peaks are ~3-20 s wide, so smaller
    // values narrow the window to the data and reduce broad-hypothesis interference.
    let assumed_fwhm_sec = std::env::var("ASSUMED_FWHM_SEC")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(36.0);
    // Claim-extent trace knobs (Change A): missed-scan tolerance and the RT half-width guard for the
    // XIC that follows a real elution beyond the ~2σ scoring window. TRACE_MAX_HALF_WIDTH_SEC is in
    // seconds; default 30 s (0.5 min).
    let trace_missed = std::env::var("TRACE_MISSED_SCANS")
        .ok()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(1);
    let trace_half_width_min = std::env::var("TRACE_MAX_HALF_WIDTH_SEC")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .map(|s| s / 60.0)
        .unwrap_or(0.5);
    // Chromatographic-persistence gate (Change A): reject features whose traced extent spans fewer
    // than this many distinct scans. Default 2 drops single-scan noise doublets; MIN_FEATURE_SCANS=1
    // disables it.
    let min_feature_scans = std::env::var("MIN_FEATURE_SCANS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(2);
    // Seed intensity floor (bounds detection cost). Lower it to reach higher coverage (the default
    // 1000 exhausts seeds ~90% ΣTIC on CA/Lumos before the coverage target is hit).
    let min_seed_intensity = std::env::var("MIN_SEED_INTENSITY")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1000.0);
    // Noise-decoy m/z shift (NOISE_DECOY_SHIFT=<frac>, default 0 = real detector). A non-zero fraction
    // (e.g. 0.5 = half-tooth) rigidly offsets the whole comb off the seed so every tooth samples noise
    // between real isotopes — generates a faithful LOW-intensity junk null for target-decoy training.
    let decoy_mz_shift_frac = std::env::var("NOISE_DECOY_SHIFT")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0);
    if decoy_mz_shift_frac != 0.0 {
        eprintln!("NOISE_DECOY_SHIFT: comb shifted {decoy_mz_shift_frac} x (¹³C/z) off the seed — noise-decoy pass");
    }
    // Score model (Change B): SCORE_MODEL=normalized selects the noise-floor-truncated normalised
    // correlation (default raw sum). NOISE_PCT is the percentile of peak intensity used as η (default 5).
    let score_model = match std::env::var("SCORE_MODEL").as_deref() {
        Ok("normalized") | Ok("normalised") => ScoreModel::NormalizedNoiseFloor,
        _ => ScoreModel::RawSum,
    };
    let noise_pct = std::env::var("NOISE_PCT")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(5.0);
    let noise_floor = if matches!(score_model, ScoreModel::NormalizedNoiseFloor) {
        estimate_noise_floor(&engine, noise_pct)
    } else {
        0.0
    };
    // Change B variant knobs: AMP_SEED=1 uses the seed intensity as the apex amplitude A (vs the
    // least-squares fit); SCORE_COSINE=1 divides additionally by ‖I‖ for a bounded cosine shape-fit.
    let score_use_seed_amplitude = matches!(std::env::var("AMP_SEED").as_deref(), Ok("1") | Ok("true"));
    let score_cosine = matches!(std::env::var("SCORE_COSINE").as_deref(), Ok("1") | Ok("true"));
    // Auto-stop at the coverage knee (opt-in; default OFF). DETECT_KNEE=1 enables it; the marginal
    // %ΣTIC-per-seed collapse is a speed/recall knob (sacrifices the low-abundance tail), not a
    // correctness fix. DETECT_KNEE_FRAC = stop-below fraction of the early-window slope (default
    // 0.02); DETECT_KNEE_WINDOW = rolling window width in seeds (default 20000); DETECT_KNEE_EPS =
    // absolute flat-tail slope floor (fraction of ΣTIC per seed, default 1e-7).
    let knee_stop_enabled = matches!(std::env::var("DETECT_KNEE").as_deref(), Ok("1") | Ok("true"));
    let knee_slope_frac = std::env::var("DETECT_KNEE_FRAC")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.02);
    let knee_window_seeds = std::env::var("DETECT_KNEE_WINDOW")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(20_000);
    let knee_abs_eps = std::env::var("DETECT_KNEE_EPS")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1e-7);
    // Auto-stop at the seed-rejection-rate threshold (opt-in; default OFF — the detector ships UNCAPPED).
    // DETECT_REJECT_STOP=1 enables it; stops once the rolling fraction of considered seeds being rejected
    // reaches DETECT_REJECT_FRAC (default 0.80) over the reject window. On the parallel paths the window is
    // sized from data per tile (DETECT_TILE2D_REJECT_WINDOW_FRAC, default 5%); on serial it is
    // DETECT_REJECT_WINDOW scored seeds (default 20000). A 2026-07-08 A/B (see trace_kernel
    // reject_stop_enabled docs) found the cap not worth the recall cost, so it stays off by default; the
    // defaults here are the best setting found, for callers who opt in for detect speed on large files.
    let reject_stop_enabled =
        matches!(std::env::var("DETECT_REJECT_STOP").as_deref(), Ok("1") | Ok("true"));
    let reject_stop_frac = std::env::var("DETECT_REJECT_FRAC")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.80);
    let reject_stop_window = std::env::var("DETECT_REJECT_WINDOW")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(20_000);
    // --- top-down preset -----------------------------------------------------------------------
    // TOPDOWN=1 flips the bottom-up defaults to a top-down-proteomics configuration: high charge
    // range (proteoforms ionize to z≈60), a long isotope comb (a ~20 kDa envelope spans ~40-50
    // significant peaks, so the 12-tooth bottom-up comb truncates below the envelope apex), a
    // slightly higher observed-isotope floor (rich envelopes make ≥3 teeth a cheap false-harmonic
    // filter), and a wider trace half-width (larger species elute broader). Every value below is a
    // *default* — the individual env vars (MIN_CHARGE/MAX_CHARGE/MAX_ISOTOPES/MIN_ISOTOPES_OBS/
    // TRACE_MAX_HALF_WIDTH_SEC) still override it, so the preset is a starting point for the sweep.
    let topdown = matches!(std::env::var("TOPDOWN").as_deref(), Ok("1") | Ok("true"));
    // Joint multi-charge (charge-ladder) detection — default ON for TOPDOWN, off otherwise;
    // DETECT_MULTICHARGE=0/1 overrides. Caps the charge range at 30 by default (per the design: a
    // proteoform's ladder is scored jointly, and z>30 is rare / mostly harmonic noise).
    let multicharge = match std::env::var("DETECT_MULTICHARGE").as_deref() {
        Ok("1") | Ok("true") => true,
        Ok("0") | Ok("false") => false,
        _ => topdown,
    };
    let min_charge = std::env::var("MIN_CHARGE")
        .ok()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(1);
    let max_charge = std::env::var("MAX_CHARGE")
        .ok()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(if multicharge { 30 } else if topdown { 60 } else { 6 });
    let min_charge_states = std::env::var("MIN_CHARGE_STATES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(2);
    // Joint cross-charge mono-offset search half-width (¹³C units) on the multi-charge path. Default 3;
    // MC_MONO_KMAX=0 disables it (mono stays at the averagine anchor) for A/B.
    let multicharge_mono_kmax = std::env::var("MC_MONO_KMAX")
        .ok()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(3);
    let max_isotopes = std::env::var("MAX_ISOTOPES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(if topdown { 60 } else { 12 });
    let min_isotopes_observed = std::env::var("MIN_ISOTOPES_OBS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(if topdown { 3 } else { 2 });
    // Top-down species elute broader; widen the trace half-width guard unless the caller set it.
    let trace_half_width_min = if std::env::var("TRACE_MAX_HALF_WIDTH_SEC").is_err() && topdown {
        1.0
    } else {
        trace_half_width_min
    };
    if topdown {
        eprintln!(
            "TOPDOWN preset: charge {min_charge}..={max_charge}, max_isotopes {max_isotopes}, \
             min_isotopes_observed {min_isotopes_observed}, trace half-width {:.0} s",
            trace_half_width_min * 60.0
        );
    }
    // Recharge (envelope-fit charge re-selection in refinement) must consider the full top-down
    // charge range, else a detected high-z feature yields an empty candidate set and is dropped.
    flashlfq_core::feature_refinement::set_recharge_max_charge(max_charge);
    // Top-down monoisotope-offset fit (averagine envelope fit over ±k ¹³C). Defaults ON for the
    // TOPDOWN preset (k=3), off otherwise; TD_MONO_FIT=<k> overrides (TD_MONO_FIT=0 disables).
    let td_mono_fit_kmax = std::env::var("TD_MONO_FIT")
        .ok()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(if topdown { 3 } else { 0 });
    flashlfq_core::feature_refinement::set_td_mono_fit_kmax(td_mono_fit_kmax);
    if td_mono_fit_kmax > 0 {
        eprintln!("TD mono-offset fit: ENABLED (±{td_mono_fit_kmax} ¹³C averagine fit, mass ≥ 3 kDa)");
    }
    // Cross-charge off-by-one bridge width in the consensus. Bottom-up default 2. TOPDOWN default 0
    // (DISABLED): the bridge merges features ~1 Da apart as one proteoform's off-by-one and votes one
    // mass, but top-down is dense with genuinely-close species — notably deamidation (+0.984 Da) sits
    // only 0.019 Da from a +1.003 Da isotope step (~1 ppm at 15 kDa) — so bridging conflates distinct
    // proteoforms onto a wrong mass. A/B on Jurkat: bridge 4→58%, 2→61%, 1→65%, 0→75% strict recall.
    let offbyone_units = std::env::var("OFFBYONE_UNITS")
        .ok()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(if topdown { 0 } else { 2 });
    flashlfq_core::feature_refinement::set_offbyone_units(offbyone_units);
    let base = TraceKernelParameters {
        ppm_tolerance: 10.0,
        min_charge,
        max_charge,
        max_isotopes,
        min_isotopes_observed,
        min_seed_intensity,
        coverage_target,
        weight_model,
        lattice_mode,
        rt_profile,
        score_model,
        noise_floor,
        score_use_seed_amplitude,
        score_cosine,
        trace_missed_scans_allowed: trace_missed,
        trace_max_half_width_minutes: trace_half_width_min,
        min_feature_scans,
        knee_stop_enabled,
        knee_window_seeds,
        knee_slope_frac,
        knee_abs_eps,
        reject_stop_enabled,
        reject_stop_window,
        reject_stop_frac,
        multicharge_enabled: multicharge,
        min_charge_states,
        multicharge_mono_kmax,
        decoy_mz_shift_frac,
        ..TraceKernelParameters::default()
    };
    if multicharge {
        eprintln!(
            "multi-charge detection: ENABLED (joint charge-ladder over z {min_charge}..={max_charge}, \
             ≥{min_charge_states} charge states, mono-offset search ±{multicharge_mono_kmax} ¹³C; serial — no tiling)"
        );
    }
    if knee_stop_enabled {
        eprintln!(
            "knee auto-stop: ENABLED (window {} seeds, slope_frac {}, abs_eps {:.0e}) — stops early at the coverage knee",
            knee_window_seeds, knee_slope_frac, knee_abs_eps
        );
    }
    if reject_stop_enabled {
        eprintln!(
            "reject-rate auto-stop: ENABLED (frac {:.2}; serial rolling window {} scored seeds; the \
             parallel paths size each tile's window from data via DETECT_TILE2D_REJECT_WINDOW_FRAC) — \
             stops when that fraction of considered seeds is being rejected",
            reject_stop_frac, reject_stop_window
        );
    }
    eprintln!(
        "score model: {:?}  (η = {:.1} @ p{:.0}, A = {}, {})",
        score_model,
        noise_floor,
        noise_pct,
        if score_use_seed_amplitude { "seed" } else { "least-squares" },
        if score_cosine { "cosine" } else { "template-norm" }
    );
    // σ from the data by default (Change A: measure FWHM via XIC half-max, ASSUMED_FWHM_SEC is the
    // fallback). FIXED_SIGMA=1 forces the assumed-FWHM path for A/B comparison.
    let data_driven_sigma = std::env::var("FIXED_SIGMA").is_err();
    let params = if data_driven_sigma {
        base.with_rt_from_index(&engine, assumed_fwhm_sec)
    } else {
        base.with_rt_from_scans(engine.scan_info(), assumed_fwhm_sec)
    };
    if data_driven_sigma {
        eprintln!(
            "σ source: data-driven (measured FWHM via XIC half-max; fallback {assumed_fwhm_sec} s) \
             → {:.2} s FWHM",
            params.rt_sigma_minutes * FWHM_TO_SIGMA * 60.0
        );
    } else {
        eprintln!("σ source: assumed FWHM {assumed_fwhm_sec} s (FIXED_SIGMA)");
    }
    eprintln!("comb weight model: {weight_model:?}");
    eprintln!(
        "detecting (charge {}..={}, {} ppm, σ_rt {:.4} min, ±{} scans, trace: {} missed / {:.0} s half-width, min {} scans, seed floor {:.0}, coverage {:.0}%) ...",
        params.min_charge,
        params.max_charge,
        params.ppm_tolerance,
        params.rt_sigma_minutes,
        params.half_window_scans,
        params.trace_missed_scans_allowed,
        params.trace_max_half_width_minutes * 60.0,
        params.min_feature_scans,
        params.min_seed_intensity,
        params.coverage_target * 100.0
    );
    // LOAD_DETECTED=<path> skips the (expensive) detect stage and loads a compact detected-feature
    // cache written by a prior run's SAVE_DETECTED — for fast iteration on the refine/resolve stages
    // (which are ~10% of runtime) without paying the ~5-min detect each time. The cache keeps only the
    // scalar fields plus the single most-intense claimed peak (the anchor), which is all the default
    // shift-apex refine path reads from `peaks`; see build_feature_slices.
    let t1 = Instant::now();
    let mut detected = if let Ok(cache) = std::env::var("LOAD_DETECTED") {
        let d = load_detected_cache(&cache);
        eprintln!("  loaded {} detected features from cache (detect SKIPPED) <- {cache}", d.len());
        d
    } else {
        detect_features(&engine, &params)
    };
    let detect_dur = t1.elapsed();
    timings.push(("detect".into(), detect_dur.as_secs_f64()));
    let detected_intensity: f64 = detected.iter().map(|f| f.summed_intensity).sum();
    eprintln!(
        "  {} features detected, explained {:.1}% of ΣTIC  ({:.1?})",
        detected.len(),
        100.0 * detected_intensity / total_intensity,
        detect_dur
    );
    if let Ok(cache) = std::env::var("SAVE_DETECTED") {
        save_detected_cache(&cache, &detected);
        eprintln!("  saved detected-feature cache ({} features) -> {cache}", detected.len());
    }
    write_detected_tsv(&detected_path, &detected);
    eprintln!("  wrote {} detected features -> {detected_path}", detected.len());

    // MZ_STRADDLE=1: measure the m/z-tiling straddle exposure for the planned 2-D parallel detector.
    // For each detected feature, emit its claimed isotope envelope's m/z span and the largest tooth on
    // each side of the apex (the most-intense claimed peak). A 2-D m/z tile border that falls inside an
    // envelope can orphan a minor tooth into a neighbouring tile, where it may seed a mis-anchored
    // feature — so the span distribution (vs the tile width, ≥ 2·reach_mz) bounds how often re-anchoring
    // is needed. Off-apex tooth intensity vs the seed floor says whether an orphaned tooth could even
    // seed. Analyse the emitted TSV in Python.
    if std::env::var("MZ_STRADDLE").is_ok() {
        let reach_mz =
            params.max_isotopes as f64 * (C13_MINUS_C12 / params.min_charge.max(1) as f64);
        let spath = sibling(out_path, "straddle");
        if let Some(mut w) = open_out(&spath) {
            writeln!(
                w,
                "apex_mz\tapex_int\tcharge\tn_peaks\tmz_min\tmz_max\tspan\t\
                 left_max_int\tleft_dist\tright_max_int\tright_dist"
            )
            .unwrap();
            for f in &detected {
                if f.peaks.is_empty() {
                    continue;
                }
                let apex = f
                    .peaks
                    .iter()
                    .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
                    .unwrap();
                let amz = apex.m() as f64;
                let (mut mzmin, mut mzmax) = (f64::INFINITY, f64::NEG_INFINITY);
                let (mut lmax, mut ldist, mut rmax, mut rdist) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
                for p in &f.peaks {
                    let (pm, pi) = (p.m() as f64, p.intensity as f64);
                    mzmin = mzmin.min(pm);
                    mzmax = mzmax.max(pm);
                    if pm < amz - 1e-6 {
                        if pi > lmax {
                            lmax = pi;
                            ldist = amz - pm;
                        }
                    } else if pm > amz + 1e-6 && pi > rmax {
                        rmax = pi;
                        rdist = pm - amz;
                    }
                }
                writeln!(
                    w,
                    "{:.5}\t{:.1}\t{}\t{}\t{:.5}\t{:.5}\t{:.5}\t{:.1}\t{:.5}\t{:.1}\t{:.5}",
                    amz, apex.intensity, f.charge, f.peaks.len(), mzmin, mzmax, mzmax - mzmin,
                    lmax, ldist, rmax, rdist
                )
                .unwrap();
            }
            let _ = w.flush();
        }
        eprintln!(
            "  MZ_STRADDLE: reach_mz={:.2} Da, seed floor={:.0}; wrote per-feature envelope spans -> {spath}",
            reach_mz, params.min_seed_intensity
        );
        return;
    }

    // Escape hatch for diagnostics: skip the (potentially intractable at huge feature counts)
    // refine + O(n^2) charge-consensus and stop after detection.
    if std::env::var("DETECT_ONLY").is_ok() {
        eprintln!("DETECT_ONLY set — skipping refine/resolve/compare after {} detected features.", detected.len());
        return;
    }

    // --- IsoDec charge re-assignment (opt-in) --------------------------------------------------
    // ISODEC_CHARGE=1: replace each feature's detector-assigned charge with the native IsoDec neural
    // predictor's call, run on the raw apex-scan peaks in IsoDec's local m/z window. The detector
    // over-calls charge on weak features (spurious high z); IsoDec's charge is more reliable, and since
    // our mono mass is derived from the apex mass *at the charge*, a corrected charge corrects the mass.
    // Only overrides when IsoDec is confident (>= ISODEC_MINPEAKS window peaks, non-zero call). Mono m/z
    // and mass are recomputed at the new charge from the feature's most-intense (anchor) peak.
    // Default ON for TOPDOWN (validated: +10pp intersection strict recall, off-by-one gap 15.7→5.8pp);
    // ISODEC_CHARGE=0 disables. Needs a full detect (uses the features' claimed peaks) and ~4× runtime.
    let isodec_charge = match std::env::var("ISODEC_CHARGE").as_deref() {
        Ok("1") | Ok("true") => true,
        Ok("0") | Ok("false") => false,
        _ => topdown && std::env::var("LOAD_DETECTED").is_err(),
    };
    if isodec_charge {
        use flashlfq_core::deconvolution::averagine_mono_from_most_intense;
        use flashlfq_core::isodec::default_model;
        let min_peaks = std::env::var("ISODEC_MINPEAKS")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(6);
        let model = default_model();
        // Feed IsoDec the FWHM-averaged composite over the footprint instead of the single apex scan.
        // Validated: the higher SNR lets IsoDec nail the charge on more features — intersection strict
        // 94.2%→98.2%, off-by-one gap 5.8→1.8pp. Default ON; ISODEC_AVERAGE=0 uses the apex scan.
        let isodec_average = !matches!(std::env::var("ISODEC_AVERAGE").as_deref(), Ok("0") | Ok("false"));
        let iso_half = {
            let fwhm_s = params.rt_sigma_minutes * FWHM_TO_SIGMA * 60.0;
            let sp_s = median_ms1_scan_spacing_minutes(engine.scan_info()) * 60.0;
            (flashlfq_core::feature_refinement::derived_avg_scans(fwhm_s, sp_s) / 2).max(1)
        };
        let iso_avg_params = flashlfq_core::spectral_averaging::SpectralAveragingParameters::default();
        if isodec_average {
            eprintln!("  ISODEC_CHARGE: averaging {}-scan composite per feature (apex ±{iso_half})", iso_half * 2 + 1);
        }
        let ti = Instant::now();
        let mut changed = 0usize;
        for f in detected.iter_mut() {
            let anchor = match f.peaks.iter().max_by(|a, b| a.intensity.total_cmp(&b.intensity)) {
                Some(p) => *p,
                None => continue,
            };
            let anchor_mz = anchor.mz as f64;
            // Feed IsoDec the RAW apex-scan peaks in this feature's own m/z footprint (the span of its
            // claimed teeth, padded) — i.e. the claimed teeth PLUS the unclaimed intervening peaks.
            // Claimed teeth alone are blind to an under-called low harmonic (a real z=10 called z=5 has
            // its extra peaks at the half-spacing offsets, never claimed); the footprint's raw peaks
            // include them so IsoDec can see the true higher charge. Bounded to the footprint so it does
            // not grab co-eluting neighbours the way a fixed wide window does. Needs full detect (the
            // LOAD_DETECTED cache keeps only the anchor peak, so the footprint would be a point).
            let (mut mzmin, mut mzmax) = (f64::INFINITY, f64::NEG_INFINITY);
            for p in &f.peaks {
                let m = p.mz as f64;
                mzmin = mzmin.min(m);
                mzmax = mzmax.max(m);
            }
            if !mzmin.is_finite() {
                continue;
            }
            let apex = (f.apex_scan_index.max(0) as usize).min(scans.len().saturating_sub(1));
            let (lo_m, hi_m) = (mzmin - 0.1, mzmax + 0.1);
            let (wmz, wint): (Vec<f64>, Vec<f32>) = if isodec_average {
                // Composite over apex ± iso_half scans in the footprint (binned, TIC-normalised).
                let lo_s = apex.saturating_sub(iso_half);
                let hi_s = (apex + iso_half).min(scans.len().saturating_sub(1));
                let mut xs: Vec<Vec<f64>> = Vec::new();
                let mut ys: Vec<Vec<f64>> = Vec::new();
                for ss in &scans[lo_s..=hi_s] {
                    let a = ss.mz.partition_point(|&m| m < lo_m);
                    let b = ss.mz.partition_point(|&m| m <= hi_m);
                    xs.push(ss.mz[a..b].to_vec());
                    ys.push(ss.intensity[a..b].to_vec());
                }
                let (cmz, cint) =
                    flashlfq_core::spectral_averaging::average_spectra(&xs, &ys, &iso_avg_params);
                (cmz, cint.iter().map(|&v| v as f32).collect())
            } else {
                let s = &scans[apex];
                let lo = s.mz.partition_point(|&m| m < lo_m);
                let hi = s.mz.partition_point(|&m| m <= hi_m);
                (
                    s.mz[lo..hi].to_vec(),
                    s.intensity[lo..hi].iter().map(|&v| v as f32).collect(),
                )
            };
            if wmz.len() < min_peaks {
                continue;
            }
            let z = model.predict_charge(&wmz, &wint);
            if z >= 1 && z <= params.max_charge && z != f.charge {
                let apex_mass = anchor_mz * z as f64 - z as f64 * flashlfq_core::isotopic_envelope::PROTON_MASS;
                let mono = averagine_mono_from_most_intense(apex_mass);
                f.charge = z;
                f.monoisotopic_mass = mono;
                f.mono_mz = mass_to_mz_f64(mono, z);
                changed += 1;
            }
        }
        eprintln!(
            "  ISODEC_CHARGE: re-assigned {changed}/{} features' charge via IsoDec (min {min_peaks} window peaks)  ({:.1?})",
            detected.len(),
            ti.elapsed()
        );
    }

    // --- refine --------------------------------------------------------------------------------
    let mut avg = flashlfq_core::spectral_averaging::SpectralAveragingParameters::default();
    // Derive the composite averaging window from the run's measured FWHM (recovered from the σ the
    // detector set) and MS1 scan spacing, instead of the fixed apex±1 floor. Only affects methods
    // that actually build a composite (shift_composite / classic); the default apex-only path
    // (shift_apex, average_spectra=false) never reads it, so this leaves the default output unchanged.
    {
        let fwhm_seconds = params.rt_sigma_minutes * FWHM_TO_SIGMA * 60.0;
        let spacing_seconds = median_ms1_scan_spacing_minutes(engine.scan_info()) * 60.0;
        avg.avg_scans =
            flashlfq_core::feature_refinement::derived_avg_scans(fwhm_seconds, spacing_seconds);
        eprintln!(
            "  refine averaging window: {} scans (FWHM {:.2} s / spacing {:.3} s)",
            avg.avg_scans, fwhm_seconds, spacing_seconds
        );
    }
    let decon = ClassicDeconvolutionParameters::new(
        params.min_charge,
        params.max_charge,
        10.0,
        3.0,
        Polarity::Positive,
    );
    // CENSOR_CLAIMED=1 removes peaks claimed by OTHER (stronger, already-assigned) features from
    // each feature's composite before deconvolution — a subtractive-decon experiment. Since the
    // detector claims greedily tallest-first, this hands weaker features a window with co-eluting
    // interferents removed.
    let censor_claimed = std::env::var("CENSOR_CLAIMED").is_ok();
    let all_claimed: HashSet<PeakKey> = if censor_claimed {
        detected
            .iter()
            .flat_map(|f| f.peaks.iter().map(|p| p.key()))
            .collect()
    } else {
        HashSet::new()
    };
    if censor_claimed {
        eprintln!(
            "  CENSOR_CLAIMED: subtracting {} claimed peaks from other features' decon windows",
            all_claimed.len()
        );
    }

    // REFINE_METHOD selects the deconvolution used to place the refined monoisotope:
    //   classic (default) | shift_composite | shift_apex  (detector-anchored FlashLFQ-style shift).
    // Default: the detector-anchored shift decon on the apex scan (REFINE_METHOD=classic to opt out
    // back to the parity-locked classic deconvolution; shift_composite selects the averaged composite).
    // TOPDOWN defaults to the multi-envelope refine (best top-down config); bottom-up stays shift_apex.
    let refine_method = std::env::var("REFINE_METHOD")
        .unwrap_or_else(|_| if topdown { "multi".into() } else { "shift_apex".into() });
    let use_shift_apex = refine_method == "shift_apex";
    let use_shift = use_shift_apex || refine_method == "shift_composite";
    // Multi-envelope refine: model the composite window as a combination of co-eluting averagine
    // envelopes (one per grid-local-maximum apex) fit jointly by NNLS. TD_MULTI_COMPONENTS caps the
    // number of envelopes (default 4).
    let use_multi = refine_method == "multi";
    let multi_components = std::env::var("TD_MULTI_COMPONENTS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(4);
    if use_multi {
        eprintln!("  refine method: multi-envelope joint fit (≤{multi_components} co-eluting averagines per window)");
    }
    // Apex-only refinement is the default, so the averaged composite is built ONLY when a method
    // actually reads it (shift_composite). On CA/Lumos 10-min the apex scan beats the averaged
    // composite: 97.4% vs 96.0% recall, 94.4% vs 92.3% charge accuracy — so we skip the composite
    // (and its per-scan slicing) in the default apex path. Revisit once data-dependent averaging
    // lands for longer gradients, where averaging may pay off again.
    let average_spectra = !use_shift_apex;
    // Charge re-selection by the envelope-fit cosine (fit + explained + completeness) defaults ON for
    // the shift methods; disable with RECHARGE=0. Recovers charge-halved features (the light-z2 class).
    let recharge = std::env::var("RECHARGE")
        .map(|v| v != "0" && !v.eq_ignore_ascii_case("false"))
        .unwrap_or(true);
    eprintln!("  refine method: {refine_method}{}", if recharge && use_shift { " + recharge (envelope-fit cosine)" } else { "" });

    // NEIGHBOR_REFINE=1: mask co-eluting neighbours' peaks (off this feature's own grid) from the
    // charge-selection fit, walk-back and score — so a low-scoring feature in a crowded window is judged
    // on the signal plausibly its own. Experiment path; builds a NeighborIndex over the detections.
    // NEIGHBOR_REFINE: mask co-eluting neighbours' peaks from the shift fit. Value picks the neighbour
    // context: "detected"/"1" = raw detections; "refined" = a two-pass build (refine once, then mask
    // against those corrected placements). NEIGHBOR_MIN_RATIO (default 5.0) = only defer to neighbours
    // at least that many times more intense — a weak feature in a strong neighbour's shadow.
    let neighbor_mode = std::env::var("NEIGHBOR_REFINE").ok();
    // "iterative": refine in descending-score order, each feature locking its corrected grid so later
    // (lower-scoring) features mask its peaks. The confident features claim their signal first.
    let iterative = use_shift && neighbor_mode.as_deref() == Some("iterative");
    let neighbor_refine = use_shift && !iterative && neighbor_mode.is_some();
    let neighbor_refined_context = neighbor_mode.as_deref() == Some("refined");
    let neighbor_min_ratio = std::env::var("NEIGHBOR_MIN_RATIO")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(5.0);
    let neighbor_idx = if neighbor_refine {
        let src = if neighbor_refined_context {
            // Pass 1: plain refine, then index the corrected placements for the masked pass below.
            let pass1: Vec<RefinedFeature> = detected
                .iter()
                .filter_map(|f| refine_feature_shift_with(&env_decon, f, &scans, &avg, 20.0, use_shift_apex, recharge, average_spectra))
                .collect();
            eprintln!(
                "  NEIGHBOR_REFINE=refined: two-pass, masking neighbours >= {neighbor_min_ratio}x (from {} refined)",
                pass1.len()
            );
            NeighborIndex::build_from_refined(&pass1, 0.05)
        } else {
            eprintln!(
                "  NEIGHBOR_REFINE: masking peaks of co-eluting detected neighbours >= {neighbor_min_ratio}x this feature's intensity"
            );
            NeighborIndex::build(&detected, 0.05)
        };
        Some(src)
    } else {
        None
    };

    let t2 = Instant::now();
    let mut refined: Vec<RefinedFeature> = Vec::with_capacity(detected.len());
    let progress_every = 1000usize;
    if iterative {
        // Only features whose fit clears this bar lock their grid into the mask — "higher-scoring"
        // is not enough; a mediocre-but-higher feature is still uncertain and would add collateral.
        let lock_min_score = std::env::var("NEIGHBOR_LOCK_MIN_SCORE")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.9);
        eprintln!("  NEIGHBOR_REFINE=iterative: score-ordered; only features with fit >= {lock_min_score} lock their grids");
        // Pass 1: plain refine to get each feature's initial score.
        let init: Vec<Option<RefinedFeature>> = detected
            .iter()
            .map(|f| refine_feature_shift_with(&env_decon, f, &scans, &avg, 20.0, use_shift_apex, recharge, average_spectra))
            .collect();
        // Process indices in descending initial score.
        let mut order: Vec<usize> = (0..detected.len()).filter(|&i| init[i].is_some()).collect();
        order.sort_by(|&a, &b| {
            init[b].as_ref().unwrap().decon_score.total_cmp(&init[a].as_ref().unwrap().decon_score)
        });
        let mut locked = LockedGrids::new(0.05);
        let mut out: Vec<Option<RefinedFeature>> = vec![None; detected.len()];
        let mut n_locked = 0usize;
        for (n, &i) in order.iter().enumerate() {
            let f = &detected[i];
            let spacing = C13_MINUS_C12 / f.charge.max(1) as f64;
            let win_min = (f.mono_mz - 1.5).max(0.0);
            let win_max = f.mono_mz + (f.num_isotopes_observed as f64 + 3.0) * spacing + 1.0;
            let mask = filter_off_own_grid(locked.positions(f.apex_rt, win_min, win_max), f);
            let r = refine_feature_shift_neighbor_with(&env_decon, f, &scans, &avg, 20.0, use_shift_apex, recharge, &mask, average_spectra)
                .or_else(|| init[i].clone());
            if let Some(rr) = &r {
                // Only confident features become mask sources for the lower-scoring ones that follow.
                if rr.decon_score >= lock_min_score {
                    let mono_mz = mass_to_mz_f64(rr.refined_monoisotopic_mass, rr.refined_charge);
                    let kmax = f.num_isotopes_observed as i32 + 2;
                    locked.add(rr.detected.apex_rt, mono_mz, rr.refined_charge, kmax);
                    n_locked += 1;
                }
            }
            out[i] = r;
            if (n + 1) % progress_every == 0 || n + 1 == order.len() {
                eprintln!("    iterative refined {}/{} ({n_locked} locked)  [{:?}]", n + 1, order.len(), t2.elapsed());
            }
        }
        refined = out.into_iter().flatten().collect();
    } else {
        for (i, f) in detected.iter().enumerate() {
            let r = if use_multi {
                refine_feature_multi(f, &scans, &avg, 20.0, multi_components)
            } else if let Some(idx) = &neighbor_idx {
                let mask = neighbor_mask_for(idx, f, neighbor_min_ratio);
                refine_feature_shift_neighbor_with(&env_decon, f, &scans, &avg, 20.0, use_shift_apex, recharge, &mask, average_spectra)
            } else if use_shift {
                refine_feature_shift_with(&env_decon, f, &scans, &avg, 20.0, use_shift_apex, recharge, average_spectra)
            } else if censor_claimed {
                refine_feature_censored(f, &scans, &avg, &decon, &all_claimed)
            } else {
                refine_feature(f, &scans, &avg, &decon)
            };
            if let Some(r) = r {
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
    }
    let refine_dur = t2.elapsed();
    timings.push(("refine".into(), refine_dur.as_secs_f64()));
    eprintln!(
        "  {} / {} features refined against averaged composites  ({:.1?})",
        refined.len(),
        detected.len(),
        refine_dur
    );

    // JOINT_FIT: post-refine joint linear-model pass over low-scoring features (see apply_joint_fit_pass).
    if use_shift && std::env::var("JOINT_FIT").is_ok() {
        let max_score = std::env::var("JOINT_FIT_MAX_SCORE")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.7);
        let min_gain = std::env::var("JOINT_FIT_MIN_GAIN")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.02);
        let tj = Instant::now();
        let moved = apply_joint_fit_pass(&mut refined, &scans, max_score, min_gain);
        timings.push(("joint fit".into(), tj.elapsed().as_secs_f64()));
        eprintln!(
            "  JOINT_FIT: joint linear-model pass moved {moved} low-scoring features (< {max_score})  ({:.1?})",
            tj.elapsed()
        );
    }

    write_refined_tsv(&refined_path, &refined);
    eprintln!("  wrote {} refined features -> {refined_path}", refined.len());

    // --- four-way decon comparator (opt-in diagnostic) -----------------------------------------
    // FOUR_WAY_DECON=1 runs the classic×shift on composite×apex comparator over every detected
    // feature, reporting how often the four monoisotope views disagree (→ candidates for advanced
    // multi-envelope decon) and writing a per-feature disagreements TSV. Gated because it roughly
    // doubles the decon cost; off by default so normal runs are unaffected.
    if std::env::var("FOUR_WAY_DECON").is_ok() {
        let t_fw = Instant::now();
        let fw_path = sibling(out_path, "disagreements");
        run_four_way(&detected, &scans, &avg, &decon, &fw_path);
        timings.push(("four-way decon".into(), t_fw.elapsed().as_secs_f64()));
    }

    // DIFF_REFINE=1 exports the features CLASSIC refine drops (no envelope) but the detector-anchored
    // SHIFT refine keeps — the completeness gain behind the +2.4% recall. Feeds make_gain_targets.py.
    if std::env::var("DIFF_REFINE").is_ok() {
        let dpath = sibling(out_path, "refinediff");
        run_refine_diff(&detected, &scans, &avg, &decon, &dpath);
    }

    // --- resolve charge-state consensus --------------------------------------------------------
    // TD_CONSENSUS=apex groups cross-charge by the robust apex (most-abundant) neutral mass and resolves
    // each proteoform's monoisotope once from the consensus apex (top-down); default is the mono-keyed
    // consensus. RT window for apex grouping is wider (co-eluting charge states share an apex RT).
    // TOPDOWN defaults to apex-mass cross-charge consensus; TD_CONSENSUS=mono forces the mono-keyed path.
    let consensus_by_apex = match std::env::var("TD_CONSENSUS").as_deref() {
        Ok("apex") => true,
        Ok("mono") => false,
        _ => topdown,
    };
    // Apex grouping can tolerate an integer-¹³C apex-isotope drift between charge states (TD_APEX_SHIFT).
    // Default 0 (tight): A/B showed shift=1 REGRESSES badly (intersection strict 84.2%→48.0%) — allowing
    // the drift over-merges genuinely-close distinct species (deamidation +0.984, off-by-N proteoforms),
    // and that cost dominates the occasional benefit of grouping a drifting charge state.
    let apex_shift = std::env::var("TD_APEX_SHIFT")
        .ok()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(0);
    let t3 = Instant::now();
    let resolved = if consensus_by_apex {
        eprintln!("  consensus: apex-mass cross-charge grouping (TD_CONSENSUS=apex, ±{apex_shift} ¹³C apex shift)");
        resolve_consensus_by_apex(&refined, 15.0, 0.3, apex_shift)
    } else {
        resolve_charge_state_consensus(&refined, 10.0, 0.1)
    };
    let consensus_dur = t3.elapsed();
    timings.push(("charge-state consensus".into(), consensus_dur.as_secs_f64()));
    eprintln!(
        "  {} peptide-level features after charge-state consensus  ({:.1?})",
        resolved.len(),
        consensus_dur
    );

    let t4 = Instant::now();
    write_tsv(out_path, &resolved);
    let write_dur = t4.elapsed();
    timings.push(("write resolved".into(), write_dur.as_secs_f64()));
    eprintln!("wrote {} -> {}", resolved.len(), out_path);

    // MSALIGN_OUT=1: also emit the resolved features as an MS1 TopFD/msDeconv `.msalign`
    // (readable by mzLib's Ms1Align reader). Opt-in; the TSV outputs above are unchanged.
    if std::env::var("MSALIGN_OUT").is_ok() {
        let msalign_path = match out_path.strip_suffix(".tsv") {
            Some(stem) => format!("{stem}.ms1.msalign"),
            None => format!("{out_path}.ms1.msalign"),
        };
        let src = std::path::Path::new(spectra_path)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(spectra_path);
        let feats = flashlfq_core::feature_export::resolved_to_ms1align_features(&resolved);
        match flashlfq_core::feature_export::write_ms1_align_file(&msalign_path, &feats, src) {
            Ok(()) => eprintln!(
                "  MSALIGN_OUT: wrote {} MS1 features -> {msalign_path}",
                feats.len()
            ),
            Err(e) => eprintln!("  WARN: could not write msalign {msalign_path}: {e}"),
        }
    }

    if let Some(ref_path) = reference_path {
        let t5 = Instant::now();
        compare_to_reference(ref_path, &resolved);
        timings.push(("compare".into(), t5.elapsed().as_secs_f64()));
    }

    // --- timing summary (chat + log file) ------------------------------------------------------
    let total = run_start.elapsed().as_secs_f64();
    let mut lines: Vec<String> = Vec::new();
    lines.push("=== timing summary ===".into());
    lines.push(format!("spectra file: {spectra_path}"));
    lines.push(format!(
        "{} MS1 scans, {} peaks, {} detected, {} refined, {} resolved features",
        scans.len(),
        n_peaks,
        detected.len(),
        refined.len(),
        resolved.len()
    ));
    for (name, secs) in &timings {
        lines.push(format!("  {name:<24} {secs:8.2} s  ({:4.1}%)", 100.0 * secs / total));
    }
    lines.push(format!("  {:<24} {total:8.2} s", "TOTAL"));

    for l in &lines {
        eprintln!("{l}");
    }
    // Append to the log file so repeated runs accumulate a history (never blocks the run).
    match std::fs::OpenOptions::new().create(true).append(true).open(&log_path) {
        Ok(mut f) => {
            for l in &lines {
                let _ = writeln!(f, "{l}");
            }
            let _ = writeln!(f);
            eprintln!("timing log appended to {log_path}");
        }
        Err(e) => eprintln!("  WARN: could not write timing log {log_path} ({e})"),
    }
}

/// Exports features that classic `refine_feature` drops (returns `None`) but detector-anchored
/// `refine_feature_shift` (apex) keeps — the completeness gain. Columns feed `make_gain_targets.py`
/// (which joins to the reference): shift-refined mono, apex RT, charge, most-abundant peak m/z.
fn run_refine_diff(
    detected: &[DetectedFeature],
    scans: &[flashlfq_core::peak_indexing::Scan],
    avg: &flashlfq_core::spectral_averaging::SpectralAveragingParameters,
    decon: &ClassicDeconvolutionParameters,
    path: &str,
) {
    let mut w = match open_out(path) {
        Some(w) => w,
        None => return,
    };
    writeln!(w, "Mono\tRT\tCharge\tPk MZ\tSummed Intensity").unwrap();
    let mut n_gain = 0usize;
    for f in detected {
        // Classic kept it → not a gain.
        if refine_feature(f, scans, avg, decon).is_some() {
            continue;
        }
        // Apex refine (default): use_apex = true, so no averaged composite is built.
        if let Some(r) = refine_feature_shift(f, scans, avg, 20.0, true, false, false) {
            n_gain += 1;
            let pk = f
                .peaks
                .iter()
                .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
                .map(|p| p.m() as f64)
                .unwrap_or(f.mono_mz);
            writeln!(
                w,
                "{:.5}\t{:.4}\t{}\t{:.5}\t{:.4e}",
                r.refined_monoisotopic_mass, f.apex_rt, f.charge, pk, f.summed_intensity
            )
            .unwrap();
        }
    }
    let _ = w.flush();
    eprintln!(
        "  DIFF_REFINE: {n_gain} features shift-kept but classic-dropped -> {path}"
    );
}

/// Runs the four-way decon comparator over every detected feature, reports the disagreement rate to
/// the chat, and writes a per-disagreeing-feature TSV (the candidates for advanced multi-envelope
/// decon). `decon` uses the pipeline's classic parameters; the shift views use a 20 ppm match tol.
fn run_four_way(
    detected: &[DetectedFeature],
    scans: &[flashlfq_core::peak_indexing::Scan],
    avg: &flashlfq_core::spectral_averaging::SpectralAveragingParameters,
    decon: &ClassicDeconvolutionParameters,
    path: &str,
) {
    let short = |v: DeconView| match v {
        DeconView::ClassicComposite => "cc",
        DeconView::ClassicApex => "ca",
        DeconView::ShiftComposite => "sc",
        DeconView::ShiftApex => "sa",
    };

    let mut w = open_out(path);
    if let Some(w) = w.as_mut() {
        writeln!(
            w,
            "Detector Mono\tCharge\tApex RT\tConsensus k\tNum Verdicts\tViews (view:k:conf)"
        )
        .unwrap();
    }

    // Full per-feature verdict table: each view's monoisotopic mass (blank if the view produced
    // nothing). Lets the within-method (composite-vs-apex) and ground-truth analysis run in Python.
    let vpath = path.replace(".disagreements.tsv", ".verdicts.tsv");
    let mut vw = open_out(&vpath);
    if let Some(vw) = vw.as_mut() {
        writeln!(vw, "Detector Mono\tCharge\tApex RT\tcc_mono\tca_mono\tsc_mono\tsa_mono").unwrap();
    }
    let mono_of = |fw: &FourWayDecon, view: DeconView| -> Option<f64> {
        fw.verdicts
            .iter()
            .find(|v| v.view == view)
            .map(|v| v.monoisotopic_mass)
    };
    let fmt = |m: Option<f64>| m.map(|x| format!("{x:.5}")).unwrap_or_default();
    // Within-method composite-vs-apex agreement (same integer ¹³C offset), counted where both views exist.
    let mut classic_agree = 0usize;
    let mut classic_disagree = 0usize;
    let mut shift_agree = 0usize;
    let mut shift_disagree = 0usize;
    let same_k = |a: f64, b: f64| ((a - b) / C13_MINUS_C12).round() == 0.0;

    // NEIGHBOR_AWARE=1 gates the shift-decon anchor: peaks belonging to a co-eluting already-detected
    // feature (its predicted isotope grid) are excluded as anchor candidates, so shift-decon can't
    // "distant grab" a stronger neighbour several isotopes away. Builds an RT-bucketed neighbour index.
    // DETECTOR_ANCHOR=1 takes precedence: anchor the shift views on the detector's own most-abundant
    // claimed peak (cannot grab any foreign peak). NEIGHBOR_AWARE=1 gates against co-eluting neighbours.
    let detector_anchor = std::env::var("DETECTOR_ANCHOR").is_ok();
    let neighbor_aware = !detector_anchor && std::env::var("NEIGHBOR_AWARE").is_ok();
    let neighbors = if neighbor_aware {
        eprintln!("  NEIGHBOR_AWARE: gating shift anchors against co-eluting detected features");
        Some(NeighborIndex::build(detected, 0.05))
    } else {
        None
    };
    if detector_anchor {
        eprintln!("  DETECTOR_ANCHOR: shift views anchor on the detector's most-abundant claimed peak");
    }
    const GRID_PPM: f64 = 15.0;

    let mut total = 0usize; // features that produced >=1 verdict
    let mut with_verdicts_hist = [0usize; 5]; // count by number of verdicts (0..=4)
    let mut unanimous = 0usize;
    let mut disagreed = 0usize;
    for (i, f) in detected.iter().enumerate() {
        let fw: FourWayDecon = if detector_anchor {
            four_way_decon_detector_anchor(f, scans, avg, decon, 20.0)
        } else {
            match &neighbors {
                Some(idx) => {
                    // Same shift window the comparator uses, to gather the neighbours' forbidden teeth.
                    let spacing = C13_MINUS_C12 / f.charge.max(1) as f64;
                    let win_min = (f.mono_mz - 1.5).max(0.0);
                    let win_max = f.mono_mz + (f.num_isotopes_observed as f64 + 3.0) * spacing + 1.0;
                    let forbidden = idx.forbidden_positions(i, win_min, win_max);
                    four_way_decon_gated(f, scans, avg, decon, 20.0, &forbidden, GRID_PPM)
                }
                None => four_way_decon(f, scans, avg, decon, 20.0),
            }
        };
        with_verdicts_hist[fw.verdicts.len().min(4)] += 1;
        if fw.verdicts.is_empty() {
            continue;
        }
        total += 1;

        // Per-feature verdict row + within-method (composite vs apex) agreement.
        let (cc, ca, sc, sa) = (
            mono_of(&fw, DeconView::ClassicComposite),
            mono_of(&fw, DeconView::ClassicApex),
            mono_of(&fw, DeconView::ShiftComposite),
            mono_of(&fw, DeconView::ShiftApex),
        );
        if let Some(vw) = vw.as_mut() {
            writeln!(
                vw,
                "{:.5}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
                fw.detector_mono,
                fw.charge,
                f.apex_rt,
                fmt(cc),
                fmt(ca),
                fmt(sc),
                fmt(sa)
            )
            .unwrap();
        }
        if let (Some(a), Some(b)) = (cc, ca) {
            if same_k(a, b) {
                classic_agree += 1;
            } else {
                classic_disagree += 1;
            }
        }
        if let (Some(a), Some(b)) = (sc, sa) {
            if same_k(a, b) {
                shift_agree += 1;
            } else {
                shift_disagree += 1;
            }
        }
        if fw.needs_advanced {
            disagreed += 1;
            if let Some(w) = w.as_mut() {
                let views: Vec<String> = fw
                    .verdicts
                    .iter()
                    .map(|v| format!("{}:{:+}:{}", short(v.view), v.offset_k, if v.confident { 1 } else { 0 }))
                    .collect();
                writeln!(
                    w,
                    "{:.5}\t{}\t{:.4}\t{}\t{}\t{}",
                    fw.detector_mono,
                    fw.charge,
                    f.apex_rt,
                    fw.consensus_k.map(|k| k.to_string()).unwrap_or_default(),
                    fw.verdicts.len(),
                    views.join(" ")
                )
                .unwrap();
            }
        } else {
            unanimous += 1;
        }
    }
    if let Some(mut w) = w {
        let _ = w.flush();
    }
    if let Some(mut vw) = vw {
        let _ = vw.flush();
        eprintln!("  wrote per-feature verdicts -> {vpath}");
    }

    let cpair = classic_agree + classic_disagree;
    let spair = shift_agree + shift_disagree;
    eprintln!("\n=== within-method composite-vs-apex agreement ===");
    eprintln!(
        "  classic (cc vs ca): {} agree / {} disagree  ({:.1}% agree of {} with both)",
        classic_agree, classic_disagree, 100.0 * classic_agree as f64 / cpair.max(1) as f64, cpair
    );
    eprintln!(
        "  shift   (sc vs sa): {} agree / {} disagree  ({:.1}% agree of {} with both)",
        shift_agree, shift_disagree, 100.0 * shift_agree as f64 / spair.max(1) as f64, spair
    );

    eprintln!("\n=== four-way decon comparator ===");
    eprintln!("  detected features: {}", detected.len());
    eprintln!(
        "  produced >=1 verdict: {} (verdict-count histogram [0..4]: {:?})",
        total, with_verdicts_hist
    );
    eprintln!(
        "  unanimous (confident placement): {} ({:.1}%)",
        unanimous,
        100.0 * unanimous as f64 / total.max(1) as f64
    );
    eprintln!(
        "  disagreed (→ advanced multi-envelope): {} ({:.1}%)  → {path}",
        disagreed,
        100.0 * disagreed as f64 / total.max(1) as f64
    );
}

/// Writes the raw detected features (pre-refinement, straight from the trace kernel) to a TSV,
/// sorted by summed intensity descending. Lets you inspect what the detector alone produced.
fn write_detected_tsv(path: &str, detected: &[DetectedFeature]) {
    let mut rows: Vec<&DetectedFeature> = detected.iter().collect();
    rows.sort_by(|a, b| b.summed_intensity.total_cmp(&a.summed_intensity));
    let mut w = match open_out(path) {
        Some(w) => w,
        None => return,
    };
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

/// Writes a **compact detected-feature cache** for fast refine/resolve iteration (see `LOAD_DETECTED`).
/// Keeps every scalar field plus only the single most-intense claimed peak (the anchor) — all the
/// default shift-apex refine path reads from `peaks`. Tab-separated, one feature per line.
fn save_detected_cache(path: &str, detected: &[DetectedFeature]) {
    let mut w = match open_out(path) {
        Some(w) => w,
        None => return,
    };
    writeln!(
        w,
        "mono_mass\tcharge\tmono_mz\tapex_scan\tapex_rt\tstart_rt\tend_rt\tsummed_int\tscore\t\
         num_iso\tanchor_mz\tanchor_int\tanchor_scan\tanchor_rt"
    )
    .unwrap();
    for d in detected {
        // Anchor = most-intense claimed peak (what refine uses); fall back to the mono m/z if empty.
        let anchor = d
            .peaks
            .iter()
            .max_by(|a, b| a.intensity.total_cmp(&b.intensity));
        let (amz, aint, ascan, art) = match anchor {
            Some(p) => (p.mz, p.intensity, p.zero_based_scan_index, p.retention_time),
            None => (d.mono_mz as f32, 0.0, d.apex_scan_index, d.apex_rt as f32),
        };
        writeln!(
            w,
            "{:.6}\t{}\t{:.6}\t{}\t{:.5}\t{:.5}\t{:.5}\t{:.6e}\t{:.6e}\t{}\t{:.6}\t{:.6e}\t{}\t{:.5}",
            d.monoisotopic_mass,
            d.charge,
            d.mono_mz,
            d.apex_scan_index,
            d.apex_rt,
            d.start_rt,
            d.end_rt,
            d.summed_intensity,
            d.score,
            d.num_isotopes_observed,
            amz,
            aint,
            ascan,
            art
        )
        .unwrap();
    }
    w.flush().unwrap();
}

/// Loads a compact detected-feature cache written by [`save_detected_cache`]. Reconstructs each
/// `DetectedFeature` with `peaks` = the single anchor peak (sufficient for the shift-apex refine path).
fn load_detected_cache(path: &str) -> Vec<DetectedFeature> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read detected cache {path}: {e}"));
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if i == 0 {
            continue; // header
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 14 {
            continue;
        }
        let anchor = IndexedMassSpectralPeak {
            mz: f[10].parse().unwrap_or(0.0),
            intensity: f[11].parse().unwrap_or(0.0),
            zero_based_scan_index: f[12].parse().unwrap_or(0),
            retention_time: f[13].parse().unwrap_or(0.0),
        };
        out.push(DetectedFeature {
            monoisotopic_mass: f[0].parse().unwrap_or(0.0),
            charge: f[1].parse().unwrap_or(1),
            mono_mz: f[2].parse().unwrap_or(0.0),
            apex_scan_index: f[3].parse().unwrap_or(0),
            apex_rt: f[4].parse().unwrap_or(0.0),
            start_rt: f[5].parse().unwrap_or(0.0),
            end_rt: f[6].parse().unwrap_or(0.0),
            summed_intensity: f[7].parse().unwrap_or(0.0),
            score: f[8].parse().unwrap_or(0.0),
            num_isotopes_observed: f[9].parse().unwrap_or(0),
            peaks: vec![anchor],
        });
    }
    out
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
    let mut w = match open_out(path) {
        Some(w) => w,
        None => return,
    };
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

/// Per-peak **ppm-error spread** across a member's apex-scan isotope envelope — an intensity-orthogonal
/// fit-quality score. For the member's detected feature (charge `z`, mono comb `mono_mz`), take the
/// peaks in the apex scan, assign each to its nearest isotope tooth `mono_mz + k·(¹³C/z)`, keep the
/// tallest peak per tooth, and return the population standard deviation (ppm) of those per-tooth mass
/// errors. A real isotope envelope shares one small calibration offset across teeth → low spread; a
/// noise coincidence's teeth scatter → high spread. Returns `None` for fewer than 2 usable teeth
/// (single-isotope features can't be scored this way — the caller treats that as "worst").
fn apex_ppm_spread(m: &RefinedFeature) -> Option<f64> {
    let f = &m.detected;
    let z = f.charge.max(1);
    let spacing = C13_MINUS_C12 / z as f64;
    if spacing <= 0.0 {
        return None;
    }
    let apex = f.apex_scan_index;
    // isotope tooth index k -> (tallest intensity so far, its ppm error)
    let mut best: HashMap<i64, (f64, f64)> = HashMap::new();
    for p in &f.peaks {
        if p.zero_based_scan_index != apex {
            continue;
        }
        let mz = p.mz as f64;
        let k = ((mz - f.mono_mz) / spacing).round();
        let expected = f.mono_mz + k * spacing;
        if expected <= 0.0 {
            continue;
        }
        let ppm = (mz - expected) / expected * 1e6;
        // Guard against a stray peak that binned to a tooth it isn't really on (the detector claimed
        // within 10 ppm; 20 ppm leaves margin without admitting garbage into the spread).
        if ppm.abs() > 20.0 {
            continue;
        }
        let e = best.entry(k as i64).or_insert((f64::NEG_INFINITY, 0.0));
        if p.intensity as f64 > e.0 {
            *e = (p.intensity as f64, ppm);
        }
    }
    if best.len() < 2 {
        return None;
    }
    let ppms: Vec<f64> = best.values().map(|v| v.1).collect();
    let mean = ppms.iter().sum::<f64>() / ppms.len() as f64;
    let var = ppms.iter().map(|p| (p - mean).powi(2)).sum::<f64>() / ppms.len() as f64;
    Some(var.sqrt())
}

/// Pearson correlation of two equal-length traces; `None` if fewer than 3 points or either trace is
/// flat (zero variance — a correlation is undefined).
fn pearson(a: &[f64], b: &[f64]) -> Option<f64> {
    let n = a.len();
    if n < 3 {
        return None;
    }
    let (ma, mb) = (a.iter().sum::<f64>() / n as f64, b.iter().sum::<f64>() / n as f64);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for i in 0..n {
        let (da, db) = (a[i] - ma, b[i] - mb);
        sab += da * db;
        saa += da * da;
        sbb += db * db;
    }
    if saa <= 0.0 || sbb <= 0.0 {
        return None;
    }
    Some(sab / (saa.sqrt() * sbb.sqrt()))
}

/// **Isotopologue co-elution correlation** — the mean pairwise Pearson correlation of the feature's
/// per-isotope XICs (intensity vs scan) over its elution window. A real feature's isotope traces rise
/// and fall together (they are one chromatographic peak sampled at ¹³C-spaced m/z), so this is near 1;
/// a noise coincidence's teeth vary independently, so it is low. Intensity-orthogonal and a cornerstone
/// signal in Dinosaur/MaxQuant-style detectors.
///
/// `top_n` restricts to the `top_n` most intense isotope teeth (0 = all): the strongest teeth have the
/// most reliable traces, so top-3/top-5 trade coverage for per-trace SNR. Returns `-2.0` (out-of-range
/// sentinel) when it cannot be computed: fewer than 2 usable teeth or fewer than 3 shared scans.
fn isotope_corr(m: &RefinedFeature, top_n: usize) -> f64 {
    const SENTINEL: f64 = -2.0;
    let f = &m.detected;
    let spacing = C13_MINUS_C12 / f.charge.max(1) as f64;
    if spacing <= 0.0 {
        return SENTINEL;
    }
    // isotope tooth index k -> (scan index -> summed intensity) = that tooth's XIC.
    let mut teeth: HashMap<i64, HashMap<i32, f64>> = HashMap::new();
    for p in &f.peaks {
        let k = (((p.mz as f64 - f.mono_mz) / spacing).round()) as i64;
        *teeth
            .entry(k)
            .or_default()
            .entry(p.zero_based_scan_index)
            .or_insert(0.0) += p.intensity as f64;
    }
    if teeth.len() < 2 {
        return SENTINEL;
    }
    // Rank teeth by total intensity; keep the top_n (0 = all).
    let mut ranked: Vec<(f64, &HashMap<i32, f64>)> =
        teeth.values().map(|xic| (xic.values().sum::<f64>(), xic)).collect();
    ranked.sort_by(|a, b| b.0.total_cmp(&a.0));
    let keep = if top_n == 0 { ranked.len() } else { top_n.min(ranked.len()) };
    let sel: Vec<&HashMap<i32, f64>> = ranked[..keep].iter().map(|(_, x)| *x).collect();
    if sel.len() < 2 {
        return SENTINEL;
    }
    // Shared scan axis (union), each trace zero-filled where a tooth had no peak in that scan.
    let scans: Vec<i32> = sel
        .iter()
        .flat_map(|x| x.keys().copied())
        .collect::<BTreeSet<i32>>()
        .into_iter()
        .collect();
    if scans.len() < 3 {
        return SENTINEL;
    }
    let traces: Vec<Vec<f64>> = sel
        .iter()
        .map(|x| scans.iter().map(|s| *x.get(s).unwrap_or(&0.0)).collect())
        .collect();
    let (mut sum, mut n) = (0.0, 0usize);
    for i in 0..traces.len() {
        for j in (i + 1)..traces.len() {
            if let Some(r) = pearson(&traces[i], &traces[j]) {
                sum += r;
                n += 1;
            }
        }
    }
    if n == 0 {
        SENTINEL
    } else {
        sum / n as f64
    }
}

/// Writes the resolved features to a human-readable TSV, sorted by summed intensity descending.
fn write_tsv(path: &str, resolved: &[ResolvedFeature]) {
    let mut rows: Vec<&ResolvedFeature> = resolved.iter().collect();
    rows.sort_by(|a, b| b.summed_intensity.total_cmp(&a.summed_intensity));

    let mut w = match open_out(path) {
        Some(w) => w,
        None => return,
    };
    // Columns lead with the observed-feature answer — the detected m/z, the RT extent, and every
    // charge state seen — then follow with the derived mass/intensity fields.
    //   `Detected m/z (primary)` is the tallest observed isotope-peak m/z of the tallest charge
    //     (the detector seed) — the direct analogue of base FlashLFQ's observed `Peak MZ`, which for
    //     heavier peptides sits ~1 ¹³C step above the monoisotope.
    //   `Per-Charge Detected m/z` lists that same observed m/z for EVERY detected charge, as
    //     `z<charge>:<m/z>` pairs (ascending by charge) — so a peptide seen at z2 and z3 shows both.
    //   `Mono m/z (primary)` is the monoisotopic-peak m/z at the primary charge (derived from the
    //     consensus neutral mass), for reference alongside the observed detected m/z.
    writeln!(
        w,
        "Detected m/z (primary)\tRT Start\tRT Apex\tRT End\tCharge States\tPer-Charge Detected m/z\t\
         Num Charge States\tPrimary Charge\tMonoisotopic Mass\tMono m/z (primary)\tSummed Intensity\t\
         Cross-Charge Support\tNum Members\tDecon Score\tMin Decon Score\tMax Num Isotopes\tPPM Spread\t\
         IsoCorr All\tIsoCorr Top5\tIsoCorr Top3"
    )
    .unwrap();
    // Most-abundant observed isotope-peak m/z of one member's detection (its tallest claimed peak).
    let member_detected_mz = |m: &RefinedFeature| -> Option<f64> {
        m.detected
            .peaks
            .iter()
            .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
            .map(|p| p.m() as f64)
    };
    for r in rows {
        // Primary member = the tallest member; its charge and its seed (tallest) peak drive the m/z.
        let primary_member = r.members.iter().max_by(|a, b| {
            a.detected
                .summed_intensity
                .total_cmp(&b.detected.summed_intensity)
        });
        let primary_charge = primary_member.map(|m| m.refined_charge).unwrap_or(0);
        let mono_mz = if primary_charge != 0 {
            mass_to_mz_f64(r.monoisotopic_mass, primary_charge)
        } else {
            0.0
        };
        let detected_mz = primary_member.and_then(member_detected_mz).unwrap_or(0.0);
        // Observed detected m/z per charge state: for each detected charge, the tallest member of
        // that charge and its most-abundant peak m/z, ascending by charge (`z2:497.2584;z3:331.8416`).
        let per_charge_mz: Vec<String> = r
            .charge_states
            .iter()
            .map(|&z| {
                let mz = r
                    .members
                    .iter()
                    .filter(|m| m.refined_charge == z)
                    .max_by(|a, b| {
                        a.detected
                            .summed_intensity
                            .total_cmp(&b.detected.summed_intensity)
                    })
                    .and_then(member_detected_mz)
                    .unwrap_or(0.0);
                format!("z{z}:{mz:.4}")
            })
            .collect();
        let charges: Vec<String> = r.charge_states.iter().map(|c| c.to_string()).collect();
        // Intensity-orthogonal fit-quality scores (see apex_ppm_spread). Decon score = the primary
        // member's envelope-fit cosine; min decon score = the weakest member's (a multi-charge feature
        // is only as trustworthy as its worst-fitting charge). Max num isotopes = envelope completeness.
        // PPM spread uses a 999 sentinel when it can't be computed (single-isotope apex) so an ascending
        // "low=better" ranking pushes those un-scorable features to the bottom, as intended.
        let decon_score = primary_member.map(|m| m.decon_score).unwrap_or(0.0);
        let min_decon = r
            .members
            .iter()
            .map(|m| m.decon_score)
            .fold(f64::INFINITY, f64::min);
        let min_decon = if min_decon.is_finite() { min_decon } else { 0.0 };
        let max_isotopes = r
            .members
            .iter()
            .map(|m| m.detected.num_isotopes_observed)
            .max()
            .unwrap_or(0);
        let ppm_spread = primary_member.and_then(apex_ppm_spread).unwrap_or(999.0);
        // Isotopologue co-elution correlation, three isotope-selection variants (see isotope_corr).
        let iso_all = primary_member.map(|m| isotope_corr(m, 0)).unwrap_or(-2.0);
        let iso_top5 = primary_member.map(|m| isotope_corr(m, 5)).unwrap_or(-2.0);
        let iso_top3 = primary_member.map(|m| isotope_corr(m, 3)).unwrap_or(-2.0);
        writeln!(
            w,
            "{:.5}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{:.5}\t{:.5}\t{:.4e}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{:.4}\t\
             {:.4}\t{:.4}\t{:.4}",
            detected_mz,
            r.start_rt,
            r.apex_rt,
            r.end_rt,
            charges.join(";"),
            per_charge_mz.join(";"),
            r.charge_states.len(),
            primary_charge,
            r.monoisotopic_mass,
            mono_mz,
            r.summed_intensity,
            r.cross_charge_support,
            r.members.len(),
            decon_score,
            min_decon,
            max_isotopes,
            ppm_spread,
            iso_all,
            iso_top5,
            iso_top3
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
