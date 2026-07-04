//! Trace-kernel feature detection — the untargeted MS1 detector.
//!
//! This is the **new** algorithm (no mzLib counterpart, so no C# golden): it locates
//! "peptide-shaped objects" directly in the raw MS1 data via a sparse 2D matched filter — an
//! **isotope comb** in the m/z dimension × a **Gaussian** in the retention-time dimension — scored
//! against the [`crate::peak_indexing::PeakIndexingEngine`]. See
//! `agent_info/Feature-Detection-Design.md` ("Trace kernel (detection)") for the rationale.
//!
//! ## What it does
//! For each seed peak (tallest first, greedy claim-as-you-go), it scores charge hypotheses
//! `z = min..=max`. Each hypothesis lays an isotope comb spaced `(C13 − C12)/z` in m/z and weighted
//! by the expected isotope envelope, times a Gaussian in RT centred on the seed's scan, and sums the
//! observed intensity the comb lands on. The best-scoring `z` wins (**cross-z non-max
//! suppression**); its peaks are claimed so overlapping harmonics (a real z=2 is a subset of z=4/z=6
//! combs) cannot re-fire. Accepted hypotheses become [`DetectedFeature`] records.
//!
//! ## Design decisions realised here
//! - **The Gaussian template *is* the shape test** — there is no separate data-vs-data correlation
//!   gate at detection time. A matched filter degrades gracefully on tailed (real) peaks; tail-aware
//!   integration bounds are `cut_peak`'s job, downstream.
//! - **The comb runs both directions** from the seed: the seed is the *most intense* peak, which for
//!   heavier masses is not the monoisotopic one, so the mono is placed at `seed − i*·spacing/z` where
//!   `i*` is the most-abundant isotope index of the envelope model.
//! - **Comb weights: closed-form Poisson is the default.** The peptide isotope envelope is
//!   `Binomial(n_C, 0.0107) ≈ Poisson(λ)`, `λ ≈ 0.00048·M` — one parameter, no table lookup
//!   ([`poisson_comb_weights`]). An averagine-table alternative is planned to benchmark against; the
//!   [`CombWeightModel`] enum is the seam for it.
//! - **Evaluate sparsely.** Only the comb's expected `(m/z, scan)` points are looked up in the index;
//!   nothing is rasterised.
//!
//! ## Not yet here (follow-ups)
//! - **Detect-then-refine**: averaging the RT window ([`crate::spectral_averaging`]) + a final
//!   [`crate::deconvolution`] pass on the composite. This module is the *detector/assembler*; the
//!   refiner is wired separately.
//! - **Charge-state consensus** across co-eluting z of the same neutral mass (needs the
//!   `DeconEnvelope` candidate-mass list). Grouping here is per-hypothesis, one charge at a time.
//! - **Averagine comb weights** (benchmark alternative to Poisson).

use std::collections::HashSet;

use crate::isotopic_envelope::{C13_MINUS_C12, PROTON_MASS};
use crate::peak_indexing::{IndexedMassSpectralPeak, PeakIndexingEngine, PeakKey, ScanInfo};
use crate::tolerance::PpmTolerance;

/// Poisson rate per dalton for the closed-form comb: `λ ≈ 0.00048·M`. This is (carbons per Da)
/// × (¹³C natural abundance) ≈ `(1 / averagineUnitMass · averageC) · 0.0107`, i.e. how many ¹³C
/// substitutions a peptide of mass `M` carries on average. The Poisson in that count *is* the
/// isotope envelope.
pub const POISSON_LAMBDA_PER_DA: f64 = 0.00048;

/// Full-width-at-half-maximum → Gaussian σ conversion factor: `FWHM = 2·√(2·ln2)·σ ≈ 2.3548·σ`.
pub const FWHM_TO_SIGMA: f64 = 2.354_820_045_030_949;

/// How the isotope comb's per-peak weights are produced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CombWeightModel {
    /// Closed-form `Poisson(λ = 0.00048·M)` over the ¹³C-substitution count. One parameter, no table
    /// lookup; the fast default.
    Poisson,
    /// The real averagine isotope envelope
    /// ([`crate::deconvolution::averagine_comb_weights`]), a table lookup. More accurate than Poisson
    /// near ~1.8 kDa where the envelope mode shifts off the monoisotope — the regime where a Poisson
    /// `i*` can misplace the monoisotope by one ¹³C unit (off-by-one).
    Averagine,
}

/// Parameters governing the trace-kernel detector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TraceKernelParameters {
    /// Lowest charge hypothesis to test (bottom-up default: 1).
    pub min_charge: i32,
    /// Highest charge hypothesis to test (bottom-up default: 6).
    pub max_charge: i32,
    /// m/z match tolerance, ppm, for looking up comb teeth in the index.
    pub ppm_tolerance: f64,
    /// Gaussian σ along retention time, in **minutes** (the same unit as `ScanInfo::retention_time`).
    pub rt_sigma_minutes: f64,
    /// **Legacy / superseded by [`Self::rt_half_window_minutes`].** Formerly bounded the window by a
    /// fixed scan count; retained for API/compat but no longer consumed by the detector (the window
    /// is now a time window). Still set by `with_rt_from_scans` for reference.
    pub half_window_scans: i32,
    /// Which comb-weight model to use.
    pub weight_model: CombWeightModel,
    /// Smallest envelope weight (relative to the tallest = 1) that still contributes a comb tooth.
    pub min_isotope_weight: f64,
    /// Hard cap on the number of comb teeth (isotopes) considered.
    pub max_isotopes: usize,
    /// Minimum number of distinct comb teeth that must be observed for a hypothesis to be accepted.
    /// Two (a doublet) is the floor — a lone peak is not a feature.
    pub min_isotopes_observed: usize,
    /// Seeds with intensity below this are not considered (noise floor). `0.0` disables the floor
    /// and seeds from every peak. Because seeds are visited intensity-descending, this also bounds
    /// runtime: the seed loop stops as soon as it drops below the floor.
    pub min_seed_intensity: f64,
    /// Stop once this fraction of the total MS1 intensity has been explained (claimed by accepted
    /// features). This is the design's primary stopping criterion — "explain most of the big
    /// signal, not every peak". Denominator = Σ of all indexed peak intensities. `1.0` (or more)
    /// disables the cap and detects until seeds are exhausted / fall below `min_seed_intensity`.
    pub coverage_target: f64,
    /// Half-width of the retention-time window (minutes) the matched filter evaluates around the
    /// seed apex. This bounds the window **in time**, not in scan count: DDA interleaves a variable
    /// number of MS2 scans between MS1 scans, so a fixed scan-count window spans wildly different
    /// times — producing over-wide features that over-claim and split one elution into several.
    /// Typically ≈ 2σ (`with_rt_from_scans` sets it there). Supersedes `half_window_scans`.
    pub rt_half_window_minutes: f64,
}

impl Default for TraceKernelParameters {
    /// Bottom-up defaults: charge 1–6, 10 ppm, Poisson weights, ≥2 observed isotopes. The RT σ and
    /// window are left at ~6 s / ±3 scans placeholders — callers should set them from the data via
    /// [`TraceKernelParameters::with_rt_from_scans`].
    fn default() -> Self {
        TraceKernelParameters {
            min_charge: 1,
            max_charge: 6,
            ppm_tolerance: 10.0,
            rt_sigma_minutes: 0.1,
            half_window_scans: 3,
            weight_model: CombWeightModel::Poisson,
            min_isotope_weight: 1e-3,
            max_isotopes: 12,
            min_isotopes_observed: 2,
            min_seed_intensity: 0.0,
            coverage_target: 1.0,
            rt_half_window_minutes: 0.5,
        }
    }
}

impl TraceKernelParameters {
    /// Derives the RT σ and scan half-window from the data, given an assumed chromatographic peak
    /// width. `assumed_fwhm_seconds` is the design's "~36 s peaks" starting assumption; the σ is
    /// `FWHM / 2.3548` and the half-window spans ±2σ in scans, using the median MS1 scan spacing.
    ///
    /// The FWHM (not the full peak width) drives the averaging window so co-eluting neighbours are
    /// not pulled into the composite; the same σ is reused as the detector's RT Gaussian width.
    pub fn with_rt_from_scans(mut self, scan_info: &[ScanInfo], assumed_fwhm_seconds: f64) -> Self {
        let fwhm_minutes = assumed_fwhm_seconds / 60.0;
        let sigma_minutes = fwhm_minutes / FWHM_TO_SIGMA;
        let spacing = median_ms1_scan_spacing_minutes(scan_info).max(f64::MIN_POSITIVE);
        self.rt_sigma_minutes = sigma_minutes;
        self.half_window_scans = ((2.0 * sigma_minutes) / spacing).round().max(1.0) as i32;
        // The matched filter is bounded in *time* (see `rt_half_window_minutes`); ±2σ covers the peak.
        self.rt_half_window_minutes = 2.0 * sigma_minutes;
        self
    }
}

/// A detected untargeted MS1 feature: one charge state's isotope envelope traced across RT.
#[derive(Debug, Clone)]
pub struct DetectedFeature {
    /// Monoisotopic neutral mass inferred from the mono comb position and charge.
    pub monoisotopic_mass: f64,
    /// Charge state (the winning hypothesis).
    pub charge: i32,
    /// m/z of the monoisotopic comb tooth (`seed − i*·spacing/z`).
    pub mono_mz: f64,
    /// Zero-based scan index of the feature's most intense claimed peak.
    pub apex_scan_index: i32,
    /// Retention time of the apex, minutes.
    pub apex_rt: f64,
    /// Earliest RT among claimed peaks, minutes.
    pub start_rt: f64,
    /// Latest RT among claimed peaks, minutes.
    pub end_rt: f64,
    /// Sum of claimed peak intensities (the feature's explained signal).
    pub summed_intensity: f64,
    /// Matched-filter response — the detector's score for this feature.
    pub score: f64,
    /// Number of distinct isotope teeth that contributed at least one observed peak.
    pub num_isotopes_observed: usize,
    /// The peaks this feature claims (across the RT window and all matched isotopes).
    pub peaks: Vec<IndexedMassSpectralPeak>,
}

/// Neutral mass from an m/z at a given charge: `|z|·mz − z·ProtonMass` (f64). Matches
/// `ClassExtensions.ToMass` / [`crate::isotopic_envelope::mz_to_mass_f32`] but keeps f64 precision.
#[inline]
fn mz_to_mass(mz: f64, charge: i32) -> f64 {
    (charge.abs() as f64) * mz - (charge as f64) * PROTON_MASS
}

/// Closed-form Poisson comb weights for a peptide of neutral mass `neutral_mass`.
///
/// Returns the isotope envelope `w[k] = e^{-λ} λ^k / k!` (`λ = 0.00048·neutral_mass`) as a
/// probability-mass vector, truncated once a weight falls below `min_weight × w_max` (past the mode)
/// or `max_isotopes` teeth are reached. The vector is **normalised so its maximum weight is 1.0**,
/// which makes the seed (the tallest observed peak) align naturally with the tallest template tooth.
pub fn poisson_comb_weights(neutral_mass: f64, min_weight: f64, max_isotopes: usize) -> Vec<f64> {
    let lambda = POISSON_LAMBDA_PER_DA * neutral_mass.max(0.0);
    let mut weights: Vec<f64> = Vec::new();
    // w[0] = e^{-λ}; w[k] = w[k-1]·λ/k. Build up to max_isotopes, tracking the max for normalisation.
    let mut w = (-lambda).exp();
    let mut max_w = w;
    weights.push(w);
    for k in 1..max_isotopes {
        w = w * lambda / (k as f64);
        weights.push(w);
        if w > max_w {
            max_w = w;
        }
        // Stop once we are past the mode (weights descending) and below the relative floor.
        if w < min_weight * max_w && w < weights[k - 1] {
            break;
        }
    }
    if max_w > 0.0 {
        for wk in weights.iter_mut() {
            *wk /= max_w;
        }
    }
    weights
}

/// Index of the most-abundant (tallest) tooth in a weight vector. Ties resolve to the lower index.
fn most_abundant_index(weights: &[f64]) -> usize {
    let mut best = 0;
    for (i, &w) in weights.iter().enumerate() {
        if w > weights[best] {
            best = i;
        }
    }
    best
}

/// Gaussian value `exp(-½ (Δ/σ)²)`.
#[inline]
fn gaussian(delta: f64, sigma: f64) -> f64 {
    let z = delta / sigma;
    (-0.5 * z * z).exp()
}

/// Median spacing between consecutive MS1 scan retention times (minutes). Returns 0.0 for < 2 scans.
///
/// `scan_info` is assumed ordered by scan (as the index builds it). Uses the standard median
/// convention (average of the two middle order statistics for an even count).
pub fn median_ms1_scan_spacing_minutes(scan_info: &[ScanInfo]) -> f64 {
    if scan_info.len() < 2 {
        return 0.0;
    }
    let mut diffs: Vec<f64> = scan_info
        .windows(2)
        .map(|w| w[1].retention_time - w[0].retention_time)
        .collect();
    diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = diffs.len();
    if n % 2 == 1 {
        diffs[n / 2]
    } else {
        (diffs[n / 2 - 1] + diffs[n / 2]) / 2.0
    }
}

/// The outcome of scoring one `(seed, charge)` hypothesis.
struct HypothesisScore {
    response: f64,
    charge: i32,
    mono_mz: f64,
    peaks: Vec<IndexedMassSpectralPeak>,
    num_isotopes_observed: usize,
}

/// Scores a single charge hypothesis for a seed peak: lays the isotope comb (anchored so the
/// most-abundant tooth sits on the seed), evaluates the RT Gaussian across the scan window, and
/// sums `weight · gaussian · observed_intensity` over every comb `(m/z, scan)` point. Peaks already
/// in `claimed` are treated as absent (this is what makes cross-feature NMS work).
fn score_hypothesis(
    engine: &PeakIndexingEngine,
    seed: &IndexedMassSpectralPeak,
    charge: i32,
    params: &TraceKernelParameters,
    ppm: &PpmTolerance,
    claimed: &HashSet<PeakKey>,
) -> HypothesisScore {
    let seed_mz = seed.m() as f64;
    let seed_mass = mz_to_mass(seed_mz, charge);
    let weights = match params.weight_model {
        CombWeightModel::Poisson => {
            poisson_comb_weights(seed_mass, params.min_isotope_weight, params.max_isotopes)
        }
        CombWeightModel::Averagine => crate::deconvolution::averagine_comb_weights(
            seed_mass,
            params.min_isotope_weight,
            params.max_isotopes,
        ),
    };
    let i_star = most_abundant_index(&weights);
    let spacing = C13_MINUS_C12 / charge as f64;
    let mono_mz = seed_mz - (i_star as f64) * spacing;

    let scan_info = engine.scan_info();
    let n_scans = scan_info.len() as i32;
    let apex = seed.zero_based_scan_index;
    let rt_apex = seed.retention_time as f64;

    let mut response = 0.0;
    let mut peaks: Vec<IndexedMassSpectralPeak> = Vec::new();
    let mut observed_isotopes: HashSet<usize> = HashSet::new();
    // Peaks already used *within this hypothesis*. For higher charges the comb spacing (1.0033/z) is
    // small, so two adjacent isotope slots can resolve to the same physical peak; without this a peak
    // would be double-counted in the response and intensity and would inflate the isotope count.
    let mut used: HashSet<PeakKey> = HashSet::new();

    // Collect the scans within the *time* window, walking outward from the apex in both directions
    // and stopping as soon as RT leaves the window. Scans are RT-ordered, so this is a bounded walk
    // that adapts to the local (uneven) MS1 spacing instead of a fixed scan count.
    let rt_win = params.rt_half_window_minutes;
    let mut window_scans: Vec<i32> = Vec::new();
    let mut s = apex;
    while s >= 0 && (scan_info[s as usize].retention_time - rt_apex).abs() <= rt_win {
        window_scans.push(s);
        s -= 1;
    }
    let mut s = apex + 1;
    while s < n_scans && (scan_info[s as usize].retention_time - rt_apex).abs() <= rt_win {
        window_scans.push(s);
        s += 1;
    }

    for &s in &window_scans {
        let rt_s = scan_info[s as usize].retention_time;
        let g = gaussian(rt_s - rt_apex, params.rt_sigma_minutes);

        for (k, &wk) in weights.iter().enumerate() {
            let expected_mz = mono_mz + (k as f64) * spacing;
            let peak = match engine.get_indexed_peak(expected_mz, s, ppm) {
                Some(p) => p,
                None => continue,
            };
            let key = peak.key();
            if claimed.contains(&key) {
                continue;
            }
            // Skip a peak already consumed by another isotope slot of this same hypothesis.
            if !used.insert(key) {
                continue;
            }
            response += wk * g * peak.intensity as f64;
            peaks.push(*peak);
            observed_isotopes.insert(k);
        }
    }

    HypothesisScore {
        response,
        charge,
        mono_mz,
        peaks,
        num_isotopes_observed: observed_isotopes.len(),
    }
}

/// Runs the trace-kernel detector over an indexed run.
///
/// Seeds are every indexed peak taken tallest-first (the busiest RT regions hold the most signal, so
/// this maximises explained-signal-per-feature — the coverage objective). For each unclaimed seed,
/// charge hypotheses `min..=max` are scored and the best-response one wins (cross-z NMS). A
/// hypothesis is accepted when it observes at least `min_isotopes_observed` teeth and has positive
/// response; its peaks are then claimed so overlapping harmonics cannot re-fire. Returns the detected
/// features in acceptance order (tallest-seed first).
pub fn detect_features(
    engine: &PeakIndexingEngine,
    params: &TraceKernelParameters,
) -> Vec<DetectedFeature> {
    let ppm = PpmTolerance::new(params.ppm_tolerance);

    // Seeds: all peaks, tallest first. Stable ordering (intensity desc) mirrors get_all_xics.
    let mut seeds = engine.all_peaks();
    seeds.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));

    // Coverage bookkeeping: denominator = Σ all peak intensities; stop once the claimed fraction
    // reaches `coverage_target`.
    let total_intensity: f64 = seeds.iter().map(|p| p.intensity as f64).sum();
    let coverage_stop = if params.coverage_target < 1.0 && total_intensity > 0.0 {
        params.coverage_target * total_intensity
    } else {
        f64::INFINITY
    };
    let mut explained_intensity = 0.0;

    let mut claimed: HashSet<PeakKey> = HashSet::new();
    let mut features: Vec<DetectedFeature> = Vec::new();

    for seed in &seeds {
        // Seeds are intensity-descending, so once we fall below the floor every remaining seed is
        // too — stop rather than continue.
        if (seed.intensity as f64) < params.min_seed_intensity {
            break;
        }
        if claimed.contains(&seed.key()) {
            continue;
        }

        // Score every charge hypothesis; keep the highest response (cross-z non-max suppression).
        let mut best: Option<HypothesisScore> = None;
        for z in params.min_charge..=params.max_charge {
            if z == 0 {
                continue;
            }
            let score = score_hypothesis(engine, seed, z, params, &ppm, &claimed);
            let better = match &best {
                None => true,
                Some(b) => score.response > b.response,
            };
            if better {
                best = Some(score);
            }
        }

        let best = match best {
            Some(b) => b,
            None => continue,
        };

        if best.num_isotopes_observed < params.min_isotopes_observed || best.response <= 0.0 {
            // Not a feature; retire this seed so we do not reconsider it.
            claimed.insert(seed.key());
            continue;
        }

        // Claim the feature's peaks so overlapping harmonics / neighbours cannot re-fire.
        for p in &best.peaks {
            claimed.insert(p.key());
        }
        let feature = build_feature(&best);
        explained_intensity += feature.summed_intensity;
        features.push(feature);

        // Stop once we have explained the target fraction of the total MS1 signal.
        if explained_intensity >= coverage_stop {
            break;
        }
    }

    features
}

/// Assembles a [`DetectedFeature`] from an accepted hypothesis (apex = tallest claimed peak;
/// RT bounds and summed intensity from the claimed peaks).
fn build_feature(hyp: &HypothesisScore) -> DetectedFeature {
    let apex = hyp
        .peaks
        .iter()
        .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
        .expect("accepted hypothesis has at least one peak");
    let start_rt = hyp
        .peaks
        .iter()
        .map(|p| p.retention_time as f64)
        .fold(f64::INFINITY, f64::min);
    let end_rt = hyp
        .peaks
        .iter()
        .map(|p| p.retention_time as f64)
        .fold(f64::NEG_INFINITY, f64::max);
    let summed_intensity = hyp.peaks.iter().map(|p| p.intensity as f64).sum();

    DetectedFeature {
        monoisotopic_mass: mz_to_mass(hyp.mono_mz, hyp.charge),
        charge: hyp.charge,
        mono_mz: hyp.mono_mz,
        apex_scan_index: apex.zero_based_scan_index,
        apex_rt: apex.retention_time as f64,
        start_rt,
        end_rt,
        summed_intensity,
        score: hyp.response,
        num_isotopes_observed: hyp.num_isotopes_observed,
        peaks: hyp.peaks.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isotopic_envelope::mass_to_mz_f64;
    use crate::peak_indexing::Scan;

    fn approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "expected {b}, got {a} (tol {tol})");
    }

    #[test]
    fn poisson_weights_light_mass_mono_is_tallest() {
        // M = 1000 → λ = 0.48, mode at k=0 (mono tallest), strictly descending.
        let w = poisson_comb_weights(1000.0, 1e-3, 12);
        assert_eq!(most_abundant_index(&w), 0);
        assert!(w[0] >= w[1] && w[1] >= w[2]);
        approx(w[0], 1.0, 1e-12); // normalised so the max is 1
    }

    #[test]
    fn poisson_weights_heavy_mass_mode_shifts_up() {
        // M = 5000 → λ = 2.4, mode at k=2 — the seed (tallest peak) is NOT the monoisotopic peak,
        // which is exactly why the comb must look below the seed.
        let w = poisson_comb_weights(5000.0, 1e-3, 20);
        assert_eq!(most_abundant_index(&w), 2);
    }

    #[test]
    fn median_scan_spacing_is_robust() {
        let scans = synthetic_envelope_scans().0;
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        // RTs are 0.1 min apart.
        approx(median_ms1_scan_spacing_minutes(engine.scan_info()), 0.1, 1e-9);
    }

    /// Builds a clean charge-2 isotope envelope (monoisotopic neutral mass 1000) eluting across
    /// 9 scans with a Gaussian RT profile (apex at scan 4). Each isotope tooth carries its Poisson
    /// weight × the apex intensity × the RT Gaussian. Returns the scans and the true mono m/z.
    fn synthetic_envelope_scans() -> (Vec<Scan>, f64) {
        let mono_mass = 1000.0;
        let charge = 2;
        let mono_mz = mass_to_mz_f64(mono_mass, charge); // ~501.007
        let spacing = C13_MINUS_C12 / charge as f64;
        let weights = poisson_comb_weights(mono_mass, 1e-4, 8);
        let apex_intensity = 1.0e7;
        let rt_sigma = 0.15;

        let n_scans = 9;
        let apex_scan = 4;
        let mut scans = Vec::new();
        for s in 0..n_scans {
            let rt = 10.0 + s as f64 * 0.1;
            let g = gaussian(rt - (10.0 + apex_scan as f64 * 0.1), rt_sigma);
            let mut mz = Vec::new();
            let mut intensity = Vec::new();
            for (k, &wk) in weights.iter().enumerate() {
                mz.push(mono_mz + k as f64 * spacing);
                intensity.push(apex_intensity * wk * g);
            }
            scans.push(Scan {
                mz,
                intensity,
                one_based_scan_number: s + 1,
                retention_time: rt,
                msn_order: 1,
            });
        }
        (scans, mono_mz)
    }

    #[test]
    fn detects_charge_two_envelope_with_correct_mass() {
        let (scans, mono_mz) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &params);
        assert!(!features.is_empty(), "should detect at least one feature");

        // The top feature (tallest seed) must be the charge-2 envelope with mono mass ~1000.
        let top = &features[0];
        assert_eq!(top.charge, 2, "charge-2 envelope must beat other z hypotheses");
        approx(top.monoisotopic_mass, 1000.0, 0.01);
        approx(top.mono_mz, mono_mz, 1e-4);
        assert!(top.num_isotopes_observed >= 2);
        assert_eq!(top.apex_scan_index, 4, "apex is the max-intensity scan");
        assert!(top.score > 0.0);
    }

    #[test]
    fn detects_charge_two_envelope_with_averagine_weights() {
        // The averagine comb-weight model should detect the same synthetic z=2 envelope with the
        // correct mass/charge, exercising the CombWeightModel::Averagine path end to end.
        let (scans, mono_mz) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            weight_model: CombWeightModel::Averagine,
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &params);
        assert!(!features.is_empty(), "averagine model should detect the envelope");
        let top = &features[0];
        assert_eq!(top.charge, 2);
        approx(top.monoisotopic_mass, 1000.0, 0.01);
        approx(top.mono_mz, mono_mz, 1e-4);
    }

    #[test]
    fn charge_two_beats_charge_one_and_four_on_response() {
        // Directly confirm the cross-z ranking on the same seed: score z=1,2,4 for the apex mono
        // peak and check z=2 wins.
        let (scans, _) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            ..TraceKernelParameters::default()
        };
        let ppm = PpmTolerance::new(params.ppm_tolerance);
        let claimed = HashSet::new();

        // Seed = the tallest peak (apex scan mono).
        let mut seeds = engine.all_peaks();
        seeds.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));
        let seed = seeds[0];

        let s1 = score_hypothesis(&engine, &seed, 1, &params, &ppm, &claimed).response;
        let s2 = score_hypothesis(&engine, &seed, 2, &params, &ppm, &claimed).response;
        let s4 = score_hypothesis(&engine, &seed, 4, &params, &ppm, &claimed).response;
        assert!(s2 > s1, "z=2 ({s2}) should beat z=1 ({s1})");
        assert!(s2 > s4, "z=2 ({s2}) should beat z=4 ({s4})");
    }

    #[test]
    fn claiming_prevents_double_detection() {
        // A single clean envelope should yield exactly one feature — its peaks get claimed, so no
        // second feature is assembled from the same signal.
        let (scans, _) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &params);
        assert_eq!(features.len(), 1, "one envelope → one feature");
    }

    #[test]
    fn with_rt_from_scans_derives_sigma_and_window() {
        let (scans, _) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        // 36 s FWHM → σ = 0.6 min / 2.3548 ≈ 0.2548 min; spacing 0.1 min → half-window ≈ round(5.1) = 5.
        let params = TraceKernelParameters::default()
            .with_rt_from_scans(engine.scan_info(), 36.0);
        approx(params.rt_sigma_minutes, (36.0 / 60.0) / FWHM_TO_SIGMA, 1e-9);
        assert_eq!(params.half_window_scans, 5);
    }
}
