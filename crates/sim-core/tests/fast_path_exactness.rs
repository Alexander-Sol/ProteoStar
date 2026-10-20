//! The full pipeline's centroiding shortcuts (`ProfileCentroider::centroid_above`: pruning
//! Gaussians far below the floor, skipping grid evaluations that cannot reach it, and widening
//! the interaction reach to match) must keep the peaks the floor keeps, at the right values. The
//! C# goldens are too sparse to exercise them, so this checks a crowded synthetic run against the
//! exact centroider and against the **true profile** — summed by brute force from the model
//! pieces, sharing nothing with the centroider.
//!
//! The exact path (the C# algorithm) is not itself the truth at the margins:
//!
//! - It finds maxima on a σ/4 grid anchored at each cluster's first Gaussian, so where the
//!   profile is nearly flat — a shoulder or wobble on a flank, a flat top — whether a maximum is
//!   reported, and where, depends on where the grid falls. Pruning moves cluster starts.
//! - It treats Gaussians > 6σ apart as independent, so near a giant it can be off by the giant's
//!   `e^-18` tail, and can report a "peak" on the giant's slope. The fast path reaches further.
//!
//! So nearly every peak must match the exact path at the contract tolerance, and every peak that
//! doesn't must be explained by the true profile: a spot too flat for the grid, a real maximum
//! the fast path has right, or an exact-path artefact. A real maximum that only the exact path
//! found fails the test, as does a sampled fast peak whose value the true profile contradicts.

use std::f64::consts::PI;

use sim_core::model::{to_mz, IsotopeEnvelope, PeakWidthModel, ProteoformModel};
use sim_core::spec::*;
use sim_core::Simulator;

#[derive(Clone, Copy)]
struct Tolerance {
    /// m/z tolerance in units of the local σ, or `None` for the contract's 1e-6 Th.
    mz_in_sigmas: Option<f64>,
    intensity_rel: f64,
}

/// Contract §5 tolerance.
const STRICT: Tolerance = Tolerance { mz_in_sigmas: None, intensity_rel: 1e-4 };
/// The exact path's own precision next to a much taller peak.
const LOOSE: Tolerance = Tolerance { mz_in_sigmas: Some(0.01), intensity_rel: 1e-3 };
/// Sanity cap on peaks outside the strict tolerance, each of which must also be explained.
const MAX_NON_STRICT_FRACTION: f64 = 1e-2;
/// A point is grid-sensitive if, on at least one side, the profile changes by less than this
/// fraction across one σ/4 grid step. An isolated Gaussian changes by 3.1 % across the step
/// beside its apex, so no clean peak qualifies.
const GRID_FLAT: f64 = 0.01;
/// Every this-many-th fast peak has its value checked against the true profile.
const TRUTH_SAMPLE_STRIDE: usize = 200;
const K: f64 = 1.2e-7;

impl Tolerance {
    fn mz_th(&self, mz: f64) -> f64 {
        self.mz_in_sigmas.map_or(1e-6, |n| n * K * mz.powf(1.5))
    }
}

fn crowded_spec(reduction: SimReduction) -> SimSpec {
    // Deterministic spread; overlapping envelopes across the whole m/z range.
    let frac = |i: usize, k: u64| ((i as u64 * k) % 1000) as f64 / 1000.0;
    let proteoforms = (0..300)
        .map(|i| SimProteoform {
            monoisotopic_mass: 5000.0 + 45_000.0 * frac(i, 7919),
            abundance: 1e8 * (1.0 + 20.0 * frac(i, 104_729)),
            rt: SimRtProfile { mu: 20.0 + 4.0 * frac(i, 1543), sigma: 0.2, tau: 0.1 },
            charge: SimChargeDistribution::Gaussian { mu_z: 8.0 + 30.0 * frac(i, 389), sigma_z: 2.5 },
            identifier: None,
        })
        .collect();

    SimSpec {
        schema_version: CURRENT_SCHEMA_VERSION,
        generator: None,
        source_run: None,
        instrument: SimInstrument::default(),
        charge_range: SimChargeRange { min: 3, max: 60 },
        peak_width_model: SimPeakWidthModel::Orbitrap { k: K },
        scan_times: vec![19.0, 21.0, 22.5, 24.0, 26.0],
        reduction: Some(reduction),
        noise: None,
        proteoforms,
        background: Vec::new(),
    }
}

/// One scan's true profile: every Gaussian of the forward model, summed out to 10σ.
struct TrueProfile {
    /// (centre, height, σ), ascending centre.
    gaussians: Vec<(f64, f64, f64)>,
    max_sigma: f64,
    width: PeakWidthModel,
}

impl TrueProfile {
    fn new(spec: &SimSpec, t: f64) -> Self {
        let width: PeakWidthModel = spec.peak_width_model.into();
        let mut gaussians = Vec::new();
        for p in &spec.proteoforms {
            let m = ProteoformModel::from(p);
            let w = m.abundance * m.rt.evaluate(t);
            if w.is_nan() || w <= 1e-18 {
                continue;
            }
            let env = IsotopeEnvelope::new(m.monoisotopic_mass);
            for z in spec.charge_range.min..=spec.charge_range.max {
                let fz = m.charge.evaluate(z);
                for i in 0..env.len() {
                    let c = to_mz(env.neutral_masses[i], z);
                    let sigma = width.sigma_at(c);
                    gaussians.push((c, w * fz * env.intensities[i] / (sigma * (2.0 * PI).sqrt()), sigma));
                }
            }
        }
        gaussians.sort_by(|a, b| a.0.total_cmp(&b.0));
        let max_sigma = gaussians.iter().fold(0.0f64, |m, g| m.max(g.2));
        Self { gaussians, max_sigma, width }
    }

    fn at(&self, x: f64) -> f64 {
        let reach = 10.0 * self.max_sigma;
        let start = self.gaussians.partition_point(|g| g.0 < x - reach);
        self.gaussians[start..]
            .iter()
            .take_while(|g| g.0 <= x + reach)
            .map(|&(c, h, sigma)| {
                let d = x - c;
                h * (-d * d / (2.0 * sigma * sigma)).exp()
            })
            .sum()
    }

    /// Whether `intensity` is the profile's value at `mz`, within the contract tolerance.
    fn agrees(&self, mz: f64, intensity: f64) -> bool {
        let truth = self.at(mz);
        (intensity - truth).abs() / truth <= STRICT.intensity_rel
    }

    /// Whether `mz` is a local maximum at a resolution far finer than the grid.
    fn is_local_max(&self, mz: f64) -> bool {
        let delta = self.width.sigma_at(mz) / 32.0;
        let peak = self.at(mz) * (1.0 + 1e-9);
        peak >= self.at(mz - delta) && peak >= self.at(mz + delta)
    }

    fn is_real_peak(&self, mz: f64, intensity: f64) -> bool {
        self.agrees(mz, intensity) && self.is_local_max(mz)
    }

    /// Whether the profile at `mz` is too flat for the σ/4 grid to place a maximum robustly.
    fn is_grid_sensitive(&self, mz: f64) -> bool {
        let sigma = self.width.sigma_at(mz);
        let peak = self.at(mz);
        [-1.0, 1.0].into_iter().any(|dir| {
            [8.0, 4.0].into_iter().all(|div| (self.at(mz + dir * sigma / div) / peak - 1.0).abs() < GRID_FLAT)
        })
    }

    fn shape(&self, mz: f64, intensity: f64) -> String {
        let sigma = self.width.sigma_at(mz);
        let v: Vec<String> =
            (-16..=16).map(|k| format!("{:.4}", self.at(mz + k as f64 * sigma / 8.0) / intensity)).collect();
        format!("profile / peak at −2σ…+2σ in σ/8 steps: {}", v.join(" "))
    }
}

/// Peaks of `of` with no counterpart in `among` within `tol`.
fn unmatched(of: &[(f64, f64)], among: &[(f64, f64)], tol: Tolerance) -> Vec<(f64, f64)> {
    of.iter()
        .copied()
        .filter(|&(mz, y)| {
            let dmz = tol.mz_th(mz);
            let j = among.partition_point(|&(m, _)| m < mz - dmz);
            !among[j..]
                .iter()
                .take_while(|&&(m, _)| m <= mz + dmz)
                .any(|&(_, v)| (v - y).abs() / y <= tol.intensity_rel)
        })
        .collect()
}

fn assert_fast_path_correct(reduction: SimReduction) {
    let spec = crowded_spec(reduction);
    let sim = Simulator::new(&spec).unwrap();

    let exact = sim.centroid_run();
    let global_max = exact.iter().flat_map(|c| c.intensity.iter().copied()).fold(0.0, f64::max);
    let floor = reduction
        .absolute_intensity_floor
        .max(global_max * reduction.relative_intensity_threshold);
    let fast = sim.simulate_run();

    let (mut kept, mut non_strict, mut sampled) = (0usize, 0usize, 0usize);
    for (s, (exact, fast)) in exact.iter().zip(&fast).enumerate() {
        let truth = TrueProfile::new(&spec, spec.scan_times[s]);
        let want: Vec<(f64, f64)> = exact
            .mz
            .iter()
            .zip(&exact.intensity)
            .filter(|&(_, &y)| y >= floor)
            .map(|(&mz, &y)| (mz, y))
            .collect();
        let got: Vec<(f64, f64)> = fast.mz.iter().copied().zip(fast.intensity.iter().copied()).collect();

        // A fast peak off the exact path must be a real peak (or at a flat spot).
        let fast_only = unmatched(&got, &want, STRICT);
        for &(mz, y) in &fast_only {
            assert!(
                truth.is_grid_sensitive(mz) || truth.is_real_peak(mz, y),
                "scan {s}: fast peak at m/z {mz} ({:.3}× floor) differs from the exact path and is \
                 not a real maximum of the true profile; {}",
                y / floor,
                truth.shape(mz, y)
            );
        }

        // An exact peak the fast path lacks at the contract tolerance must either be present
        // within the exact path's own precision (that fast counterpart was checked above), or be
        // disowned by the true profile: an exact-path artefact, or at a flat spot.
        let exact_only = unmatched(&want, &got, STRICT);
        let loosely_missing = unmatched(&want, &got, LOOSE);
        for &(mz, y) in &exact_only {
            let explained = !loosely_missing.contains(&(mz, y))
                || truth.is_grid_sensitive(mz)
                || !truth.is_real_peak(mz, y);
            assert!(
                explained,
                "scan {s}: real maximum at m/z {mz} ({:.3}× floor) is missing from the fast path; {}",
                y / floor,
                truth.shape(mz, y)
            );
        }

        // Spot-check fast values against the truth directly.
        for &(mz, y) in got.iter().step_by(TRUTH_SAMPLE_STRIDE) {
            assert!(
                truth.agrees(mz, y) || truth.is_grid_sensitive(mz),
                "scan {s}: fast peak at m/z {mz} has intensity {y}, true profile {}",
                truth.at(mz)
            );
            sampled += 1;
        }

        non_strict += fast_only.len() + exact_only.len();
        kept += want.len();
    }

    let total: usize = exact.iter().map(|c| c.mz.len()).sum();
    println!(
        "kept {kept} of {total} exact centroids at floor {floor:e}; {non_strict} differ from the \
         exact path beyond the contract tolerance, all explained by the true profile; \
         {sampled} sampled fast peaks agree with it"
    );
    assert!(kept > 1000 && kept < total, "the floor should matter in this test");
    assert!(
        non_strict as f64 <= MAX_NON_STRICT_FRACTION * kept as f64,
        "{non_strict} of {kept} peaks differ from the exact path"
    );
}

#[test]
fn fast_path_correct_under_relative_floor() {
    assert_fast_path_correct(SimReduction::default());
}

/// An absolute floor at the default noise model's detection limit — the regime the live path
/// runs in when noise is on.
#[test]
fn fast_path_correct_under_noise_level_floor() {
    assert_fast_path_correct(SimReduction {
        relative_intensity_threshold: 0.0,
        absolute_intensity_floor: 1.94 * 1120.0,
    });
}
