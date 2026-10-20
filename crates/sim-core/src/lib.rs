//! Pure-Rust port of the mzLib `TopDownSimulator` MS1 simulation (forward model, centroider,
//! reduction, noise), driven by a SimSpec v1 file.
//!
//! See `agent_info/TopDown-Simulator-Rust-Port-Plan.md` (Phase 1) and
//! `agent_info/TopDown-SimSpec-Contract-v1.md`. Two granularities, matching the workbench's
//! interaction model:
//!
//! - [`simulate_scan`] — one scan, the live path while parameters are being edited.
//! - [`simulate_run`] — every scan, the commit path that gets indexed as a dataset.
//!
//! Both run the C# `Simulator.PrepareMs1` pipeline: centroid the summed profile, threshold
//! against the detection floor, then (if the spec has noise) jitter the signal, add the sampled
//! noise floor and collapse unresolved peaks. [`Simulator::centroid_scan`] stops after the first
//! step, exactly — it is what `Simulator.SimulateCentroided` and the `golden-v1` golden produce,
//! bit for bit.
//!
//! The full pipeline is faster than the C# rather than identical to it: it skips Gaussians far
//! below the detection floor before centroiding (survivors agree with C# within the contract
//! tolerance), and its noise uses native random streams (agrees with C# statistically).

// `!(x > 0.0)` is the C# guard idiom this crate ports, and it is deliberate: unlike `x <= 0.0`
// it also rejects NaN, which is what keeps the port's NaN behaviour identical to the C#.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

pub mod averagine;
pub mod centroider;
pub mod model;
pub mod noise;
pub mod spec;

use std::fmt;
use std::sync::OnceLock;

use flashlfq_core::peak_indexing::Scan;
use rayon::prelude::*;

use crate::centroider::{Centroids, Parallelism, ProfileCentroider};
use crate::model::{PeakWidthModel, ProteoformModel};
use crate::noise::{NoiseFloorModel, PeakJitterModel};
use crate::spec::{SimReduction, SimSpec};

#[derive(Debug)]
pub enum SimError {
    Io(std::io::Error),
    Json(serde_json::Error),
    UnsupportedSchema(u32),
    InvalidSpec(String),
    ScanOutOfRange { index: usize, count: usize },
}

impl fmt::Display for SimError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::Json(e) => write!(f, "invalid SimSpec JSON: {e}"),
            Self::UnsupportedSchema(v) => write!(
                f,
                "SimSpec schema v{v} is not supported; this build reads v{}",
                spec::CURRENT_SCHEMA_VERSION
            ),
            Self::InvalidSpec(msg) => write!(f, "invalid SimSpec: {msg}"),
            Self::ScanOutOfRange { index, count } => {
                write!(f, "scan index {index} is out of range for a {count}-scan run")
            }
        }
    }
}

impl std::error::Error for SimError {}

/// The live path: one fully simulated scan (reduced, and with noise if the spec has it).
pub fn simulate_scan(spec: &SimSpec, scan_index: usize) -> Result<Scan, SimError> {
    Simulator::new(spec)?.simulate_scan(scan_index)
}

/// The commit path: every scan of the run, in scan order.
pub fn simulate_run(spec: &SimSpec) -> Result<Vec<Scan>, SimError> {
    Ok(Simulator::new(spec)?.simulate_run())
}

struct NoiseConfig {
    seed: u64,
    floor: NoiseFloorModel,
    jitter: Option<PeakJitterModel>,
}

/// A spec decoded into its models, ready to simulate any scan of the run.
pub struct Simulator {
    centroider: ProfileCentroider,
    scan_times: Vec<f64>,
    width: PeakWidthModel,
    reduction: SimReduction,
    noise: Option<NoiseConfig>,
    /// Brightest unreduced centroid over the whole run; the relative reduction floor is a
    /// fraction of it. Only needed when the spec has no noise (noise brings its own floor).
    global_max: OnceLock<f64>,
    /// No scan is reduced against a lower floor than this, so the full pipeline's centroiding
    /// can skip work below it (see [`ProfileCentroider::centroid_above`]).
    floor_lower_bound: f64,
}

impl Simulator {
    pub fn new(spec: &SimSpec) -> Result<Self, SimError> {
        spec.validate()?;

        // Background is rendered exactly like the identified proteoforms; it only differs in
        // being left out of the feature truth.
        let models: Vec<ProteoformModel> =
            spec.proteoforms.iter().chain(&spec.background).map(ProteoformModel::from).collect();
        let width = PeakWidthModel::from(spec.peak_width_model);
        let centroider = ProfileCentroider::new(
            models,
            spec.charge_range.min,
            spec.charge_range.max,
            width,
        )?;

        let noise = spec.noise.as_ref().map(|n| NoiseConfig {
            seed: n.seed as u64,
            floor: NoiseFloorModel::from(&n.floor),
            jitter: PeakJitterModel::from_spec(&n.jitter),
        });
        let reduction = spec.reduction.unwrap_or_default();

        // A lower bound on the floor every scan will be reduced against, so the shortcuts in
        // centroiding never drop anything the floor would have kept.
        let floor_lower_bound = match &noise {
            Some(n) => n.floor.min_detection_floor(),
            None => {
                // The relative floor needs the run's brightest centroid, which is only known after
                // centroiding. The tallest single Gaussian is a lower bound on it and cheap.
                let tallest = spec
                    .scan_times
                    .par_iter()
                    .map(|&t| centroider.max_component_height(t))
                    .reduce(|| 0.0, f64::max);
                reduction
                    .absolute_intensity_floor
                    .max(tallest * reduction.relative_intensity_threshold)
            }
        };

        Ok(Self {
            centroider,
            scan_times: spec.scan_times.clone(),
            width,
            reduction,
            noise,
            global_max: OnceLock::new(),
            floor_lower_bound,
        })
    }

    pub fn scan_count(&self) -> usize {
        self.scan_times.len()
    }

    /// The noise-free, unreduced centroids of one scan (C# `SimulateCentroided`).
    pub fn centroid_scan(&self, scan_index: usize) -> Result<Centroids, SimError> {
        let t = self.scan_time(scan_index)?;
        Ok(self.centroider.centroid(t, Parallelism::AcrossClusters))
    }

    /// [`Self::centroid_scan`] for every scan, in parallel.
    pub fn centroid_run(&self) -> Vec<Centroids> {
        self.scan_times
            .par_iter()
            .map(|&t| self.centroider.centroid(t, Parallelism::Serial))
            .collect()
    }

    /// One scan through the full pipeline.
    ///
    /// Without noise, the reduction floor is relative to the brightest peak of the *whole run*,
    /// so the first call centroids every scan to find it (cached for later calls on this
    /// instance). With noise the floor is the noise model's and no other scan is touched.
    pub fn simulate_scan(&self, scan_index: usize) -> Result<Scan, SimError> {
        let signal = self.centroider.centroid_above(
            self.scan_time(scan_index)?,
            self.floor_lower_bound,
            Parallelism::AcrossClusters,
        );
        let global_max = match self.noise {
            Some(_) => 0.0,
            None => *self
                .global_max
                .get_or_init(|| global_max_intensity(&self.pruned_centroid_run())),
        };
        Ok(self.finish_scan(scan_index, signal, global_max))
    }

    /// Every scan through the full pipeline, in parallel. Output is independent of scheduling:
    /// each scan's noise comes from its own seeded stream.
    pub fn simulate_run(&self) -> Vec<Scan> {
        let signal = self.pruned_centroid_run();
        let global_max = *self.global_max.get_or_init(|| global_max_intensity(&signal));
        signal
            .into_par_iter()
            .enumerate()
            .map(|(s, centroids)| self.finish_scan(s, centroids, global_max))
            .collect()
    }

    fn pruned_centroid_run(&self) -> Vec<Centroids> {
        self.scan_times
            .par_iter()
            .map(|&t| self.centroider.centroid_above(t, self.floor_lower_bound, Parallelism::Serial))
            .collect()
    }

    fn scan_time(&self, scan_index: usize) -> Result<f64, SimError> {
        self.scan_times
            .get(scan_index)
            .copied()
            .ok_or(SimError::ScanOutOfRange { index: scan_index, count: self.scan_times.len() })
    }

    /// Reduce, then inject noise — in that order, because the noise-free relative floor sits far
    /// above the noise amplitude and would erase any noise added before it.
    fn finish_scan(&self, scan_index: usize, signal: Centroids, global_max: f64) -> Scan {
        let reduced = match &self.noise {
            // Noise switches reduction to the noise model's own detection limit, so signal and
            // noise are thresholded under one rule.
            Some(n) => threshold(signal, |mz| n.floor.detection_floor_at(mz)),
            None => {
                let r = self.reduction;
                let floor = r
                    .absolute_intensity_floor
                    .max(global_max * r.relative_intensity_threshold);
                threshold(signal, |_| floor)
            }
        };

        let spectrum = match &self.noise {
            Some(n) => noise::inject_scan(
                &reduced,
                scan_index,
                n.seed,
                &n.floor,
                n.jitter.as_ref(),
                &self.width,
            ),
            None => reduced,
        };

        Scan {
            mz: spectrum.mz,
            intensity: spectrum.intensity,
            one_based_scan_number: scan_index as i32 + 1,
            retention_time: self.scan_times[scan_index],
            msn_order: 1,
        }
    }
}

fn global_max_intensity(scans: &[Centroids]) -> f64 {
    scans
        .iter()
        .flat_map(|c| c.intensity.iter().copied())
        .fold(0.0, |max, y| if y > max { y } else { max })
}

/// Drops every peak below the floor at its m/z.
fn threshold(c: Centroids, floor_at: impl Fn(f64) -> f64) -> Centroids {
    let (mz, intensity) = c
        .mz
        .into_iter()
        .zip(c.intensity)
        .filter(|&(mz, y)| !(y < floor_at(mz)))
        .unzip();
    Centroids { mz, intensity }
}
