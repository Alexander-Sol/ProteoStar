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
//! Both build a [`Simulator`] and discard it. A workbench keeps one instead and brings it to each
//! edited spec with [`Simulator::update`], which rebuilds only what the edit invalidated.
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
pub mod tic;

use std::fmt;
use std::sync::Mutex;

use flashlfq_core::peak_indexing::Scan;
use rayon::prelude::*;

use crate::centroider::{Centroids, EnvelopeLayout, Parallelism, ProfileCentroider, ScanMax};
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

impl NoiseConfig {
    fn from_spec(spec: &SimSpec) -> Option<Self> {
        spec.noise.as_ref().map(|n| Self {
            seed: n.seed as u64,
            floor: NoiseFloorModel::from(&n.floor),
            jitter: PeakJitterModel::from_spec(&n.jitter),
        })
    }
}

/// The pre-noise, pre-reduction centroids of the scan simulated last, and the prune floor they
/// were found above. A noise or reduction edit that keeps the prune floor reuses them.
struct ViewedSignal {
    scan_index: usize,
    prune_floor: f64,
    signal: Centroids,
}

/// A spec decoded into its models, ready to simulate any scan of the run — and to be brought to
/// an edited spec with [`Self::update`], which recomputes only what the edit invalidated.
///
/// Every output is a pure function of the current spec: an updated simulator produces exactly
/// what [`Self::new`] on the same spec would. The caches below only skip work.
pub struct Simulator {
    spec: SimSpec,
    centroider: ProfileCentroider,
    reduction: SimReduction,
    noise: Option<NoiseConfig>,
    /// Per scan, its brightest centroid; `None` until computed or after an edit invalidated it.
    /// The noise-free relative floor is a fraction of their maximum. They depend only on the
    /// signal model, so noise and reduction edits keep them.
    scan_maxima: Mutex<Vec<Option<ScanMax>>>,
    viewed: Mutex<Option<ViewedSignal>>,
}

impl Simulator {
    pub fn new(spec: &SimSpec) -> Result<Self, SimError> {
        spec.validate()?;
        let centroider = ProfileCentroider::new(
            signal_models(spec),
            spec.charge_range.min,
            spec.charge_range.max,
            PeakWidthModel::from(spec.peak_width_model),
        )?;
        Ok(Self {
            centroider,
            reduction: spec.reduction.unwrap_or_default(),
            noise: NoiseConfig::from_spec(spec),
            scan_maxima: Mutex::new(vec![None; spec.scan_times.len()]),
            viewed: Mutex::new(None),
            spec: spec.clone(),
        })
    }

    /// The spec this simulator currently simulates.
    pub fn spec(&self) -> &SimSpec {
        &self.spec
    }

    /// Brings the simulator to `spec`, keeping every envelope and cached scan maximum the
    /// difference leaves valid (see `agent_info/TopDown-Simulator-Live-Session-Plan.md`). On
    /// error the simulator is unchanged.
    pub fn update(&mut self, spec: &SimSpec) -> Result<(), SimError> {
        if *spec == self.spec {
            return Ok(());
        }
        spec.validate()?;

        let layout = EnvelopeLayout {
            min_charge: spec.charge_range.min,
            max_charge: spec.charge_range.max,
            width: PeakWidthModel::from(spec.peak_width_model),
        };
        let change = self.centroider.update(signal_models(spec), layout)?;

        let maxima = self.scan_maxima.get_mut().unwrap_or_else(|e| e.into_inner());
        if spec.scan_times != self.spec.scan_times || change.everything {
            *maxima = vec![None; spec.scan_times.len()];
        } else if !change.is_empty() {
            let times = &spec.scan_times;
            maxima.par_iter_mut().enumerate().for_each(|(s, m)| {
                if let Some(max) = m {
                    if !max.survives(change.tallest_changed_component(times[s])) {
                        *m = None;
                    }
                }
            });
        }

        let viewed = self.viewed.get_mut().unwrap_or_else(|e| e.into_inner());
        if spec.scan_times != self.spec.scan_times || !change.is_empty() {
            *viewed = None;
        }

        self.reduction = spec.reduction.unwrap_or_default();
        self.noise = NoiseConfig::from_spec(spec);
        self.spec = spec.clone();
        Ok(())
    }

    pub fn scan_count(&self) -> usize {
        self.spec.scan_times.len()
    }

    /// The noise-free, unreduced centroids of one scan (C# `SimulateCentroided`).
    pub fn centroid_scan(&self, scan_index: usize) -> Result<Centroids, SimError> {
        let t = self.scan_time(scan_index)?;
        Ok(self.centroider.centroid(t, Parallelism::AcrossClusters))
    }

    /// [`Self::centroid_scan`] for every scan, in parallel.
    pub fn centroid_run(&self) -> Vec<Centroids> {
        self.spec
            .scan_times
            .par_iter()
            .map(|&t| self.centroider.centroid(t, Parallelism::Serial))
            .collect()
    }

    /// One scan through the full pipeline.
    ///
    /// Without noise, the reduction floor is relative to the brightest centroid of the *whole
    /// run*, so this needs every scan's maximum — cached, so only the scans no earlier call (or
    /// edit) has left valid are computed. With noise the floor is the noise model's and no other
    /// scan is touched.
    pub fn simulate_scan(&self, scan_index: usize) -> Result<Scan, SimError> {
        let t = self.scan_time(scan_index)?;
        let floor = self.reduction_floor();
        let prune_floor = prune_floor(floor);

        let mut viewed = self.viewed.lock().unwrap_or_else(|e| e.into_inner());
        let signal = match &*viewed {
            Some(v) if v.scan_index == scan_index && v.prune_floor.to_bits() == prune_floor.to_bits() => {
                v.signal.clone()
            }
            _ => {
                let signal = self.centroider.centroid_above(t, prune_floor, Parallelism::AcrossClusters);
                *viewed = Some(ViewedSignal { scan_index, prune_floor, signal: signal.clone() });
                signal
            }
        };
        drop(viewed);

        Ok(self.finish_scan(scan_index, signal, floor))
    }

    /// Every scan through the full pipeline, in parallel. Output is independent of scheduling:
    /// each scan's noise comes from its own seeded stream.
    pub fn simulate_run(&self) -> Vec<Scan> {
        let floor = self.reduction_floor();
        let prune_floor = prune_floor(floor);
        self.spec
            .scan_times
            .par_iter()
            .enumerate()
            .map(|(s, &t)| {
                let signal = self.centroider.centroid_above(t, prune_floor, Parallelism::Serial);
                self.finish_scan(s, signal, floor)
            })
            .collect()
    }

    /// The brightest centroid of the run, computing only the scan maxima not already cached.
    pub fn global_max(&self) -> f64 {
        let mut maxima = self.scan_maxima.lock().unwrap_or_else(|e| e.into_inner());
        let times = &self.spec.scan_times;
        maxima.par_iter_mut().enumerate().for_each(|(s, m)| {
            if m.is_none() {
                *m = Some(self.centroider.scan_max(times[s], Parallelism::Serial));
            }
        });
        maxima
            .iter()
            .map(|m| m.expect("every scan maximum was just filled").brightest)
            .fold(0.0, |max, y| if y > max { y } else { max })
    }

    /// How many scans' maxima are cached (diagnostic: what the last edit left valid).
    pub fn cached_scan_maxima(&self) -> usize {
        self.scan_maxima.lock().unwrap_or_else(|e| e.into_inner()).iter().flatten().count()
    }

    /// The floor signal is reduced against: flat without noise, the noise model's m/z-dependent
    /// detection limit (represented by its minimum) with it.
    fn reduction_floor(&self) -> ReductionFloor {
        match &self.noise {
            Some(n) => ReductionFloor::Noise(n.floor.min_detection_floor()),
            None => {
                let r = self.reduction;
                ReductionFloor::Flat(
                    r.absolute_intensity_floor.max(self.global_max() * r.relative_intensity_threshold),
                )
            }
        }
    }

    fn scan_time(&self, scan_index: usize) -> Result<f64, SimError> {
        self.spec
            .scan_times
            .get(scan_index)
            .copied()
            .ok_or(SimError::ScanOutOfRange { index: scan_index, count: self.spec.scan_times.len() })
    }

    /// Reduce, then inject noise — in that order, because the noise-free relative floor sits far
    /// above the noise amplitude and would erase any noise added before it.
    fn finish_scan(&self, scan_index: usize, signal: Centroids, floor: ReductionFloor) -> Scan {
        let spectrum = match (&self.noise, floor) {
            // Noise switches reduction to the noise model's own detection limit, so signal and
            // noise are thresholded under one rule.
            (Some(n), _) => noise::inject_scan(
                &threshold(signal, |mz| n.floor.detection_floor_at(mz)),
                scan_index,
                n.seed,
                &n.floor,
                n.jitter.as_ref(),
                &self.centroider.width(),
            ),
            (None, ReductionFloor::Flat(floor)) => threshold(signal, |_| floor),
            (None, ReductionFloor::Noise(_)) => unreachable!("a noise floor implies noise"),
        };

        Scan {
            mz: spectrum.mz,
            intensity: spectrum.intensity,
            one_based_scan_number: scan_index as i32 + 1,
            retention_time: self.spec.scan_times[scan_index],
            msn_order: 1,
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum ReductionFloor {
    /// The noise model's lowest detection limit anywhere in m/z.
    Noise(f64),
    /// One floor for every peak (no noise).
    Flat(f64),
}

/// The floor centroiding may prune against: a lower bound on every floor the scan is reduced
/// against, so the shortcuts in [`ProfileCentroider::centroid_above`] never drop anything the
/// floor would have kept. Rounded down to a power of two, so small noise and reduction edits
/// keep the same prune floor and the viewed scan's cached signal stays valid.
fn prune_floor(floor: ReductionFloor) -> f64 {
    let (ReductionFloor::Noise(f) | ReductionFloor::Flat(f)) = floor;
    if f > 0.0 && f.is_finite() {
        // Clearing the mantissa of a positive normal double rounds it down to a power of two;
        // a subnormal becomes 0, which only disables pruning.
        f64::from_bits(f.to_bits() & 0xFFF0_0000_0000_0000)
    } else {
        0.0
    }
}

/// The proteoforms the spectra are rendered from. Background is rendered exactly like the
/// identified proteoforms; it only differs in being left out of the feature truth.
fn signal_models(spec: &SimSpec) -> Vec<ProteoformModel> {
    spec.proteoforms.iter().chain(&spec.background).map(ProteoformModel::from).collect()
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
