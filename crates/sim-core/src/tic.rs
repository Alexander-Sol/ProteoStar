//! The run's TIC without rendering its spectra: what the workbench's chromatogram view redraws on
//! every edit.
//!
//! A scan's TIC is the sum of its centroids after reduction. A proteoform's centroids, simulated
//! alone, scale exactly with its weight `A_p·g_p(t)` (its profile is linear in it), so they are
//! computed once, at unit weight, per (mass, charge distribution, layout). A scan's share is then
//! the weight times the summed unit heights of the centroids that clear the floor at that weight:
//! one binary search over the centroids sorted by the weight at which each clears it. An abundance
//! or RT edit costs one profile evaluation and one search per scan; nothing is re-centroided.
//!
//! Exact for proteoforms whose peaks do not overlap. Where two proteoforms' peaks merge into one
//! centroid the full simulation reports the merged maximum, which this counts as the sum of both.
//! Jitter dropping a peak below the floor is not modelled; otherwise jitter is mean-preserving.
//! Without noise, the floor is relative to the run's brightest centroid, which this takes as the
//! brightest isolated-proteoform centroid. The noise floor's share is exact: the same seeded draws
//! the full pipeline makes, summed.

use std::collections::HashMap;
use std::sync::Arc;

use rayon::prelude::*;
use serde::Serialize;

use crate::centroider::{Centroids, Parallelism, ProfileCentroider};
use crate::model::{EmgProfile, PeakWidthModel, ProteoformModel};
use crate::noise::{relative_noise_level_at, stream_for_scan, NoiseFloorModel, NOISE_STREAM};
use crate::spec::{SimNoise, SimReduction, SimSpec};
use crate::{signal_models, SimError};

/// A selected proteoform's XIC covers its centroids at least this fraction of its tallest: the
/// peaks a real XIC can be extracted at without being mostly other species. At 0.1 a fitted
/// Jurkat proteoform's real XIC picked up co-eluting species' peaks a third as tall as its own.
pub const XIC_TARGET_FRACTION: f64 = 0.3;

/// Same cut-off the forward model applies to a species' weight in a scan.
const MINIMUM_CONTRIBUTION_WEIGHT: f64 = 1e-18;

#[derive(Debug, Clone, Default, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RunTic {
    /// Per scan, the summed centroids of every proteoform (and background) that clear the floor.
    pub signal: Vec<f64>,
    /// Per scan, the summed noise peaks; zeros without noise.
    pub noise: Vec<f64>,
    /// The XIC of the proteoform asked for, if any.
    pub selected: Option<ProteoformXic>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ProteoformXic {
    /// Per scan, its summed centroids at [`Self::target_mzs`] that clear the floor.
    pub intensity: Vec<f64>,
    /// Its centroids at least [`XIC_TARGET_FRACTION`] of its tallest, ascending m/z: where to
    /// extract the matching XIC from a real run.
    pub target_mzs: Vec<f64>,
}

/// The detection floor at m/z `x` is `scale · shape(x)`. The shape fixes the order in which a
/// proteoform's centroids drop out as its weight falls, so only a shape change re-sorts them.
#[derive(Debug, Clone, Copy, PartialEq)]
enum FloorShape {
    /// The noise model's m/z-dependent detection limit.
    NoiseLevel,
    /// One floor everywhere (no noise).
    Flat,
}

impl FloorShape {
    fn at(self, mz: f64) -> f64 {
        match self {
            Self::NoiseLevel => relative_noise_level_at(mz),
            Self::Flat => 1.0,
        }
    }
}

/// A set of unit-weight centroids sorted by the weight at which each clears the floor.
#[derive(Debug, Clone, Default)]
struct Survival {
    /// Ascending: centroid `k` clears the floor at weight `w` when `w / scale ≥ thresholds[k]`.
    thresholds: Vec<f64>,
    /// `cumulative[k]` is the summed unit height of the first `k` centroids.
    cumulative: Vec<f64>,
}

impl Survival {
    fn new(unit: &Centroids, shape: FloorShape, keep: impl Fn(usize) -> bool) -> Self {
        let mut pairs: Vec<(f64, f64)> = (0..unit.mz.len())
            .filter(|&k| keep(k) && unit.intensity[k] > 0.0)
            .map(|k| (shape.at(unit.mz[k]) / unit.intensity[k], unit.intensity[k]))
            .collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut cumulative = Vec::with_capacity(pairs.len() + 1);
        let mut sum = 0.0;
        cumulative.push(sum);
        for &(_, h) in &pairs {
            sum += h;
            cumulative.push(sum);
        }
        Self { thresholds: pairs.into_iter().map(|p| p.0).collect(), cumulative }
    }

    /// The summed centroids at `weight` against a floor of `scale · shape`.
    fn intensity(&self, weight: f64, scale: f64) -> f64 {
        if !(weight > MINIMUM_CONTRIBUTION_WEIGHT) {
            return 0.0;
        }
        let x = if scale > 0.0 { weight / scale } else { f64::INFINITY };
        weight * self.cumulative[self.thresholds.partition_point(|&c| !(c > x))]
    }
}

/// One proteoform's unit-weight centroids: everything about it that an abundance or RT edit keeps.
#[derive(Debug)]
struct Envelope {
    unit: Centroids,
    survival: Survival,
    tallest: f64,
}

#[derive(Debug, Clone)]
struct Species {
    model: ProteoformModel,
    envelope: Arc<Envelope>,
    /// `A_p·g_p(t)` per scan.
    weights: Arc<Vec<f64>>,
}

/// What an envelope depends on besides the run layout: the mass and the charge distribution.
type EnvelopeKey = (u64, u64, u64);

fn envelope_key(m: &ProteoformModel) -> EnvelopeKey {
    (m.monoisotopic_mass.to_bits(), m.charge.mu_z.to_bits(), m.charge.sigma_z.to_bits())
}

/// What a species' per-scan weights depend on besides the scan times.
type WeightKey = (u64, u64, u64, u64);

fn weight_key(m: &ProteoformModel) -> WeightKey {
    let EmgProfile { mu, sigma, tau } = m.rt;
    (m.abundance.to_bits(), mu.to_bits(), sigma.to_bits(), tau.to_bits())
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Layout {
    min_charge: i32,
    max_charge: i32,
    width: PeakWidthModel,
    shape: FloorShape,
}

/// A spec's TIC, kept current across edits by [`Self::update`], which redoes only the proteoforms
/// an edit touched (and the noise floor's share only when the noise changed).
pub struct TicModel {
    spec: SimSpec,
    layout: Layout,
    species: Vec<Species>,
    noise_tic: Vec<f64>,
}

impl TicModel {
    pub fn new(spec: &SimSpec) -> Result<Self, SimError> {
        spec.validate()?;
        let layout = layout_of(spec);
        let species = build_species(spec, &layout, &[], &[])?;
        Ok(Self { noise_tic: noise_tic(spec), layout, species, spec: spec.clone() })
    }

    pub fn spec(&self) -> &SimSpec {
        &self.spec
    }

    /// Brings the model to `spec`. On error it is unchanged.
    pub fn update(&mut self, spec: &SimSpec) -> Result<(), SimError> {
        if *spec == self.spec {
            return Ok(());
        }
        spec.validate()?;
        let layout = layout_of(spec);
        let times_kept = spec.scan_times == self.spec.scan_times;
        let species = build_species(
            spec,
            &layout,
            if layout == self.layout { &self.species } else { &[] },
            if times_kept { &self.species } else { &[] },
        )?;
        if !times_kept || spec.noise != self.spec.noise {
            self.noise_tic = noise_tic(spec);
        }
        self.species = species;
        self.layout = layout;
        self.spec = spec.clone();
        Ok(())
    }

    /// The TIC per scan, and the XIC of proteoform `selected` (an index into `proteoforms`).
    pub fn trace(&self, selected: Option<usize>) -> RunTic {
        let n = self.spec.scan_times.len();
        let scale = self.floor_scale();
        let mut signal = vec![0.0; n];
        for s in &self.species {
            for (total, &w) in signal.iter_mut().zip(s.weights.iter()) {
                *total += s.envelope.survival.intensity(w, scale);
            }
        }

        // Background species follow the proteoforms in `species`; only a proteoform is selectable.
        let selected = selected.filter(|&i| i < self.spec.proteoforms.len()).map(|i| {
            let s = &self.species[i];
            let e = &s.envelope;
            let keep = |k: usize| e.unit.intensity[k] >= XIC_TARGET_FRACTION * e.tallest;
            let targets = Survival::new(&e.unit, self.layout.shape, keep);
            ProteoformXic {
                intensity: s.weights.iter().map(|&w| targets.intensity(w, scale)).collect(),
                target_mzs: (0..e.unit.mz.len()).filter(|&k| keep(k)).map(|k| e.unit.mz[k]).collect(),
            }
        });

        RunTic { signal, noise: self.noise_tic.clone(), selected }
    }

    /// The `scale` of the floor's `scale · shape`.
    fn floor_scale(&self) -> f64 {
        match &self.spec.noise {
            Some(n) => NoiseFloorModel::from(&n.floor).detection_floor_scale(),
            None => {
                let r: SimReduction = self.spec.reduction.unwrap_or_default();
                let brightest = self
                    .species
                    .iter()
                    .flat_map(|s| s.weights.iter().map(move |&w| w * s.envelope.tallest))
                    .fold(0.0f64, f64::max);
                r.absolute_intensity_floor.max(brightest * r.relative_intensity_threshold)
            }
        }
    }
}

fn layout_of(spec: &SimSpec) -> Layout {
    Layout {
        min_charge: spec.charge_range.min,
        max_charge: spec.charge_range.max,
        width: PeakWidthModel::from(spec.peak_width_model),
        shape: if spec.noise.is_some() { FloorShape::NoiseLevel } else { FloorShape::Flat },
    }
}

/// Every species of `spec`, reusing the envelopes in `same_layout` (built under this layout) and
/// the weights in `same_times` (on these scan times) whose inputs are unchanged. Matched by value,
/// not position, so adding, removing or duplicating a proteoform builds only that one.
fn build_species(
    spec: &SimSpec,
    layout: &Layout,
    same_layout: &[Species],
    same_times: &[Species],
) -> Result<Vec<Species>, SimError> {
    let envelopes: HashMap<EnvelopeKey, &Arc<Envelope>> =
        same_layout.iter().map(|s| (envelope_key(&s.model), &s.envelope)).collect();
    let weights: HashMap<WeightKey, &Arc<Vec<f64>>> =
        same_times.iter().map(|s| (weight_key(&s.model), &s.weights)).collect();

    signal_models(spec)
        .into_par_iter()
        .map(|model| {
            let envelope = match envelopes.get(&envelope_key(&model)) {
                Some(&e) => e.clone(),
                None => Arc::new(build_envelope(&model, layout)?),
            };
            let weights = match weights.get(&weight_key(&model)) {
                Some(&w) => w.clone(),
                None => Arc::new(weights_for(&model, &spec.scan_times)),
            };
            Ok(Species { model, envelope, weights })
        })
        .collect()
}

fn weights_for(model: &ProteoformModel, times: &[f64]) -> Vec<f64> {
    times
        .iter()
        .map(|&t| {
            let g = model.rt.evaluate(t);
            if g > 0.0 { model.abundance * g } else { 0.0 }
        })
        .collect()
}

/// The proteoform's centroids alone at unit weight: simulated at an RT where `A·g(t)` is 1.
fn build_envelope(model: &ProteoformModel, layout: &Layout) -> Result<Envelope, SimError> {
    let unit_rt = EmgProfile { mu: 0.0, sigma: 1.0, tau: 0.0 };
    let peak = unit_rt.evaluate(0.0);
    let alone = ProteoformModel { abundance: 1.0 / peak, rt: unit_rt, ..model.clone() };
    let centroider =
        ProfileCentroider::new(vec![alone], layout.min_charge, layout.max_charge, layout.width)?;
    let unit = centroider.centroid(0.0, Parallelism::Serial);
    let survival = Survival::new(&unit, layout.shape, |_| true);
    let tallest = unit.intensity.iter().copied().fold(0.0f64, f64::max);
    Ok(Envelope { unit, survival, tallest })
}

/// The noise floor's share of each scan's TIC: the draws the full pipeline makes for that scan.
fn noise_tic(spec: &SimSpec) -> Vec<f64> {
    let n = spec.scan_times.len();
    match &spec.noise {
        None => vec![0.0; n],
        Some(SimNoise { seed, floor, .. }) => {
            let model = NoiseFloorModel::from(floor);
            let seed = *seed as u64;
            (0..n)
                .into_par_iter()
                .map(|s| model.sample_scan_tic(&mut stream_for_scan(seed, s as u32, NOISE_STREAM)))
                .collect()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit(mz: &[f64], h: &[f64]) -> Centroids {
        Centroids { mz: mz.to_vec(), intensity: h.to_vec() }
    }

    #[test]
    fn survival_keeps_the_centroids_that_clear_the_floor() {
        let s = Survival::new(&unit(&[700.0, 800.0, 900.0], &[1.0, 4.0, 2.0]), FloorShape::Flat, |_| true);
        // Floor 3 at weight 1: only the height-4 centroid clears it.
        assert_eq!(s.intensity(1.0, 3.0), 4.0);
        // At weight 2 the heights are 2, 8, 4: all but the first clear 3.
        assert_eq!(s.intensity(2.0, 3.0), 12.0);
        // Exactly on the floor counts, as in the full pipeline's threshold.
        assert_eq!(s.intensity(1.5, 3.0), 1.5 * 6.0);
        assert_eq!(s.intensity(10.0, 3.0), 70.0);
        assert_eq!(s.intensity(0.0, 3.0), 0.0);
    }
}
