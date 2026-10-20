//! Forward model → centroid spectra. Port of `Simulation/ProfileCentroider.cs`.
//!
//! One centroid per local maximum of the summed profile, at that maximum's height — what an
//! instrument reports. Per scan, the active (proteoform, charge, isotopologue) Gaussians are
//! sorted by centre and split into clusters wherever neighbours are more than
//! [`INTERACTION_WINDOW_IN_SIGMAS`] σ apart; each cluster is walked on a σ/4 grid and every grid
//! maximum is Newton-polished on the analytic profile.

use std::f64::consts::PI;

use rayon::prelude::*;

use crate::model::{to_mz, IsotopeEnvelope, PeakWidthModel, ProteoformModel};
use crate::SimError;

/// Gaussians further apart than this cannot shift each other's maxima measurably.
const INTERACTION_WINDOW_IN_SIGMAS: f64 = 6.0;
/// Two resolvable maxima are at least ~2σ apart, so σ/4 cannot skip one.
const GRID_POINTS_PER_SIGMA: f64 = 4.0;
/// Same cut-off the C# forward model applies to a species' weight in a scan.
const MINIMUM_CONTRIBUTION_WEIGHT: f64 = 1e-18;
const MAX_NEWTON_ITERATIONS: usize = 20;

/// [`ProfileCentroider::centroid_above`]: Gaussians shorter than this fraction of the floor are
/// dropped. At 1e-6 even a peak sitting right on the floor moves by far less than the contract's
/// 1e-4 relative tolerance.
pub const PRUNE_FRACTION_OF_FLOOR: f64 = 1e-6;

/// [`ProfileCentroider::centroid_above`]: the grid walk skips stretches whose reachable Gaussians
/// sum to less than this fraction of the floor. Below 1 so that a maximum just above the floor,
/// whose grid neighbours sit a few percent lower, is never skipped.
pub const SKIP_FRACTION_OF_FLOOR: f64 = 0.5;

/// One (proteoform, charge) envelope's Gaussians.
#[derive(Debug, Clone)]
struct Envelope {
    charge_weight: f64,
    centers: Vec<f64>,
    unit_heights: Vec<f64>,
    inv_2_sigma2: Vec<f64>,
    sigmas: Vec<f64>,
    /// Largest entry of `unit_heights`.
    max_unit_height: f64,
}

#[derive(Debug, Clone, Copy)]
struct Component {
    center: f64,
    height: f64,
    inv_2_sigma2: f64,
    sigma: f64,
}

/// How to spend threads on one scan. The output is identical either way.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Parallelism {
    /// Walk the scan's clusters in parallel — for a single scan, where it is the only work.
    AcrossClusters,
    /// One thread — for a scan inside a run that is already parallel across scans, where nested
    /// parallelism costs more than it gains.
    Serial,
}

/// A centroided spectrum, ascending m/z.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Centroids {
    pub mz: Vec<f64>,
    pub intensity: Vec<f64>,
}

/// The charge states envelopes are built over, and the peak width that shapes every Gaussian —
/// the run-level inputs every envelope depends on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnvelopeLayout {
    pub min_charge: i32,
    pub max_charge: i32,
    pub width: PeakWidthModel,
}

/// One proteoform's share of the precomputed state.
#[derive(Debug, Clone)]
struct Species {
    model: ProteoformModel,
    /// Depends only on the mass; kept so that a layout change does not redo averagine.
    kernel: IsotopeEnvelope,
    /// Ascending charge.
    envelopes: Vec<Envelope>,
    /// The tallest component this proteoform can contribute per unit of `A_p·g_p(t)`: max over
    /// its envelopes of `f(z) · unit height`. Lets it be skipped whole when even its tallest
    /// Gaussian is below the pruning threshold.
    max_unit_component: f64,
}

impl Species {
    fn build(model: ProteoformModel, kernel: IsotopeEnvelope, layout: &EnvelopeLayout) -> Result<Self, SimError> {
        let envelopes = envelopes_for(&model, &kernel, layout)?;
        let max_unit_component =
            envelopes.iter().fold(0.0f64, |m, e| m.max(e.charge_weight * e.max_unit_height));
        Ok(Self { model, kernel, envelopes, max_unit_component })
    }

    /// The tallest Gaussian this proteoform puts in a scan at `scan_time`, with the same cut-off
    /// [`ProfileCentroider::max_component_height`] applies.
    fn tallest_component(&self, scan_time: f64) -> f64 {
        tallest_component(&self.model, self.max_unit_component, scan_time)
    }

    fn reach(&self) -> SpeciesReach {
        SpeciesReach { model: self.model.clone(), max_unit_component: self.max_unit_component }
    }
}

/// Just enough of a [`Species`] to bound what it adds to a scan.
#[derive(Debug, Clone)]
struct SpeciesReach {
    model: ProteoformModel,
    max_unit_component: f64,
}

impl SpeciesReach {
    fn tallest_component(&self, scan_time: f64) -> f64 {
        tallest_component(&self.model, self.max_unit_component, scan_time)
    }
}

#[inline]
fn tallest_component(model: &ProteoformModel, max_unit_component: f64, scan_time: f64) -> f64 {
    let w = rt_weight(model, scan_time);
    if w > MINIMUM_CONTRIBUTION_WEIGHT { w * max_unit_component } else { 0.0 }
}

/// `A_p · g_p(t)`.
#[inline]
fn rt_weight(p: &ProteoformModel, scan_time: f64) -> f64 {
    let rt = p.rt.evaluate(scan_time);
    if rt > 0.0 { p.abundance * rt } else { 0.0 }
}

/// What an [`ProfileCentroider::update`] changed, for invalidating anything cached per scan.
#[derive(Debug, Clone, Default)]
pub struct SignalChange {
    /// Every envelope was rebuilt: nothing computed from the old signal can be kept.
    pub everything: bool,
    /// The proteoforms the update replaced, removed or added, as they were before and are after.
    before: Vec<SpeciesReach>,
    after: Vec<SpeciesReach>,
}

impl SignalChange {
    pub fn is_empty(&self) -> bool {
        !self.everything && self.before.is_empty() && self.after.is_empty()
    }

    /// The tallest Gaussian any changed proteoform contributed before, or contributes after, at
    /// `scan_time`. Meaningless when [`Self::everything`] is set.
    pub fn tallest_changed_component(&self, scan_time: f64) -> f64 {
        self.before
            .iter()
            .chain(&self.after)
            .map(|s| s.tallest_component(scan_time))
            .fold(0.0, f64::max)
    }
}

/// The brightest centroid of one scan, and the tallest single Gaussian under it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanMax {
    pub brightest: f64,
    pub tallest_component: f64,
}

impl ScanMax {
    /// Whether this scan's value is unaffected by a change whose tallest changed component in
    /// the scan is `changed`. Below the prune threshold [`ProfileCentroider::scan_max`] used, the
    /// changed proteoforms were and are dropped from the scan whole, so its component set — and
    /// hence its tallest component and brightest centroid — is the same either way.
    pub fn survives(&self, changed: f64) -> bool {
        changed < PRUNE_FRACTION_OF_FLOOR * self.tallest_component
    }
}

/// Precomputed envelopes for a set of proteoforms. Read-only while simulating, so scans can be
/// centroided in parallel against one instance; [`Self::update`] edits it between simulations.
#[derive(Debug, Clone)]
pub struct ProfileCentroider {
    species: Vec<Species>,
    layout: EnvelopeLayout,
}

impl ProfileCentroider {
    pub fn new(
        proteoforms: Vec<ProteoformModel>,
        min_charge: i32,
        max_charge: i32,
        width: PeakWidthModel,
    ) -> Result<Self, SimError> {
        let layout = EnvelopeLayout { min_charge, max_charge, width };
        check_layout(&layout)?;
        Ok(Self { species: build_species(proteoforms, &layout)?, layout })
    }

    pub fn width(&self) -> PeakWidthModel {
        self.layout.width
    }

    /// Brings the centroider to `proteoforms` under `layout`, rebuilding only the envelopes the
    /// difference touches. The result is identical to [`Self::new`] on the same inputs. On error
    /// nothing has changed.
    pub fn update(
        &mut self,
        proteoforms: Vec<ProteoformModel>,
        layout: EnvelopeLayout,
    ) -> Result<SignalChange, SimError> {
        check_layout(&layout)?;

        if layout != self.layout {
            // Every envelope depends on the layout; only the averagine kernels carry over.
            let kernels: Vec<IsotopeEnvelope> = proteoforms
                .par_iter()
                .enumerate()
                .map(|(i, m)| match self.species.get(i) {
                    Some(s) if s.model.monoisotopic_mass == m.monoisotopic_mass => s.kernel.clone(),
                    _ => IsotopeEnvelope::new(m.monoisotopic_mass),
                })
                .collect();
            let species = proteoforms
                .into_par_iter()
                .zip(kernels)
                .map(|(m, k)| Species::build(m, k, &layout))
                .collect::<Result<_, _>>()?;
            self.species = species;
            self.layout = layout;
            return Ok(SignalChange { everything: true, ..Default::default() });
        }

        // Unchanged proteoforms at both ends; the middle of the old list becomes the middle of
        // the new one.
        let old = &self.species;
        let same = |s: &Species, m: &ProteoformModel| s.model == *m;
        let prefix = old.iter().zip(&proteoforms).take_while(|(s, m)| same(s, m)).count();
        let max_suffix = old.len().min(proteoforms.len()) - prefix;
        let suffix = old
            .iter()
            .rev()
            .zip(proteoforms.iter().rev())
            .take(max_suffix)
            .take_while(|(s, m)| same(s, m))
            .count();
        let old_middle = prefix..old.len() - suffix;
        let new_middle_len = proteoforms.len() - suffix - prefix;
        let paired = old_middle.len() == new_middle_len;

        // Build every replacement before committing any, so a failure leaves `self` untouched.
        let replacements: Vec<Option<Species>> = proteoforms[prefix..prefix + new_middle_len]
            .par_iter()
            .enumerate()
            .map(|(i, m)| {
                let previous = paired.then(|| &old[prefix + i]);
                match previous {
                    Some(s) if s.model == *m => Ok(None),
                    // Envelopes depend on mass and charge distribution only.
                    Some(s) if s.model.monoisotopic_mass == m.monoisotopic_mass
                        && s.model.charge == m.charge =>
                    {
                        Ok(Some(Species { model: m.clone(), ..s.clone() }))
                    }
                    Some(s) if s.model.monoisotopic_mass == m.monoisotopic_mass => {
                        Species::build(m.clone(), s.kernel.clone(), &self.layout).map(Some)
                    }
                    _ => Species::build(m.clone(), IsotopeEnvelope::new(m.monoisotopic_mass), &self.layout)
                        .map(Some),
                }
            })
            .collect::<Result<_, _>>()?;

        let mut change = SignalChange::default();
        if paired {
            for (i, r) in replacements.into_iter().enumerate() {
                if let Some(new) = r {
                    change.after.push(new.reach());
                    change.before.push(std::mem::replace(&mut self.species[prefix + i], new).reach());
                }
            }
        } else {
            let added: Vec<Species> =
                replacements.into_iter().map(|r| r.expect("an unpaired proteoform is always built")).collect();
            change.after = added.iter().map(Species::reach).collect();
            change.before = self.species.splice(old_middle, added).map(|s| s.reach()).collect();
        }
        Ok(change)
    }

    /// The centroids of the summed profile at one retention time, ascending m/z. Exact: every
    /// component contributes, as in C# `ProfileCentroider`.
    pub fn centroid(&self, scan_time: f64, parallelism: Parallelism) -> Centroids {
        self.centroid_above(scan_time, 0.0, parallelism)
    }

    /// The centroids that can survive a detection floor of at least `floor_lower_bound`, found
    /// without doing the work for the ones that cannot. 0 is [`Self::centroid`].
    ///
    /// The forward model renders every isotopologue of every charge of every proteoform whose
    /// EMG is above 1e-18 — through the EMG tails, nearly all of them — and the floor later
    /// discards most of the centroids. Two savings, both relative to the floor:
    ///
    /// - Gaussians shorter than [`PRUNE_FRACTION_OF_FLOOR`] of it are dropped. Each can move a
    ///   surviving peak by at most its own height, far inside the contract's 1e-4 tolerance.
    /// - Grid points where every Gaussian within reach sums to less than
    ///   [`SKIP_FRACTION_OF_FLOOR`] of it are not evaluated, since the profile there cannot be
    ///   that tall.
    ///
    /// The interaction reach widens to match. The exact path treats Gaussians more than 6σ apart
    /// as independent, which holds when tiny Gaussians chain the run into large clusters but not
    /// once pruning has removed them: a giant's tail 6σ out is still `e^-18` of the giant, which
    /// can exceed a nearby small peak. So clusters split, and evaluation truncates, only where
    /// every Gaussian that could cross has fallen below the pruning threshold — judged locally,
    /// so that one giant does not widen the reach everywhere in the scan.
    ///
    /// Centroids shorter than the skip threshold may or may not be reported; the floor removes
    /// them either way.
    pub fn centroid_above(
        &self,
        scan_time: f64,
        floor_lower_bound: f64,
        parallelism: Parallelism,
    ) -> Centroids {
        let prune_below = PRUNE_FRACTION_OF_FLOOR * floor_lower_bound;
        let components = self.active_components(scan_time, prune_below);
        let skip_below = SKIP_FRACTION_OF_FLOOR * floor_lower_bound;

        // Distance, in σ, at which a Gaussian of height `h` falls to `prune_below`.
        let reach_for = |h: f64| -> f64 {
            if prune_below > 0.0 && h > prune_below {
                INTERACTION_WINDOW_IN_SIGMAS.max((2.0 * (h / prune_below).ln()).sqrt())
            } else {
                INTERACTION_WINDOW_IN_SIGMAS
            }
        };
        // Beyond this nothing in the scan matters, so it bounds every local search.
        let global_reach = reach_for(components.iter().fold(0.0f64, |m, c| m.max(c.height)));
        // σ grows with m/z, so a Gaussian across a gap can be a little wider than either
        // neighbour; over a reach it is far less than this margin.
        const SIGMA_MARGIN: f64 = 1.01;

        let mut clusters = Vec::new();
        let (mut start, mut cluster_tallest) = (0, 0.0f64);
        for k in 1..=components.len() {
            cluster_tallest = cluster_tallest.max(components[k - 1].height);
            let split = k == components.len() || {
                let (left, right) = (&components[k - 1], &components[k]);
                let gap = right.center - left.center;
                let neighbour_sigma = left.sigma.max(right.sigma);
                let sigma = neighbour_sigma * SIGMA_MARGIN;
                if gap <= INTERACTION_WINDOW_IN_SIGMAS * neighbour_sigma {
                    false
                } else if prune_below == 0.0 {
                    // The exact path: C#'s rule, nothing more.
                    true
                } else if gap > global_reach * sigma {
                    true
                } else {
                    // The tallest Gaussian close enough to the gap to reach across it.
                    let d = global_reach * sigma;
                    let a = components.partition_point(|c| c.center < left.center - d);
                    let b = components.partition_point(|c| c.center <= right.center + d);
                    let nearby = components[a..b].iter().fold(0.0f64, |m, c| m.max(c.height));
                    gap > reach_for(nearby) * sigma
                }
            };
            if split {
                clusters.push((start..k, reach_for(cluster_tallest)));
                start = k;
                cluster_tallest = 0.0;
            }
        }

        // Clusters are independent (nothing carries from one to the next), so one scan's
        // clusters can be walked in parallel; concatenating in order keeps m/z ascending.
        let walk = |(range, reach): (std::ops::Range<usize>, f64)| {
            let mut part = Centroids::default();
            self.pick_cluster(&components[range], reach, skip_below, &mut part);
            part
        };
        let parts: Vec<Centroids> = match parallelism {
            Parallelism::AcrossClusters => clusters.into_par_iter().with_min_len(64).map(walk).collect(),
            Parallelism::Serial => clusters.into_iter().map(walk).collect(),
        };

        let total = parts.iter().map(|p| p.mz.len()).sum();
        let mut out = Centroids { mz: Vec::with_capacity(total), intensity: Vec::with_capacity(total) };
        for part in parts {
            out.mz.extend(part.mz);
            out.intensity.extend(part.intensity);
        }
        out
    }

    /// The tallest single Gaussian at this retention time. A lower bound on the tallest centroid,
    /// since the summed profile is at least as tall as any one of its (non-negative) terms.
    pub fn max_component_height(&self, scan_time: f64) -> f64 {
        self.species.iter().map(|s| s.tallest_component(scan_time)).fold(0.0, f64::max)
    }

    /// The brightest centroid at this retention time, found with the scan's tallest Gaussian `T`
    /// as the floor bound: the brightest centroid is at least `T`, so [`Self::centroid_above`]
    /// always reports it, while everything below `PRUNE_FRACTION_OF_FLOOR · T` is pruned and only
    /// the neighbourhoods of the tallest Gaussians are walked.
    ///
    /// Depends on the signal model alone (not on noise or reduction), so it can be cached per
    /// scan across those edits; [`ScanMax::survives`] says when a signal edit leaves it valid.
    pub fn scan_max(&self, scan_time: f64, parallelism: Parallelism) -> ScanMax {
        let tallest_component = self.max_component_height(scan_time);
        let brightest = if tallest_component > 0.0 {
            self.centroid_above(scan_time, tallest_component, parallelism)
                .intensity
                .iter()
                .fold(0.0, |max, &y| if y > max { y } else { max })
        } else {
            0.0
        };
        ScanMax { brightest, tallest_component }
    }

    fn active_components(&self, scan_time: f64, min_height: f64) -> Vec<Component> {
        let mut components = Vec::new();
        for s in &self.species {
            let mut weight = rt_weight(&s.model, scan_time);
            if min_height > 0.0 && weight * s.max_unit_component < min_height {
                weight = 0.0;
            }
            if weight <= MINIMUM_CONTRIBUTION_WEIGHT {
                continue;
            }

            for e in &s.envelopes {
                let weight = weight * e.charge_weight;
                if weight <= MINIMUM_CONTRIBUTION_WEIGHT || weight * e.max_unit_height < min_height {
                    continue;
                }

                for i in 0..e.centers.len() {
                    let height = weight * e.unit_heights[i];
                    if height > 0.0 && height >= min_height {
                        components.push(Component {
                            center: e.centers[i],
                            height,
                            inv_2_sigma2: e.inv_2_sigma2[i],
                            sigma: e.sigmas[i],
                        });
                    }
                }
            }
        }

        components.sort_by(|a, b| a.center.total_cmp(&b.center));
        components
    }

    /// Finds the maxima of the profile made by one cluster of components, evaluating each point
    /// from the components within `reach` × the cluster's widest σ. With `skip_below` > 0, only
    /// maxima at least that tall are guaranteed (see [`Self::centroid_above`]).
    fn pick_cluster(&self, cluster: &[Component], reach: f64, skip_below: f64, out: &mut Centroids) {
        // A lone Gaussian peaks exactly at its centre; no grid needed.
        if cluster.len() == 1 {
            if cluster[0].height >= skip_below {
                out.mz.push(cluster[0].center);
                out.intensity.push(cluster[0].height);
            }
            return;
        }

        let max_sigma = cluster.iter().fold(0.0f64, |m, c| m.max(c.sigma));
        let window = reach * max_sigma;

        let first = cluster[0];
        let last = cluster[cluster.len() - 1];
        let lo = first.center - first.sigma;
        let hi = last.center + last.sigma;

        // cluster[leave..enter] is exactly the set `evaluate` would find by binary search at x
        // (centres within x ± window); x only increases, so it is maintained by sliding.
        // `window_sum`, the sum of their heights, bounds the profile at x from above.
        let skipping = skip_below > 0.0;
        let (mut enter, mut leave, mut window_sum) = (0usize, 0usize, 0.0f64);

        let (mut prev2, mut prev1, mut prev_x) = (f64::NAN, f64::NAN, f64::NAN);
        let mut last_accepted_mz = f64::NEG_INFINITY;
        let mut x = lo;
        loop {
            while enter < cluster.len() && cluster[enter].center <= x + window {
                window_sum += cluster[enter].height;
                enter += 1;
            }
            while leave < enter && cluster[leave].center < x - window {
                window_sum -= cluster[leave].height;
                leave += 1;
            }

            let y = if skipping && window_sum < skip_below {
                // Below anything reportable; 0 keeps the maximum test below correct for every
                // maximum that is. The walk still takes this grid step rather than jumping
                // ahead: where Newton cannot improve on a grid maximum (a flat top or shoulder)
                // the grid point *is* the reported m/z, so the grid must stay aligned with the
                // exact path's.
                0.0
            } else {
                profile_at(&cluster[leave..enter], x)
            };

            // prev1 is a maximum when it rises from prev2 and does not fall short of y.
            if !prev2.is_nan() && prev1 > prev2 && prev1 >= y {
                let step = self.layout.width.sigma_at(prev_x) / GRID_POINTS_PER_SIGMA;
                let (peak_mz, peak_height) = polish(cluster, prev_x, window, step);
                if peak_height < skip_below {
                    // Possibly an artefact of a skipped neighbour read as 0, and below the floor
                    // regardless.
                } else if peak_mz - last_accepted_mz > 1e-9 * peak_mz {
                    out.mz.push(peak_mz);
                    out.intensity.push(peak_height);
                    last_accepted_mz = peak_mz;
                } else if let Some(last) = out.intensity.last_mut() {
                    if peak_height > *last {
                        *last = peak_height;
                    }
                }
            }

            if x >= hi {
                break;
            }
            prev2 = prev1;
            prev1 = y;
            prev_x = x;
            x = hi.min(x + self.layout.width.sigma_at(x) / GRID_POINTS_PER_SIGMA);
        }
    }
}

fn check_layout(layout: &EnvelopeLayout) -> Result<(), SimError> {
    if layout.min_charge < 1 || layout.max_charge < layout.min_charge {
        return Err(SimError::InvalidSpec("charge range must satisfy 1 <= min <= max".into()));
    }
    Ok(())
}

/// Proteoforms are independent, so they are built in parallel; collecting in order keeps the
/// result identical to a serial build.
fn build_species(proteoforms: Vec<ProteoformModel>, layout: &EnvelopeLayout) -> Result<Vec<Species>, SimError> {
    proteoforms
        .into_par_iter()
        .map(|m| {
            let kernel = IsotopeEnvelope::new(m.monoisotopic_mass);
            Species::build(m, kernel, layout)
        })
        .collect()
}

/// One proteoform's (charge) envelopes, ascending charge.
fn envelopes_for(
    model: &ProteoformModel,
    kernel: &IsotopeEnvelope,
    layout: &EnvelopeLayout,
) -> Result<Vec<Envelope>, SimError> {
    let EnvelopeLayout { min_charge, max_charge, width } = *layout;
    let n = kernel.len();
    let mut envelopes = Vec::new();

    for z in min_charge..=max_charge {
        let fz = model.charge.evaluate(z);
        if !(fz > 0.0) || n == 0 {
            continue;
        }

        let mut centers = Vec::with_capacity(n);
        let mut heights = Vec::with_capacity(n);
        let mut inv_2_sigma2 = Vec::with_capacity(n);
        let mut sigmas = Vec::with_capacity(n);
        for i in 0..n {
            let center = to_mz(kernel.neutral_masses[i], z);
            let sigma = width.sigma_at(center);
            if !(sigma > 0.0) || !sigma.is_finite() {
                return Err(SimError::InvalidSpec(format!(
                    "σ_m must be finite and positive; got {sigma} at m/z {center}"
                )));
            }
            centers.push(center);
            sigmas.push(sigma);
            // Unit-area Gaussians weighted by the normalized isotopologue abundance.
            heights.push(kernel.intensities[i] / (sigma * (2.0 * PI).sqrt()));
            inv_2_sigma2.push(1.0 / (2.0 * sigma * sigma));
        }

        let max_unit_height = heights.iter().copied().fold(0.0f64, f64::max);
        envelopes.push(Envelope {
            charge_weight: fz,
            centers,
            unit_heights: heights,
            inv_2_sigma2,
            sigmas,
            max_unit_height,
        });
    }

    Ok(envelopes)
}

/// C# `Math.Clamp`.
#[inline]
fn clamp(v: f64, min: f64, max: f64) -> f64 {
    if v < min {
        min
    } else if v > max {
        max
    } else {
        v
    }
}

/// Newton's method on S'(x) = 0 from a grid maximum, each step held inside one grid cell so it
/// cannot jump to a neighbouring maximum.
fn polish(cluster: &[Component], x0: f64, window: f64, step: f64) -> (f64, f64) {
    let mut x = x0;
    for _ in 0..MAX_NEWTON_ITERATIONS {
        let (_, d1, d2) = evaluate(cluster, x, window);
        if !(d2 < 0.0) {
            break;
        }

        let delta = clamp(-d1 / d2, -step, step);
        let next = clamp(x + delta, x0 - step, x0 + step);
        if (next - x).abs() < 1e-12 * x {
            x = next;
            break;
        }
        x = next;
    }

    (x, evaluate(cluster, x, window).0)
}

/// The profile at `x` from exactly the components in `within_window`, summed in the same order
/// as [`evaluate`] — so bit-identical to its first output, without the derivatives.
#[inline]
fn profile_at(within_window: &[Component], x: f64) -> f64 {
    let mut s = 0.0;
    for c in within_window {
        let d = x - c.center;
        s += c.height * (-d * d * c.inv_2_sigma2).exp();
    }
    s
}

/// The profile and its first two derivatives at `x`.
fn evaluate(cluster: &[Component], x: f64, window: f64) -> (f64, f64, f64) {
    let start = cluster.partition_point(|c| c.center < x - window);
    let (mut s, mut d1, mut d2) = (0.0, 0.0, 0.0);
    for c in &cluster[start..] {
        if !(c.center <= x + window) {
            break;
        }
        let d = x - c.center;
        let g = c.height * (-d * d * c.inv_2_sigma2).exp();
        let two_q = 2.0 * c.inv_2_sigma2;
        s += g;
        d1 -= g * two_q * d;
        d2 += g * (two_q * two_q * d * d - two_q);
    }
    (s, d1, d2)
}
