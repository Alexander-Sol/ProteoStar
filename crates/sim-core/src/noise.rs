//! The noise subsystem: sampled FT noise floor, per-peak jitter, and unresolved-peak collapse.
//!
//! Models ported from `Noise/NoiseFloorModel.cs`, `SignalToNoiseDistribution.cs`,
//! `PeakJitterModel.cs`, `NoiseRelativeFloor.cs` and `NoiseInjector.cs` — same distributions and
//! calibration tables, but **not** the same random streams: this uses a native xoshiro256++
//! generator, so output matches C# statistically rather than draw for draw.

use std::f64::consts::PI;

use crate::centroider::Centroids;
use crate::model::PeakWidthModel;
use crate::spec::{SimNoiseFloor, SimNoiseJitter};

#[inline]
fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// xoshiro256++ with a cached second Box–Muller variate.
#[derive(Clone)]
pub struct Rng {
    s: [u64; 4],
    spare_normal: Option<f64>,
}

impl Rng {
    pub fn from_seed(seed: u64) -> Self {
        let mut z = seed;
        let mut s = [0u64; 4];
        for word in &mut s {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            *word = splitmix64(z);
        }
        Self { s, spare_normal: None }
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let s = &mut self.s;
        let result = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// Uniform in [0, 1).
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Standard normal, Box–Muller, both variates used.
    #[inline]
    pub fn normal(&mut self) -> f64 {
        if let Some(z) = self.spare_normal.take() {
            return z;
        }
        // 1 - u is in (0, 1], keeping the log finite.
        let r = (-2.0 * (1.0 - self.next_f64()).ln()).sqrt();
        let (sin, cos) = (2.0 * PI * self.next_f64()).sin_cos();
        self.spare_normal = Some(r * sin);
        r * cos
    }

    /// Exponential with mean 1.
    #[inline]
    fn exp1(&mut self) -> f64 {
        -(1.0 - self.next_f64()).ln()
    }
}

/// Stream ids, so jitter and noise draw from independent sequences: enabling jitter must not
/// change the noise floor a seed produces.
pub const NOISE_STREAM: u32 = 0;
pub const JITTER_STREAM: u32 = 1;

/// An independent stream per (seed, scan, stream), so output does not depend on which scans are
/// simulated or in what order.
pub fn stream_for_scan(seed: u64, scan_index: u32, stream: u32) -> Rng {
    let key = ((scan_index as u64) << 1) | stream as u64;
    Rng::from_seed(splitmix64(seed ^ splitmix64(key)))
}

/// Knuth's product method below λ=30, a rounded normal approximation above.
fn sample_poisson(rng: &mut Rng, lambda: f64) -> usize {
    if lambda < 30.0 {
        let limit = (-lambda).exp();
        let mut product = rng.next_f64();
        let mut count = 0;
        while product > limit {
            count += 1;
            product *= rng.next_f64();
        }
        return count;
    }
    (lambda + lambda.sqrt() * rng.normal()).round().max(0.0) as usize
}

const BIN_WIDTH: f64 = 100.0;
const FIRST_BIN_START: f64 = 600.0;
const FIRST_BIN_CENTRE: f64 = FIRST_BIN_START + 0.5 * BIN_WIDTH;

/// Measured noise peaks per scan in each 100-Th bin from m/z 600 (Jurkat quiet window).
const JURKAT_PEAKS_PER_SCAN_PER_BIN: [f64; 14] = [
    3527.0, 4303.0, 3504.0, 2519.0, 1337.0, 745.0, 383.0, 235.0, 141.0, 87.5, 50.8, 29.0, 16.4,
    8.2,
];

/// Noise amplitude in each bin relative to the first (rep2 fract7 wash).
const RELATIVE_NOISE_LEVEL_PER_BIN: [f64; 14] = [
    1.00, 1.25, 1.51, 1.76, 2.22, 2.30, 2.47, 2.44, 2.51, 2.48, 2.47, 2.50, 2.41, 2.41,
];

/// Noise amplitude at `mz` relative to the reference m/z: linear between bin centres, flat outside
/// the calibrated range.
#[inline]
pub fn relative_noise_level_at(mz: f64) -> f64 {
    let table = &RELATIVE_NOISE_LEVEL_PER_BIN;
    let last = table.len() - 1;
    // Bin centres are evenly spaced, so the bracketing pair is found by division.
    let x = (mz - FIRST_BIN_CENTRE) / BIN_WIDTH;
    if !(x > 0.0) {
        table[0]
    } else if x >= last as f64 {
        table[last]
    } else {
        let i = x as usize;
        let frac = x - i as f64;
        table[i] + frac * (table[i + 1] - table[i])
    }
}

/// Whether a peak is a forward-model isotopologue or an FT noise peak. Collapse treats them
/// differently.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeakOrigin {
    Signal,
    FtNoise,
}

#[derive(Debug, Clone, Copy)]
pub struct SimulatedPeak {
    pub mz: f64,
    pub intensity: f64,
    pub origin: PeakOrigin,
}

/// The unstructured FT noise floor. Port of `NoiseFloorModel` (the scalar-knob form; the
/// per-scan-conditioned tables are not v1).
#[derive(Debug, Clone)]
pub struct NoiseFloorModel {
    noise_level_at_reference_mz: f64,
    density_scale: f64,
    level_dispersion: f64,
    sn_threshold: f64,
    sn_mean_excess: f64,
}

impl From<&SimNoiseFloor> for NoiseFloorModel {
    fn from(f: &SimNoiseFloor) -> Self {
        Self {
            noise_level_at_reference_mz: f.noise_level_at_reference_mz,
            density_scale: f.density_scale,
            level_dispersion: f.level_dispersion,
            sn_threshold: f.signal_to_noise.threshold,
            sn_mean_excess: f.signal_to_noise.mean_excess,
        }
    }
}

impl NoiseFloorModel {
    /// Noise amplitude at `mz`: linear between bin centres, flat outside the calibrated range.
    #[inline]
    pub fn noise_level_at(&self, mz: f64) -> f64 {
        self.noise_level_at_reference_mz * relative_noise_level_at(mz)
    }

    /// The scale of [`Self::detection_floor_at`]: the floor at `mz` is this times
    /// [`relative_noise_level_at`]`(mz)`.
    pub fn detection_floor_scale(&self) -> f64 {
        self.sn_threshold * self.noise_level_at_reference_mz
    }

    /// The summed intensity of the noise peaks [`Self::sample_scan`] draws from `rng`: the noise
    /// floor's share of that scan's TIC (collapse merges peaks but keeps their sum).
    /// Makes exactly [`Self::sample_scan`]'s draws, without storing or sorting the peaks.
    pub fn sample_scan_tic(&self, rng: &mut Rng) -> f64 {
        let mut tic = 0.0;
        for (bin, &peaks) in JURKAT_PEAKS_PER_SCAN_PER_BIN.iter().enumerate() {
            let expected = peaks * self.density_scale;
            if expected <= 0.0 {
                continue;
            }
            let bin_start = FIRST_BIN_START + bin as f64 * BIN_WIDTH;
            for _ in 0..sample_poisson(rng, expected) {
                let mz = bin_start + rng.next_f64() * BIN_WIDTH;
                let mut level = self.noise_level_at(mz);
                if self.level_dispersion > 0.0 {
                    level *= (self.level_dispersion * rng.normal()).exp();
                }
                tic += level * (self.sn_threshold + self.sn_mean_excess * rng.exp1());
            }
        }
        tic
    }

    /// The detection limit a peak must clear: S/N threshold × the local noise level.
    #[inline]
    pub fn detection_floor_at(&self, mz: f64) -> f64 {
        self.sn_threshold * self.noise_level_at(mz)
    }

    /// The lowest [`Self::detection_floor_at`] anywhere in m/z.
    pub fn min_detection_floor(&self) -> f64 {
        let min_relative = RELATIVE_NOISE_LEVEL_PER_BIN.iter().copied().fold(f64::INFINITY, f64::min);
        self.sn_threshold * self.noise_level_at_reference_mz * min_relative
    }

    /// Draws one scan's noise peaks into `into`, ascending m/z: Poisson count per bin, uniform
    /// position within the bin, S/N = threshold + Exp(meanExcess) times a lognormal level.
    pub fn sample_scan(&self, rng: &mut Rng, into: &mut Vec<SimulatedPeak>) {
        for (bin, &peaks) in JURKAT_PEAKS_PER_SCAN_PER_BIN.iter().enumerate() {
            let expected = peaks * self.density_scale;
            if expected <= 0.0 {
                continue;
            }

            let count = sample_poisson(rng, expected);
            if count == 0 {
                continue;
            }

            let bin_start = FIRST_BIN_START + bin as f64 * BIN_WIDTH;
            let first = into.len();
            into.reserve(count);
            for _ in 0..count {
                let mz = bin_start + rng.next_f64() * BIN_WIDTH;
                let mut level = self.noise_level_at(mz);
                if self.level_dispersion > 0.0 {
                    level *= (self.level_dispersion * rng.normal()).exp();
                }
                let signal_to_noise = self.sn_threshold + self.sn_mean_excess * rng.exp1();
                into.push(SimulatedPeak {
                    mz,
                    intensity: level * signal_to_noise,
                    origin: PeakOrigin::FtNoise,
                });
            }

            // Bins are visited in ascending m/z and never overlap, so sorting each slice keeps
            // the whole list ascending.
            into[first..].sort_unstable_by(|a, b| a.mz.total_cmp(&b.mz));
        }
    }
}

/// Scan-to-scan measurement error on real peaks. Port of `PeakJitterModel`.
#[derive(Debug, Clone, Copy)]
pub struct PeakJitterModel {
    pub mz_ppm_at_unit_signal_to_noise: f64,
    pub mz_ppm_floor: f64,
    pub common_mode_mz_ppm: f64,
    pub log_intensity_sigma_at_unit_signal_to_noise: f64,
    pub log_intensity_sigma_floor: f64,
    pub common_mode_log_intensity_sigma: f64,
    pub drop_peaks_falling_below_floor: bool,
}

impl PeakJitterModel {
    /// The jitter a spec asks for, or `None` when it is disabled or would change nothing.
    pub fn from_spec(j: &SimNoiseJitter) -> Option<Self> {
        if !j.enabled {
            return None;
        }
        let m = Self {
            mz_ppm_at_unit_signal_to_noise: j.mz_ppm_at_unit_signal_to_noise,
            mz_ppm_floor: j.mz_ppm_floor,
            common_mode_mz_ppm: j.common_mode_mz_ppm,
            log_intensity_sigma_at_unit_signal_to_noise: j.log_intensity_sigma_at_unit_signal_to_noise,
            log_intensity_sigma_floor: j.log_intensity_sigma_floor,
            common_mode_log_intensity_sigma: j.common_mode_log_intensity_sigma,
            drop_peaks_falling_below_floor: j.drop_peaks_falling_below_floor,
        };
        m.is_active().then_some(m)
    }

    fn is_active(&self) -> bool {
        self.mz_ppm_at_unit_signal_to_noise > 0.0
            || self.mz_ppm_floor > 0.0
            || self.common_mode_mz_ppm > 0.0
            || self.log_intensity_sigma_at_unit_signal_to_noise > 0.0
            || self.log_intensity_sigma_floor > 0.0
            || self.common_mode_log_intensity_sigma > 0.0
    }

    /// Per-peak m/z σ in ppm: `sqrt(a²/(S/N) + b²)`.
    fn mz_ppm_sigma(&self, signal_to_noise: f64) -> f64 {
        let a = self.mz_ppm_at_unit_signal_to_noise;
        let b = self.mz_ppm_floor;
        if !(signal_to_noise > 0.0) || !signal_to_noise.is_finite() {
            return b;
        }
        (a * a / signal_to_noise + b * b).sqrt()
    }

    /// σ of log intensity: `sqrt(c²/(S/N) + d²)`.
    fn log_intensity_sigma(&self, signal_to_noise: f64) -> f64 {
        let c = self.log_intensity_sigma_at_unit_signal_to_noise;
        let d = self.log_intensity_sigma_floor;
        if !(signal_to_noise > 0.0) || !signal_to_noise.is_finite() {
            return d;
        }
        (c * c / signal_to_noise + d * d).sqrt()
    }
}

/// Mean-preserving multiplicative factor `exp(N(−σ²/2, σ))`.
#[inline]
fn intensity_factor(rng: &mut Rng, log_sigma: f64) -> f64 {
    if !(log_sigma > 0.0) {
        return 1.0;
    }
    (-0.5 * log_sigma * log_sigma + log_sigma * rng.normal()).exp()
}

/// Two centroids closer than this many σ are reported as one.
const MERGE_WITHIN_SIGMAS: f64 = 1.0;

/// Applies measurement error to one reduced scan: jitter the signal, add the sampled floor,
/// collapse what the instrument could not resolve.
pub fn inject_scan(
    signal: &Centroids,
    scan_index: usize,
    seed: u64,
    model: &NoiseFloorModel,
    jitter: Option<&PeakJitterModel>,
    width: &PeakWidthModel,
) -> Centroids {
    let mut signal_peaks = Vec::with_capacity(signal.mz.len());
    match jitter {
        None => signal_peaks.extend(signal.mz.iter().zip(&signal.intensity).map(
            |(&mz, &intensity)| SimulatedPeak { mz, intensity, origin: PeakOrigin::Signal },
        )),
        Some(j) => {
            let mut rng = stream_for_scan(seed, scan_index as u32, JITTER_STREAM);
            add_jittered_signal(signal, &mut rng, &mut signal_peaks, model, j);
            // Jitter moves peaks by a few ppm, which can swap close neighbours.
            signal_peaks.sort_unstable_by(|a, b| a.mz.total_cmp(&b.mz));
        }
    }

    let mut noise_peaks = Vec::new();
    let mut rng = stream_for_scan(seed, scan_index as u32, NOISE_STREAM);
    model.sample_scan(&mut rng, &mut noise_peaks);

    collapse(&merge_sorted(&signal_peaks, &noise_peaks), width)
}

/// Merges two m/z-ascending lists.
fn merge_sorted(a: &[SimulatedPeak], b: &[SimulatedPeak]) -> Vec<SimulatedPeak> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i].mz <= b[j].mz {
            out.push(a[i]);
            i += 1;
        } else {
            out.push(b[j]);
            j += 1;
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

fn add_jittered_signal(
    signal: &Centroids,
    rng: &mut Rng,
    into: &mut Vec<SimulatedPeak>,
    model: &NoiseFloorModel,
    jitter: &PeakJitterModel,
) {
    // Common-mode terms are drawn once per scan and shared by every peak in it.
    let common_mode_ppm = if jitter.common_mode_mz_ppm > 0.0 {
        jitter.common_mode_mz_ppm * rng.normal()
    } else {
        0.0
    };
    let common_mode_intensity = intensity_factor(rng, jitter.common_mode_log_intensity_sigma);

    for (&mz, &intensity) in signal.mz.iter().zip(&signal.intensity) {
        let noise_level = model.noise_level_at(mz);
        let signal_to_noise =
            if noise_level > 0.0 { intensity / noise_level } else { f64::INFINITY };

        let ppm = common_mode_ppm + jitter.mz_ppm_sigma(signal_to_noise) * rng.normal();
        let jittered_mz = mz * (1.0 + ppm * 1e-6);

        let jittered_intensity = intensity
            * common_mode_intensity
            * intensity_factor(rng, jitter.log_intensity_sigma(signal_to_noise));

        if jitter.drop_peaks_falling_below_floor
            && jittered_intensity < model.detection_floor_at(jittered_mz)
        {
            continue;
        }

        into.push(SimulatedPeak { mz: jittered_mz, intensity: jittered_intensity, origin: PeakOrigin::Signal });
    }
}

/// Collapses groups within [`MERGE_WITHIN_SIGMAS`] of their first member into one centroid at
/// the intensity-weighted m/z. Pure-noise groups are left alone: their density was measured on
/// already-reported centroids, so merging them again would double-count.
fn collapse(sorted: &[SimulatedPeak], width: &PeakWidthModel) -> Centroids {
    let mut out = Centroids {
        mz: Vec::with_capacity(sorted.len()),
        intensity: Vec::with_capacity(sorted.len()),
    };

    let mut i = 0;
    while i < sorted.len() {
        let group_start = sorted[i].mz;
        let tolerance = MERGE_WITHIN_SIGMAS * width.sigma_at(group_start);

        let mut j = i + 1;
        while j < sorted.len() && sorted[j].mz - group_start <= tolerance {
            j += 1;
        }

        let group = &sorted[i..j];
        if group.len() == 1 || !group.iter().any(|p| p.origin == PeakOrigin::Signal) {
            for p in group {
                out.mz.push(p.mz);
                out.intensity.push(p.intensity);
            }
            i = j;
            continue;
        }

        let (mut weighted_mz, mut summed) = (0.0, 0.0);
        for p in group {
            weighted_mz += p.mz * p.intensity;
            summed += p.intensity;
        }
        out.mz.push(if summed > 0.0 { weighted_mz / summed } else { group_start });
        out.intensity.push(summed);
        i = j;
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::SimNoiseFloor;

    fn mean_sd(v: &[f64]) -> (f64, f64) {
        let n = v.len() as f64;
        let mean = v.iter().sum::<f64>() / n;
        let var = v.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1.0);
        (mean, var.sqrt())
    }

    #[test]
    fn noise_level_interpolates_between_bin_centres_and_clamps() {
        let m = NoiseFloorModel::from(&SimNoiseFloor { noise_level_at_reference_mz: 1.0, ..Default::default() });
        assert_eq!(m.noise_level_at(100.0), 1.00);
        assert_eq!(m.noise_level_at(650.0), 1.00);
        assert!((m.noise_level_at(700.0) - 1.125).abs() < 1e-12);
        assert!((m.noise_level_at(750.0) - 1.25).abs() < 1e-12);
        assert_eq!(m.noise_level_at(1950.0), 2.41);
        assert_eq!(m.noise_level_at(5000.0), 2.41);
        assert_eq!(m.min_detection_floor(), 1.94);
    }

    #[test]
    fn rng_moments() {
        let mut rng = Rng::from_seed(42);
        let u: Vec<f64> = (0..200_000).map(|_| rng.next_f64()).collect();
        assert!(u.iter().all(|x| (0.0..1.0).contains(x)));
        let (mu, sd) = mean_sd(&u);
        assert!((mu - 0.5).abs() < 0.005 && (sd - (1.0f64 / 12.0).sqrt()).abs() < 0.005);

        let z: Vec<f64> = (0..200_000).map(|_| rng.normal()).collect();
        let (mu, sd) = mean_sd(&z);
        assert!(mu.abs() < 0.01 && (sd - 1.0).abs() < 0.01, "normal mean {mu} sd {sd}");
    }

    #[test]
    fn streams_are_deterministic_and_distinct() {
        let draw = |seed, scan, stream| stream_for_scan(seed, scan, stream).next_u64();
        assert_eq!(draw(7, 3, 0), draw(7, 3, 0));
        assert_ne!(draw(7, 3, 0), draw(7, 4, 0));
        assert_ne!(draw(7, 3, 0), draw(7, 3, 1));
        assert_ne!(draw(7, 3, 0), draw(8, 3, 0));
    }

    #[test]
    fn noise_floor_has_calibrated_density_and_signal_to_noise() {
        let floor = SimNoiseFloor { level_dispersion: 0.0, ..Default::default() };
        let model = NoiseFloorModel::from(&floor);
        let expected: f64 = JURKAT_PEAKS_PER_SCAN_PER_BIN.iter().sum();

        let mut counts = Vec::new();
        let mut sn = Vec::new();
        for s in 0..40 {
            let mut peaks = Vec::new();
            model.sample_scan(&mut stream_for_scan(1, s, NOISE_STREAM), &mut peaks);
            assert!(peaks.windows(2).all(|w| w[0].mz <= w[1].mz));
            counts.push(peaks.len() as f64);
            sn.extend(peaks.iter().map(|p| p.intensity / model.noise_level_at(p.mz)));
        }
        let (mean_count, _) = mean_sd(&counts);
        assert!((mean_count / expected - 1.0).abs() < 0.01, "{mean_count} vs {expected}");

        // S/N = threshold + Exp(meanExcess): mean threshold + meanExcess, min ≥ threshold.
        let (mean_sn, _) = mean_sd(&sn);
        assert!((mean_sn - (1.94 + 1.165)).abs() < 0.01, "mean S/N {mean_sn}");
        assert!(sn.iter().all(|&x| x >= 1.94 - 1e-9));
    }

    #[test]
    fn scan_tic_sums_the_sampled_peaks() {
        let model = NoiseFloorModel::from(&SimNoiseFloor::default());
        for s in 0..5 {
            let mut peaks = Vec::new();
            model.sample_scan(&mut stream_for_scan(9, s, NOISE_STREAM), &mut peaks);
            let summed: f64 = peaks.iter().map(|p| p.intensity).sum();
            let tic = model.sample_scan_tic(&mut stream_for_scan(9, s, NOISE_STREAM));
            assert!((tic / summed - 1.0).abs() < 1e-12, "scan {s}: {tic} vs {summed}");
        }
    }

    #[test]
    fn jitter_follows_its_laws() {
        // One bright peak, jitter only (no noise floor), many scans.
        let model = NoiseFloorModel::from(&SimNoiseFloor { density_scale: 0.0, ..Default::default() });
        let jitter = PeakJitterModel::from_spec(&SimNoiseJitter::default()).unwrap();
        let width = PeakWidthModel::Orbitrap { k: 1.2e-7 };
        let (mz0, y0) = (1000.0, 1e7);
        let signal = Centroids { mz: vec![mz0], intensity: vec![y0] };

        let (mut ppm, mut log_y, mut ys) = (Vec::new(), Vec::new(), Vec::new());
        for s in 0..20_000 {
            let out = inject_scan(&signal, s, 3, &model, Some(&jitter), &width);
            assert_eq!(out.mz.len(), 1);
            ppm.push((out.mz[0] / mz0 - 1.0) * 1e6);
            log_y.push((out.intensity[0] / y0).ln());
            ys.push(out.intensity[0] / y0);
        }

        let sn = y0 / model.noise_level_at(mz0);
        let want_ppm = (jitter.mz_ppm_sigma(sn).powi(2) + jitter.common_mode_mz_ppm.powi(2)).sqrt();
        let (_, sd_ppm) = mean_sd(&ppm);
        assert!((sd_ppm / want_ppm - 1.0).abs() < 0.03, "ppm sd {sd_ppm} vs {want_ppm}");

        let (_, sd_log) = mean_sd(&log_y);
        let want_log = jitter.log_intensity_sigma(sn);
        assert!((sd_log / want_log - 1.0).abs() < 0.03, "log sd {sd_log} vs {want_log}");

        // Mean-preserving, not median-preserving.
        let (mean_y, _) = mean_sd(&ys);
        assert!((mean_y - 1.0).abs() < 0.01, "mean intensity factor {mean_y}");
    }

    #[test]
    fn collapse_merges_only_groups_with_signal() {
        let width = PeakWidthModel::Constant { sigma: 0.01 };
        let peaks = [
            SimulatedPeak { mz: 700.000, intensity: 1.0, origin: PeakOrigin::FtNoise },
            SimulatedPeak { mz: 700.005, intensity: 1.0, origin: PeakOrigin::FtNoise },
            SimulatedPeak { mz: 800.000, intensity: 3.0, origin: PeakOrigin::Signal },
            SimulatedPeak { mz: 800.008, intensity: 1.0, origin: PeakOrigin::FtNoise },
        ];
        let out = collapse(&peaks, &width);
        assert_eq!(out.mz.len(), 3);
        assert_eq!(out.intensity[2], 4.0);
        assert!((out.mz[2] - 800.002).abs() < 1e-12);
    }
}
