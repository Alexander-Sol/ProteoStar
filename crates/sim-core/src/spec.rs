//! SimSpec v1 — the simulator input contract (`agent_info/TopDown-SimSpec-Contract-v1.md`).
//!
//! Serde mirror of the C# `SimSpec` records (`mzLib/TopDownSimulator/Simulation/SimulationSpec.cs`).
//! camelCase field names; optional blocks are `Option`s. The same type is the in-app save format,
//! so it serializes back out in the same shape it was read.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::SimError;

/// The schema version this build reads. A spec with any other version is refused.
pub const CURRENT_SCHEMA_VERSION: u32 = 1;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimSpec {
    pub schema_version: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub generator: Option<SimGenerator>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_run: Option<SimSourceRun>,
    pub instrument: SimInstrument,
    pub charge_range: SimChargeRange,
    pub peak_width_model: SimPeakWidthModel,
    /// MS1 retention-time axis, minutes, one entry per scan.
    pub scan_times: Vec<f64>,
    /// `None` means the default reduction (relative 1e-4, absolute 0), as C# `PrepareMs1` applies.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reduction: Option<SimReduction>,
    /// `None` is a clean analytical simulation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub noise: Option<SimNoise>,
    pub proteoforms: Vec<SimProteoform>,
    /// Rendered into the spectra but excluded from feature truth.
    #[serde(default)]
    pub background: Vec<SimProteoform>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimGenerator {
    pub tool: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mz_lib_commit: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub written_utc: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimSourceRun {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimInstrument {
    pub polarity: String,
    pub analyzer: String,
    pub ms_level: i32,
}

impl Default for SimInstrument {
    fn default() -> Self {
        Self { polarity: "Positive".into(), analyzer: "Orbitrap".into(), ms_level: 1 }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SimChargeRange {
    pub min: i32,
    pub max: i32,
}

/// `σ_m(mz)`: `Orbitrap` is `k·(m/z)^1.5`, `Constant` is a fixed σ.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum SimPeakWidthModel {
    Orbitrap { k: f64 },
    Constant { sigma: f64 },
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimReduction {
    pub relative_intensity_threshold: f64,
    pub absolute_intensity_floor: f64,
}

impl Default for SimReduction {
    fn default() -> Self {
        Self { relative_intensity_threshold: 1e-4, absolute_intensity_floor: 0.0 }
    }
}

/// The live-editable noise configuration (contract §2a).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimNoise {
    /// Truncated to its low 32 bits, as the C# `int noiseSeed` the simulator takes.
    pub seed: i64,
    pub floor: SimNoiseFloor,
    pub jitter: SimNoiseJitter,
}

impl SimNoise {
    /// The measured Orbitrap floor + jitter defaults, enabled.
    pub fn orbitrap_defaults(seed: i64) -> Self {
        Self { seed, floor: SimNoiseFloor::default(), jitter: SimNoiseJitter::default() }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimNoiseFloor {
    pub noise_level_at_reference_mz: f64,
    pub density_scale: f64,
    pub level_dispersion: f64,
    pub signal_to_noise: SimSignalToNoise,
}

impl Default for SimNoiseFloor {
    fn default() -> Self {
        Self {
            noise_level_at_reference_mz: 1120.0,
            density_scale: 1.0,
            level_dispersion: 0.25,
            signal_to_noise: SimSignalToNoise::default(),
        }
    }
}

/// S/N = threshold + Exp(meanExcess).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimSignalToNoise {
    pub threshold: f64,
    pub mean_excess: f64,
}

impl Default for SimSignalToNoise {
    /// C# `SignalToNoiseDistribution.OrbitrapCentroids`.
    fn default() -> Self {
        Self { threshold: 1.94, mean_excess: 1.165 }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimNoiseJitter {
    pub enabled: bool,
    pub mz_ppm_at_unit_signal_to_noise: f64,
    pub mz_ppm_floor: f64,
    pub common_mode_mz_ppm: f64,
    pub log_intensity_sigma_at_unit_signal_to_noise: f64,
    pub log_intensity_sigma_floor: f64,
    pub common_mode_log_intensity_sigma: f64,
    pub drop_peaks_falling_below_floor: bool,
}

impl Default for SimNoiseJitter {
    /// C# `PeakJitterModel.Orbitrap`, enabled.
    fn default() -> Self {
        Self {
            enabled: true,
            mz_ppm_at_unit_signal_to_noise: 9.06,
            mz_ppm_floor: 0.81,
            common_mode_mz_ppm: 1.887,
            log_intensity_sigma_at_unit_signal_to_noise: 1.026,
            log_intensity_sigma_floor: 0.165,
            common_mode_log_intensity_sigma: 0.0,
            drop_peaks_falling_below_floor: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimProteoform {
    /// Neutral monoisotopic mass, Da.
    pub monoisotopic_mass: f64,
    /// `A_p`.
    pub abundance: f64,
    pub rt: SimRtProfile,
    pub charge: SimChargeDistribution,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub identifier: Option<String>,
}

/// Unit-area EMG elution profile; minutes.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SimRtProfile {
    pub mu: f64,
    pub sigma: f64,
    pub tau: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum SimChargeDistribution {
    #[serde(rename_all = "camelCase")]
    Gaussian { mu_z: f64, sigma_z: f64 },
}

impl SimSpec {
    pub fn from_json(json: &str) -> Result<Self, SimError> {
        // Peek at the version first so a future-schema file reports as such rather than as
        // whatever field happened to fail to deserialize.
        #[derive(Deserialize)]
        #[serde(rename_all = "camelCase")]
        struct Version {
            schema_version: u32,
        }
        let version: Version = serde_json::from_str(json).map_err(SimError::Json)?;
        if version.schema_version != CURRENT_SCHEMA_VERSION {
            return Err(SimError::UnsupportedSchema(version.schema_version));
        }

        let spec: SimSpec = serde_json::from_str(json).map_err(SimError::Json)?;
        spec.validate()?;
        Ok(spec)
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self, SimError> {
        let text = std::fs::read_to_string(path).map_err(SimError::Io)?;
        Self::from_json(&text)
    }

    pub fn to_json(&self) -> Result<String, SimError> {
        serde_json::to_string_pretty(self).map_err(SimError::Json)
    }

    pub fn write(&self, path: impl AsRef<Path>) -> Result<(), SimError> {
        std::fs::write(path, self.to_json()?).map_err(SimError::Io)
    }

    /// The invariants the C# model constructors enforce, checked up front so a spec edited in
    /// the UI fails with a message rather than a NaN spectrum.
    pub fn validate(&self) -> Result<(), SimError> {
        let invalid = |msg: String| Err(SimError::InvalidSpec(msg));

        let SimChargeRange { min, max } = self.charge_range;
        if min < 1 || max < min {
            return invalid(format!("charge range must satisfy 1 <= min <= max; got [{min}, {max}]"));
        }

        match self.peak_width_model {
            SimPeakWidthModel::Orbitrap { k } if !(k > 0.0 && k.is_finite()) => {
                return invalid(format!("Orbitrap k must be finite and positive; got {k}"));
            }
            SimPeakWidthModel::Constant { sigma } if !(sigma > 0.0 && sigma.is_finite()) => {
                return invalid(format!("constant sigma must be finite and positive; got {sigma}"));
            }
            _ => {}
        }

        if self.scan_times.iter().any(|t| !t.is_finite()) {
            return invalid("scan times must be finite".into());
        }
        if self.scan_times.windows(2).any(|w| !(w[1] > w[0])) {
            return invalid("scan times must be strictly ascending".into());
        }

        if let Some(r) = self.reduction {
            if !(r.relative_intensity_threshold >= 0.0) {
                return invalid("relative intensity threshold must be non-negative".into());
            }
            if !(r.absolute_intensity_floor >= 0.0) {
                return invalid("absolute intensity floor must be non-negative".into());
            }
        }

        if let Some(n) = self.noise {
            let f = n.floor;
            if !(f.noise_level_at_reference_mz > 0.0 && f.noise_level_at_reference_mz.is_finite()) {
                return invalid("noise level must be finite and positive".into());
            }
            if !(f.density_scale >= 0.0 && f.density_scale.is_finite()) {
                return invalid("density scale must be finite and non-negative".into());
            }
            if !(f.level_dispersion >= 0.0 && f.level_dispersion.is_finite()) {
                return invalid("level dispersion must be finite and non-negative".into());
            }
            let sn = f.signal_to_noise;
            if !(sn.threshold > 0.0 && sn.threshold.is_finite()) {
                return invalid("S/N threshold must be finite and positive".into());
            }
            if !(sn.mean_excess > 0.0 && sn.mean_excess.is_finite()) {
                return invalid("S/N mean excess must be finite and positive".into());
            }
        }

        for p in self.proteoforms.iter().chain(&self.background) {
            let name = p.identifier.as_deref().unwrap_or("<unnamed>");
            if !(p.monoisotopic_mass > 0.0 && p.monoisotopic_mass.is_finite()) {
                return invalid(format!("{name}: monoisotopic mass must be finite and positive"));
            }
            if !p.abundance.is_finite() {
                return invalid(format!("{name}: abundance must be finite"));
            }
            if !(p.rt.sigma > 0.0) {
                return invalid(format!("{name}: RT sigma must be positive"));
            }
            let SimChargeDistribution::Gaussian { sigma_z, .. } = p.charge;
            if !(sigma_z > 0.0) {
                return invalid(format!("{name}: charge sigma must be positive"));
            }
        }

        Ok(())
    }
}
