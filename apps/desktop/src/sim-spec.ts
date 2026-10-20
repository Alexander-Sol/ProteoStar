// SimSpec v1 — the simulator input contract (agent_info/TopDown-SimSpec-Contract-v1.md).
//
// Mirror of the Rust `sim_core::spec::SimSpec` and the C# `SimSpec`
// (mzLib/TopDownSimulator/Simulation/SimulationSpec.cs). camelCase field names cross the Tauri
// boundary to `open_simulated_dataset` / `simulate_scan` verbatim. This is also the in-app save
// format — the workbench holds a live `SimSpec` the parameter panel edits.

export const SIM_SCHEMA_VERSION = 1;

export interface SimSpec {
  schemaVersion: number;
  generator?: SimGenerator;
  sourceRun?: SimSourceRun;
  instrument: SimInstrument;
  chargeRange: SimChargeRange;
  peakWidthModel: SimPeakWidthModel;
  /** MS1 retention-time axis, minutes, strictly ascending, one entry per scan. */
  scanTimes: number[];
  /** Omit for the default reduction (relative 1e-4, absolute 0). */
  reduction?: SimReduction;
  /** Omit for a clean analytical simulation (no floor, no jitter). */
  noise?: SimNoise;
  proteoforms: SimProteoform[];
  /** Rendered into the spectra but excluded from feature truth. */
  background?: SimProteoform[];
}

export interface SimGenerator {
  tool: string;
  mzLibCommit?: string;
  writtenUtc?: string;
}

export interface SimSourceRun {
  fileName?: string;
  description?: string;
}

export interface SimInstrument {
  polarity: string;
  analyzer: string;
  msLevel: number;
}

export interface SimChargeRange {
  min: number;
  max: number;
}

/** `σ_m(mz)`: `Orbitrap` is `k·(m/z)^1.5`; `Constant` is a fixed σ. */
export type SimPeakWidthModel =
  | { kind: "Orbitrap"; k: number }
  | { kind: "Constant"; sigma: number };

export interface SimReduction {
  relativeIntensityThreshold: number;
  absoluteIntensityFloor: number;
}

/** The live-editable noise configuration (contract §2a). */
export interface SimNoise {
  seed: number;
  floor: SimNoiseFloor;
  jitter: SimNoiseJitter;
}

export interface SimNoiseFloor {
  noiseLevelAtReferenceMz: number;
  densityScale: number;
  levelDispersion: number;
  signalToNoise: SimSignalToNoise;
}

/** S/N = threshold + Exp(meanExcess). */
export interface SimSignalToNoise {
  threshold: number;
  meanExcess: number;
}

export interface SimNoiseJitter {
  enabled: boolean;
  mzPpmAtUnitSignalToNoise: number;
  mzPpmFloor: number;
  commonModeMzPpm: number;
  logIntensitySigmaAtUnitSignalToNoise: number;
  logIntensitySigmaFloor: number;
  commonModeLogIntensitySigma: number;
  dropPeaksFallingBelowFloor: boolean;
}

export interface SimProteoform {
  /** Neutral monoisotopic mass, Da. */
  monoisotopicMass: number;
  abundance: number;
  rt: SimRtProfile;
  charge: SimChargeDistribution;
  identifier?: string;
}

/** Unit-area EMG elution profile; minutes. */
export interface SimRtProfile {
  mu: number;
  sigma: number;
  tau: number;
}

export type SimChargeDistribution = { kind: "Gaussian"; muZ: number; sigmaZ: number };

/** The measured Orbitrap floor + jitter defaults, enabled — mirrors `SimNoise::orbitrap_defaults`. */
export function orbitrapNoiseDefaults(seed = 0): SimNoise {
  return {
    seed,
    floor: {
      noiseLevelAtReferenceMz: 1120.0,
      densityScale: 1.0,
      levelDispersion: 0.25,
      signalToNoise: { threshold: 1.94, meanExcess: 1.165 }
    },
    jitter: {
      enabled: true,
      mzPpmAtUnitSignalToNoise: 9.06,
      mzPpmFloor: 0.81,
      commonModeMzPpm: 1.887,
      logIntensitySigmaAtUnitSignalToNoise: 1.026,
      logIntensitySigmaFloor: 0.165,
      commonModeLogIntensitySigma: 0.0,
      dropPeaksFallingBelowFloor: true
    }
  };
}
