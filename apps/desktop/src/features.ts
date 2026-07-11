// Feature-finding overlay: load the runner's resolved-feature TSV and derive the
// isotope-envelope m/z grids used to overlay predicted peaks on a spectrum.

import { invoke, Channel } from "@tauri-apps/api/core";

import type { Feature, ProgressEvent } from "./contract";

const PROTON_MASS = 1.0072764668;
// ¹³C − ¹²C mass difference — the isotope spacing in Da (per charge, divide by z).
const C13_C12 = 1.0033548;

/** Load resolved features from a TSV file produced by `detect_features_tsv`. */
export function loadFeatures(path: string): Promise<Feature[]> {
  return invoke<Feature[]>("load_features", { path });
}

export interface DetectOptions {
  maxCharge?: number;
  minSeedIntensity?: number;
}

/**
 * Run feature finding in-process on an open dataset (handle), streaming
 * detect/refine/resolve progress. Uses the top-down pipeline with a charge cap
 * (default 25 — the guardrail that avoids the full-range detect crash).
 */
export function runFeatureDetection(
  handle: number,
  options: DetectOptions = {},
  onProgress?: (p: ProgressEvent) => void
): Promise<Feature[]> {
  const channel = new Channel<ProgressEvent>();
  if (onProgress) channel.onmessage = onProgress;
  return invoke<Feature[]>("run_feature_detection", { handle, options, onProgress: channel });
}

/** Monoisotopic neutral mass → m/z at a charge state. */
export function massToMz(mass: number, charge: number): number {
  return (mass + charge * PROTON_MASS) / Math.max(1, charge);
}

/**
 * Predicted isotope-peak m/z positions for one charge state of a feature:
 * `monoMz + k · (C13_C12 / z)` for k = 0 … count-1. These are the ticks the
 * spectrum overlay draws so you can see whether the claimed envelope lines up
 * with real peaks (and where the monoisotope was placed).
 */
export function isotopeGrid(mass: number, charge: number, count = 8): number[] {
  const monoMz = massToMz(mass, charge);
  const spacing = C13_C12 / Math.max(1, charge);
  return Array.from({ length: count }, (_, k) => monoMz + k * spacing);
}

/** True if `rt` falls within the feature's traced elution extent. */
export function featureElutesAt(feature: Feature, rt: number, pad = 0): boolean {
  return rt >= feature.rtStart - pad && rt <= feature.rtEnd + pad;
}
