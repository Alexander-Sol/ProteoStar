// MS1 peak annotation: pick the most prominent peaks in view and infer each one's charge from the
// isotope spacing of its neighbors, so the spectrum can label them "m/z · z".
//
// Charge inference is a local heuristic — it looks for a run of isotope peaks spaced by
// ~1.00235/z Da around a target peak. It intentionally only assigns a charge when it finds a
// contiguous isotope run (>= 2 matched siblings), and returns null otherwise rather than guessing.

import type { PeakAnnotation, SpectrumPlotPeak } from "@msbrowser/plot-adapter";

/** Δm/z between adjacent isotope peaks of a +1 ion (≈ one neutron); the spacing for charge z is /z. */
export const ISOTOPE_SPACING = 1.00235;

export interface PeakLabel extends PeakAnnotation {
  /** Inferred charge, or null when no confident isotope spacing was found. */
  charge: number | null;
}

/** True if any m/z in the ascending-sorted `sortedMz` lies within `ppm` of `mz` (binary search). */
export function hasPeakNear(sortedMz: readonly number[], mz: number, ppm: number): boolean {
  const n = sortedMz.length;
  if (n === 0) return false;
  const tol = (mz * ppm) / 1e6;
  const loTarget = mz - tol;
  let lo = 0;
  let hi = n;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (sortedMz[mid] < loTarget) lo = mid + 1;
    else hi = mid;
  }
  return lo < n && sortedMz[lo] <= mz + tol;
}

/** Index of the entry in ascending `sortedMz` closest to `mz`. */
export function nearestIndex(sortedMz: readonly number[], mz: number): number {
  const n = sortedMz.length;
  if (n === 0) return -1;
  let lo = 0;
  let hi = n - 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (sortedMz[mid] < mz) lo = mid + 1;
    else hi = mid;
  }
  // `lo` is the first entry >= mz; compare it with its predecessor.
  if (lo > 0 && Math.abs(sortedMz[lo - 1] - mz) <= Math.abs(sortedMz[lo] - mz)) return lo - 1;
  return lo;
}

/** Count contiguous isotope siblings around `targetMz` at spacing `s`, stopping at the first gap in
 *  each direction (isotope envelopes are contiguous). */
function contiguousMatches(
  sortedMz: readonly number[],
  targetMz: number,
  s: number,
  ppm: number,
  maxSteps: number
): number {
  let count = 0;
  for (const dir of [1, -1] as const) {
    for (let k = 1; k <= maxSteps; k++) {
      if (hasPeakNear(sortedMz, targetMz + dir * k * s, ppm)) count++;
      else break;
    }
  }
  return count;
}

/**
 * Infer the charge of the peak at `sortedMz[targetIndex]` from its isotope spacing. Tries each
 * candidate charge, scoring by how many contiguous isotope siblings sit at ~1.00235/z spacing.
 * Returns the best-supported charge (>= 2 siblings), preferring the lower charge on ties (more
 * conservative), or null when no candidate is well supported.
 */
export function estimateCharge(
  sortedMz: readonly number[],
  targetIndex: number,
  opts: { maxCharge?: number; ppm?: number; maxSteps?: number } = {}
): number | null {
  const { maxCharge = 8, ppm = 15, maxSteps = 4 } = opts;
  if (targetIndex < 0 || targetIndex >= sortedMz.length) return null;
  const targetMz = sortedMz[targetIndex];

  let bestZ: number | null = null;
  let bestScore = 0;
  for (let z = 1; z <= maxCharge; z++) {
    const score = contiguousMatches(sortedMz, targetMz, ISOTOPE_SPACING / z, ppm, maxSteps);
    // Strictly-greater keeps the lowest charge on a tie (loop ascends z).
    if (score > bestScore) {
      bestScore = score;
      bestZ = z;
    }
  }
  return bestScore >= 2 ? bestZ : null;
}

/** Label text for a peak: `647.65 · z2`, or just the m/z when charge is unknown. */
export function peakLabelText(mz: number, charge: number | null): string {
  return charge ? `${mz.toFixed(2)} · z${charge}` : mz.toFixed(2);
}

/**
 * Choose the most prominent peaks within the current x-window and build their m/z + charge labels.
 * Peaks are taken tallest-first, de-duplicated so two labels don't stack on the same peak, capped at
 * `maxLabels`.
 */
export function computePeakLabels(
  peaks: readonly SpectrumPlotPeak[],
  opts: {
    xMin?: number | null;
    xMax?: number | null;
    maxLabels?: number;
    maxCharge?: number;
    ppm?: number;
    minSeparation?: number;
  } = {}
): PeakLabel[] {
  const {
    xMin = null,
    xMax = null,
    maxLabels = 6,
    maxCharge = 8,
    ppm = 15,
    minSeparation = 0.5
  } = opts;

  const inView =
    xMin !== null && xMax !== null ? peaks.filter((p) => p.mz >= xMin && p.mz <= xMax) : peaks;
  if (inView.length === 0) return [];

  const sortedMz = [...inView].map((p) => p.mz).sort((a, b) => a - b);
  const byIntensity = [...inView].sort((a, b) => b.intensity - a.intensity);

  const labels: PeakLabel[] = [];
  for (const p of byIntensity) {
    if (labels.length >= maxLabels) break;
    if (labels.some((l) => Math.abs(l.mz - p.mz) < minSeparation)) continue; // don't stack labels
    const charge = estimateCharge(sortedMz, nearestIndex(sortedMz, p.mz), { maxCharge, ppm });
    labels.push({ mz: p.mz, intensity: p.intensity, charge, text: peakLabelText(p.mz, charge) });
  }
  return labels;
}
