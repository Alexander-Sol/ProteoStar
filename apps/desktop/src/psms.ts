// PSM overlay: load a MetaMorpheus `.psmtsv`, link each PSM to a detected feature
// (mass + RT + charge join, done here where both lists live), and extract the
// per-isotope XICs the selector displays for a chosen PSM.

import { invoke } from "@tauri-apps/api/core";

import type { DatasetProvider, Feature, Psm, PsmLink, XicPoint } from "./contract";
import { isotopeGrid } from "./features";

/** Load PSMs from a MetaMorpheus `.psmtsv` file. */
export function loadPsms(path: string): Promise<Psm[]> {
  return invoke<Psm[]>("load_psms", { path });
}

export interface LinkOptions {
  /** Mass match tolerance in ppm (default 10). */
  massPpm?: number;
  /** RT padding (min) added either side of the feature's elution window (default 0.2). */
  rtPad?: number;
  /** Require the PSM's charge to be one of the feature's resolved charge states (default true). */
  requireCharge?: boolean;
}

/**
 * Link each PSM to at most one detected feature. A match requires the PSM's monoisotopic
 * mass within `massPpm` of the feature's, its MS2 RT inside the feature's padded elution
 * window (skipped when the psmtsv gives no RT), and — by default — its charge among the
 * feature's charge states. On multiple candidates the smallest mass error wins. PSMs that
 * match nothing get a null link (they stay in the list, using their theoretical m/z).
 *
 * Features are searched via a mass-sorted index + binary search, so this stays fast even on
 * the 100k+ feature lists a top-down run produces.
 */
export function linkPsms(
  psms: readonly Psm[],
  features: readonly Feature[],
  opts: LinkOptions = {}
): PsmLink[] {
  const massPpm = opts.massPpm ?? 10;
  const rtPad = opts.rtPad ?? 0.2;
  const requireCharge = opts.requireCharge ?? true;

  // Feature indices ordered by ascending monoisotopic mass, with a parallel mass array
  // for binary search. Built once and reused across every PSM.
  const order = features.map((_, i) => i).sort(
    (a, b) => features[a].monoisotopicMass - features[b].monoisotopicMass
  );
  const masses = order.map((i) => features[i].monoisotopicMass);

  // First position in `masses` whose value is ≥ target.
  const lowerBound = (target: number): number => {
    let lo = 0;
    let hi = masses.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (masses[mid] < target) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  };

  const unlinked: PsmLink = { featureIndex: null, massPpmError: null, rtDelta: null };

  return psms.map((psm) => {
    const mass = psm.monoisotopicMass;
    if (!(mass > 0)) return unlinked;
    const tol = (mass * massPpm) / 1e6;

    let best: PsmLink | null = null;
    for (let k = lowerBound(mass - tol); k < masses.length && masses[k] <= mass + tol; k++) {
      const fi = order[k];
      const f = features[fi];
      if (requireCharge && !f.chargeStates.includes(psm.precursorCharge)) continue;

      let rtDelta: number | null = null;
      if (psm.ms2RetentionTime >= 0) {
        if (
          psm.ms2RetentionTime < f.rtStart - rtPad ||
          psm.ms2RetentionTime > f.rtEnd + rtPad
        ) {
          continue;
        }
        rtDelta = psm.ms2RetentionTime - f.rtApex;
      }

      const ppm = (Math.abs(mass - f.monoisotopicMass) / f.monoisotopicMass) * 1e6;
      if (best === null || ppm < best.massPpmError!) {
        best = { featureIndex: fi, massPpmError: ppm, rtDelta };
      }
    }
    return best ?? unlinked;
  });
}

export interface IsotopeXicOptions {
  /** How many isotope peaks to trace (mono, +1, …). Default 3. */
  numIsotopes?: number;
  /** Half-width of each isotope's m/z extraction window, in ppm. Default 15. */
  halfWindowPpm?: number;
  rtRange?: { min: number; max: number };
  maxPoints?: number;
}

/**
 * Extract one XIC per isotope for a species at (`mass`, `charge`). Each isotope's m/z comes
 * from the averagine ladder and is summed over a ±`halfWindowPpm` window per scan via the
 * existing `getRangeXic` path. Returns one `XicPoint[]` per isotope, index 0 = monoisotopic.
 * All traces share the same scan set (same RT window + decimation), so they align by index —
 * the summed-envelope view just adds them position-wise.
 */
export async function fetchIsotopeXics(
  provider: DatasetProvider,
  mass: number,
  charge: number,
  opts: IsotopeXicOptions = {}
): Promise<XicPoint[][]> {
  const num = opts.numIsotopes ?? 3;
  const halfPpm = opts.halfWindowPpm ?? 15;
  const grid = isotopeGrid(mass, charge, num);
  return Promise.all(
    grid.map((mz) => {
      const half = (mz * halfPpm) / 1e6;
      return provider
        .getRangeXic(mz - half, mz + half, { rtRange: opts.rtRange, maxPoints: opts.maxPoints })
        .then((pts) => [...pts]);
    })
  );
}

/** Sum aligned per-isotope XICs into a single envelope trace (position-wise; all share scans). */
export function sumXics(xics: readonly XicPoint[][]): XicPoint[] {
  if (xics.length === 0) return [];
  const base = xics[0];
  return base.map((pt, i) => ({
    scanIndex: pt.scanIndex,
    retentionTime: pt.retentionTime,
    intensity: xics.reduce((s, x) => s + (x[i]?.intensity ?? 0), 0)
  }));
}
