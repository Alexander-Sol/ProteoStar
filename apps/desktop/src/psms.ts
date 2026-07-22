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
  /** How many (most abundant) isotope peaks to trace. Default 3. */
  numIsotopes?: number;
  /** Half-width of each isotope's m/z extraction window, in ppm. Default 15. */
  halfWindowPpm?: number;
  rtRange?: { min: number; max: number };
  maxPoints?: number;
}

/** One species' isotope XICs plus the isotope indices they correspond to (for labelling). */
export interface IsotopeXics {
  /** Isotope indices traced, ascending (0 = monoisotopic). e.g. `[0,1,2]` for a small peptide,
   *  `[6,7,8]` for a large proteoform whose envelope apex sits well above the monoisotope. */
  indices: number[];
  /** One XIC per entry in `indices`, aligned by position (all share the same scan set). */
  traces: XicPoint[][];
}

/**
 * Approximate relative isotopologue abundances for a neutral `mass`, as a Poisson distribution in
 * the number of heavy isotopes. ¹³C dominates a peptide/protein envelope: averagine has ≈4.9384 C
 * per ~111.05 Da residue and ¹³C occurs at ≈1.07%, giving a Poisson mean λ ≈ mass · 4.76e-4. Good
 * enough to pick *which* isotopologues are the tallest (the only use here); it ignores N/O/S/H
 * heavy isotopes, which only slightly broaden the real envelope. Values are unnormalised.
 */
export function isotopeAbundances(mass: number, count: number): number[] {
  const lambda = Math.max(0, mass * 4.757e-4);
  const out = new Array<number>(count);
  let term = Math.exp(-lambda); // Poisson pmf at k = 0
  for (let k = 0; k < count; k++) {
    out[k] = term;
    term = (term * lambda) / (k + 1); // pmf(k+1) = pmf(k) · λ/(k+1)
  }
  return out;
}

/** Indices of the `n` most abundant isotopologues for `mass` (see [`isotopeAbundances`]), returned
 *  in ascending index order. For a small peptide this is `[0,…,n-1]`; for a large proteoform it is
 *  the `n` peaks straddling the envelope apex (which sits above the monoisotope). */
export function topIsotopeIndices(mass: number, n: number): number[] {
  const lambda = Math.max(0, mass * 4.757e-4);
  const count = Math.max(n, Math.ceil(lambda) + 4); // search a window that covers the apex
  return isotopeAbundances(mass, count)
    .map((a, i) => [a, i] as const)
    .sort((x, y) => y[0] - x[0])
    .slice(0, n)
    .map(([, i]) => i)
    .sort((a, b) => a - b);
}

/**
 * Extract one XIC for each of the `numIsotopes` most abundant isotopologues of a species at
 * (`mass`, `charge`). Each isotope's m/z comes from the averagine ladder ([`isotopeGrid`]) and is
 * summed over a ±`halfWindowPpm` window per scan via the existing `getRangeXic` path. Returns the
 * traced isotope indices and one `XicPoint[]` per index (aligned by position; all share the same
 * scan set, so the summed-envelope view just adds them position-wise).
 */
export async function fetchIsotopeXics(
  provider: DatasetProvider,
  mass: number,
  charge: number,
  opts: IsotopeXicOptions = {}
): Promise<IsotopeXics> {
  const num = opts.numIsotopes ?? 3;
  const halfPpm = opts.halfWindowPpm ?? 15;
  const indices = topIsotopeIndices(mass, num);
  const maxIdx = indices.length > 0 ? indices[indices.length - 1] : 0;
  const grid = isotopeGrid(mass, charge, maxIdx + 1);
  const traces = await Promise.all(
    indices.map((k) => {
      const mz = grid[k];
      const half = (mz * halfPpm) / 1e6;
      return provider
        .getRangeXic(mz - half, mz + half, { rtRange: opts.rtRange, maxPoints: opts.maxPoints })
        .then((pts) => [...pts]);
    })
  );
  return { indices, traces };
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
