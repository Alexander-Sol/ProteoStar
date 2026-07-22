import { describe, expect, it } from "vitest";

import {
  ISOTOPE_SPACING,
  computePeakLabels,
  estimateCharge,
  nearestIndex,
  peakLabelText
} from "./annotate";

// Build an isotope envelope of `n` peaks at charge `z` starting at `baseMz`, spaced 1.00235/z.
function envelope(baseMz: number, z: number, n: number, intensity = 1000) {
  return Array.from({ length: n }, (_, k) => ({
    mz: baseMz + (k * ISOTOPE_SPACING) / z,
    intensity: intensity * Math.pow(0.7, k)
  }));
}

describe("estimateCharge", () => {
  it("infers z=2 from a half-Da isotope run", () => {
    const peaks = envelope(500, 2, 4);
    const mz = peaks.map((p) => p.mz).sort((a, b) => a - b);
    expect(estimateCharge(mz, nearestIndex(mz, 500))).toBe(2);
  });

  it("infers z=1 from a one-Da isotope run", () => {
    const peaks = envelope(800, 1, 3);
    const mz = peaks.map((p) => p.mz).sort((a, b) => a - b);
    expect(estimateCharge(mz, nearestIndex(mz, 800))).toBe(1);
  });

  it("infers z=3 from a third-Da isotope run", () => {
    const peaks = envelope(650, 3, 4);
    const mz = peaks.map((p) => p.mz).sort((a, b) => a - b);
    expect(estimateCharge(mz, nearestIndex(mz, 650))).toBe(3);
  });

  it("returns null for an isolated peak with no siblings", () => {
    const mz = [400, 1000, 1600];
    expect(estimateCharge(mz, nearestIndex(mz, 1000))).toBeNull();
  });
});

describe("computePeakLabels", () => {
  it("labels the prominent peaks with m/z and inferred charge", () => {
    const peaks = [...envelope(500, 2, 4, 5000), ...envelope(800, 1, 3, 3000)];
    const labels = computePeakLabels(peaks, { maxLabels: 6 });

    const mono2 = labels.find((l) => Math.abs(l.mz - 500) < 1e-6);
    const mono1 = labels.find((l) => Math.abs(l.mz - 800) < 1e-6);
    expect(mono2?.charge).toBe(2);
    expect(mono2?.text).toBe("500.00<br>z2");
    expect(mono1?.charge).toBe(1);
  });

  it("respects the x-window and the label cap", () => {
    const peaks = [...envelope(500, 2, 4, 5000), ...envelope(800, 1, 3, 3000)];
    const labels = computePeakLabels(peaks, { xMin: 490, xMax: 520, maxLabels: 2 });
    expect(labels.length).toBe(2);
    expect(labels.every((l) => l.mz >= 490 && l.mz <= 520)).toBe(true);
  });

  it("does not stack two labels on the same peak", () => {
    const peaks = envelope(500, 2, 4, 5000);
    const labels = computePeakLabels(peaks, { minSeparation: 0.5 });
    const mzs = labels.map((l) => l.mz);
    for (let i = 0; i < mzs.length; i++)
      for (let j = i + 1; j < mzs.length; j++) expect(Math.abs(mzs[i] - mzs[j])).toBeGreaterThan(0.5);
  });
});

describe("peakLabelText", () => {
  it("includes charge when known and omits it otherwise", () => {
    expect(peakLabelText(647.653, 2)).toBe("647.65<br>z2");
    expect(peakLabelText(647.653, null)).toBe("647.65");
  });
});
