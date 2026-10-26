import { describe, expect, it } from "vitest";

import {
  createDefaultViewport,
  fitSpectrumTracesYRange,
  fitSpectrumYRange,
  hasPersistedY,
  resolveSpectrumTracesYRange,
  resolveSpectrumYRange
} from "../src/viewport";
import type { SpectrumPlotPeak } from "../src/types";

const peaks: SpectrumPlotPeak[] = [
  { mz: 400, intensity: 100 },
  { mz: 500, intensity: 1000 },
  { mz: 600, intensity: 300 }
];

describe("createDefaultViewport", () => {
  it("is fully unbounded (auto-fit both axes)", () => {
    expect(createDefaultViewport()).toEqual({ xMin: null, xMax: null, yMin: null, yMax: null });
  });
});

describe("fitSpectrumYRange", () => {
  it("fits to the tallest peak with 5% headroom", () => {
    expect(fitSpectrumYRange(peaks, null, null)).toEqual([0, 1050]);
  });

  it("restricts the fit to peaks inside the x-window", () => {
    // Only the 600 peak (intensity 300) is in [560, 620]; the tall 500 peak is excluded.
    expect(fitSpectrumYRange(peaks, 560, 620)?.[1]).toBeCloseTo(315);
  });

  it("returns undefined when no peaks are in view", () => {
    expect(fitSpectrumYRange(peaks, 700, 800)).toBeUndefined();
  });
});

describe("hasPersistedY", () => {
  it("is true only when both y bounds are set", () => {
    expect(hasPersistedY({ xMin: null, xMax: null, yMin: 0, yMax: 500 })).toBe(true);
    expect(hasPersistedY({ xMin: null, xMax: null, yMin: null, yMax: 500 })).toBe(false);
    expect(hasPersistedY(createDefaultViewport())).toBe(false);
  });
});

describe("resolveSpectrumYRange", () => {
  it("uses the persisted y-range when present (stable envelope height)", () => {
    const vp = { xMin: 490, xMax: 510, yMin: 0, yMax: 1234 };
    expect(resolveSpectrumYRange(vp, peaks)).toEqual([0, 1234]);
  });

  it("auto-fits when no persisted y-range is set", () => {
    const vp = createDefaultViewport();
    expect(resolveSpectrumYRange(vp, peaks)).toEqual([0, 1050]);
  });
});

describe("fitSpectrumTracesYRange", () => {
  const sim: SpectrumPlotPeak[] = [
    { mz: 450, intensity: 2000 },
    { mz: 600, intensity: 50 }
  ];
  const real = { slotIndex: 0 as const, peaks, color: "#000" };
  const mirrored = { slotIndex: 1 as const, peaks: sim, color: "#f00", mirrored: true };

  it("matches fitSpectrumYRange without a mirrored trace", () => {
    expect(fitSpectrumTracesYRange([real], null, null)).toEqual([0, 1050]);
  });

  it("extends below zero by the tallest mirrored peak", () => {
    expect(fitSpectrumTracesYRange([real, mirrored], null, null)).toEqual([-2100, 1050]);
  });

  it("windows each direction independently", () => {
    const [lo, hi] = fitSpectrumTracesYRange([real, mirrored], 560, 620)!;
    expect(lo).toBeCloseTo(-52.5);
    expect(hi).toBeCloseTo(315);
  });

  it("fits a mirrored trace alone", () => {
    expect(fitSpectrumTracesYRange([mirrored], null, null)).toEqual([-2100, 0]);
  });

  it("resolves to a persisted y-range when present", () => {
    const vp = { xMin: null, xMax: null, yMin: -10, yMax: 10 };
    expect(resolveSpectrumTracesYRange(vp, [real, mirrored])).toEqual([-10, 10]);
  });
});
