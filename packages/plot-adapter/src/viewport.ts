import type { PlotViewport, SpectrumPlotPeak } from "./types";

export function createDefaultViewport(): PlotViewport {
  return {
    xMin: null,
    xMax: null,
    yMin: null,
    yMax: null
  };
}

/** True when the viewport carries a usable persisted y-range (both bounds set). */
export function hasPersistedY(viewport: PlotViewport): boolean {
  return (
    viewport.yMin !== null &&
    viewport.yMin !== undefined &&
    viewport.yMax !== null &&
    viewport.yMax !== undefined
  );
}

/**
 * Auto-fit y-range for a spectrum: `[0, 1.05 * max intensity]` over the peaks inside the current
 * x-window (or all peaks when x is unbounded). Returns `undefined` when there are no peaks in view,
 * so the caller can let Plotly autorange.
 */
export function fitSpectrumYRange(
  peaks: readonly SpectrumPlotPeak[],
  xMin: number | null,
  xMax: number | null
): [number, number] | undefined {
  const inView =
    xMin !== null && xMax !== null
      ? peaks.filter((p) => p.mz >= xMin && p.mz <= xMax)
      : peaks;
  if (inView.length === 0) return undefined;
  let maxY = 0;
  for (const p of inView) if (p.intensity > maxY) maxY = p.intensity;
  return [0, maxY * 1.05];
}

/**
 * Resolve the y-range a spectrum should render with: the viewport's persisted y-range when present
 * (holds the envelope height stable across scan steps), otherwise the auto-fit to the visible
 * x-window. This is the single source of truth for the "stop the y-axis continually rescaling"
 * behavior — kept pure so it's unit-testable without Plotly.
 */
export function resolveSpectrumYRange(
  viewport: PlotViewport,
  peaks: readonly SpectrumPlotPeak[]
): [number, number] | undefined {
  if (hasPersistedY(viewport)) {
    return [viewport.yMin as number, viewport.yMax as number];
  }
  return fitSpectrumYRange(peaks, viewport.xMin, viewport.xMax);
}
