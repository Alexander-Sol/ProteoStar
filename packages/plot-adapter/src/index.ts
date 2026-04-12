export const PLOT_ADAPTER_NAME = "plot-adapter";

export interface PlotViewport {
  xMin: number | null;
  xMax: number | null;
}

export function createDefaultViewport(): PlotViewport {
  return {
    xMin: null,
    xMax: null
  };
}
