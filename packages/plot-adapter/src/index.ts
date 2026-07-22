export const PLOT_ADAPTER_NAME = "plot-adapter";

export { TicPlot, SpectrumPlot } from "./plots";
export {
  createDefaultViewport,
  fitSpectrumYRange,
  hasPersistedY,
  resolveSpectrumYRange
} from "./viewport";
export type {
  EnvelopeLine,
  FeatureMarker,
  NumericRange,
  PeakAnnotation,
  PlotViewport,
  RtRegion,
  SlotIndex,
  SpectrumPlotEvent,
  SpectrumPlotPeak,
  SpectrumPlotProps,
  SpectrumPlotTrace,
  TicPlotEvent,
  TicPlotPoint,
  TicPlotProps,
  TicPlotTrace
} from "./types";
