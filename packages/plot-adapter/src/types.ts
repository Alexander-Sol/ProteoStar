export interface PlotViewport {
  xMin: number | null;
  xMax: number | null;
}

export interface NumericRange {
  min: number;
  max: number;
}

export interface TicPlotPoint {
  scanIndex: number;
  retentionTime: number;
  intensity: number;
}

export interface SpectrumPlotPeak {
  mz: number;
  intensity: number;
}

export type TicPlotEvent =
  | { type: "point-click"; point: TicPlotPoint }
  | { type: "point-hover"; point: TicPlotPoint | null }
  | { type: "range-select"; range: NumericRange };

export type SpectrumPlotEvent =
  | { type: "point-hover"; peak: SpectrumPlotPeak | null }
  | { type: "range-select"; range: NumericRange };

export interface TicPlotProps {
  points: readonly TicPlotPoint[];
  viewport: PlotViewport;
  selectedScanIndex: number | null;
  rangeSelectionEnabled: boolean;
  onEvent(event: TicPlotEvent): void;
}

export interface SpectrumPlotProps {
  peaks: readonly SpectrumPlotPeak[];
  viewport: PlotViewport;
  rangeSelectionEnabled: boolean;
  onEvent(event: SpectrumPlotEvent): void;
}
