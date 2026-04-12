export type SlotIndex = 0 | 1;

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

export interface TicPlotTrace {
  slotIndex: SlotIndex;
  points: readonly TicPlotPoint[];
  selectedScanIndex: number | null;
  color: string;
}

export interface SpectrumPlotTrace {
  slotIndex: SlotIndex;
  peaks: readonly SpectrumPlotPeak[];
  color: string;
}

export type TicPlotEvent =
  | { type: "point-click"; slotIndex: SlotIndex; point: TicPlotPoint }
  | { type: "point-hover"; point: TicPlotPoint | null }
  | { type: "range-select"; range: NumericRange };

export type SpectrumPlotEvent =
  | { type: "point-hover"; peak: SpectrumPlotPeak | null }
  | { type: "range-select"; range: NumericRange };

export interface TicPlotProps {
  traces: readonly TicPlotTrace[];
  viewport: PlotViewport;
  rangeSelectionEnabled: boolean;
  onEvent(event: TicPlotEvent): void;
}

export interface SpectrumPlotProps {
  traces: readonly SpectrumPlotTrace[];
  viewport: PlotViewport;
  rangeSelectionEnabled: boolean;
  onEvent(event: SpectrumPlotEvent): void;
}
