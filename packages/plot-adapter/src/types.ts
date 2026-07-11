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

// A feature-finding result rendered on the TIC as a marker at its apex RT (the
// "feature rug"). Clicking it selects the feature.
export interface FeatureMarker {
  featureIndex: number;
  retentionTime: number;
  color: string;
  label: string;
}

// A shaded RT band (e.g. the selected feature's traced elution extent).
export interface RtRegion {
  min: number;
  max: number;
  color: string;
}

// A predicted isotope-peak m/z line drawn over a spectrum.
export interface EnvelopeLine {
  mz: number;
  color: string;
  label?: string;
}

export type TicPlotEvent =
  | { type: "area-click"; retentionTime: number }
  | { type: "feature-click"; featureIndex: number }
  | { type: "point-hover"; point: TicPlotPoint | null }
  | { type: "range-select"; range: NumericRange };

export type SpectrumPlotEvent =
  | { type: "point-hover"; peak: SpectrumPlotPeak | null }
  | { type: "range-select"; range: NumericRange };

export interface TicPlotProps {
  traces: readonly TicPlotTrace[];
  viewport: PlotViewport;
  rangeSelectionEnabled: boolean;
  /** Feature apex markers overlaid along the baseline (optional). */
  featureRug?: readonly FeatureMarker[];
  /** Shaded RT band(s) — typically the selected feature's elution extent (optional). */
  regions?: readonly RtRegion[];
  onEvent(event: TicPlotEvent): void;
}

export interface SpectrumPlotProps {
  traces: readonly SpectrumPlotTrace[];
  viewport: PlotViewport;
  rangeSelectionEnabled: boolean;
  /** Predicted isotope-peak m/z lines overlaid on the spectrum (optional). */
  envelope?: readonly EnvelopeLine[];
  onEvent(event: SpectrumPlotEvent): void;
}
