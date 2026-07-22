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
  /** When any trace sets this, the y-axis auto-range is driven by those traces alone (others
   *  keep their true absolute height but are free to clip out of view). A summed-XIC overlay
   *  sets it so the plot zooms to the XIC's abundance while the much taller TIC runs off the top. */
  yScale?: boolean;
  /** Hover-series name for this trace (defaults to the plot's `yTitle`). */
  label?: string;
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
  /** True when an observed peak sits at this comb position in the displayed scan — drawn as a
   *  solid, opaque line versus a faint dotted line for a missing (predicted-only) tooth. */
  detected?: boolean;
}

export type TicPlotEvent =
  | { type: "area-click"; retentionTime: number }
  | { type: "feature-click"; featureIndex: number }
  | { type: "point-hover"; point: TicPlotPoint | null }
  | { type: "range-select"; range: NumericRange };

export type SpectrumPlotEvent =
  | { type: "point-hover"; peak: SpectrumPlotPeak | null }
  | { type: "peak-click"; peak: SpectrumPlotPeak }
  | { type: "range-select"; range: NumericRange };

export interface TicPlotProps {
  traces: readonly TicPlotTrace[];
  viewport: PlotViewport;
  rangeSelectionEnabled: boolean;
  /** Feature apex markers overlaid along the baseline (optional). */
  featureRug?: readonly FeatureMarker[];
  /** Shaded RT band(s) — typically the selected feature's elution extent (optional). */
  regions?: readonly RtRegion[];
  /** y-axis label + hover series name (default "TIC"). Set to "XIC" for chromatogram reuse. */
  yTitle?: string;
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
