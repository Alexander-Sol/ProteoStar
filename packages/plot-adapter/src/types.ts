export type SlotIndex = 0 | 1;

export interface PlotViewport {
  xMin: number | null;
  xMax: number | null;
  /**
   * Optional persisted y-axis bounds. When both are set the plot honors them instead of auto-fitting
   * the y-axis to the visible peaks — this is what holds an isotope envelope at a stable height while
   * stepping through neighboring scans (otherwise the y-axis refits per scan and the envelope shrinks
   * as a taller peak enters view). `null`/absent means "auto-fit y to the visible x-window".
   */
  yMin?: number | null;
  yMax?: number | null;
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

// A text label anchored at a spectrum peak (m/z + inferred charge for prominent MS1 peaks).
export interface PeakAnnotation {
  /** m/z of the peak the label points at (x position). */
  mz: number;
  /** Peak intensity (y position the label sits above). */
  intensity: number;
  /** Label text, e.g. `647.65 · z2`. */
  text: string;
}

// A feature-membership highlight anchored to an observed spectrum peak: a marker drawn at the peak's
// (m/z, intensity) apex in the owning feature's colour, instead of a full-height comb line. Used for
// the feature overlays (single selected feature + all features eluting at the scan).
export interface SpectrumPeakHighlight {
  mz: number;
  intensity: number;
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
  /** Plotly `uirevision`. While this value is unchanged, the user's interactive zoom/pan is
   *  preserved across re-renders (so e.g. clicking a point to load a spectrum doesn't snap the
   *  view back to the prop range). Change it to intentionally re-apply the range props (reframe).
   *  Leave undefined to keep the legacy behaviour (every render re-applies the range). */
  uirevision?: string | number;
  onEvent(event: TicPlotEvent): void;
}

export interface SpectrumPlotProps {
  traces: readonly SpectrumPlotTrace[];
  viewport: PlotViewport;
  rangeSelectionEnabled: boolean;
  /** Predicted isotope-peak m/z lines overlaid on the spectrum (optional). */
  envelope?: readonly EnvelopeLine[];
  /** Feature-membership markers drawn on matched peaks (optional) — the feature overlays use these
   *  instead of full-height comb lines. */
  highlights?: readonly SpectrumPeakHighlight[];
  /** Text labels on prominent peaks (m/z + inferred charge for MS1). Optional. */
  annotations?: readonly PeakAnnotation[];
  /** Plotly `uirevision`: while unchanged, the user's interactive zoom/pan is preserved across
   *  re-renders (e.g. stepping scans with the arrow keys). Change it to re-apply the range props. */
  uirevision?: string | number;
  onEvent(event: SpectrumPlotEvent): void;
}
