"use client";

import { type ReactElement } from "react";
import type { Layout, PlotData, Shape } from "plotly.js";
import Plot from "react-plotly.js";

import type {
  EnvelopeLine,
  FeatureMarker,
  NumericRange,
  PeakAnnotation,
  PlotViewport,
  RtRegion,
  SlotIndex,
  SpectrumPeakHighlight,
  SpectrumPlotProps,
  TicPlotPoint,
  TicPlotProps,
  TicPlotTrace
} from "./types";
import { resolveSpectrumYRange } from "./viewport";

interface PlotPointEvent {
  points?: Array<{ customdata?: unknown }>;
}

interface PlotSelectionEvent {
  range?: {
    x?: unknown;
  };
}

type TicCustomData = TicPlotPoint & { slotIndex: SlotIndex };

// Margins must match the layout margin values below.
const TIC_MARGIN = { l: 56, r: 18, t: 20, b: 44 };

export function TicPlot(props: TicPlotProps): ReactElement {
  const { traces, viewport, rangeSelectionEnabled, featureRug, regions, onEvent } = props;
  const yTitle = props.yTitle ?? "TIC";

  const data: PlotData[] = traces.flatMap((trace) => buildTicTraceData(trace, yTitle));
  if (featureRug && featureRug.length > 0) {
    data.push(buildFeatureRugData(featureRug));
  }

  const layout: Partial<Layout> = {
    autosize: true,
    margin: TIC_MARGIN,
    dragmode: rangeSelectionEnabled ? "select" : "zoom",
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "#ffffff",
    font: { color: "#24364d", family: "Inter, Arial, sans-serif" },
    xaxis: {
      title: { text: "Retention time (min)" },
      range: toPlotlyRange(viewport),
      gridcolor: "#dfe7f2",
      zeroline: false
    },
    yaxis: {
      title: { text: yTitle },
      range: visibleYRange(traces, viewport),
      gridcolor: "#dfe7f2",
      zeroline: false
    },
    // Preserve interactive zoom/pan across re-renders while this is stable; the caller bumps it
    // to intentionally re-apply the range props (reframe). Undefined → legacy re-apply-every-render.
    uirevision: props.uirevision,
    shapes: buildRegionShapes(regions),
    showlegend: false,
    hovermode: "closest"
  };

  return (
    <div style={{ width: "100%", height: "100%" }}>
      <Plot
        data={data}
        layout={layout}
        config={{
          responsive: true,
          displaylogo: false,
          modeBarButtonsToRemove: [
            "lasso2d",
            "zoomIn2d",
            "zoomOut2d",
            "autoScale2d",
            "toggleSpikelines"
          ]
        }}
        style={{ width: "100%", height: "100%" }}
        useResizeHandler
        onClick={(event: PlotPointEvent) => {
          // Plotly's own click event (reliable even under the pan drag-layer, which
          // swallows DOM clicks): a feature-rug point carries `featureIndex`, a TIC
          // line point carries `retentionTime`.
          const cd = readCustomData<{ featureIndex?: number; retentionTime?: number }>(
            event?.points?.[0]?.customdata
          );
          if (!cd) return;
          if (typeof cd.featureIndex === "number") {
            onEvent({ type: "feature-click", featureIndex: cd.featureIndex });
          } else if (typeof cd.retentionTime === "number") {
            onEvent({ type: "area-click", retentionTime: cd.retentionTime });
          }
        }}
        onHover={(event: PlotPointEvent) => {
          const raw = readCustomData<TicCustomData>(event?.points?.[0]?.customdata);
          const point =
            raw && !("featureIndex" in raw)
              ? (({ slotIndex: _s, ...p }) => p)(raw) as TicPlotPoint
              : null;
          onEvent({ type: "point-hover", point });
        }}
        onUnhover={() => {
          onEvent({ type: "point-hover", point: null });
        }}
        onSelected={(event: PlotSelectionEvent) => {
          const range = readRange(event?.range?.x);
          if (range) {
            onEvent({ type: "range-select", range });
          }
        }}
      />
    </div>
  );
}

function buildTicTraceData(trace: TicPlotTrace, yTitle: string): PlotData[] {
  const hoverName = trace.label ?? yTitle;
  const lineTrace = {
    type: "scattergl",
    // Invisible (transparent) markers over the line give Plotly a reliable per-point
    // click/hover hit-target — pure "lines" click hit-testing is unreliable, and the
    // pan drag-layer swallows raw DOM clicks. The line itself is drawn via `line`.
    mode: "lines+markers",
    x: trace.points.map((p) => p.retentionTime),
    y: trace.points.map((p) => p.intensity),
    customdata: trace.points.map((p) => ({
      ...p,
      slotIndex: trace.slotIndex
    })) as unknown as PlotData["customdata"],
    line: { color: trace.color, width: 2 },
    marker: { size: 8, color: "rgba(0,0,0,0)" },
    hovertemplate: `RT %{x:.3f} min<br>${hoverName} %{y:.0f}<extra></extra>`
  } as unknown as PlotData;

  return [lineTrace];
}

/** One markers trace placing a feature glyph at each apex RT along the baseline. */
function buildFeatureRugData(markers: readonly FeatureMarker[]): PlotData {
  return {
    type: "scattergl",
    mode: "markers",
    x: markers.map((m) => m.retentionTime),
    y: markers.map(() => 0),
    customdata: markers.map((m) => ({
      featureIndex: m.featureIndex,
      label: m.label
    })) as unknown as PlotData["customdata"],
    marker: {
      symbol: "triangle-up",
      size: 9,
      color: markers.map((m) => m.color),
      line: { width: 0.5, color: "#1b2a3d" }
    },
    hovertemplate: "%{customdata.label}<extra></extra>"
  } as unknown as PlotData;
}

/** Shaded vertical RT bands for the given regions (drawn behind the traces). */
function buildRegionShapes(regions: readonly RtRegion[] | undefined): Partial<Shape>[] {
  if (!regions || regions.length === 0) return [];
  return regions.map((r) => ({
    type: "rect" as const,
    xref: "x" as const,
    yref: "paper" as const,
    x0: r.min,
    x1: r.max,
    y0: 0,
    y1: 1,
    fillcolor: r.color,
    line: { width: 0 },
    layer: "below" as const
  }));
}

/** Vertical lines at predicted isotope-peak m/z positions over a spectrum. A `detected` tooth (an
 *  observed peak sits at that position in the displayed scan) is drawn solid and opaque; a
 *  predicted-only tooth is a faint dotted line. */
function buildEnvelopeShapes(envelope: readonly EnvelopeLine[] | undefined): Partial<Shape>[] {
  if (!envelope || envelope.length === 0) return [];
  return envelope.map((e) => ({
    type: "line" as const,
    xref: "x" as const,
    yref: "paper" as const,
    x0: e.mz,
    x1: e.mz,
    y0: 0,
    y1: 1,
    line: {
      color: e.color,
      width: e.detected ? 1.8 : 1,
      dash: e.detected ? ("solid" as const) : ("dot" as const)
    },
    opacity: e.detected ? 0.95 : 0.4,
    layer: "below" as const
  }));
}

// m/z half-width of a feature-highlight rectangle — a touch wider than the 0.001 peak bar so the
// colour shows around each peak. Fixed in data units, so it's most visible at inspection zoom.
const HIGHLIGHT_WIDTH = 0.06;

/** Feature-membership highlights: a semi-transparent colored rectangle behind each matched peak (at
 *  the peak's true m/z, as tall as the peak), coloured by the owning feature — drawn behind the peaks
 *  so the real peak stays visible on top. One `bar` trace (cheap even for many features). The peak
 *  `{mz, intensity}` rides in customdata so a click still resolves to a peak (e.g. a walkthrough seed). */
function buildHighlightTrace(
  highlights: readonly SpectrumPeakHighlight[] | undefined
): PlotData {
  const hs = highlights ?? [];
  return {
    type: "bar",
    x: hs.map((h) => h.mz),
    y: hs.map((h) => h.intensity),
    width: HIGHLIGHT_WIDTH,
    marker: { color: hs.map((h) => h.color), line: { width: 0 } },
    opacity: 0.5,
    customdata: hs.map((h) => ({
      mz: h.mz,
      intensity: h.intensity
    })) as unknown as PlotData["customdata"],
    hoverinfo: "skip"
  } as unknown as PlotData;
}

/** Text labels above prominent peaks (m/z, with the inferred charge on a line beneath). Drawn
 *  arrow-less and horizontal (parallel to the x-axis), centered over and sitting just above the
 *  peak apex. */
function buildPeakAnnotations(
  annotations: readonly PeakAnnotation[] | undefined
): Partial<Layout>["annotations"] {
  if (!annotations || annotations.length === 0) return [];
  return annotations.map((a) => ({
    x: a.mz,
    y: a.intensity,
    text: a.text,
    showarrow: false,
    xanchor: "center" as const,
    yanchor: "bottom" as const,
    textangle: "0",
    align: "center" as const,
    font: { size: 10, color: "#24364d", family: "Inter, Arial, sans-serif" }
  }));
}

export function SpectrumPlot(props: SpectrumPlotProps): ReactElement {
  const { traces, viewport, rangeSelectionEnabled, envelope, highlights, annotations, onEvent } =
    props;
  const allPeaks = traces.flatMap((t) => t.peaks);

  const peakTraces: PlotData[] = traces.flatMap((trace) => {
    // The 0.001-m/z-wide bars are the visual, but far too thin to hover or click. Overlay an
    // invisible wide-marker scatter at each peak apex to give Plotly a reliable hover/click
    // hit-target (mirrors the TIC's line+markers approach). The bar itself skips hover so the
    // tooltip and click both resolve to the marker trace's customdata.
    const bar = {
      type: "bar",
      x: trace.peaks.map((p) => p.mz),
      y: trace.peaks.map((p) => p.intensity),
      marker: { color: trace.color },
      width: 0.001,
      hoverinfo: "skip"
    } as unknown as PlotData;
    const hit = {
      type: "scattergl",
      mode: "markers",
      x: trace.peaks.map((p) => p.mz),
      y: trace.peaks.map((p) => p.intensity),
      customdata: trace.peaks.map((p) => ({ ...p })) as unknown as PlotData["customdata"],
      marker: { size: 12, color: "rgba(0,0,0,0)" },
      hovertemplate: "m/z %{x:.4f}<br>Intensity %{y:.0f}<extra></extra>"
    } as unknown as PlotData;
    return [bar, hit];
  });
  // Highlight trace FIRST so the colored backing rectangles render BEHIND the peaks. Always present
  // (empty when there's nothing to mark) so the trace COUNT stays constant across scan steps — a
  // varying count disrupts Plotly's `uirevision`, which holds the user's zoom fixed while stepping.
  const data: PlotData[] = [buildHighlightTrace(highlights), ...peakTraces];

  const layout: Partial<Layout> = {
    autosize: true,
    margin: { l: 56, r: 18, t: 20, b: 44 },
    dragmode: rangeSelectionEnabled ? "select" : "pan",
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "#ffffff",
    font: { color: "#24364d", family: "Inter, Arial, sans-serif" },
    xaxis: {
      title: { text: "m/z" },
      // No explicit `autorange`: supplying it every render fights `uirevision`, which is what keeps a
      // user's manual zoom fixed across scan steps. A reframe bumps uirevision and re-fits from
      // `range: undefined` on its own.
      range: toPlotlyRange(viewport),
      gridcolor: "#dfe7f2",
      zeroline: false
    },
    yaxis: {
      title: { text: "Intensity" },
      // Honor a persisted y-range when present (holds envelope height stable across scan steps),
      // else auto-fit to the visible x-window.
      range: resolveSpectrumYRange(viewport, allPeaks),
      gridcolor: "#dfe7f2",
      zeroline: false
    },
    // Preserve interactive zoom/pan across re-renders (e.g. arrow-key scan stepping) while stable;
    // the caller bumps it to intentionally re-apply the range props (reframe).
    uirevision: props.uirevision,
    shapes: buildEnvelopeShapes(envelope),
    annotations: buildPeakAnnotations(annotations),
    showlegend: false,
    hovermode: "closest",
    barmode: "overlay"
  };

  return (
    <Plot
      data={data}
      layout={layout}
      config={{
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: [
          "lasso2d",
          "zoomIn2d",
          "zoomOut2d",
          "autoScale2d",
          "toggleSpikelines"
        ]
      }}
      style={{ width: "100%", height: "100%" }}
      useResizeHandler
      onClick={(event: PlotPointEvent) => {
        // A clicked bar carries its peak `{mz, intensity}` in customdata — the seed pick for the
        // feature-finding walkthrough.
        const peak = readCustomData<SpectrumPlotProps["traces"][number]["peaks"][number]>(
          event?.points?.[0]?.customdata
        );
        if (peak) onEvent({ type: "peak-click", peak });
      }}
      onHover={(event: PlotPointEvent) => {
        const peak = readCustomData<SpectrumPlotProps["traces"][number]["peaks"][number]>(
          event?.points?.[0]?.customdata
        );
        onEvent({ type: "point-hover", peak: peak ?? null });
      }}
      onUnhover={() => {
        onEvent({ type: "point-hover", peak: null });
      }}
      onSelected={(event: PlotSelectionEvent) => {
        const range = readRange(event?.range?.x);
        if (range) {
          onEvent({ type: "range-select", range });
        }
      }}
      onRelayout={(event: Record<string, unknown>) => {
        // Report the new visible m/z window so overlays can react to zoom/pan. Autorange reset
        // (double-click) sends null; an explicit zoom sends the [min, max] bounds.
        if (event["xaxis.autorange"]) {
          onEvent({ type: "xrange-change", range: null });
          return;
        }
        const min = event["xaxis.range[0]"];
        const max = event["xaxis.range[1]"];
        if (typeof min === "number" && typeof max === "number") {
          onEvent({ type: "xrange-change", range: { min, max } });
        }
      }}
    />
  );
}

/** Compute the y-axis range from points visible within the current x viewport. When any trace
 *  opts in via `yScale`, only those traces drive the range — so a summed-XIC overlay zooms the
 *  axis to its own abundance and the much taller TIC (drawn on the same absolute scale) simply
 *  runs off the top of the view. */
function visibleYRange(
  traces: readonly TicPlotTrace[],
  viewport: PlotViewport
): [number, number] | undefined {
  const scalers = traces.some((t) => t.yScale) ? traces.filter((t) => t.yScale) : traces;
  const { xMin, xMax } = viewport;
  const points = scalers.flatMap((t) =>
    xMin !== null && xMax !== null
      ? t.points.filter((p) => p.retentionTime >= xMin && p.retentionTime <= xMax)
      : t.points
  );
  if (points.length === 0) return undefined;
  const maxY = Math.max(...points.map((p) => p.intensity));
  return [0, maxY * 1.05];
}

function toPlotlyRange(viewport: { xMin: number | null; xMax: number | null }):
  | [number, number]
  | undefined {
  return viewport.xMin !== null && viewport.xMax !== null
    ? [viewport.xMin, viewport.xMax]
    : undefined;
}

function readRange(values: unknown): NumericRange | null {
  if (!Array.isArray(values) || values.length !== 2) {
    return null;
  }

  const [left, right] = values;
  if (typeof left !== "number" || typeof right !== "number") {
    return null;
  }

  return left <= right ? { min: left, max: right } : { min: right, max: left };
}

function readCustomData<T>(value: unknown): T | null {
  if (!value || typeof value !== "object") {
    return null;
  }

  return value as T;
}
