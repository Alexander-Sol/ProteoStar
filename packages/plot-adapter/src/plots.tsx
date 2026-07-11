"use client";

import { type ReactElement } from "react";
import type { Layout, PlotData, Shape } from "plotly.js";
import Plot from "react-plotly.js";

import type {
  EnvelopeLine,
  FeatureMarker,
  NumericRange,
  PlotViewport,
  RtRegion,
  SlotIndex,
  SpectrumPlotProps,
  SpectrumPlotTrace,
  TicPlotPoint,
  TicPlotProps,
  TicPlotTrace
} from "./types";

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

  const data: PlotData[] = traces.flatMap((trace) => buildTicTraceData(trace));
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
      title: { text: "TIC" },
      range: visibleYRange(traces, viewport),
      gridcolor: "#dfe7f2",
      zeroline: false
    },
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

function buildTicTraceData(trace: TicPlotTrace): PlotData[] {
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
    hovertemplate: "RT %{x:.3f} min<br>TIC %{y:.0f}<extra></extra>"
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

/** Vertical lines at predicted isotope-peak m/z positions over a spectrum. */
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
    line: { color: e.color, width: 1, dash: "dot" as const },
    layer: "below" as const
  }));
}

export function SpectrumPlot(props: SpectrumPlotProps): ReactElement {
  const { traces, viewport, rangeSelectionEnabled, envelope, onEvent } = props;

  const data: PlotData[] = traces.map((trace) => ({
    type: "bar",
    x: trace.peaks.map((p) => p.mz),
    y: trace.peaks.map((p) => p.intensity),
    customdata: trace.peaks.map((p) => ({ ...p })) as unknown as PlotData["customdata"],
    marker: { color: trace.color },
    width: 0.001,
    hovertemplate: "m/z %{x:.4f}<br>Intensity %{y:.0f}<extra></extra>"
  } as unknown as PlotData));

  const layout: Partial<Layout> = {
    autosize: true,
    margin: { l: 56, r: 18, t: 20, b: 44 },
    dragmode: rangeSelectionEnabled ? "select" : "pan",
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "#ffffff",
    font: { color: "#24364d", family: "Inter, Arial, sans-serif" },
    xaxis: {
      title: { text: "m/z" },
      range: toPlotlyRange(viewport),
      gridcolor: "#dfe7f2",
      zeroline: false
    },
    yaxis: {
      title: { text: "Intensity" },
      range: visibleSpectrumYRange(traces, viewport),
      gridcolor: "#dfe7f2",
      zeroline: false
    },
    shapes: buildEnvelopeShapes(envelope),
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
    />
  );
}

/** Compute the y-axis range from points visible within the current x viewport. */
function visibleYRange(
  traces: readonly TicPlotTrace[],
  viewport: PlotViewport
): [number, number] | undefined {
  const { xMin, xMax } = viewport;
  const points = traces.flatMap((t) =>
    xMin !== null && xMax !== null
      ? t.points.filter((p) => p.retentionTime >= xMin && p.retentionTime <= xMax)
      : t.points
  );
  if (points.length === 0) return undefined;
  const maxY = Math.max(...points.map((p) => p.intensity));
  return [0, maxY * 1.05];
}

/** Compute the y-axis range from peaks visible within the current x viewport. */
function visibleSpectrumYRange(
  traces: readonly SpectrumPlotTrace[],
  viewport: PlotViewport
): [number, number] | undefined {
  const { xMin, xMax } = viewport;
  const peaks = traces.flatMap((t) =>
    xMin !== null && xMax !== null
      ? t.peaks.filter((p) => p.mz >= xMin && p.mz <= xMax)
      : t.peaks
  );
  if (peaks.length === 0) return undefined;
  const maxY = Math.max(...peaks.map((p) => p.intensity));
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
