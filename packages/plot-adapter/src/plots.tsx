"use client";

import { useRef, type ReactElement } from "react";
import type { Layout, PlotData } from "plotly.js";
import Plot from "react-plotly.js";

import type {
  NumericRange,
  SlotIndex,
  SpectrumPlotProps,
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
// A mousedown–mouseup that moves more than this (px) is a drag, not a click.
const DRAG_THRESHOLD_PX = 4;

export function TicPlot(props: TicPlotProps): ReactElement {
  const { traces, viewport, rangeSelectionEnabled, onEvent } = props;

  const containerRef = useRef<HTMLDivElement>(null);
  const mousedownPos = useRef<{ x: number; y: number } | null>(null);
  const wasDrag = useRef(false);

  const data: PlotData[] = traces.flatMap((trace) => buildTicTraceData(trace));

  const layout: Partial<Layout> = {
    autosize: true,
    margin: TIC_MARGIN,
    dragmode: rangeSelectionEnabled ? "select" : "pan",
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
      gridcolor: "#dfe7f2",
      zeroline: false
    },
    showlegend: false,
    hovermode: "closest"
  };

  // Compute the effective x range so we can convert pixel → RT.
  // When viewport is unset Plotly auto-ranges to the data extents.
  function effectiveXRange(): { min: number; max: number } | null {
    if (viewport.xMin !== null && viewport.xMax !== null) {
      return { min: viewport.xMin, max: viewport.xMax };
    }
    const rts = traces.flatMap((t) => t.points.map((p) => p.retentionTime));
    if (rts.length === 0) return null;
    return { min: Math.min(...rts), max: Math.max(...rts) };
  }

  function handleMouseDown(e: React.MouseEvent) {
    mousedownPos.current = { x: e.clientX, y: e.clientY };
    wasDrag.current = false;
  }

  function handleMouseMove(e: React.MouseEvent) {
    if (!mousedownPos.current) return;
    const dx = Math.abs(e.clientX - mousedownPos.current.x);
    const dy = Math.abs(e.clientY - mousedownPos.current.y);
    if (dx > DRAG_THRESHOLD_PX || dy > DRAG_THRESHOLD_PX) {
      wasDrag.current = true;
    }
  }

  function handleClick(e: React.MouseEvent<HTMLDivElement>) {
    if (wasDrag.current) return;
    if (!containerRef.current) return;

    const rect = containerRef.current.getBoundingClientRect();
    const plotLeft = TIC_MARGIN.l;
    const plotRight = rect.width - TIC_MARGIN.r;
    const x = e.clientX - rect.left;

    // Ignore clicks in the axis / margin areas.
    if (x < plotLeft || x > plotRight) return;

    const range = effectiveXRange();
    if (!range) return;

    const fraction = (x - plotLeft) / (plotRight - plotLeft);
    const retentionTime = range.min + fraction * (range.max - range.min);
    onEvent({ type: "area-click", retentionTime });
  }

  return (
    <div
      ref={containerRef}
      style={{ width: "100%", height: "100%" }}
      onMouseDown={handleMouseDown}
      onMouseMove={handleMouseMove}
      onClick={handleClick}
    >
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
          const raw = readCustomData<TicCustomData>(event?.points?.[0]?.customdata);
          const point = raw ? (({ slotIndex: _s, ...p }) => p)(raw) as TicPlotPoint : null;
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
    mode: "lines",
    x: trace.points.map((p) => p.retentionTime),
    y: trace.points.map((p) => p.intensity),
    customdata: trace.points.map((p) => ({
      ...p,
      slotIndex: trace.slotIndex
    })) as unknown as PlotData["customdata"],
    line: { color: trace.color, width: 2 },
    hovertemplate: "RT %{x:.3f} min<br>TIC %{y:.0f}<extra></extra>"
  } as unknown as PlotData;

  return [lineTrace];
}

export function SpectrumPlot(props: SpectrumPlotProps): ReactElement {
  const { traces, viewport, rangeSelectionEnabled, onEvent } = props;

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
      gridcolor: "#dfe7f2",
      zeroline: false
    },
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
