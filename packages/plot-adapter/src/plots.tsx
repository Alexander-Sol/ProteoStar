"use client";

import type { ReactElement } from "react";
import type { Layout, PlotData } from "plotly.js";
import Plot from "react-plotly.js";

import type {
  NumericRange,
  SpectrumPlotProps,
  TicPlotPoint,
  TicPlotProps
} from "./types";

interface PlotPointEvent {
  points?: Array<{ customdata?: unknown }>;
}

interface PlotSelectionEvent {
  range?: {
    x?: unknown;
  };
}

export function TicPlot(props: TicPlotProps): ReactElement {
  const { points, selectedScanIndex, viewport, rangeSelectionEnabled, onEvent } = props;
  const data: PlotData[] = [
    {
      type: "scattergl",
      mode: "lines+markers",
      x: points.map((point) => point.retentionTime),
      y: points.map((point) => point.intensity),
      customdata: points.map((point) => ({ ...point })) as unknown as PlotData["customdata"],
      line: { color: "#2b6cb0", width: 2 },
      marker: {
        color: points.map((point) =>
          point.scanIndex === selectedScanIndex ? "#d9485f" : "#2b6cb0"
        ),
        size: points.map((point) => (point.scanIndex === selectedScanIndex ? 8 : 6))
      },
      hovertemplate: "RT %{x:.3f} min<br>TIC %{y:.0f}<extra></extra>"
    } as unknown as PlotData
  ];

  const layout: Partial<Layout> = {
    autosize: true,
    margin: { l: 56, r: 18, t: 20, b: 44 },
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
        const point = readCustomData<TicPlotPoint>(event?.points?.[0]?.customdata);
        if (point) {
          onEvent({ type: "point-click", point });
        }
      }}
      onHover={(event: PlotPointEvent) => {
        const point = readCustomData<TicPlotPoint>(event?.points?.[0]?.customdata);
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
  );
}

export function SpectrumPlot(props: SpectrumPlotProps): ReactElement {
  const { peaks, viewport, rangeSelectionEnabled, onEvent } = props;
  const data: PlotData[] = [
    {
      type: "bar",
      x: peaks.map((peak) => peak.mz),
      y: peaks.map((peak) => peak.intensity),
      customdata: peaks.map((peak) => ({ ...peak })) as unknown as PlotData["customdata"],
      marker: { color: "#4c51bf" },
      hovertemplate: "m/z %{x:.4f}<br>Intensity %{y:.0f}<extra></extra>"
    } as unknown as PlotData
  ];

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
    hovermode: "closest"
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
        const peak = readCustomData<SpectrumPlotProps["peaks"][number]>(
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
