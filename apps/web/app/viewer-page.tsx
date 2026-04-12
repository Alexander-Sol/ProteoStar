"use client";

import {
  createImspDatasetProvider,
  type DatasetProvider,
  type ScanSummary,
  type Spectrum,
  type TicPoint
} from "@msbrowser/imsp-core";
import {
  type PlotViewport,
  type SpectrumPlotPeak,
  type TicPlotPoint
} from "@msbrowser/plot-adapter";
import {
  FileOpenButton,
  MetricReadout,
  Panel,
  PanelActionButton,
  PanelHeader,
  StatusBanner,
  ViewerShell,
  WorkspaceBadge
} from "@msbrowser/ui";
import {
  createViewerStore,
  type NumericRange,
  type ViewerDataset
} from "@msbrowser/viewer-state";
import { useEffect, useMemo, useState, useSyncExternalStore, startTransition } from "react";

import { ClientSpectrumPlot, ClientTicPlot } from "./viewer-plots";

const viewerStore = createViewerStore();

interface LoadedDataset {
  fileName: string;
  provider: DatasetProvider;
  metadata: Awaited<ReturnType<DatasetProvider["getMetadata"]>>;
  scanSummaries: readonly ScanSummary[];
  ticTrace: readonly TicPoint[];
}

export function ViewerPage() {
  const viewerState = useViewerState();
  const [loadedDataset, setLoadedDataset] = useState<LoadedDataset | null>(null);
  const [selectedSpectrum, setSelectedSpectrum] = useState<Spectrum | null>(null);
  const [hoveredTicPoint, setHoveredTicPoint] = useState<TicPlotPoint | null>(null);
  const [hoveredSpectrumPeak, setHoveredSpectrumPeak] =
    useState<SpectrumPlotPeak | null>(null);

  const ticViewport = toViewport(viewerState.ticPanel.range);
  const spectrumViewport = toViewport(viewerState.spectrumPanel.range);

  const selectedScanSummary =
    loadedDataset && viewerState.selectedScanIndex !== null
      ? loadedDataset.scanSummaries[viewerState.selectedScanIndex] ?? null
      : null;

  useEffect(() => {
    const selectedScanIndex = viewerState.selectedScanIndex;
    if (!loadedDataset || selectedScanIndex === null) {
      setSelectedSpectrum(null);
      return;
    }

    let cancelled = false;
    void loadedDataset.provider.getSpectrumForScan(selectedScanIndex).then((spectrum) => {
      if (!cancelled) {
        setSelectedSpectrum(spectrum);
      }
    });

    return () => {
      cancelled = true;
    };
  }, [loadedDataset, viewerState.selectedScanIndex]);

  const ticPoints = useMemo<readonly TicPlotPoint[]>(() => {
    return (loadedDataset?.ticTrace ?? []).map((point) => ({
      scanIndex: point.scanIndex,
      retentionTime: point.retentionTime,
      intensity: point.tic
    }));
  }, [loadedDataset]);

  const spectrumPeaks = useMemo<readonly SpectrumPlotPeak[]>(() => {
    return (selectedSpectrum?.peaks ?? []).map((peak) => ({
      mz: peak.mz,
      intensity: peak.intensity
    }));
  }, [selectedSpectrum]);

  return (
    <ViewerShell
      title="Interactive IMSP Viewer"
      subtitle="Open a local `.imsp` file, inspect the TIC in the top panel, and click a scan to reconstruct its spectrum in the lower panel."
      toolbar={
        <>
          <FileOpenButton onSelect={(file) => void handleFileOpen(file, setLoadedDataset)} />
          <WorkspaceBadge
            label={loadedDataset ? loadedDataset.fileName : viewerState.dataset.status}
          />
        </>
      }
    >
      {viewerState.dataset.status === "loading" ? (
        <StatusBanner tone="info">Reading the selected IMSP file…</StatusBanner>
      ) : null}
      {viewerState.dataset.status === "error" ? (
        <StatusBanner tone="error">{viewerState.dataset.errorMessage}</StatusBanner>
      ) : null}
      {!loadedDataset && viewerState.dataset.status !== "loading" ? (
        <StatusBanner tone="muted">
          Use the file button above to open an IMSP fixture or your own local file.
        </StatusBanner>
      ) : null}

      <Panel
        header={
          <PanelHeader
            title="Total Ion Chromatogram"
            subtitle="Click a point to select the nearest scan. Pin the panel to enable drag-to-zoom on retention time."
            readouts={
              <>
                <MetricReadout
                  label="Selected Scan"
                  value={selectedScanSummary?.oneBasedScanNumber ?? "None"}
                />
                <MetricReadout
                  label="Hover RT"
                  value={formatNumber(hoveredTicPoint?.retentionTime, 3, "min")}
                />
                <MetricReadout
                  label="View"
                  value={formatRange(viewerState.ticPanel.range, 3, "min")}
                />
              </>
            }
            actions={
              <>
                <PanelActionButton
                  onClick={() =>
                    viewerStore
                      .getState()
                      .dispatch({ type: "panel/toggle-pinned", panelId: "tic" })
                  }
                  pressed={viewerState.ticPanel.pinned}
                >
                  Pin Zoom
                </PanelActionButton>
                <PanelActionButton
                  onClick={() =>
                    viewerStore.getState().dispatch({ type: "panel/reset", panelId: "tic" })
                  }
                >
                  Reset
                </PanelActionButton>
              </>
            }
          />
        }
      >
        {loadedDataset ? (
          <ClientTicPlot
            onEvent={(event) => {
              if (event.type === "point-click") {
                viewerStore.getState().dispatch({
                  type: "selection/select-nearest-scan",
                  retentionTime: event.point.retentionTime
                });
                return;
              }

              if (event.type === "point-hover") {
                setHoveredTicPoint(event.point);
                return;
              }

              viewerStore.getState().dispatch({
                type: "panel/zoom",
                panelId: "tic",
                range: event.range
              });
            }}
            points={ticPoints}
            rangeSelectionEnabled={viewerState.ticPanel.pinned}
            selectedScanIndex={viewerState.selectedScanIndex}
            viewport={ticViewport}
          />
        ) : (
          <StatusBanner tone="muted">
            The TIC will appear here after you load an IMSP dataset.
          </StatusBanner>
        )}
      </Panel>

      <Panel
        header={
          <PanelHeader
            title="Mass Spectrum"
            subtitle="The selected scan is reconstructed on demand from the IMSP dataset. Pin the panel to enable drag-to-zoom on m/z."
            readouts={
              <>
                <MetricReadout
                  label="Retention Time"
                  value={formatNumber(selectedScanSummary?.retentionTime, 3, "min")}
                />
                <MetricReadout
                  label="Hover m/z"
                  value={formatNumber(hoveredSpectrumPeak?.mz, 4)}
                />
                <MetricReadout
                  label="View"
                  value={formatRange(viewerState.spectrumPanel.range, 4)}
                />
              </>
            }
            actions={
              <>
                <PanelActionButton
                  onClick={() =>
                    viewerStore
                      .getState()
                      .dispatch({ type: "panel/toggle-pinned", panelId: "spectrum" })
                  }
                  pressed={viewerState.spectrumPanel.pinned}
                >
                  Pin Zoom
                </PanelActionButton>
                <PanelActionButton
                  onClick={() =>
                    viewerStore
                      .getState()
                      .dispatch({ type: "panel/reset", panelId: "spectrum" })
                  }
                >
                  Reset
                </PanelActionButton>
              </>
            }
          />
        }
      >
        {selectedSpectrum ? (
          <ClientSpectrumPlot
            onEvent={(event) => {
              if (event.type === "point-hover") {
                setHoveredSpectrumPeak(event.peak);
                return;
              }

              viewerStore.getState().dispatch({
                type: "panel/zoom",
                panelId: "spectrum",
                range: event.range
              });
            }}
            peaks={spectrumPeaks}
            rangeSelectionEnabled={viewerState.spectrumPanel.pinned}
            viewport={spectrumViewport}
          />
        ) : (
          <StatusBanner tone="muted">
            Select a scan in the TIC to reconstruct and display its spectrum.
          </StatusBanner>
        )}
      </Panel>
    </ViewerShell>
  );
}

function useViewerState() {
  return useSyncExternalStore(
    viewerStore.subscribe,
    viewerStore.getState,
    viewerStore.getState
  );
}

async function handleFileOpen(
  file: File,
  setLoadedDataset: (dataset: LoadedDataset | null) => void
): Promise<void> {
  viewerStore.getState().dispatch({ type: "dataset/load-started" });

  try {
    const buffer = await file.arrayBuffer();
    const provider = createImspDatasetProvider(buffer);
    const [metadata, scanSummaries, ticTrace] = await Promise.all([
      provider.getMetadata(),
      provider.getScanSummaries(),
      provider.getTicTrace()
    ]);

    const viewerDataset: ViewerDataset = {
      metadata: {
        retentionTimeRange: metadata.retentionTimeRange,
        mzRange: metadata.mzRange,
        scanCount: metadata.scanCount
      },
      scanSummaries: scanSummaries.map((scan) => ({
        scanIndex: scan.scanIndex,
        oneBasedScanNumber: scan.oneBasedScanNumber,
        retentionTime: scan.retentionTime,
        tic: scan.tic
      }))
    };

    startTransition(() => {
      setLoadedDataset({
        fileName: file.name,
        provider,
        metadata,
        scanSummaries,
        ticTrace
      });
      viewerStore.getState().dispatch({
        type: "dataset/load-succeeded",
        dataset: viewerDataset
      });
    });
  } catch (error) {
    const message = error instanceof Error ? error.message : "Failed to open IMSP file";
    startTransition(() => {
      setLoadedDataset(null);
      viewerStore
        .getState()
        .dispatch({ type: "dataset/load-failed", errorMessage: message });
    });
  }
}

function toViewport(range: NumericRange | null): PlotViewport {
  return {
    xMin: range?.min ?? null,
    xMax: range?.max ?? null
  };
}

function formatNumber(
  value: number | undefined,
  decimals: number,
  suffix?: string
): string {
  if (value === undefined) {
    return "-";
  }

  return `${value.toFixed(decimals)}${suffix ? ` ${suffix}` : ""}`;
}

function formatRange(
  range: NumericRange | null,
  decimals: number,
  suffix?: string
): string {
  if (!range) {
    return "Full range";
  }

  return `${formatNumber(range.min, decimals, suffix)} to ${formatNumber(range.max, decimals, suffix)}`;
}
