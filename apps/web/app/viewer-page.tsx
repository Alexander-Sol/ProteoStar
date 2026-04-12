"use client";

import {
  type Spectrum,
} from "@msbrowser/imsp-core";
import {
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
  type ViewerStore
} from "@msbrowser/viewer-state";
import React, { useEffect, useState, useSyncExternalStore, startTransition } from "react";

import {
  type LoadedDataset,
  formatNumber,
  formatRange,
  loadViewerDataset,
  toSpectrumPlotPeaks,
  toTicPlotPoints,
  toViewport
} from "./viewer-controller";
import { ClientSpectrumPlot, ClientTicPlot } from "./viewer-plots";

const DEFAULT_SPECTRUM_RANGE = { min: 200, max: 1200 };

export function ViewerPage() {
  const [viewerStore] = useState<ViewerStore>(() => createViewerStore());
  const viewerState = useViewerState(viewerStore);
  const [loadedDataset, setLoadedDataset] = useState<LoadedDataset | null>(null);
  const [selectedSpectrum, setSelectedSpectrum] = useState<Spectrum | null>(null);
  const [hoveredTicPoint, setHoveredTicPoint] = useState<TicPlotPoint | null>(null);
  const [hoveredSpectrumPeak, setHoveredSpectrumPeak] =
    useState<SpectrumPlotPeak | null>(null);

  const ticViewport = toViewport(viewerState.ticPanel.range);
  const spectrumRange = viewerState.spectrumPanel.range ?? DEFAULT_SPECTRUM_RANGE;
  const spectrumViewport = toViewport(spectrumRange);

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

  useEffect(() => {
    return () => {
      loadedDataset?.dispose?.();
    };
  }, [loadedDataset]);

  const ticPoints = toTicPlotPoints(loadedDataset?.ticTrace ?? []);
  const spectrumPeaks = toSpectrumPlotPeaks(selectedSpectrum);

  return (
    <ViewerShell
      title="IMSP Viewer"
      subtitle="Open a local `.imsp` file to inspect the TIC and the selected scan spectrum."
      toolbar={
        <>
          <FileOpenButton
            onSelect={(file) =>
              void handleFileOpen(file, viewerStore, setLoadedDataset, setHoveredTicPoint, setHoveredSpectrumPeak)
            }
          />
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
            subtitle="Click to select a scan. Pin to drag-zoom retention time."
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
                viewerStore.getState().dispatch({ type: "panel/reset", panelId: "spectrum" });
                setHoveredSpectrumPeak(null);
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
            subtitle="The selected scan is reconstructed on demand. Pin to drag-zoom m/z."
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
                  value={formatRange(spectrumRange, 4)}
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

function useViewerState(viewerStore: ViewerStore) {
  return useSyncExternalStore(
    viewerStore.subscribe,
    viewerStore.getState,
    viewerStore.getState
  );
}

async function handleFileOpen(
  file: File,
  viewerStore: ViewerStore,
  setLoadedDataset: (dataset: LoadedDataset | null) => void,
  setHoveredTicPoint: (point: TicPlotPoint | null) => void,
  setHoveredSpectrumPeak: (peak: SpectrumPlotPeak | null) => void
): Promise<void> {
  viewerStore.getState().dispatch({ type: "dataset/load-started" });
  setHoveredTicPoint(null);
  setHoveredSpectrumPeak(null);

  try {
    const { loadedDataset, viewerDataset } = await loadViewerDataset(file);

    startTransition(() => {
      setLoadedDataset(loadedDataset);
      viewerStore.getState().dispatch({
        type: "dataset/load-succeeded",
        dataset: viewerDataset
      });
      if (viewerDataset.scanSummaries.length > 0) {
        viewerStore.getState().dispatch({ type: "selection/set-scan", scanIndex: 0 });
      }
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
