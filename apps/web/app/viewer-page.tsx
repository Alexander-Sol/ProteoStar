"use client";

import {
  type Spectrum,
} from "@msbrowser/imsp-core";
import {
  type SlotIndex,
  type SpectrumPlotTrace,
  type TicPlotTrace,
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
import { SLOT_COLORS } from "./slot-colors";

const DEFAULT_SPECTRUM_RANGE = { min: 200, max: 1200 };

type DatasetPair = [LoadedDataset | null, LoadedDataset | null];
type SpectrumPair = [Spectrum | null, Spectrum | null];

export function ViewerPage() {
  const [viewerStore] = useState<ViewerStore>(() => createViewerStore());
  const viewerState = useViewerState(viewerStore);
  const [loadedDatasets, setLoadedDatasets] = useState<DatasetPair>([null, null]);
  const [selectedSpectra, setSelectedSpectra] = useState<SpectrumPair>([null, null]);
  const [hoveredTicPoint, setHoveredTicPoint] = useState<TicPlotPoint | null>(null);
  const [hoveredSpectrumPeak, setHoveredSpectrumPeak] = useState<
    SpectrumPlotTrace["peaks"][number] | null
  >(null);

  const ticViewport = toViewport(viewerState.ticPanel.range);
  const spectrumRange = viewerState.spectrumPanel.range ?? DEFAULT_SPECTRUM_RANGE;
  const spectrumViewport = toViewport(spectrumRange);

  const slot0 = viewerState.datasetSlots[0];
  const slot1 = viewerState.datasetSlots[1];

  const scanSummary0 =
    loadedDatasets[0] && slot0.selectedScanIndex !== null
      ? loadedDatasets[0].scanSummaries[slot0.selectedScanIndex] ?? null
      : null;
  const scanSummary1 =
    loadedDatasets[1] && slot1.selectedScanIndex !== null
      ? loadedDatasets[1].scanSummaries[slot1.selectedScanIndex] ?? null
      : null;

  // Load spectrum for slot 0 when its selected scan changes.
  useEffect(() => {
    const idx = slot0.selectedScanIndex;
    if (!loadedDatasets[0] || idx === null) {
      setSelectedSpectra((prev) => [null, prev[1]]);
      return;
    }
    let cancelled = false;
    void loadedDatasets[0].provider.getSpectrumForScan(idx).then((spectrum) => {
      if (!cancelled) setSelectedSpectra((prev) => [spectrum, prev[1]]);
    });
    return () => { cancelled = true; };
  }, [loadedDatasets[0], slot0.selectedScanIndex]);

  // Load spectrum for slot 1 when its selected scan changes.
  useEffect(() => {
    const idx = slot1.selectedScanIndex;
    if (!loadedDatasets[1] || idx === null) {
      setSelectedSpectra((prev) => [prev[0], null]);
      return;
    }
    let cancelled = false;
    void loadedDatasets[1].provider.getSpectrumForScan(idx).then((spectrum) => {
      if (!cancelled) setSelectedSpectra((prev) => [prev[0], spectrum]);
    });
    return () => { cancelled = true; };
  }, [loadedDatasets[1], slot1.selectedScanIndex]);

  // Dispose providers on unmount or dataset replacement.
  useEffect(() => {
    return () => { loadedDatasets[0]?.dispose?.(); };
  }, [loadedDatasets[0]]);
  useEffect(() => {
    return () => { loadedDatasets[1]?.dispose?.(); };
  }, [loadedDatasets[1]]);

  // Build TIC traces — one per loaded dataset.
  const ticTraces: TicPlotTrace[] = [];
  if (loadedDatasets[0]) {
    ticTraces.push({
      slotIndex: 0,
      points: toTicPlotPoints(loadedDatasets[0].ticTrace),
      selectedScanIndex: slot0.selectedScanIndex,
      color: SLOT_COLORS[0]
    });
  }
  if (loadedDatasets[1]) {
    ticTraces.push({
      slotIndex: 1,
      points: toTicPlotPoints(loadedDatasets[1].ticTrace),
      selectedScanIndex: slot1.selectedScanIndex,
      color: SLOT_COLORS[1]
    });
  }

  // Build spectrum traces — one per loaded spectrum.
  const spectrumTraces: SpectrumPlotTrace[] = [];
  if (selectedSpectra[0]) {
    spectrumTraces.push({
      slotIndex: 0,
      peaks: toSpectrumPlotPeaks(selectedSpectra[0]),
      color: SLOT_COLORS[0]
    });
  }
  if (selectedSpectra[1]) {
    spectrumTraces.push({
      slotIndex: 1,
      peaks: toSpectrumPlotPeaks(selectedSpectra[1]),
      color: SLOT_COLORS[1]
    });
  }

  const neitherLoaded =
    slot0.load.status !== "loading" &&
    slot0.load.status !== "ready" &&
    slot1.load.status !== "loading" &&
    slot1.load.status !== "ready";

  return (
    <ViewerShell
      title="IMSP Viewer"
      subtitle="Open up to two local `.imsp` files to compare their TIC and spectra."
      toolbar={
        <>
          <FileOpenButton
            color={SLOT_COLORS[0]}
            onSelect={(file) =>
              void handleFileOpen(file, 0, viewerStore, setLoadedDatasets, setHoveredTicPoint, setHoveredSpectrumPeak)
            }
          />
          <WorkspaceBadge
            label={loadedDatasets[0] ? loadedDatasets[0].fileName : slot0.load.status}
          />
          <FileOpenButton
            label="Open Dataset B (.imsp)"
            color={SLOT_COLORS[1]}
            onSelect={(file) =>
              void handleFileOpen(file, 1, viewerStore, setLoadedDatasets, setHoveredTicPoint, setHoveredSpectrumPeak)
            }
          />
          <WorkspaceBadge
            label={loadedDatasets[1] ? loadedDatasets[1].fileName : slot1.load.status}
          />
        </>
      }
    >
      {[slot0, slot1].map((slot, i) =>
        slot.load.status === "loading" ? (
          <StatusBanner tone="info" key={i}>
            Dataset {i + 1}: Reading file…
          </StatusBanner>
        ) : slot.load.status === "error" ? (
          <StatusBanner tone="error" key={i}>
            Dataset {i + 1}: {slot.load.errorMessage}
          </StatusBanner>
        ) : null
      )}
      {neitherLoaded ? (
        <StatusBanner tone="muted">
          Use the file buttons above to open one or two IMSP datasets.
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
                  label="A: Scan"
                  value={scanSummary0?.oneBasedScanNumber ?? "—"}
                />
                <MetricReadout
                  label="B: Scan"
                  value={scanSummary1?.oneBasedScanNumber ?? "—"}
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
        {ticTraces.length > 0 ? (
          <ClientTicPlot
            onEvent={(event) => {
              if (event.type === "point-click") {
                const { retentionTime } = event.point;
                const currentState = viewerStore.getState();
                (([0, 1]) as SlotIndex[]).forEach((slotIndex) => {
                  if (currentState.datasetSlots[slotIndex].load.status === "ready") {
                    viewerStore.getState().dispatch({
                      type: "selection/select-nearest-scan",
                      slotIndex,
                      retentionTime
                    });
                  }
                });
                if (!viewerStore.getState().spectrumPanel.pinned) {
                  viewerStore.getState().dispatch({ type: "panel/reset", panelId: "spectrum" });
                }
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
            traces={ticTraces}
            rangeSelectionEnabled={viewerState.ticPanel.pinned}
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
                  label="A: RT"
                  value={formatNumber(scanSummary0?.retentionTime, 3, "min")}
                />
                <MetricReadout
                  label="B: RT"
                  value={formatNumber(scanSummary1?.retentionTime, 3, "min")}
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
        {spectrumTraces.length > 0 ? (
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
            traces={spectrumTraces}
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
  slotIndex: SlotIndex,
  viewerStore: ViewerStore,
  setLoadedDatasets: React.Dispatch<React.SetStateAction<DatasetPair>>,
  setHoveredTicPoint: (point: TicPlotPoint | null) => void,
  setHoveredSpectrumPeak: (peak: null) => void
): Promise<void> {
  viewerStore.getState().dispatch({ type: "dataset/load-started", slotIndex });
  setHoveredTicPoint(null);
  setHoveredSpectrumPeak(null);

  try {
    const { loadedDataset, viewerDataset } = await loadViewerDataset(file);

    startTransition(() => {
      setLoadedDatasets((prev) => {
        const next: DatasetPair = [prev[0], prev[1]];
        next[slotIndex] = loadedDataset;
        return next;
      });
      viewerStore.getState().dispatch({
        type: "dataset/load-succeeded",
        slotIndex,
        dataset: viewerDataset
      });
      if (viewerDataset.scanSummaries.length > 0) {
        viewerStore.getState().dispatch({
          type: "selection/set-scan",
          slotIndex,
          scanIndex: 0
        });
      }
    });
  } catch (error) {
    const message = error instanceof Error ? error.message : "Failed to open IMSP file";
    startTransition(() => {
      setLoadedDatasets((prev) => {
        const next: DatasetPair = [prev[0], prev[1]];
        next[slotIndex] = null;
        return next;
      });
      viewerStore
        .getState()
        .dispatch({ type: "dataset/load-failed", slotIndex, errorMessage: message });
    });
  }
}
