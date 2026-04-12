import { ViewerStateError } from "./errors";
import type {
  NumericRange,
  PanelId,
  PanelState,
  ViewerCommand,
  ViewerScanSummary,
  ViewerState
} from "./types";

export function reduceViewerState(
  state: ViewerState,
  command: ViewerCommand
): ViewerState {
  switch (command.type) {
    case "dataset/load-started":
      return {
        ...state,
        dataset: {
          status: "loading",
          dataset: null,
          errorMessage: null
        },
        selectedScanIndex: null,
        ticPanel: createDefaultPanelState("tic"),
        spectrumPanel: createDefaultPanelState("spectrum")
      };

    case "dataset/load-succeeded":
      return {
        ...state,
        dataset: {
          status: "ready",
          dataset: command.dataset,
          errorMessage: null
        },
        selectedScanIndex: null,
        ticPanel: createDefaultPanelState("tic"),
        spectrumPanel: createDefaultPanelState("spectrum")
      };

    case "dataset/load-failed":
      return {
        ...state,
        dataset: {
          status: "error",
          dataset: null,
          errorMessage: command.errorMessage
        },
        selectedScanIndex: null,
        ticPanel: createDefaultPanelState("tic"),
        spectrumPanel: createDefaultPanelState("spectrum")
      };

    case "selection/set-scan":
      if (command.scanIndex === null) {
        return {
          ...state,
          selectedScanIndex: null
        };
      }

      assertDatasetReady(state);
      assertScanIndexExists(state, command.scanIndex);

      return {
        ...state,
        selectedScanIndex: command.scanIndex
      };

    case "selection/select-nearest-scan": {
      assertDatasetReady(state);
      const nearestScan = findNearestScan(
        state.dataset.dataset.scanSummaries,
        command.retentionTime
      );

      return {
        ...state,
        selectedScanIndex: nearestScan?.scanIndex ?? null
      };
    }

    case "panel/set-active":
      return {
        ...state,
        activePanel: command.panelId
      };

    case "panel/set-pinned":
      return updatePanel(state, command.panelId, {
        ...getPanel(state, command.panelId),
        pinned: command.pinned
      });

    case "panel/toggle-pinned": {
      const panel = getPanel(state, command.panelId);
      return updatePanel(state, command.panelId, {
        ...panel,
        pinned: !panel.pinned
      });
    }

    case "panel/zoom": {
      const panel = getPanel(state, command.panelId);
      if (!panel.pinned) {
        return state;
      }

      return updatePanel(state, command.panelId, {
        ...panel,
        range: normalizeRange(command.range)
      });
    }

    case "panel/reset":
      return updatePanel(state, command.panelId, {
        ...getPanel(state, command.panelId),
        range: null
      });
  }
}

function findNearestScan(
  scanSummaries: readonly ViewerScanSummary[],
  retentionTime: number
): ViewerScanSummary | null {
  let nearestScan: ViewerScanSummary | null = null;
  let smallestDistance = Number.POSITIVE_INFINITY;

  for (const scan of scanSummaries) {
    const distance = Math.abs(scan.retentionTime - retentionTime);
    if (distance < smallestDistance) {
      nearestScan = scan;
      smallestDistance = distance;
    }
  }

  return nearestScan;
}

function assertDatasetReady(state: ViewerState): asserts state is ViewerState & {
  dataset: { status: "ready"; dataset: NonNullable<ViewerState["dataset"]["dataset"]> };
} {
  if (state.dataset.status !== "ready" || !state.dataset.dataset) {
    throw new ViewerStateError(
      "DATASET_NOT_READY",
      "A loaded dataset is required for this viewer command"
    );
  }
}

function assertScanIndexExists(state: ViewerState, scanIndex: number): void {
  const scan = state.dataset.dataset?.scanSummaries[scanIndex];
  if (!scan) {
    throw new ViewerStateError(
      "SCAN_INDEX_OUT_OF_RANGE",
      `Scan index ${scanIndex} is outside the valid range`
    );
  }
}

function getPanel(state: ViewerState, panelId: PanelId): PanelState {
  return panelId === "tic" ? state.ticPanel : state.spectrumPanel;
}

function updatePanel(
  state: ViewerState,
  panelId: PanelId,
  panelState: PanelState
): ViewerState {
  return panelId === "tic"
    ? { ...state, ticPanel: panelState }
    : { ...state, spectrumPanel: panelState };
}

function normalizeRange(range: NumericRange): NumericRange {
  return range.min <= range.max
    ? range
    : {
        min: range.max,
        max: range.min
      };
}

function createDefaultPanelState(panelId: PanelId): PanelState {
  return panelId === "spectrum"
    ? { pinned: true, range: null }
    : { pinned: false, range: null };
}
