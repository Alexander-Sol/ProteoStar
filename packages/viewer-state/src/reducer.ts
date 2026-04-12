import { ViewerStateError } from "./errors";
import type {
  DatasetSlotState,
  NumericRange,
  PanelId,
  PanelState,
  SlotIndex,
  ViewerCommand,
  ViewerScanSummary,
  ViewerState
} from "./types";

export function reduceViewerState(
  state: ViewerState,
  command: ViewerCommand
): ViewerState {
  switch (command.type) {
    case "dataset/load-started": {
      const { slotIndex } = command;
      const otherIndex: SlotIndex = slotIndex === 0 ? 1 : 0;
      const shouldResetPanels =
        slotIndex === 0 && state.datasetSlots[otherIndex].load.status === "idle";
      const newSlot: DatasetSlotState = {
        load: { status: "loading", dataset: null, errorMessage: null },
        selectedScanIndex: null
      };
      let next = updateSlot(state, slotIndex, newSlot);
      if (shouldResetPanels) {
        next = {
          ...next,
          ticPanel: createDefaultPanelState("tic"),
          spectrumPanel: createDefaultPanelState("spectrum")
        };
      }
      return next;
    }

    case "dataset/load-succeeded": {
      const { slotIndex, dataset } = command;
      return updateSlot(state, slotIndex, {
        load: { status: "ready", dataset, errorMessage: null },
        selectedScanIndex: null
      });
    }

    case "dataset/load-failed": {
      const { slotIndex, errorMessage } = command;
      return updateSlot(state, slotIndex, {
        load: { status: "error", dataset: null, errorMessage },
        selectedScanIndex: null
      });
    }

    case "selection/set-scan": {
      const { slotIndex } = command;
      if (command.scanIndex === null) {
        return updateSlot(state, slotIndex, {
          ...getSlot(state, slotIndex),
          selectedScanIndex: null
        });
      }

      const slot = getSlot(state, slotIndex);
      assertSlotDatasetReady(slot);
      assertScanIndexExists(slot, command.scanIndex);

      return updateSlot(state, slotIndex, {
        ...slot,
        selectedScanIndex: command.scanIndex
      });
    }

    case "selection/select-nearest-scan": {
      const { slotIndex } = command;
      const slot = getSlot(state, slotIndex);
      assertSlotDatasetReady(slot);
      const nearestScan = findNearestScan(
        slot.load.dataset.scanSummaries,
        command.retentionTime
      );

      return updateSlot(state, slotIndex, {
        ...slot,
        selectedScanIndex: nearestScan?.scanIndex ?? null
      });
    }

    case "panel/set-active":
      return {
        ...state,
        activePanel: command.panelId
      };

    case "panel/set-pinned": {
      const nextPinned = command.pinned;
      let next = updatePanel(state, command.panelId, {
        ...getPanel(state, command.panelId),
        pinned: nextPinned
      });
      if (nextPinned) {
        const otherId = getOtherPanelId(command.panelId);
        next = updatePanel(next, otherId, { ...getPanel(next, otherId), pinned: false });
      }
      return next;
    }

    case "panel/toggle-pinned": {
      const panel = getPanel(state, command.panelId);
      const nextPinned = !panel.pinned;
      let next = updatePanel(state, command.panelId, { ...panel, pinned: nextPinned });
      if (nextPinned) {
        const otherId = getOtherPanelId(command.panelId);
        next = updatePanel(next, otherId, { ...getPanel(next, otherId), pinned: false });
      }
      return next;
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

// --- Slot helpers ---

function getSlot(state: ViewerState, slotIndex: SlotIndex): DatasetSlotState {
  return state.datasetSlots[slotIndex];
}

function updateSlot(
  state: ViewerState,
  slotIndex: SlotIndex,
  slot: DatasetSlotState
): ViewerState {
  const slots: [DatasetSlotState, DatasetSlotState] = [
    state.datasetSlots[0],
    state.datasetSlots[1]
  ];
  slots[slotIndex] = slot;
  return { ...state, datasetSlots: slots };
}

// --- Panel helpers ---

function getPanel(state: ViewerState, panelId: PanelId): PanelState {
  return panelId === "tic" ? state.ticPanel : state.spectrumPanel;
}

function getOtherPanelId(panelId: PanelId): PanelId {
  return panelId === "tic" ? "spectrum" : "tic";
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

// --- Scan helpers ---

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

function assertSlotDatasetReady(
  slot: DatasetSlotState
): asserts slot is DatasetSlotState & {
  load: { status: "ready"; dataset: NonNullable<DatasetSlotState["load"]["dataset"]> };
} {
  if (slot.load.status !== "ready" || !slot.load.dataset) {
    throw new ViewerStateError(
      "DATASET_NOT_READY",
      "A loaded dataset is required for this viewer command"
    );
  }
}

function assertScanIndexExists(slot: DatasetSlotState, scanIndex: number): void {
  const scan = slot.load.dataset?.scanSummaries[scanIndex];
  if (!scan) {
    throw new ViewerStateError(
      "SCAN_INDEX_OUT_OF_RANGE",
      `Scan index ${scanIndex} is outside the valid range`
    );
  }
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
