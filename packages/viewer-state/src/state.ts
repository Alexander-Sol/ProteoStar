import type { DatasetSlotState, PanelState, ViewerState } from "./types";

const DEFAULT_TIC_PANEL_STATE: PanelState = {
  pinned: false,
  range: null
};

const DEFAULT_SPECTRUM_PANEL_STATE: PanelState = {
  pinned: true,
  range: null
};

function createIdleSlot(): DatasetSlotState {
  return {
    load: { status: "idle", dataset: null, errorMessage: null },
    selectedScanIndex: null
  };
}

export function createInitialViewerState(): ViewerState {
  return {
    activePanel: "tic",
    datasetSlots: [createIdleSlot(), createIdleSlot()],
    ticPanel: { ...DEFAULT_TIC_PANEL_STATE },
    spectrumPanel: { ...DEFAULT_SPECTRUM_PANEL_STATE }
  };
}
