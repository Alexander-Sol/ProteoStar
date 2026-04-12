import type { PanelState, ViewerState } from "./types";

const DEFAULT_PANEL_STATE: PanelState = {
  pinned: false,
  range: null
};

export function createInitialViewerState(): ViewerState {
  return {
    activePanel: "tic",
    dataset: {
      status: "idle",
      dataset: null,
      errorMessage: null
    },
    selectedScanIndex: null,
    ticPanel: { ...DEFAULT_PANEL_STATE },
    spectrumPanel: { ...DEFAULT_PANEL_STATE }
  };
}
