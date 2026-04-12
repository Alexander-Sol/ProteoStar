export type PanelId = "tic" | "spectrum";

export interface ViewerState {
  activePanel: PanelId;
  pinnedPanel: PanelId | null;
}

export function createInitialViewerState(): ViewerState {
  return {
    activePanel: "tic",
    pinnedPanel: null
  };
}
