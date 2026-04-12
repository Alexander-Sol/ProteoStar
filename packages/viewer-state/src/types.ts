export type PanelId = "tic" | "spectrum";

export interface NumericRange {
  min: number;
  max: number;
}

export interface ViewerDatasetMetadata {
  retentionTimeRange: NumericRange | null;
  mzRange: NumericRange | null;
  scanCount: number;
}

export interface ViewerScanSummary {
  scanIndex: number;
  oneBasedScanNumber: number;
  retentionTime: number;
  tic: number;
}

export interface ViewerDataset {
  metadata: ViewerDatasetMetadata;
  scanSummaries: readonly ViewerScanSummary[];
}

export type DatasetStatus = "idle" | "loading" | "ready" | "error";

export interface DatasetLoadState {
  status: DatasetStatus;
  dataset: ViewerDataset | null;
  errorMessage: string | null;
}

export interface PanelState {
  pinned: boolean;
  range: NumericRange | null;
}

export interface ViewerState {
  activePanel: PanelId;
  dataset: DatasetLoadState;
  selectedScanIndex: number | null;
  ticPanel: PanelState;
  spectrumPanel: PanelState;
}

export interface ViewerStoreState extends ViewerState {
  dispatch(command: ViewerCommand): void;
}

export interface ViewerStore {
  getState(): ViewerStoreState;
  setState(
    partial:
      | Partial<ViewerStoreState>
      | ((state: ViewerStoreState) => Partial<ViewerStoreState>),
    replace?: boolean
  ): void;
  subscribe(listener: () => void): () => void;
}

export type ViewerCommand =
  | { type: "dataset/load-started" }
  | { type: "dataset/load-succeeded"; dataset: ViewerDataset }
  | { type: "dataset/load-failed"; errorMessage: string }
  | { type: "selection/set-scan"; scanIndex: number | null }
  | { type: "selection/select-nearest-scan"; retentionTime: number }
  | { type: "panel/set-active"; panelId: PanelId }
  | { type: "panel/set-pinned"; panelId: PanelId; pinned: boolean }
  | { type: "panel/toggle-pinned"; panelId: PanelId }
  | { type: "panel/zoom"; panelId: PanelId; range: NumericRange }
  | { type: "panel/reset"; panelId: PanelId };
