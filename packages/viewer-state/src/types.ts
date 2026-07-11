export type PanelId = "tic" | "spectrum";
export type SlotIndex = 0 | 1;

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

export interface DatasetSlotState {
  load: DatasetLoadState;
  selectedScanIndex: number | null;
}

export interface PanelState {
  pinned: boolean;
  range: NumericRange | null;
}

export interface ViewerState {
  activePanel: PanelId;
  datasetSlots: [DatasetSlotState, DatasetSlotState];
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
  | { type: "dataset/load-started"; slotIndex: SlotIndex }
  | { type: "dataset/load-succeeded"; slotIndex: SlotIndex; dataset: ViewerDataset }
  | { type: "dataset/load-failed"; slotIndex: SlotIndex; errorMessage: string }
  | { type: "selection/set-scan"; slotIndex: SlotIndex; scanIndex: number | null }
  | { type: "selection/select-nearest-scan"; slotIndex: SlotIndex; retentionTime: number }
  | { type: "panel/set-active"; panelId: PanelId }
  | { type: "panel/set-pinned"; panelId: PanelId; pinned: boolean }
  | { type: "panel/toggle-pinned"; panelId: PanelId }
  | { type: "panel/zoom"; panelId: PanelId; range: NumericRange }
  | { type: "panel/reset"; panelId: PanelId };
