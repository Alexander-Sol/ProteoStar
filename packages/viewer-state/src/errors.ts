export type ViewerStateErrorCode =
  | "DATASET_NOT_READY"
  | "SCAN_INDEX_OUT_OF_RANGE";

export class ViewerStateError extends Error {
  constructor(
    public readonly code: ViewerStateErrorCode,
    message: string
  ) {
    super(message);
    this.name = "ViewerStateError";
  }
}
