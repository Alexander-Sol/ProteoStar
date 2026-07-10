// MsViewer IPC contract — TypeScript surface (self-contained).
//
// Mirrors `MsViewer_IPC_Contract.md` §2 (JSON schemas) and §4 (DatasetProvider).
// The whole UI depends only on `DatasetProvider`; the Tauri adapter implements it
// and unit tests can supply a fake. Do NOT import these from imsp-core — this file
// is the single source of truth for the new viewer.

export interface NumericRange {
  min: number;
  max: number;
}

export interface DatasetMetadata {
  fileName: string;
  format: "mzml" | "thermo_raw";
  scanCount: number;
  ms1ScanCount: number;
  msLevelsPresent: number[];
  retentionTimeRange: NumericRange | null;
  mzRange: NumericRange | null;
}

export interface Precursor {
  mz: number;
  scanIndex: number; // -1 if unknown
  isolationLow: number;
  isolationHigh: number;
}

export interface ScanSummary {
  scanIndex: number;
  oneBasedScanNumber: number;
  retentionTime: number;
  tic: number;
  msLevel: number;
  precursor: Precursor | null;
}

export interface TicPoint {
  scanIndex: number;
  retentionTime: number;
  intensity: number;
}

export interface XicPoint {
  scanIndex: number;
  retentionTime: number;
  intensity: number;
}

export interface SpectrumPeak {
  mz: number;
  intensity: number;
}

export interface Spectrum {
  scanIndex: number;
  oneBasedScanNumber: number;
  retentionTime: number;
  msLevel: number;
  precursor: Precursor | null;
  peaks: readonly SpectrumPeak[];
}

// Streamed on the `open_dataset` channel (contract §2 ProgressEvent).
export interface ProgressEvent {
  phase: "reading" | "indexing" | "done";
  scansDone: number;
  scansTotal: number;
}

// JSON payload returned by `open_dataset` (contract §2 OpenResult).
export interface OpenResult {
  handle: number;
  metadata: DatasetMetadata;
}

// Rejected commands (contract §2 ViewerError).
export interface ViewerError {
  code:
    | "FILE_NOT_FOUND"
    | "UNSUPPORTED_FORMAT"
    | "READ_ERROR"
    | "INDEX_BUILD_FAILED"
    | "HANDLE_NOT_FOUND"
    | "SCAN_OUT_OF_RANGE"
    | "EMPTY_INDEX"
    | "THERMO_RUNTIME_MISSING"
    | "INTERNAL";
  message: string;
}

export interface DatasetProvider {
  getMetadata(): Promise<DatasetMetadata>;
  getScanSummaries(): Promise<readonly ScanSummary[]>;
  getNearestScan(retentionTime: number, msLevel?: number): Promise<ScanSummary | null>;
  getTicTrace(opts?: {
    rtRange?: NumericRange;
    msLevel?: number;
    maxPoints?: number;
  }): Promise<readonly TicPoint[]>;
  getRangeXic(
    mzLow: number,
    mzHigh: number,
    opts?: { rtRange?: NumericRange; maxPoints?: number }
  ): Promise<readonly XicPoint[]>;
  getSpectrum(
    scanIndex: number,
    opts?: { mzRange?: NumericRange; maxPeaks?: number }
  ): Promise<Spectrum>;
  getMs2ForPrecursor(mz: number, ms1ScanIndex: number): Promise<readonly ScanSummary[]>;
}
