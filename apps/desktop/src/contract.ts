// MsViewer IPC contract — TypeScript surface (self-contained).
//
// Mirrors `MsViewer_IPC_Contract.md` §2 (JSON schemas) and §4 (DatasetProvider).
// The whole UI depends only on `DatasetProvider`; the Tauri adapter implements it
// and unit tests can supply a fake. This file is the single source of truth for
// the viewer's IPC types.

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
  /** Set to an error message if the background peak-index build failed; null otherwise. */
  indexError: string | null;
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

// Streamed on the `open_dataset` channel (contract §2 ProgressEvent) and on the
// `run_feature_detection` channel (detect/refine/resolve phases).
export interface ProgressEvent {
  phase: "reading" | "indexing" | "detecting" | "refining" | "resolving" | "done";
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
    | "PSM_READ"
    | "FEATURE_READ"
    | "FEATURE_PARSE"
    | "INTERNAL";
  message: string;
}

// ---------------------------------------------------------------- features
// Resolved feature-finding output (one row of the runner's resolved TSV),
// returned by the `load_features` command. Independent of any open dataset —
// features are loaded from a TSV file and overlaid on the raw data.

export interface PerChargeMz {
  charge: number;
  mz: number;
}

export interface Feature {
  detectedMz: number;
  rtStart: number;
  rtApex: number;
  rtEnd: number;
  chargeStates: number[];
  perChargeMz: PerChargeMz[];
  primaryCharge: number;
  monoisotopicMass: number;
  monoMz: number;
  summedIntensity: number;
  crossChargeSupport: number;
  numMembers: number;
}

// ---------------------------------------------------------------------- PSMs
// A MetaMorpheus PSM (one `.psmtsv` row), returned by the `load_psms` command.
// Independent of any open dataset; linked to a detected feature client-side.

export interface Psm {
  fullSequence: string;
  monoisotopicMass: number;
  precursorCharge: number;
  /** Theoretical precursor m/z (mass + charge) — the XIC fallback for an unlinked PSM. */
  precursorMz: number;
  /** MS2 retention time (min); -1 when the psmtsv omits it. */
  ms2RetentionTime: number;
  /** One-based MS2 scan number; -1 when the psmtsv omits it. Used to pull the identified spectrum. */
  ms2ScanNumber: number;
  qValue: number;
  score: number;
  fileName: string;
  isDecoy: boolean;
}

// The result of joining one PSM to the detected-feature list (mass + RT + charge).
export interface PsmLink {
  /** Index into the `features` array, or null when the PSM matched no feature. */
  featureIndex: number | null;
  /** |Δmass| in ppm against the linked feature, or null when unlinked. */
  massPpmError: number | null;
  /** MS2 RT − feature apex RT (min), or null when unlinked or RT is unavailable. */
  rtDelta: number | null;
}

// ---------------------------------------------------------- ladder walkthrough
// Diagnostics for the top-down charge-state-ladder fit of one manually chosen
// seed + anchoring charge, returned by the `score_seed_ladder` command.

export interface LadderTooth {
  isotopeIndex: number;
  /** Predicted comb m/z: monoMz + k·spacing. */
  expectedMz: number;
  /** Averagine weight (tallest tooth = 1.0) — the comb tooth's relative height. */
  weight: number;
  /** Observed peak m/z at the seed's apex scan within tolerance, else null. */
  observedMz: number | null;
  observedIntensity: number | null;
  /** True if the scorer credited this tooth to this charge (matched, not deduped away). */
  credited: boolean;
}

export interface LadderCharge {
  charge: number;
  monoMz: number;
  spacing: number;
  /** Cross-window matched-filter response over credited teeth. */
  response: number;
  numIsotopesObserved: number;
  /** True if this charge cleared the isotope floor and contributes to the mass response. */
  retained: boolean;
  teeth: LadderTooth[];
}

export interface SeedLadder {
  seedMz: number;
  seedScanIndex: number;
  seedRt: number;
  zSeed: number;
  /** Teeth found by the cheap apex spacing screen. */
  screenTeeth: number;
  screenPassed: boolean;
  /** Most-abundant averagine tooth index the seed is assumed to occupy at zSeed. */
  iStar: number;
  monoMz: number;
  monoMass: number;
  massInRange: boolean;
  /** Monoisotope after the joint cross-charge cosine offset search. */
  refinedMonoMass: number;
  totalResponse: number;
  numChargeStates: number;
  accepted: boolean;
  windowScanCount: number;
  maxCharge: number;
  charges: LadderCharge[];
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
  /** On-demand MS1 spectrum nearest a retention time; works before the peak index is built. */
  getSpectrumAtRt(
    retentionTime: number,
    opts?: { mzRange?: NumericRange; maxPeaks?: number }
  ): Promise<Spectrum>;
  /** On-demand MSn spectrum by one-based scan number (a PSM's `Scan Number`); no peak index needed. */
  getMs2Spectrum(
    scanNumber: number,
    opts?: { mzRange?: NumericRange; maxPeaks?: number }
  ): Promise<Spectrum>;
  getMs2ForPrecursor(mz: number, ms1ScanIndex: number): Promise<readonly ScanSummary[]>;
}
