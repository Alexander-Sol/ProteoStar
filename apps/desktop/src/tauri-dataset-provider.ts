// Transport adapter: fulfills the `DatasetProvider` contract by calling the Rust
// core over Tauri IPC. Numeric-array queries return Apache Arrow IPC streams
// (decoded with `tableFromIPC`); the rest are plain JSON. See IPC contract §5.

import { invoke, Channel } from "@tauri-apps/api/core";
import { tableFromIPC, type Table } from "apache-arrow";

import type {
  DatasetMetadata,
  DatasetProvider,
  NumericRange,
  OpenResult,
  Precursor,
  ProgressEvent,
  ScanSummary,
  Spectrum,
  SpectrumPeak,
  TicPoint
} from "./contract";
import type { SimSpec } from "./sim-spec";

const toTable = (buf: ArrayBuffer): Table => tableFromIPC(new Uint8Array(buf));

export async function openDataset(
  path: string,
  onProgress?: (p: ProgressEvent) => void
): Promise<{ handle: number; provider: DatasetProvider }> {
  const channel = new Channel<ProgressEvent>();
  if (onProgress) channel.onmessage = onProgress;

  const { handle } = await invoke<OpenResult>("open_dataset", {
    path,
    onProgress: channel
  });

  return { handle, provider: providerForHandle(handle) };
}

/**
 * Commit path (workbench interaction model): simulate a whole run in the backend and open it as a
 * dataset. The returned provider is the same `DatasetProvider` as a real run's, so the viewer
 * shows a simulated run in a slot with no other changes. One-off: the workbench uses
 * {@link SimSession} (live preview + commit on one long-lived simulator) instead.
 */
export async function openSimulatedDataset(
  spec: SimSpec
): Promise<{ handle: number; provider: DatasetProvider }> {
  const { handle } = await invoke<OpenResult>("open_simulated_dataset", { spec });
  return { handle, provider: providerForHandle(handle) };
}

/**
 * Live path (workbench interaction model): one fully simulated scan, decoded as a {@link Spectrum}.
 * Stateless — it builds a simulator from scratch per call — so for repeated edits prefer
 * {@link SimSession.scan}, which recomputes only what each edit changed.
 */
export async function simulateScan(
  spec: SimSpec,
  scanIndex: number,
  opts: { mzRange?: NumericRange; maxPeaks?: number } = {}
): Promise<Spectrum> {
  const buf = await invoke<ArrayBuffer>("simulate_scan", {
    spec,
    scanIndex,
    mzMin: opts.mzRange?.min,
    mzMax: opts.mzRange?.max,
    maxPeaks: opts.maxPeaks
  });
  return decodeSpectrum(toTable(buf));
}

/** Release a dataset handle's backend state (e.g. a superseded committed simulation). */
export function closeDataset(handle: number): Promise<void> {
  return invoke("close_dataset", { handle });
}

/** Load a `.simspec.json`; the backend checks its schema version and validates it. */
export function readSimSpec(path: string): Promise<SimSpec> {
  return invoke<SimSpec>("read_sim_spec", { path });
}

export function writeSimSpec(path: string, spec: SimSpec): Promise<void> {
  return invoke("write_sim_spec", { path, spec });
}

/** The demo build's bundled inputs to open at launch; `null` in a normal build. */
export interface DemoDataset {
  mzml: string;
  psms: string;
  simSpec: string;
}

export function demoDataset(): Promise<DemoDataset | null> {
  return invoke<DemoDataset | null>("demo_dataset");
}

type ScanOpts = { mzRange?: NumericRange; maxPeaks?: number };

interface PendingScan {
  spec: SimSpec;
  scanIndex: number;
  opts: ScanOpts;
  resolve: (s: Spectrum | null) => void;
  reject: (e: unknown) => void;
}

/**
 * A live simulator session: one long-lived backend `Simulator` that each edit updates in place,
 * recomputing only what the edit invalidated. The panel always passes its whole current
 * `SimSpec`; the backend diffs it against the last one.
 */
export class SimSession {
  private inFlight = false;
  private pending: PendingScan | null = null;

  private constructor(readonly id: number) {}

  static async open(spec: SimSpec): Promise<SimSession> {
    return new SimSession(await invoke<number>("open_sim_session", { spec }));
  }

  /**
   * Live path, latest-wins: at most one request is in flight. A call made while one is running
   * replaces any earlier waiting call, which resolves to `null` (superseded) — so a slider drag
   * renders the newest position as soon as the backend is free, with no debounce constant.
   */
  scan(spec: SimSpec, scanIndex: number, opts: ScanOpts = {}): Promise<Spectrum | null> {
    return new Promise((resolve, reject) => {
      this.pending?.resolve(null);
      this.pending = { spec, scanIndex, opts, resolve, reject };
      if (!this.inFlight) void this.drain();
    });
  }

  private async drain(): Promise<void> {
    this.inFlight = true;
    while (this.pending) {
      const { spec, scanIndex, opts, resolve, reject } = this.pending;
      this.pending = null;
      try {
        const buf = await invoke<ArrayBuffer>("sim_session_scan", {
          session: this.id,
          spec,
          scanIndex,
          mzMin: opts.mzRange?.min,
          mzMax: opts.mzRange?.max,
          maxPeaks: opts.maxPeaks
        });
        resolve(decodeSpectrum(toTable(buf)));
      } catch (e) {
        reject(e);
      }
    }
    this.inFlight = false;
  }

  /** Commit path (Accept): simulate the whole run on this session and open it as a dataset. */
  async commit(spec: SimSpec): Promise<{ handle: number; provider: DatasetProvider }> {
    const { handle } = await invoke<OpenResult>("commit_sim_session", { session: this.id, spec });
    return { handle, provider: providerForHandle(handle) };
  }

  close(): Promise<void> {
    this.pending?.resolve(null);
    this.pending = null;
    return invoke("close_sim_session", { session: this.id });
  }
}

/** A simulated run's TIC, per scan of its spec (sim-core `RunTic`). */
export interface RunTic {
  /** The proteoforms' summed centroids above the floor. */
  signal: number[];
  /** The noise floor's summed peaks; zeros without noise. */
  noise: number[];
  /** The selected proteoform's XIC, when one was asked for. */
  selected: { intensity: number[]; targetMzs: number[] } | null;
}

interface PendingTic {
  spec: SimSpec;
  selected: number | null;
  resolve: (t: RunTic | null) => void;
  reject: (e: unknown) => void;
}

/**
 * A live TIC session: the chromatogram view's counterpart of {@link SimSession}. One long-lived
 * backend `TicModel` returns the whole run's TIC without simulating its spectra, redoing only the
 * proteoforms each edit touched. Latest-wins, like {@link SimSession.scan}.
 */
export class TicSession {
  private inFlight = false;
  private pending: PendingTic | null = null;

  private constructor(readonly id: number) {}

  static async open(spec: SimSpec): Promise<TicSession> {
    return new TicSession(await invoke<number>("open_tic_session", { spec }));
  }

  trace(spec: SimSpec, selected: number | null): Promise<RunTic | null> {
    return new Promise((resolve, reject) => {
      this.pending?.resolve(null);
      this.pending = { spec, selected, resolve, reject };
      if (!this.inFlight) void this.drain();
    });
  }

  private async drain(): Promise<void> {
    this.inFlight = true;
    while (this.pending) {
      const { spec, selected, resolve, reject } = this.pending;
      this.pending = null;
      try {
        resolve(
          await invoke<RunTic>("tic_session_trace", {
            session: this.id,
            spec,
            selected: selected ?? undefined
          })
        );
      } catch (e) {
        reject(e);
      }
    }
    this.inFlight = false;
  }

  close(): Promise<void> {
    this.pending?.resolve(null);
    this.pending = null;
    return invoke("close_tic_session", { session: this.id });
  }
}

/**
 * Per MS1 scan of an open run, the summed tallest peak within `ppm` of each of `mzs`: the real-run
 * XIC at a simulated proteoform's centroids.
 */
export async function getTargetsXic(handle: number, mzs: readonly number[], ppm: number): Promise<TicPoint[]> {
  const buf = await invoke<ArrayBuffer>("get_targets_xic", { handle, mzs, ppm });
  return decodeTrace(toTable(buf));
}

/** Build a `DatasetProvider` over Tauri IPC for an already-open dataset `handle`. */
function providerForHandle(handle: number): DatasetProvider {
  const arrow = (cmd: string, args: object): Promise<Table> =>
    invoke<ArrayBuffer>(cmd, { handle, ...args }).then(toTable);

  const provider: DatasetProvider = {
    // Re-fetch from the backend (not the cached open-time `metadata`): the open result is
    // provisional (RT range only, ms1ScanCount 0) and gets refined when indexing finishes.
    getMetadata: () => invoke<DatasetMetadata>("get_metadata", { handle }),

    getScanSummaries: async () => decodeScanSummaries(await arrow("get_scan_summaries", {})),

    getNearestScan: (retentionTime, msLevel) =>
      invoke("get_nearest_scan", { handle, retentionTime, msLevel }),

    getTicTrace: async (o = {}) =>
      decodeTrace(
        await arrow("get_tic_trace", {
          rtMin: o.rtRange?.min,
          rtMax: o.rtRange?.max,
          msLevel: o.msLevel ?? 1,
          maxPoints: o.maxPoints
        })
      ),

    getRangeXic: async (mzLow, mzHigh, o = {}) =>
      decodeTrace(
        await arrow("get_range_xic", {
          mzLow,
          mzHigh,
          rtMin: o.rtRange?.min,
          rtMax: o.rtRange?.max,
          maxPoints: o.maxPoints
        })
      ),

    getSpectrum: async (scanIndex, o = {}) =>
      decodeSpectrum(
        await arrow("get_spectrum", {
          scanIndex,
          mzMin: o.mzRange?.min,
          mzMax: o.mzRange?.max,
          maxPeaks: o.maxPeaks
        })
      ),

    getSpectrumAtRt: async (retentionTime, o = {}) =>
      decodeSpectrum(
        await arrow("get_spectrum_at_rt", {
          retentionTime,
          mzMin: o.mzRange?.min,
          mzMax: o.mzRange?.max,
          maxPeaks: o.maxPeaks
        })
      ),

    getMs2Spectrum: async (scanNumber, o = {}) =>
      decodeSpectrum(
        await arrow("get_ms2_spectrum", {
          scanNumber,
          mzMin: o.mzRange?.min,
          mzMax: o.mzRange?.max,
          maxPeaks: o.maxPeaks
        })
      ),

    getMs2ForPrecursor: (mz, ms1ScanIndex) =>
      invoke("get_ms2_for_precursor", { handle, mz, ms1ScanIndex })
  };

  return provider;
}

// `trace` schema — shared by get_tic_trace and get_range_xic.
function decodeTrace(t: Table): TicPoint[] {
  const rt = t.getChild("retentionTime")!.toArray() as Float32Array;
  const int = t.getChild("intensity")!.toArray() as Float32Array;
  const si = t.getChild("scanIndex")!.toArray() as Uint32Array;
  const out: TicPoint[] = new Array(rt.length);
  for (let i = 0; i < rt.length; i++) {
    out[i] = { scanIndex: si[i], retentionTime: rt[i], intensity: int[i] };
  }
  return out;
}

function decodeSpectrum(t: Table): Spectrum {
  const mz = t.getChild("mz")!.toArray() as Float64Array;
  const int = t.getChild("intensity")!.toArray() as Float32Array;
  const peaks: SpectrumPeak[] = new Array(mz.length);
  for (let i = 0; i < mz.length; i++) peaks[i] = { mz: mz[i], intensity: int[i] };

  const m = t.schema.metadata; // Map<string,string>
  const precMz = m.get("precursorMz");
  return {
    scanIndex: Number(m.get("scanIndex")),
    oneBasedScanNumber: Number(m.get("oneBasedScanNumber")),
    retentionTime: Number(m.get("retentionTime")),
    msLevel: Number(m.get("msLevel")),
    precursor:
      precMz === undefined
        ? null
        : {
            mz: Number(precMz),
            scanIndex: Number(m.get("precursorScanIndex") ?? -1),
            isolationLow: Number(m.get("isolationLow") ?? NaN),
            isolationHigh: Number(m.get("isolationHigh") ?? NaN)
          },
    peaks
  };
}

// `scan_summaries` schema — 8 named columns (§3). Precursor fields are nullable;
// a null `precursorMz` marks an MS1 row (precursor: null).
function decodeScanSummaries(t: Table): ScanSummary[] {
  const oneBased = t.getChild("oneBasedScanNumber")!.toArray() as Uint32Array;
  const rt = t.getChild("retentionTime")!.toArray() as Float32Array;
  const tic = t.getChild("tic")!.toArray() as Float32Array;
  const msLevel = t.getChild("msLevel")!.toArray() as Uint8Array;
  const precursorMz = t.getChild("precursorMz")!;
  const precursorScanIndex = t.getChild("precursorScanIndex")!;
  const isolationLow = t.getChild("isolationLow")!;
  const isolationHigh = t.getChild("isolationHigh")!;

  const out: ScanSummary[] = new Array(oneBased.length);
  for (let i = 0; i < oneBased.length; i++) {
    const pMz = precursorMz.get(i);
    const precursor: Precursor | null =
      pMz === null || pMz === undefined
        ? null
        : {
            mz: Number(pMz),
            scanIndex: Number(precursorScanIndex.get(i) ?? -1),
            isolationLow: Number(isolationLow.get(i) ?? NaN),
            isolationHigh: Number(isolationHigh.get(i) ?? NaN)
          };
    out[i] = {
      scanIndex: i,
      oneBasedScanNumber: oneBased[i],
      retentionTime: rt[i],
      tic: tic[i],
      msLevel: msLevel[i],
      precursor
    };
  }
  return out;
}
