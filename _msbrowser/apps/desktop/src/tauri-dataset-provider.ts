// Transport adapter: fulfills the `DatasetProvider` contract by calling the Rust
// core over Tauri IPC. Numeric-array queries return Apache Arrow IPC streams
// (decoded with `tableFromIPC`); the rest are plain JSON. See IPC contract §5.

import { invoke, Channel } from "@tauri-apps/api/core";
import { tableFromIPC, type Table } from "apache-arrow";

import type {
  DatasetProvider,
  OpenResult,
  Precursor,
  ProgressEvent,
  ScanSummary,
  Spectrum,
  SpectrumPeak,
  TicPoint
} from "./contract";

const toTable = (buf: ArrayBuffer): Table => tableFromIPC(new Uint8Array(buf));

export async function openDataset(
  path: string,
  onProgress?: (p: ProgressEvent) => void
): Promise<{ handle: number; provider: DatasetProvider }> {
  const channel = new Channel<ProgressEvent>();
  if (onProgress) channel.onmessage = onProgress;

  const { handle, metadata } = await invoke<OpenResult>("open_dataset", {
    path,
    onProgress: channel
  });

  const arrow = (cmd: string, args: object): Promise<Table> =>
    invoke<ArrayBuffer>(cmd, { handle, ...args }).then(toTable);

  const provider: DatasetProvider = {
    getMetadata: async () => metadata,

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

    getMs2ForPrecursor: (mz, ms1ScanIndex) =>
      invoke("get_ms2_for_precursor", { handle, mz, ms1ScanIndex })
  };

  return { handle, provider };
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
