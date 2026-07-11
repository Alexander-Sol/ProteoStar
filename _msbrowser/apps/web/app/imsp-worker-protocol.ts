import type {
  DatasetMetadata,
  ImspPeak,
  ScanSummary,
  Spectrum,
  TicPoint
} from "@msbrowser/imsp-core";

export type ImspWorkerRequest =
  | {
      id: number;
      type: "initialize";
      buffer: ArrayBuffer;
    }
  | {
      id: number;
      type: "getMetadata";
    }
  | {
      id: number;
      type: "getScanSummaries";
    }
  | {
      id: number;
      type: "getNearestScan";
      retentionTime: number;
    }
  | {
      id: number;
      type: "getSpectrumForScan";
      scanIndex: number;
    }
  | {
      id: number;
      type: "getPeaksInMzRange";
      mzMin: number;
      mzMax: number;
    }
  | {
      id: number;
      type: "getTicTrace";
      rtRange?: { min: number; max: number };
    };

export type ImspWorkerResult =
  | { type: "initialize"; result: true }
  | { type: "getMetadata"; result: DatasetMetadata }
  | { type: "getScanSummaries"; result: readonly ScanSummary[] }
  | { type: "getNearestScan"; result: ScanSummary | null }
  | { type: "getSpectrumForScan"; result: Spectrum }
  | { type: "getPeaksInMzRange"; result: readonly ImspPeak[] }
  | { type: "getTicTrace"; result: readonly TicPoint[] };

export type ImspWorkerResponse =
  | {
      id: number;
      ok: true;
      type: ImspWorkerResult["type"];
      result: ImspWorkerResult["result"];
    }
  | {
      id: number;
      ok: false;
      error: {
        name: string;
        message: string;
      };
    };
