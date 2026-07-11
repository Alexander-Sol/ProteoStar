import {
  type DatasetProvider,
  type ScanSummary,
  type Spectrum,
  type TicPoint
} from "@msbrowser/imsp-core";
import type { PlotViewport, SpectrumPlotPeak, TicPlotPoint } from "@msbrowser/plot-adapter";
import type { NumericRange, ViewerDataset } from "@msbrowser/viewer-state";

import { createBrowserDatasetProvider } from "./worker-dataset-provider";

export interface LoadedDataset {
  fileName: string;
  provider: DatasetProvider;
  dispose?: () => void;
  metadata: Awaited<ReturnType<DatasetProvider["getMetadata"]>>;
  scanSummaries: readonly ScanSummary[];
  ticTrace: readonly TicPoint[];
}

export interface LoadedViewerDataset {
  loadedDataset: LoadedDataset;
  viewerDataset: ViewerDataset;
}

export async function loadViewerDataset(file: File): Promise<LoadedViewerDataset> {
  const buffer = await file.arrayBuffer();
  const { provider, dispose } = await createBrowserDatasetProvider(buffer);
  const [metadata, scanSummaries, ticTrace] = await Promise.all([
    provider.getMetadata(),
    provider.getScanSummaries(),
    provider.getTicTrace()
  ]);

  return {
    loadedDataset: {
      fileName: file.name,
      provider,
      dispose,
      metadata,
      scanSummaries,
      ticTrace
    },
    viewerDataset: {
      metadata: {
        retentionTimeRange: metadata.retentionTimeRange,
        mzRange: metadata.mzRange,
        scanCount: metadata.scanCount
      },
      scanSummaries: scanSummaries.map((scan) => ({
        scanIndex: scan.scanIndex,
        oneBasedScanNumber: scan.oneBasedScanNumber,
        retentionTime: scan.retentionTime,
        tic: scan.tic
      }))
    }
  };
}

export function toTicPlotPoints(ticTrace: readonly TicPoint[]): readonly TicPlotPoint[] {
  return ticTrace.map((point) => ({
    scanIndex: point.scanIndex,
    retentionTime: point.retentionTime,
    intensity: point.tic
  }));
}

export function toSpectrumPlotPeaks(
  spectrum: Spectrum | null
): readonly SpectrumPlotPeak[] {
  return (spectrum?.peaks ?? []).map((peak) => ({
    mz: peak.mz,
    intensity: peak.intensity
  }));
}

export function toViewport(range: NumericRange | null): PlotViewport {
  return {
    xMin: range?.min ?? null,
    xMax: range?.max ?? null
  };
}

export function formatNumber(
  value: number | undefined,
  decimals: number,
  suffix?: string
): string {
  if (value === undefined) {
    return "-";
  }

  return `${value.toFixed(decimals)}${suffix ? ` ${suffix}` : ""}`;
}

export function formatRange(
  range: NumericRange | null,
  decimals: number,
  suffix?: string
): string {
  if (!range) {
    return "Full range";
  }

  return `${formatNumber(range.min, decimals, suffix)} to ${formatNumber(range.max, decimals, suffix)}`;
}
