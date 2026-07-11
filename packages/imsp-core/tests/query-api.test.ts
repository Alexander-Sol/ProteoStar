import { readFileSync } from "node:fs";
import { join } from "node:path";
import { describe, expect, it } from "vitest";

import {
  ImspParseError,
  createImspDatasetProvider
} from "../src/index";

function loadFixture(name: string): ArrayBuffer {
  const fixturePath = join(process.cwd(), "fixtures", name);
  const file = readFileSync(fixturePath);
  return file.buffer.slice(file.byteOffset, file.byteOffset + file.byteLength);
}

describe("IMSP dataset query API", () => {
  it("returns dataset metadata and scan summaries", async () => {
    const dataset = createImspDatasetProvider(loadFixture("tiny-known.imsp"));

    await expect(dataset.getMetadata()).resolves.toEqual({
      version: 1,
      binsPerDalton: 100,
      scanCount: 3,
      nonEmptyBinCount: 6,
      totalPeakCount: 9,
      retentionTimeRange: { min: 1, max: 3 },
      mzRange: { min: 500.1234, max: 1250.3456 }
    });

    await expect(dataset.getScanSummaries()).resolves.toEqual([
      { scanIndex: 0, oneBasedScanNumber: 1, retentionTime: 1, tic: 225000 },
      { scanIndex: 1, oneBasedScanNumber: 2, retentionTime: 2, tic: 250000 },
      { scanIndex: 2, oneBasedScanNumber: 3, retentionTime: 3, tic: 210000 }
    ]);
  });

  it("finds the nearest scan by retention time", async () => {
    const dataset = createImspDatasetProvider(loadFixture("tiny-known.imsp"));

    await expect(dataset.getNearestScan(0.2)).resolves.toEqual({
      scanIndex: 0,
      oneBasedScanNumber: 1,
      retentionTime: 1,
      tic: 225000
    });

    await expect(dataset.getNearestScan(1.51)).resolves.toEqual({
      scanIndex: 1,
      oneBasedScanNumber: 2,
      retentionTime: 2,
      tic: 250000
    });

    await expect(dataset.getNearestScan(2.6)).resolves.toEqual({
      scanIndex: 2,
      oneBasedScanNumber: 3,
      retentionTime: 3,
      tic: 210000
    });
  });

  it("reconstructs spectra per scan and reuses cached results", async () => {
    const dataset = createImspDatasetProvider(loadFixture("tiny-known.imsp"));

    const spectrum = await dataset.getSpectrumForScan(1);

    expect(spectrum).toEqual({
      scanIndex: 1,
      oneBasedScanNumber: 2,
      retentionTime: 2,
      peaks: [
        { mz: 500.1234, intensity: 60000 },
        { mz: 750.5678, intensity: 110000 },
        { mz: 1250.3456, intensity: 80000 }
      ]
    });

    await expect(dataset.getSpectrumForScan(1)).resolves.toBe(spectrum);
  });

  it("supports m/z range and TIC trace queries", async () => {
    const dataset = createImspDatasetProvider(loadFixture("tiny-known.imsp"));

    await expect(dataset.getPeaksInMzRange(750.56, 1000.91)).resolves.toEqual([
      { mz: 750.5678, intensity: 100000, scanIndex: 0 },
      { mz: 750.5678, intensity: 110000, scanIndex: 1 },
      { mz: 850.6789, intensity: 95000, scanIndex: 2 },
      { mz: 1000.9012, intensity: 75000, scanIndex: 0 },
      { mz: 1000.9012, intensity: 70000, scanIndex: 2 }
    ]);

    await expect(dataset.getTicTrace()).resolves.toEqual([
      { scanIndex: 0, retentionTime: 1, tic: 225000 },
      { scanIndex: 1, retentionTime: 2, tic: 250000 },
      { scanIndex: 2, retentionTime: 3, tic: 210000 }
    ]);

    await expect(dataset.getTicTrace({ min: 1, max: 2 })).resolves.toEqual([
      { scanIndex: 0, retentionTime: 1, tic: 225000 },
      { scanIndex: 1, retentionTime: 2, tic: 250000 }
    ]);
  });

  it("fails cleanly for an invalid scan index", async () => {
    const dataset = createImspDatasetProvider(loadFixture("tiny-known.imsp"));

    await expect(dataset.getSpectrumForScan(3)).rejects.toBeInstanceOf(ImspParseError);

    try {
      await dataset.getSpectrumForScan(3);
    } catch (error) {
      expect(error).toBeInstanceOf(ImspParseError);
      expect((error as ImspParseError).code).toBe("SCAN_INDEX_OUT_OF_RANGE");
    }
  });
});
