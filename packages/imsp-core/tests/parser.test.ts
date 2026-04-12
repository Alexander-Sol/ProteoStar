import { readFileSync } from "node:fs";
import { join } from "node:path";
import { describe, expect, it } from "vitest";

import {
  IMSP_BINS_PER_DALTON,
  IMSP_MAGIC,
  IMSP_VERSION,
  ImspParseError,
  parseImsp
} from "../src/index";

function loadFixture(name: string): ArrayBuffer {
  const fixturePath = join(process.cwd(), "fixtures", name);
  const file = readFileSync(fixturePath);
  return file.buffer.slice(file.byteOffset, file.byteOffset + file.byteLength);
}

describe("parseImsp", () => {
  it("parses tiny-known.imsp correctly", () => {
    const imsp = parseImsp(loadFixture("tiny-known.imsp"));

    expect(imsp.header).toEqual({
      version: IMSP_VERSION,
      binsPerDalton: IMSP_BINS_PER_DALTON,
      nonEmptyBinCount: 6,
      totalPeakCount: 9,
      scanCount: 3
    });

    expect(imsp.scans).toEqual([
      { oneBasedScanNumber: 1, retentionTime: 1, tic: 225000 },
      { oneBasedScanNumber: 2, retentionTime: 2, tic: 250000 },
      { oneBasedScanNumber: 3, retentionTime: 3, tic: 210000 }
    ]);

    expect(imsp.bins).toEqual([
      { binIndex: 50012, peakOffset: 0, peakCount: 2 },
      { binIndex: 60023, peakOffset: 2, peakCount: 1 },
      { binIndex: 75057, peakOffset: 3, peakCount: 2 },
      { binIndex: 85068, peakOffset: 5, peakCount: 1 },
      { binIndex: 100090, peakOffset: 6, peakCount: 2 },
      { binIndex: 125035, peakOffset: 8, peakCount: 1 }
    ]);

    expect(imsp.readPeak(3)).toEqual({
      mz: 750.5678,
      intensity: 100000,
      scanIndex: 0
    });

    expect(imsp.getPeaksForBin(50012)).toEqual([
      { mz: 500.1234, intensity: 50000, scanIndex: 0 },
      { mz: 500.1234, intensity: 60000, scanIndex: 1 }
    ]);

    expect(imsp.peaksInMzRange(750.56, 1000.91)).toEqual([
      { mz: 750.5678, intensity: 100000, scanIndex: 0 },
      { mz: 750.5678, intensity: 110000, scanIndex: 1 },
      { mz: 850.6789, intensity: 95000, scanIndex: 2 },
      { mz: 1000.9012, intensity: 75000, scanIndex: 0 },
      { mz: 1000.9012, intensity: 70000, scanIndex: 2 }
    ]);
  });

  it("fails cleanly for truncated input", () => {
    const bytes = new Uint8Array(loadFixture("tiny-known.imsp")).slice(0, 251);

    expect(() => parseImsp(bytes.buffer)).toThrowError(
      new ImspParseError(
        "INVALID_SECTION_BOUNDS",
        "Computed IMSP sections extend beyond the provided buffer"
      )
    );
  });

  it("fails cleanly for invalid magic", () => {
    const bytes = new Uint8Array(loadFixture("tiny-known.imsp"));
    bytes.set([0x42, 0x41, 0x44, 0x21], 0);

    expect(() => parseImsp(bytes.buffer)).toThrowError(ImspParseError);

    try {
      parseImsp(bytes.buffer);
    } catch (error) {
      expect(error).toBeInstanceOf(ImspParseError);
      expect((error as ImspParseError).code).toBe("INVALID_MAGIC");
      expect((error as ImspParseError).message).toContain(IMSP_MAGIC);
    }
  });

  it("fails when a peak index is out of range", () => {
    const imsp = parseImsp(loadFixture("tiny-known.imsp"));

    expect(() => imsp.readPeak(9)).toThrowError(ImspParseError);
  });
});
