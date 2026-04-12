export const IMSP_MAGIC = "IMSP";
export const IMSP_VERSION = 1;
export const IMSP_BINS_PER_DALTON = 100;

const HEADER_BYTES = 24;
const SCAN_BYTES = 16;
const BIN_BYTES = 12;
const PEAK_BYTES = 12;

export interface ImspHeader {
  version: number;
  binsPerDalton: number;
  nonEmptyBinCount: number;
  totalPeakCount: number;
  scanCount: number;
}

export interface ImspScan {
  oneBasedScanNumber: number;
  retentionTime: number;
  tic: number;
}

export interface ImspBinEntry {
  binIndex: number;
  peakOffset: number;
  peakCount: number;
}

export interface ImspPeak {
  mz: number;
  intensity: number;
  scanIndex: number;
}

export type ImspParseErrorCode =
  | "BUFFER_TOO_SMALL"
  | "INVALID_MAGIC"
  | "UNSUPPORTED_VERSION"
  | "UNSUPPORTED_BINS_PER_DALTON"
  | "INVALID_SECTION_BOUNDS"
  | "INVALID_FILE_LENGTH"
  | "INVALID_SCAN_ORDER"
  | "INVALID_SCAN_VALUE"
  | "INVALID_BIN_ORDER"
  | "INVALID_BIN_RANGE"
  | "INVALID_PEAK_VALUE"
  | "INVALID_PEAK_SCAN_INDEX"
  | "PEAK_BIN_MISMATCH"
  | "PEAK_INDEX_OUT_OF_RANGE"
  | "BIN_NOT_FOUND";

export class ImspParseError extends Error {
  constructor(
    public readonly code: ImspParseErrorCode,
    message: string
  ) {
    super(message);
    this.name = "ImspParseError";
  }
}

export interface ImspFile {
  header: ImspHeader;
  scans: readonly ImspScan[];
  bins: readonly ImspBinEntry[];
  readPeak(peakIndex: number): ImspPeak;
  getPeaksForBin(binIndex: number): readonly ImspPeak[];
  peaksInMzRange(mzMin: number, mzMax: number): readonly ImspPeak[];
}

export function parseImsp(buffer: ArrayBuffer): ImspFile {
  if (buffer.byteLength < HEADER_BYTES) {
    throw new ImspParseError(
      "BUFFER_TOO_SMALL",
      `IMSP buffer is smaller than the ${HEADER_BYTES}-byte header`
    );
  }

  const view = new DataView(buffer);
  const magic = String.fromCharCode(
    view.getUint8(0),
    view.getUint8(1),
    view.getUint8(2),
    view.getUint8(3)
  );

  if (magic !== IMSP_MAGIC) {
    throw new ImspParseError(
      "INVALID_MAGIC",
      `Expected magic ${IMSP_MAGIC} but found ${magic}`
    );
  }

  const version = view.getUint32(4, true);
  if (version !== IMSP_VERSION) {
    throw new ImspParseError(
      "UNSUPPORTED_VERSION",
      `Unsupported IMSP version ${version}`
    );
  }

  const header: ImspHeader = {
    version,
    binsPerDalton: view.getUint32(8, true),
    nonEmptyBinCount: view.getUint32(12, true),
    totalPeakCount: view.getUint32(16, true),
    scanCount: view.getUint32(20, true)
  };

  if (header.binsPerDalton !== IMSP_BINS_PER_DALTON) {
    throw new ImspParseError(
      "UNSUPPORTED_BINS_PER_DALTON",
      `Unsupported bins-per-dalton value ${header.binsPerDalton}`
    );
  }

  const scanTableBytes = multiplyChecked(header.scanCount, SCAN_BYTES);
  const binDirectoryBytes = multiplyChecked(header.nonEmptyBinCount, BIN_BYTES);
  const peakArrayBytes = multiplyChecked(header.totalPeakCount, PEAK_BYTES);

  const scanTableStart = HEADER_BYTES;
  const binDirectoryStart = addChecked(scanTableStart, scanTableBytes);
  const peakArrayStart = addChecked(binDirectoryStart, binDirectoryBytes);
  const expectedByteLength = addChecked(peakArrayStart, peakArrayBytes);

  if (expectedByteLength < peakArrayStart || expectedByteLength > buffer.byteLength) {
    throw new ImspParseError(
      "INVALID_SECTION_BOUNDS",
      "Computed IMSP sections extend beyond the provided buffer"
    );
  }

  if (expectedByteLength !== buffer.byteLength) {
    throw new ImspParseError(
      "INVALID_FILE_LENGTH",
      `Expected IMSP file length ${expectedByteLength} bytes but received ${buffer.byteLength}`
    );
  }

  const scans: ImspScan[] = [];
  let previousScanNumber = 0;
  for (let scanIndex = 0; scanIndex < header.scanCount; scanIndex += 1) {
    const offset = scanTableStart + scanIndex * SCAN_BYTES;
    const scan: ImspScan = {
      oneBasedScanNumber: view.getUint32(offset, true),
      retentionTime: view.getFloat64(offset + 4, true),
      tic: view.getFloat32(offset + 12, true)
    };

    if (
      scan.oneBasedScanNumber < 1 ||
      !Number.isFinite(scan.retentionTime) ||
      scan.retentionTime < 0 ||
      !Number.isFinite(scan.tic) ||
      scan.tic < 0
    ) {
      throw new ImspParseError(
        "INVALID_SCAN_VALUE",
        `Invalid scan entry at index ${scanIndex}`
      );
    }

    if (scan.oneBasedScanNumber <= previousScanNumber) {
      throw new ImspParseError(
        "INVALID_SCAN_ORDER",
        "Scan table must be ordered by ascending scan number"
      );
    }

    previousScanNumber = scan.oneBasedScanNumber;
    scans.push(scan);
  }

  const bins: ImspBinEntry[] = [];
  let previousBinIndex = -1;
  let expectedPeakOffset = 0;
  for (let binOffsetIndex = 0; binOffsetIndex < header.nonEmptyBinCount; binOffsetIndex += 1) {
    const offset = binDirectoryStart + binOffsetIndex * BIN_BYTES;
    const bin: ImspBinEntry = {
      binIndex: view.getUint32(offset, true),
      peakOffset: view.getUint32(offset + 4, true),
      peakCount: view.getUint32(offset + 8, true)
    };

    if (bin.peakCount < 1) {
      throw new ImspParseError(
        "INVALID_BIN_RANGE",
        `Bin at index ${binOffsetIndex} must contain at least one peak`
      );
    }

    if (bin.binIndex <= previousBinIndex) {
      throw new ImspParseError(
        "INVALID_BIN_ORDER",
        "Bin directory must be ordered by ascending bin index"
      );
    }

    if (bin.peakOffset !== expectedPeakOffset) {
      throw new ImspParseError(
        "INVALID_BIN_RANGE",
        "Bin directory must describe a contiguous peak array"
      );
    }

    expectedPeakOffset = addChecked(bin.peakOffset, bin.peakCount);
    if (expectedPeakOffset > header.totalPeakCount) {
      throw new ImspParseError(
        "INVALID_BIN_RANGE",
        `Bin at index ${binOffsetIndex} exceeds the peak array bounds`
      );
    }

    previousBinIndex = bin.binIndex;
    bins.push(bin);
  }

  if (expectedPeakOffset !== header.totalPeakCount) {
    throw new ImspParseError(
      "INVALID_BIN_RANGE",
      "Bin directory does not cover the full peak array"
    );
  }

  function readPeak(peakIndex: number): ImspPeak {
    if (!Number.isInteger(peakIndex) || peakIndex < 0 || peakIndex >= header.totalPeakCount) {
      throw new ImspParseError(
        "PEAK_INDEX_OUT_OF_RANGE",
        `Peak index ${peakIndex} is outside the valid range`
      );
    }

    const offset = peakArrayStart + peakIndex * PEAK_BYTES;
    const peak: ImspPeak = {
      mz: view.getUint32(offset, true) / 10000,
      intensity: view.getFloat32(offset + 4, true),
      scanIndex: view.getUint32(offset + 8, true)
    };

    if (!Number.isFinite(peak.mz) || peak.mz <= 0 || !Number.isFinite(peak.intensity) || peak.intensity <= 0) {
      throw new ImspParseError(
        "INVALID_PEAK_VALUE",
        `Invalid peak entry at index ${peakIndex}`
      );
    }

    if (peak.scanIndex >= header.scanCount) {
      throw new ImspParseError(
        "INVALID_PEAK_SCAN_INDEX",
        `Peak at index ${peakIndex} references scan ${peak.scanIndex}`
      );
    }

    return peak;
  }

  for (const bin of bins) {
    for (let peakIndex = bin.peakOffset; peakIndex < bin.peakOffset + bin.peakCount; peakIndex += 1) {
      const peak = readPeak(peakIndex);
      if (Math.round(peak.mz * header.binsPerDalton) !== bin.binIndex) {
        throw new ImspParseError(
          "PEAK_BIN_MISMATCH",
          `Peak ${peakIndex} does not belong to bin ${bin.binIndex}`
        );
      }
    }
  }

  return {
    header,
    scans,
    bins,
    readPeak,
    getPeaksForBin(binIndex: number): readonly ImspPeak[] {
      const bin = bins.find((candidate) => candidate.binIndex === binIndex);
      if (!bin) {
        throw new ImspParseError("BIN_NOT_FOUND", `Bin ${binIndex} was not found`);
      }

      return collectPeaks(bin.peakOffset, bin.peakCount, readPeak);
    },
    peaksInMzRange(mzMin: number, mzMax: number): readonly ImspPeak[] {
      const lowerMz = Math.min(mzMin, mzMax);
      const upperMz = Math.max(mzMin, mzMax);
      const lowerBin = Math.floor(lowerMz * header.binsPerDalton);
      const upperBin = Math.ceil(upperMz * header.binsPerDalton);
      const peaks: ImspPeak[] = [];

      for (const bin of bins) {
        if (bin.binIndex < lowerBin) {
          continue;
        }

        if (bin.binIndex > upperBin) {
          break;
        }

        for (let peakIndex = bin.peakOffset; peakIndex < bin.peakOffset + bin.peakCount; peakIndex += 1) {
          const peak = readPeak(peakIndex);
          if (peak.mz >= lowerMz && peak.mz <= upperMz) {
            peaks.push(peak);
          }
        }
      }

      return peaks;
    }
  };
}

function collectPeaks(
  startIndex: number,
  peakCount: number,
  readPeak: (peakIndex: number) => ImspPeak
): ImspPeak[] {
  const peaks: ImspPeak[] = [];

  for (let index = 0; index < peakCount; index += 1) {
    peaks.push(readPeak(startIndex + index));
  }

  return peaks;
}

function addChecked(left: number, right: number): number {
  const result = left + right;

  if (!Number.isSafeInteger(result)) {
    throw new ImspParseError(
      "INVALID_SECTION_BOUNDS",
      "IMSP section offsets exceed the supported integer range"
    );
  }

  return result;
}

function multiplyChecked(left: number, right: number): number {
  const result = left * right;

  if (!Number.isSafeInteger(result)) {
    throw new ImspParseError(
      "INVALID_SECTION_BOUNDS",
      "IMSP section lengths exceed the supported integer range"
    );
  }

  return result;
}
