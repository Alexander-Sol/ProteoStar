# IMSP Peak Index File Format (`.ind`)

## Overview

An `.ind` file contains all MS1 peaks from a single mass-spec run, organised into
a sparse m/z bin index.  The binning strategy is identical to that used by the
FlashLFQ `PeakIndexingEngine`: every peak is placed in a bin whose index is

    BinIndex = round(mz × BinsPerDalton)

where `BinsPerDalton = 100`.  Peaks in bin 1000 therefore cover m/z values that
round to 10.00 Th (roughly 9.995–10.005 Th).

All numeric values are **little-endian**.

---

## File Layout

```
[ Section 1 ] File Header          –  24 bytes  (fixed)
[ Section 2 ] Scan Table           –  S × 12 bytes
[ Section 3 ] Bin Directory        –  N × 12 bytes
[ Section 4 ] Peak Array           –  T × 12 bytes
```

where:
- **S** = number of MS1 scans  (`ScanCount` in the header)
- **N** = number of non-empty bins  (`NonEmptyBinCount` in the header)
- **T** = total number of peaks  (`TotalPeakCount` in the header)

---

## Section 1 — File Header (offset 0, 24 bytes)

| Offset | Bytes | Type    | Field            | Description                          |
|-------:|------:|---------|------------------|--------------------------------------|
|      0 |     4 | ASCII   | Magic            | Always the four bytes `I M S P`      |
|      4 |     4 | uint32  | Version          | Format version; currently `1`        |
|      8 |     4 | uint32  | BinsPerDalton    | Bin resolution; currently `100`      |
|     12 |     4 | uint32  | NonEmptyBinCount | Number of non-empty bins (N)         |
|     16 |     4 | uint32  | TotalPeakCount   | Total peaks across all bins (T)      |
|     20 |     4 | uint32  | ScanCount        | Number of MS1 scans (S)              |

---

## Section 2 — Scan Table (offset 24, S × 16 bytes)

One entry per MS1 scan, ordered by ascending scan number.
The zero-based position of an entry is the `ScanIndex` referenced in the peak array.

| Relative offset | Bytes | Type    | Field               | Description                                              |
|----------------:|------:|---------|---------------------|----------------------------------------------------------|
|               0 |     4 | uint32   | OneBasedScanNumber  | Instrument scan number                                   |
|               4 |     8 | float64 | RetentionTime       | Retention time (minutes)                                 |
|              12 |     4 | float32 | Tic                 | Sum of filtered peak intensities; use for TIC chromatogram |

---

## Section 3 — Bin Directory (offset `24 + S×12`, N × 12 bytes)

One entry per non-empty bin, ordered by ascending `BinIndex`.
Use this section for fast m/z-range queries without scanning the full peak array.

| Relative offset | Bytes | Type    | Field       | Description                                          |
|----------------:|------:|---------|-------------|------------------------------------------------------|
|               0 |     4 | int32   | BinIndex    | `round(representativeMz × BinsPerDalton)`            |
|               4 |     4 | uint32  | PeakOffset  | 0-based index of the first peak for this bin in §4   |
|               8 |     4 | uint32  | PeakCount   | Number of peaks belonging to this bin                |

---

## Section 4 — Peak Array (offset `24 + S×16 + N×12`, T × 12 bytes)

Peaks are stored contiguously, bin-by-bin in the same order as the bin directory.
Each record is exactly **12 bytes**, so the byte offset of peak `i` within this
section is `i × 12`.

Retention time is not stored per-peak; look it up via `scans[ZeroBasedScanIndex].retentionTime`.

| Relative offset | Bytes | Type    | Field              | Notes                                                           |
|----------------:|------:|---------|--------------------|------------------------------------------------------------------|
|               0 |     4 | uint32  | MzTenThousandths   | `round(mz × 10000)` — recover with `/ 10000`; 4 d.p. max       |
|               4 |     4 | float32 | Intensity          | Raw detector intensity (always positive; covers up to ~3.4e38)  |
|               8 |     4 | uint32  | ZeroBasedScanIndex | Index into the Scan Table (Section 2)                           |

---

## TypeScript Reader (DataView)

```typescript
const HEADER_BYTES = 24;
const SCAN_BYTES   = 16;
const BIN_BYTES    = 12;
const PEAK_BYTES   = 12;

interface ImspHeader {
  binsPerDalton:    number;
  nonEmptyBinCount: number;
  totalPeakCount:   number;
  scanCount:        number;
}

interface ImspScan {
  oneBasedScanNumber: number;
  retentionTime:      number;
  tic:                number;
}

interface ImspBinEntry {
  binIndex:   number;
  peakOffset: number;
  peakCount:  number;
}

interface ImspPeak {
  mz:        number;  // Th, 4 d.p.
  intensity: number;
  scanIndex: number;
}

function parseImsp(buffer: ArrayBuffer) {
  const v = new DataView(buffer);

  // ── Header ───────────────────────────────────────────────────────────────────
  const magic = String.fromCharCode(
    v.getUint8(0), v.getUint8(1), v.getUint8(2), v.getUint8(3)
  );
  if (magic !== 'IMSP') throw new Error('Not an IMSP file');

  const header: ImspHeader = {
    binsPerDalton:    v.getUint32( 8, true),
    nonEmptyBinCount: v.getUint32(12, true),
    totalPeakCount:   v.getUint32(16, true),
    scanCount:        v.getUint32(20, true),
  };

  // ── Derived offsets ───────────────────────────────────────────────────────────
  const scanTableStart = HEADER_BYTES;
  const binDirStart    = scanTableStart + header.scanCount        * SCAN_BYTES;
  const peakArrayStart = binDirStart    + header.nonEmptyBinCount * BIN_BYTES;

  // ── Scan Table ────────────────────────────────────────────────────────────────
  const scans: ImspScan[] = [];
  for (let i = 0; i < header.scanCount; i++) {
    const off = scanTableStart + i * SCAN_BYTES;
    scans.push({
      oneBasedScanNumber: v.getInt32  (off,      true),
      retentionTime:      v.getFloat64(off +  4, true),
      tic:                v.getFloat32(off + 12, true),
    });
  }

  // ── Bin Directory ─────────────────────────────────────────────────────────────
  const bins: ImspBinEntry[] = [];
  for (let i = 0; i < header.nonEmptyBinCount; i++) {
    const off = binDirStart + i * BIN_BYTES;
    bins.push({
      binIndex:   v.getInt32 (off,     true),
      peakOffset: v.getUint32(off + 4, true),
      peakCount:  v.getUint32(off + 8, true),
    });
  }

  // ── Peak accessor ─────────────────────────────────────────────────────────────
  function readPeak(peakIndex: number): ImspPeak {
    const off = peakArrayStart + peakIndex * PEAK_BYTES;
    return {
      mz:        v.getUint32 (off,     true) / 10000,  // ten-thousandths Th → Th
      intensity: v.getFloat32(off + 4, true),
      scanIndex: v.getUint32 (off + 8, true),
    };
  }

  // ── Convenience: peaks for an m/z range ──────────────────────────────────────
  function peaksInMzRange(mzMin: number, mzMax: number): ImspPeak[] {
    const binMin = Math.floor(mzMin * header.binsPerDalton);
    const binMax = Math.ceil (mzMax * header.binsPerDalton);
    const result: ImspPeak[] = [];
    for (const bin of bins) {
      if (bin.binIndex < binMin) continue;
      if (bin.binIndex > binMax) break;
      for (let k = 0; k < bin.peakCount; k++)
        result.push(readPeak(bin.peakOffset + k));
    }
    return result;
  }

  return { header, scans, bins, readPeak, peaksInMzRange };
}
```

---

## Visualisation Notes

The natural 2-D view is a TIC (total ion chromatogram) plot of retention time (x-axis) vs. the summed intensity of all peaks in each scan (y-axis).  For a more detailed view, plot retention time vs. m/z, with point color representing intensity. This is referred to as an XIC, or extracted ion chromatogram, plot:

| Axis  | Source field      | Typical range           |
|-------|-------------------|-------------------------|
| X     | `RetentionTime`   | 0 – ~120 min            |
| Y     | `Mz`              | ~200 – ~2000 Th         |
| Color | `Intensity`       | log₁₀ scale recommended |

For large files, down-sample by binning peaks into a 2-D grid before rendering.
The bin directory lets you restrict reads to an m/z window without touching the
rest of the file — combine with a retention-time filter on the returned peaks for
a sub-region zoom.
