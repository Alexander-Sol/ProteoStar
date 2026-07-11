# Task 003: Build Dataset Query API

## Goal

Expose viewer-friendly dataset operations above the raw parser.

## Scope

1. Implement `getMetadata()`.
2. Implement `getScanSummaries()`.
3. Implement `getNearestScan(rt)`.
4. Implement `getSpectrumForScan(scanIndex)`.
5. Implement `getPeaksInMzRange(mzMin, mzMax)`.
6. Implement `getTicTrace(rtRange?)`.

## Deliverables

1. Dataset provider interface
2. Local `.imsp` implementation of that interface
3. Basic spectrum caching strategy if needed

## Notes

Design this as the future contract surface for the C# backend.
