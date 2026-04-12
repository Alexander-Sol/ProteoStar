# Task 002: Build IMSP Parser

## Goal

Implement a validated `.imsp` parser in `packages/imsp-core`.

## Scope

1. Parse header, scan table, bin directory, and peak array accessors.
2. Validate magic bytes, version, and section boundaries.
3. Expose typed structures for scans, bins, and peak access.
4. Support parsing from `ArrayBuffer`.

## Deliverables

1. `parseImsp(buffer)` entry point
2. Stable domain types
3. Error handling for malformed files

## Tests

1. `tiny-known.imsp` parses correctly
2. Truncated input fails cleanly
3. Invalid magic fails cleanly

## Notes

Keep the parser framework-independent. UI code should not own binary parsing.
