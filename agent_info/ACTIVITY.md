# MsBrowser Activity

## Purpose

This file tracks active work, near-term priorities, decisions, and progress notes.

## Current Status

- Workspace scaffold completed
- `imsp-core` parser implemented for v1 format with fixture-backed tests
- Dataset query API implemented above the raw parser with fixture-backed tests
- Viewer-state reducer and Zustand store implemented with semantic commands and interaction tests
- Initial web UI implemented with local file loading and Plotly-backed TIC/spectrum panels
- `.imsp` format documented and fixtures created
- Git repository initialized with LFS for large files
- Root workspace, packages, and initial web app wired together

## Current Priorities

1. Scaffold the Bun workspace and package structure.
2. Implement and test the `.imsp` parser in `imsp-core`.
3. Build the dataset query API for TIC and per-scan spectrum retrieval.
4. Create viewer state for selection, pinning, and zoom behavior.
5. Build the initial web UI with Plotly-backed TIC and spectrum panels.

## Decisions Log

### Confirmed

- Canonical file extension: `.imsp`
- Initial hosting target: Vercel
- Runtime and package manager: Bun
- Web framework: Next.js
- State management: Zustand
- Plotting library for v1: Plotly
- Theme direction for v1: light theme only
- Initial data access mode: local file loading in browser

### Deferred

- Advanced cross-panel interactions
- Tauri packaging
- C# backend integration details
- Plotly replacement evaluation

## Progress Log

### 2026-04-11

- Reviewed `.imsp` format and QualBrowser screenshots.
- Defined high-level architecture and implementation roadmap.
- Confirmed need for worker-safe query boundaries and plot abstraction.
- Established requested planning structure under `agent_info/`.
- Began task `001` by scaffolding the Bun workspace, package manifests, shared TypeScript config, and initial Vitest coverage.
- Completed task `001` with a bootable Next.js app and shared package layout.
- Completed task `002` by implementing `parseImsp(buffer)` with header, scan table, bin directory, and peak access validation in `packages/imsp-core`.
- Added parser tests covering `tiny-known.imsp`, invalid magic, truncated input, and peak-index bounds.
- Completed task `003` by adding a dataset-provider interface and local `.imsp` implementation for metadata, scan summaries, nearest-scan lookup, spectrum reconstruction, m/z-range queries, and TIC trace queries.
- Added query API tests covering expected `tiny-known.imsp` results and local spectrum caching behavior.
- Completed task `004` by implementing viewer-state dataset loading, scan selection, independent TIC and spectrum pin state, semantic zoom/reset commands, and a Zustand-backed store in focused files.
- Added viewer-state tests covering initialization, nearest-scan selection, pin-gated zoom behavior, panel reset isolation, and store dispatch.
- Completed task `005` by wiring the Next.js app to local `.imsp` file loading, `imsp-core` dataset queries, `viewer-state` semantic commands, and Plotly-backed TIC and spectrum panels through the plot adapter.
- Added shared UI primitives for panel layout, file open controls, status banners, and readouts, and updated the plot adapter to emit semantic plot events while keeping Plotly client-only in the web app.

## Risks To Track

1. Spectrum-by-scan access may become a bottleneck if reconstruction is too expensive.
2. Plotly may become limiting for large interactive traces.
3. Backend and frontend may drift if dataset-query contracts are not formalized early.
4. Fixture coverage is good for now, but a true multi-megabyte medium file would improve performance validation.

## Next Checkpoint

The next meaningful checkpoint is hardening the UI workflow with integration tests, improved empty/loading/error handling, and refinement of spectrum and TIC interaction details.
