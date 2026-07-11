# MsBrowser Strategy

## Purpose

Build a TypeScript application for interactive visualization of mass spectrometry data stored in the custom `.imsp` format. The initial deployment target is a Vercel-hosted web app using Bun for development. The long-term target is a Tauri desktop application backed by a future C# service layer.

The first usable product should support the core QualBrowser-like workflow:

1. Open a local `.imsp` file.
2. Render a TIC in the top panel.
3. Select a point in the TIC to show the corresponding mass spectrum in the bottom panel.
4. Pin the TIC panel and drag-select to zoom the TIC.
5. Pin the spectrum panel and drag-select to zoom the spectrum.
6. Reset zoom state independently per panel.

## Product Principles

1. Prioritize scientific correctness over visual novelty.
2. Preserve the proven QualBrowser workflow while modernizing presentation.
3. Keep the app responsive on realistic files by isolating expensive work.
4. Design around stable query interfaces so the frontend can later switch from local parsing to a C# backend without major UI rewrites.
5. Avoid premature custom plotting. Use Plotly first behind an abstraction boundary.

## Architectural Vision

The application should be split into clear layers:

1. Platform shell
   - Web app now.
   - Tauri shell later.

2. Viewer application layer
   - Coordinates user interactions, panel state, selection state, and commands.
   - Has no knowledge of binary layout details.

3. Dataset query layer
   - Exposes dataset operations such as `getScanSummaries`, `getNearestScan`, `getSpectrumForScan`, and `getPeaksInMzRange`.
   - Can be implemented by a local TypeScript parser now and a C# backend later.

4. IMSP core layer
   - Parses binary `.imsp` files.
   - Validates headers and section boundaries.
   - Builds efficient query structures.

5. Plot adapter layer
   - Maps app-level plotting needs to Plotly.
   - Prevents Plotly-specific events and data structures from leaking into domain and viewer logic.

## Target Repository Shape

```text
apps/
  web/
packages/
  imsp-core/
  viewer-state/
  plot-adapter/
  ui/
  shared/
fixtures/
LargeFiles/
agent_info/
  STRATEGY.md
  ACTIVITY.md
  TASKS/
```

## Technical Choices

1. Runtime and package manager: Bun
2. Web framework: Next.js
3. UI framework: React
4. State management: Zustand
5. Plotting: Plotly via an adapter layer
6. Testing: Vitest
7. Future backend integration: C# service with dataset-query parity

## Data Strategy

The `.imsp` format is optimized for:

1. Fast TIC rendering from the scan table
2. Fast m/z-window lookups using the bin directory

The format is less directly optimized for repeated spectrum-by-scan rendering, so the frontend should build a scan-oriented query path on top of the raw format.

Recommended approach:

1. Parse the file into validated structures.
2. Expose raw accessors for scans, bins, and peaks.
3. Build higher-level query methods for viewer use.
4. Add caching for repeated spectrum reconstruction.
5. Move parsing and heavy queries into a worker before the UI becomes coupled to synchronous access.

## UI Strategy

The visual direction should be modern scientific software rather than a direct QualBrowser clone.

UI characteristics:

1. Light theme only for now.
2. Desktop-first layout.
3. Clean top app bar.
4. Two vertically stacked plot panels.
5. Compact panel headers with pin and reset actions.
6. Clear hover readouts and view state indicators.
7. Resizable panels after the base workflow is stable.

## Interaction Strategy

Define behavior in terms of semantic commands rather than chart-library callbacks.

Core v1 rules:

1. Clicking the TIC selects the nearest scan and updates the spectrum.
2. Drag-selecting while the TIC is pinned zooms the TIC RT range.
3. Drag-selecting while the spectrum is pinned zooms the spectrum m/z range.
4. Resetting a panel restores its default viewport.

Advanced cross-panel interactions are explicitly out of scope for v1 but should be enabled by the state model.

## Performance Strategy

Performance work should be staged rather than overbuilt.

Stage 1:

1. Parse fixtures and validate correctness.
2. Render TIC directly from scan summaries.
3. Reconstruct selected spectra on demand.

Stage 2:

1. Move parsing into a worker.
2. Cache reconstructed spectra.
3. Reduce unnecessary rerenders in the viewer state and plot adapters.

Stage 3:

1. Add large-file profiling.
2. Add downsampling or visible-range rendering where needed.
3. Re-evaluate Plotly if performance becomes limiting.

## Testing Strategy

Tests should be organized around correctness, interaction behavior, and scalability.

1. Parser tests
   - magic bytes
   - version handling
   - offset validation
   - entry decoding

2. Query tests
   - TIC generation
   - nearest scan lookup
   - spectrum reconstruction
   - m/z range extraction

3. Viewer-state tests
   - scan selection behavior
   - pin and zoom behavior
   - reset behavior

4. Integration tests
   - local file load to TIC render to spectrum update

5. Performance smoke tests
   - medium and large fixtures

## Backend Compatibility Strategy

The future C# backend should implement the same conceptual dataset operations used by the frontend. The frontend should never depend on implementation details that only make sense for local `DataView` parsing.

This means the primary frontend contract should eventually look like a dataset provider interface, with local file parsing as only one provider.

## Delivery Plan

1. Foundation
   - scaffold workspace and packages
   - add tooling

2. IMSP core
   - parser
   - validations
   - query API
   - tests

3. Viewer state
   - panel state
   - pin/zoom rules
   - semantic commands

4. Web UI
   - file open flow
   - TIC and spectrum panels
   - Plotly adapter

5. Worker and performance pass
   - async parsing and query boundary
   - caching

6. Hardening
   - polish
   - integration tests
   - deployability review

## Success Criteria For v1

1. A user can open a local `.imsp` file in the browser.
2. The TIC renders correctly from scan data.
3. Selecting the TIC updates the spectrum reliably.
4. Pinned zoom behavior works for both panels.
5. The app remains responsive on at least the current small and medium fixtures.
6. The codebase has a clean boundary between parsing, state, plotting, and UI.
