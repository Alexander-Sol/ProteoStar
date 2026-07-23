# Spectrum zoom persistence — two overlapping mechanisms

**Status:** tech-debt note / design caveat. Written 2026-07-22 after the
`worktree-playwright-e2e` → `main` merge (commit `b8db911`).

**Update 2026-07-23 (branch `fix-raw-tic-flat`):** the fragility predicted below actually bit.
After commit `282c876` wired `onRelayout` → `setSpectrumXView` (zoom-reactive highlights), an
interactive box-zoom triggered a React re-render on the zoom itself — and Plotly's `uirevision` did
**not** reliably preserve the manual x-zoom across that re-render (or across the subsequent
new-spectrum data swap). Symptoms: the zoom held once, was lost on the next spectrum selection, then
stopped working entirely. **Fixed by implementing Option A for the x-axis** (below): the
`xrange-change` handler in `App.tsx` now writes the interactive x-range into `spectrumViewport`
(without bumping `uiRev`), so the `range` prop — not `uirevision` — is the source of truth that
re-applies the zoom deterministically across re-renders. A double-click autorange reset
(`xrange-change` with `range === null`) calls `reframeSpectrum(createDefaultViewport())`. The
interactive-zoom coverage gap flagged in §"Why it's fragile" (3) is now closed by
`apps/desktop/e2e/tauri/spectrum-interactive-zoom.spec.ts`. The y-axis still auto-fits to the
visible x-window (consistent with the programmatic-zoom paths); only x is captured.

## Summary

The mass-spectrum panel keeps its zoom "fixed" across scan changes through **two independent
mechanisms that were built on separate branches and now coexist**:

1. **Plotly `uirevision`** (from the `visualize-results` line) — preserves the user's *interactive*
   zoom/pan (mouse box-zoom via the modebar, drag-pan) across re-renders.
2. **A React `viewport` + y-range selector** (from the feature-detection/walkthrough line) —
   `PlotViewport.{xMin,xMax,yMin,yMax}` drives the Plotly `range` prop, and freezes the auto-fit
   y-range so an isotope envelope doesn't shrink while arrow-stepping scans.

They compose correctly today, but there is **no single source of truth for "the current zoom,"** the
interaction is subtle, and the two paths have different test coverage. This note documents how they
fit together, why it's fragile, and options for consolidation.

## The two mechanisms

### 1. `uirevision` (interactive zoom)

- `App.tsx` holds `spectrumUiRev` (a counter) and `reframeSpectrum(vp)`, which sets
  `spectrumViewport` **and** bumps `spectrumUiRev`.
- `spectrumUiRev` is passed to `SpectrumPlot` as `uirevision`, and `plots.tsx` sets
  `layout.uirevision = props.uirevision`.
- Plotly semantics: **while `uirevision` is unchanged**, the user's interactive zoom/pan is preserved
  and the supplied `range` prop is ignored for any axis the user has touched. **When `uirevision`
  changes** (a reframe), Plotly discards the interactive state and re-applies the figure as given.
- So a *manual* zoom lives **only inside Plotly** — it is never written back into `spectrumViewport`.
- Reframe callers (each bumps `uiRev`): feature selection, PSM selection, walkthrough ladder seed
  (`computeLadder`), charge-row focus (`focusChargeView`), and the spectrum panel's "Reset zoom" button.

### 2. `viewport` + `resolveSpectrumYRange` (programmatic zoom + y-freeze)

- `PlotViewport` carries optional `yMin/yMax` (`packages/plot-adapter/src/types.ts`).
- `packages/plot-adapter/src/viewport.ts`: `resolveSpectrumYRange(viewport, peaks)` returns the
  persisted y-range when both bounds are set (`hasPersistedY`), else the auto-fit to the visible
  x-window (`fitSpectrumYRange`). `plots.tsx` feeds this into `yaxis.range`; `xaxis.range` comes from
  `toPlotlyRange(viewport)`.
- `App.tsx` mutates the viewport **without** bumping `uiRev` on the "keep the view" paths:
  - `stepScan` freezes y (`setSpectrumViewport(v => ({...v, yMin, yMax}))`) on the first arrow step
    while x-zoomed, so the envelope height stays stable across scans.
  - `handleAreaClick` clears the frozen y (`yMin/yMax → null`) so a fresh TIC pick auto-fits.

## How they compose today

| User action | `uiRev` | `viewport` | Result |
| --- | --- | --- | --- |
| Manual box-zoom (modebar), then arrow-step scans | unchanged | unchanged (null) | `uirevision` holds the manual zoom (x+y) |
| Feature / PSM / ladder select | bumped | set to the target window | `range` prop applies the programmatic zoom |
| Arrow-step while x-zoomed via a reframe | unchanged | x kept; y frozen by `stepScan` | `range` prop holds x; frozen y holds envelope height |
| Fresh TIC click (`handleAreaClick`) | unchanged | y cleared | new scan auto-fits y; x kept |
| "Reset zoom" | bumped (`createDefaultViewport()`) | all null | `uiRev` bump + `range: undefined` → Plotly autoranges |

The rule that makes it work: **every programmatic reframe passes a fresh viewport (no stale y) and
bumps `uiRev`; every "keep the view" path leaves `uiRev` alone.**

## Why it's fragile

1. **No single source of truth.** Interactive zoom lives in Plotly; programmatic zoom lives in
   `viewport`. Reasoning about "what is the current zoom?" means reasoning about both.

2. **Rendering attributes silently break one path but not the other.** Concrete incident: the merge
   added an explicit `autorange` to the spectrum axes to make "Reset zoom" re-fit. That is applied on
   *every* render, which overrode `uirevision` and made a manual zoom snap back on each scan step —
   while every programmatic path still worked. It was reverted in `b8db911`. The lesson: any change to
   `xaxis`/`yaxis` attributes can fight `uirevision` for the interactive path.

3. **The interactive path has no automated coverage.** The e2e (`Task 4` in
   `apps/desktop/e2e/tauri/msviewer-tasks.spec.ts`) creates its x-zoom via the walkthrough ladder — a
   **programmatic** reframe — and asserts the y-range holds across a step and "Reset zoom" widens. It
   never performs an actual interactive zoom, so it stayed green even with the `autorange` regression
   present. The interactive-zoom-persistence behaviour reached a human before any test noticed.

4. **Partial redundancy.** `stepScan`'s y-freeze and `uirevision` overlap: if the user has manually
   zoomed (x+y), `uirevision` already holds y across steps, so the freeze is redundant there. The
   freeze earns its keep only for the *programmatic* x-zoom case (no interactive y-state), where the
   `range` prop is what applies. This overlap is a maintenance hazard.

## Recommendations

Short term:
- Add an e2e that performs a real interactive zoom (drive a Plotly relayout / modebar box-zoom on the
  spectrum graph div), steps a scan, and asserts the range is preserved — closing the coverage gap in
  §"Why it's fragile" (3). Without it, the `uirevision` path is untested.
- Treat any edit to the spectrum `xaxis`/`yaxis` layout as able to break interactive persistence; test
  both paths after such edits.

Longer term — pick one source of truth:
- **Option A: capture interactive zoom into `viewport`.** Wire Plotly `onRelayout` to write the user's
  x/y range back into `spectrumViewport`, then drop `uirevision` and let the `range` prop be the sole
  driver. One state, one code path; `stepScan`'s y-freeze becomes unnecessary (the persisted range
  covers it). This is the cleaner end state but a larger change.
- **Option B: lean fully on `uirevision`.** Let Plotly own zoom persistence and reduce `viewport` to
  the intentional-reframe targets only. Risk: the auto-fit-y-stability that `resolveSpectrumYRange`
  provides for the programmatic x-zoom case would need re-expressing.

Option A is preferred: it removes the dual-source-of-truth entirely.

## Code references

- `apps/desktop/src/App.tsx` — `spectrumViewport`, `spectrumUiRev`, `reframeSpectrum`, `stepScan`
  (y-freeze), `handleAreaClick` (y-clear), reframe callers (`selectFeature`, PSM select,
  `computeLadder`, `focusChargeView`, "Reset zoom" button).
- `packages/plot-adapter/src/viewport.ts` — `resolveSpectrumYRange`, `hasPersistedY`,
  `fitSpectrumYRange`, `createDefaultViewport`.
- `packages/plot-adapter/src/plots.tsx` — `SpectrumPlot` layout: `xaxis.range`/`yaxis.range` and
  `layout.uirevision`.
- `packages/plot-adapter/src/types.ts` — `PlotViewport` (`yMin`/`yMax`), `SpectrumPlotProps`
  (`uirevision`).
- `apps/desktop/e2e/tauri/msviewer-tasks.spec.ts` — `Task 4` (programmatic-zoom coverage only).
