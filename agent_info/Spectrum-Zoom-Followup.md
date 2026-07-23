# Follow-up: finish consolidating spectrum zoom onto one source of truth

**Status:** proposed follow-up. Written 2026-07-23 on branch `fix-raw-tic-flat`, after fixing the
interactive-zoom regression (see `Spectrum-Zoom-Persistence.md` for the full history and the bug
that was just fixed).

## Context

The spectrum panel's zoom persistence still rests on **two overlapping mechanisms** — a React
`spectrumViewport` (drives the Plotly `range` prop) and Plotly's own `uirevision` (preserves
interactive state across re-renders). The just-landed fix implemented "Option A" from
`Spectrum-Zoom-Persistence.md`, but **only for the x-axis**: an interactive zoom/pan now writes its
x-range into `spectrumViewport` via the `xrange-change` handler in `App.tsx`, so the `range` prop —
not `uirevision` — deterministically re-applies the zoom across re-renders (new-spectrum select,
scan step). That fixed the reported bug (zoom lost on the next selection, then stopped working).

What remains is that the y-axis and the `uirevision` machinery are still in the loop, so there is
**not yet a single source of truth for "the current zoom."** This note proposes finishing the job.

## Why finish it

- **Residual fragility.** `uirevision` still governs the y-axis interactive state. As
  `Spectrum-Zoom-Persistence.md` §"Why it's fragile" documents, any future edit to the spectrum
  `xaxis`/`yaxis` layout can silently break one path while the other keeps working — exactly the
  class of bug that just bit us. One source of truth removes that whole failure mode.
- **Redundant y-freeze.** `stepScan`'s manual y-freeze (`setSpectrumViewport(v => ({...v, yMin, yMax}))`)
  overlaps with what a captured interactive y-range would provide. Two code paths express "hold the
  y-range across scans"; they can drift.
- **Interactive y-zoom is lossy today.** Because only x is captured, a user's box-zoom y-bounds are
  not persisted — y auto-fits to the visible x-window instead. That is consistent with the
  programmatic-zoom paths and is arguably the desired behavior for spectra, but it is a *behavior*
  decision that should be made explicitly, not left as a side effect of partial capture.

## Proposed work

1. **Capture y as well as x.** Extend the plot's relayout reporting so `onRelayout` emits the
   interactive **y-range** alongside x (today `plots.tsx` reports only `xrange-change`). Write both
   into `spectrumViewport` in the same handler. Decide deliberately whether interactive y-bounds are
   honored or whether y always auto-fits to the x-window (keep the current auto-fit if that is the
   intended UX, but make it a documented choice).
2. **Drop `uirevision` as a persistence mechanism.** Once the viewport is the sole driver of the
   `range` prop, `spectrumUiRev` no longer needs to gate interactive persistence. Reframes become
   plain viewport writes; the counter can be removed or repurposed only if some Plotly quirk still
   needs a forced re-apply (verify with the e2e below).
3. **Remove the now-redundant `stepScan` y-freeze** if step 1 makes the persisted y-range cover the
   "envelope height stays stable across scans" behavior it exists for.
4. **Collapse `spectrumXView` into `spectrumViewport`.** `spectrumXView` (added for zoom-reactive
   highlights) now duplicates the viewport's x-bounds. With the viewport authoritative, derive the
   highlight window from `spectrumViewport.{xMin,xMax}` and delete the separate state.

## Testing

`apps/desktop/e2e/tauri/spectrum-interactive-zoom.spec.ts` already covers the interactive-zoom path
(zoom holds → persists across a new-spectrum select → second zoom works → autorange reset widens).
Keep it green throughout the refactor; extend it to assert interactive **y** persistence (or the
deliberate auto-fit behavior) once step 1 lands. Also keep `msviewer-tasks.spec.ts` Task 4
(programmatic zoom + y-freeze + button reset) green — if the y-freeze is removed in step 3, update
that task to reflect the new single-source-of-truth behavior rather than the freeze.

## Related upstream follow-up (from the `.raw` TIC fix on this branch)

The flat-`.raw`-TIC bug fixed on this branch was rooted in mzdata: `ThermoRawReader` calls
`update_summaries()`, which recomputes each spectrum's `total ion current` cvParam from the *loaded*
peaks and **overwrites** it — so at `DetailLevel::MetadataOnly` (no peaks loaded) it stamps the TIC
to `0`. We worked around it by not attempting the metadata TIC for `.raw` (see
`crates/flashlfq-core/src/peak_indexing.rs::read_ms1_tic_metadata`). Worth filing upstream against
mzdata so a metadata-level read preserves the instrument-reported TIC instead of zeroing it; if
fixed, the `.raw` fast path could serve a smooth MS1 TIC immediately like mzML, without waiting for
the index.
