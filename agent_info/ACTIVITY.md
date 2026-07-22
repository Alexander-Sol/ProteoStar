# MsViewer overnight activity

Autonomous overnight run: implement and Playwright-test four MsViewer desktop-app tasks pulled from
`TODO.md` (MsViewer section). Update the per-task **Status** and the **Progress log** as work lands.

Source of truth for the task descriptions: `TODO.md`. When a task here is finished and verified,
also check it off there.

## Scope / ground rules

- Frontend lives in `apps/desktop/src/App.tsx` + workspace packages `packages/ui` (shell/panels) and
  `packages/plot-adapter` (Plotly TIC/Spectrum plots).
- Test harness: `apps/desktop/e2e` (see its `README.md`). Two modes:
  - **browser** — headless Chromium, IPC mocked. Fast, but the app only calls the backend after a
    dataset is open, so spectrum-dependent checks are limited here.
  - **tauri** — real webview + Rust backend; opens real data (`SmallCalibratible_Yeast.mzML`). This is
    where the label/zoom behavior is actually exercised. Reference: `e2e/tauri/open-and-spectrum.spec.ts`.
- Prefer: unit tests (vitest) for pure logic + tauri-mode Playwright for the visual/data behavior.
- Keep changes minimal and in-style; don't disturb the detection engine or the walkthrough feature.

## Tasks

### 1. No-scroll layout — spectrum plot fully visible
Everything (TIC + full spectrum) must be visible at once; no page scroll, no cut-off spectrum bottom.
- Target: `packages/ui/src/layout.tsx` (`ViewerShell`/`Panel` grid tracks), possibly a global
  `html,body,#root` height/margin reset (currently none in `index.html`/`main.tsx`).
- Likely cause: grid rows use `1fr` (= `minmax(auto,1fr)`) so a tall Plotly min-size can push content
  past the viewport; switch load-bearing tracks to `minmax(0,1fr)` and zero body margin / full-height root.
- Acceptance: with a spectrum loaded, `document.documentElement.scrollHeight <= innerHeight` (no vertical
  page scroll) and the spectrum panel's plot bottom (x-axis label) is within the viewport.
- **Status:** VERIFIED (typecheck + unit + tauri e2e). ✅
  - `apps/desktop/index.html` — added a `<style>`: `html,body{margin:0}` + `#root{height:100vh;overflow:hidden}`
    (the default 8px body margin was pushing the 100vh shell past the viewport → the small scroll).
  - `packages/ui/src/layout.tsx` — shell/content/panel grid tracks `1fr` → `minmax(0,1fr)`, `boxSizing:border-box`
    on the shell, so a tall Plotly min-height can't force overflow.
  - **e2e:** tauri `Task 1` (real data): `document.documentElement.scrollHeight (783) == innerHeight (783)`, no page
    scroll, and the spectrum graph div's bottom (773.8) is within the viewport. Fast browser-mode guard added to
    `e2e/browser/smoke.spec.ts` ("the shell does not scroll vertically": scrollHeight 720 == innerHeight 720).

### 2. Label m/z + charge for prominent MS1 features
When viewing an MS1 spectrum, annotate the most prominent peaks with m/z and inferred charge.
- Target: `packages/plot-adapter` (add optional peak annotations to `SpectrumPlot` via Plotly
  `annotations`) + `App.tsx` (compute the labels for the current MS1 spectrum) + a small charge
  estimator (isotope spacing Δm/z ≈ 1.00235/z among nearby peaks).
- Acceptance: for the reference scan, the top-N peaks carry visible `m/z … · z…` labels; charge
  estimator unit-tested on a synthetic isotope envelope.
- **Status:** VERIFIED (unit + tauri e2e). ✅
  - `apps/desktop/src/annotate.ts` — `estimateCharge` (isotope-spacing search, Δ≈1.00235/z),
    `computePeakLabels` (top-N prominent peaks in the x-window, deduped), `peakLabelText`.
  - `apps/desktop/src/annotate.test.tsx` — 8 passing tests (z=1/2/3 envelopes, null for lone peak, label text).
  - `packages/plot-adapter` — new `PeakAnnotation` type + `SpectrumPlot annotations` prop rendered as Plotly
    `annotations` (vertical labels above peaks, `buildPeakAnnotations` in `plots.tsx`).
  - `apps/desktop/src/App.tsx` — `spectrumAnnotations` useMemo (MS1 only, suppressed in walkthrough) → SpectrumPlot.
  - **e2e:** tauri `Task 2` (real yeast apex scan) reads `gd.layout.annotations` and gets 6 labels with real
    m/z + charge: `["356.19 · z1","433.74 · z2","356.69 · z2","711.38 · z1","546.32 · z1","466.77 · z2"]` — every
    label carries an m/z and at least one an inferred charge.

### 3. Resizable side panel (a.k.a. "resize the PSM panel")
**Caveat — there is no PSM panel yet.** The right side currently hosts the Features / Walkthrough
**drawers** (fixed 360px, `drawerStyle` in `App.tsx`). The PSM panel, when built, will be a right-side
drawer just like these. So this task = make the right-side drawer width user-resizable via a draggable
divider; the future PSM panel reuses the same mechanism.
- Target: `App.tsx` drawer(s) — add a drag handle on the left edge, width in state, min/max clamp; keep
  `ViewerShell rightInset` in sync so plots reflow.
- Acceptance: dragging the handle changes the drawer width (browser/tauri Playwright drag check);
  width clamps to sane bounds.
- **Status:** VERIFIED (tauri e2e). ✅ (Caveat below still stands — flag for the user.)
  - `apps/desktop/src/App.tsx` — new `ResizableDrawer` frame with a `col-resize` handle on the left edge
    (`data-testid="drawer"` / `"drawer-resize"`); `drawerWidth` state (default 360, clamp 280–760);
    `LadderDrawer`/`FeatureDrawer` now render inside it; `ViewerShell rightInset` follows `drawerWidth`.
  - **e2e:** tauri `Task 3` opens the walkthrough drawer, then a synthesized pointer drag on the handle moves the
    width to a controlled 520px and clamps at both ends (drag past-right → 760, drag past-left → 280). Start width
    isn't asserted (the running dev app persists `drawerWidth` between runs), so the test proves the drag+clamp
    behavior independent of start state.
  - **CAVEAT:** no PSM panel exists yet — implemented on the existing right-side drawer as the reusable
    mechanism the PSM panel will use. Confirm with user this is the intended interpretation.

### 4. Persist spectrum zoom — stop the y-axis continually rescaling
The x-zoom already persists across scan steps, but the y-axis is recomputed from each scan's visible
peaks every render (`visibleSpectrumYRange` in `plots.tsx`), so a user's y-zoom is overridden and an
isotope envelope of interest shrinks as a taller peak enters view on neighboring scans.
- Target: `plot-adapter` viewport (extend `PlotViewport` with optional y bounds), `plots.tsx` (honor a
  set y-range instead of always refitting; capture y from Plotly relayout), `App.tsx` (persist across
  scan steps; auto-fit only on fresh load / reset zoom).
- Acceptance: with a set viewport, stepping scans keeps the y-range fixed (envelope height stable);
  "Reset zoom" refits. Unit-test the "use persisted y when present, else fit" selector.
- **Status:** VERIFIED (unit + tauri e2e). ✅
  - `packages/plot-adapter` — `PlotViewport` gains optional `yMin/yMax`; new pure helpers in `viewport.ts`
    (`fitSpectrumYRange`, `hasPersistedY`, `resolveSpectrumYRange`); `SpectrumPlot` y-axis now uses
    `resolveSpectrumYRange` (honors persisted y, else auto-fits). Removed the old `visibleSpectrumYRange`.
  - `packages/plot-adapter/tests/viewport.test.ts` — 7 passing tests.
  - `apps/desktop/src/App.tsx` — `stepScan` freezes y (from the current fit) on the first arrow step while
    x-zoomed, so the envelope stays a stable height across scans; `handleAreaClick` clears frozen y (fresh
    pick auto-fits); "Reset zoom" (via `createDefaultViewport`) clears y too.
  - **e2e:** tauri `Task 4` seeds the walkthrough ladder on the base peak (x-zooms to a 4.0-m/z comb window),
    reads `y=[0,67787701.8]`, presses ArrowRight to step scans, and re-reads the SAME `y=[0,67787701.8]`
    (envelope height held). "Reset zoom" then widens the x-window from 4.0 → 1361.4 m/z (viewport refit). The
    real x-zoom is driven through the walkthrough seed path since that's the app's only non-feature x-zoom entry.
  - NOTE: this addresses the reported y-rescale complaint. It does NOT add capturing the user's *manual*
    mouse y-zoom via Plotly `onRelayout` — a possible enhancement, left for the follow-up.

## Progress log

### 2026-07-22
- Deleted the stale MsBrowser-era `ACTIVITY.md`; recreated it around these four tasks.
- Reconnaissance done: mapped `App.tsx`, `packages/ui/layout.tsx`, `packages/plot-adapter/plots.tsx`,
  and the e2e harness. Recorded the key findings above (no-scroll cause, y-refit cause, no PSM panel).
- Implemented all four tasks (details under each task's Status). Static verification green:
  `bun run typecheck` → exit 0; `bun run test` → 15/15 unit tests pass (8 charge/label, 7 viewport).
- Stopped here for handoff after static verification only.

### 2026-07-22 (runtime/e2e verification — all four tasks now VERIFIED ✅)
- Brought up the real app (`bun run e2e:app`, listening on tcp://127.0.0.1:6274) and drove all four tasks
  against the **real webview + Rust backend** on `SmallCalibratible_Yeast.mzML`.
- Added `apps/desktop/e2e/tauri/msviewer-tasks.spec.ts` (4 tests, one per task) plus reusable helpers in
  `e2e/tauri.helpers.ts` (`readOverflow`, `readSpectrumAnnotations`, `readSpectrum{X,Y}Range`,
  `readSpectrumScan`, `clickButtonByText`, `seedLadderAt`, `readDrawerWidth`, `dragDrawerResizeTo`,
  `pressArrow`, `isIndexing`). Added a fast no-scroll guard to `e2e/browser/smoke.spec.ts`.
- **Results (all green):** tauri e2e **9/9** (4 new task tests + reference `open-and-spectrum` + smoke);
  browser e2e **3/3**; `bun run typecheck` exit 0; `bun run test` **15/15**; e2e tsconfig typecheck exit 0.
  Per-task evidence recorded under each task's Status above.
- Fixes made while verifying: drawer width reads ~1px over style width (border box) → assert with ±2px
  slack; Task 1 measures the spectrum graph div's bottom (guaranteed by `waitForSpectrum`) instead of the
  SVG `.xtitle` text (which paints a beat later after a remount) and polls until it lays out; Task 3 is
  start-state-agnostic because the long-running dev app persists `drawerWidth` across runs.

## Handoff to next agent

**State:** DONE. All four tasks are code-complete AND runtime-verified. `bun run typecheck` exit 0,
`bun run test` 15/15, browser e2e 3/3, tauri e2e 9/9. The four MsViewer items are checked off in `TODO.md`.

**If you re-run the tauri e2e:** the real app must be up first — `cd apps/desktop && bun run e2e:app`, wait
for `listening on tcp://127.0.0.1:6274`, then in another shell `$env:PW_NO_WEBSERVER=1; bun run e2e:tauri`
(or filter to the new file: `bun x playwright test --config e2e/playwright.config.ts --project=tauri
msviewer-tasks`). On Windows kill leftovers on 1420/6274 if a run wedges (see project `CLAUDE.md`).

**Still open / flag to the user:**
- **Task 3 deviation (no PSM panel).** Implemented on the existing right-side drawer as the reusable
  resize mechanism the future PSM panel will use. Confirm this is the intended interpretation.
- Task 4 intentionally *reverses* the old "reset spectrum viewport on new scan" behavior (expected, not a
  regression); it does not yet capture the user's *manual* mouse y-zoom via Plotly `onRelayout` (follow-up).
- The real app persists UI state (e.g. `drawerWidth`) across e2e runs since it's one long-lived dev
  process; tests that depend on a pristine default should not assume it.

## Open decisions / risks
- Charge inference from a single MS1 spectrum is heuristic (isotope-spacing search). Keep it simple and
  only label peaks where a confident spacing is found; don't over-claim charge.
- Tauri-mode e2e needs a Rust build with `--features e2e-testing` and the app running on TCP 6274 — the
  first build is slow. If it proves too flaky overnight, fall back to vitest unit tests + browser-mode
  layout checks and document what wasn't run.
- Task 3 deviates from the literal wording (no PSM panel exists) — implemented on the existing drawer as
  the reusable mechanism. Flag for user confirmation.
