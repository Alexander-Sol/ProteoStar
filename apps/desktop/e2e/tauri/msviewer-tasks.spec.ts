import { expect } from "@playwright/test";
import { test } from "../tauri.fixtures";
import {
  dataPath,
  openDatasetByPath,
  waitForTic,
  readTicApex,
  clickTicAtRt,
  waitForSpectrum,
  readSpectrumBasePeak,
  readOverflow,
  readSpectrumAnnotations,
  readSpectrumYRange,
  readSpectrumXRange,
  readSpectrumScan,
  clickButtonByText,
  seedLadderAt,
  readDrawerWidth,
  dragDrawerResizeTo,
  pressArrow,
  isIndexing,
} from "../tauri.helpers";
import type { TauriPage } from "@srsholmes/tauri-playwright";

// Runtime verification of the four MsViewer tasks against the REAL webview + Rust backend
// (requires the app running via `bun run e2e:app`). Each task's pure logic is already unit-tested;
// these specs prove the integrated behavior in the running app.
//   1. No-scroll layout — spectrum fully visible, no page scroll.
//   2. MS1 peak labels — prominent peaks annotated with m/z (+ inferred charge).
//   3. Resizable drawer — drag the handle to change width, clamped to [280, 760].
//   4. Persisted spectrum y-zoom — stepping scans keeps the y-range fixed; Reset zoom refits.
const YEAST = "SmallCalibratible_Yeast.mzML";

/** Toggle the walkthrough off if it is currently on (keeps tests order-independent). */
async function ensureWalkthroughOff(page: TauriPage): Promise<void> {
  const body = await page.innerText("body");
  if (body.includes("Walkthrough: on")) await clickButtonByText(page, "Walkthrough: on");
}

/** Poll a boolean predicate from Node across evaluates (each evaluate is capped at 30s). */
async function poll(fn: () => Promise<boolean>, timeoutMs: number, label: string): Promise<void> {
  const start = Date.now();
  while (Date.now() - start < timeoutMs) {
    if (await fn()) return;
    await new Promise((r) => setTimeout(r, 250));
  }
  throw new Error(`timed out after ${timeoutMs}ms waiting for: ${label}`);
}

test.describe.configure({ mode: "serial" });

test.beforeEach(async ({ tauriPage }) => {
  test.setTimeout(120_000);
  await openDatasetByPath(tauriPage, dataPath(YEAST));
  await waitForTic(tauriPage);
  const apex = await readTicApex(tauriPage);
  await clickTicAtRt(tauriPage, apex.apexRt);
  await waitForSpectrum(tauriPage);
  await ensureWalkthroughOff(tauriPage);
});

// ---------------------------------------------------------------- Task 1: no-scroll layout
test("Task 1 — the page does not scroll and the spectrum's x-axis is within the viewport", async ({
  tauriPage,
}) => {
  // Poll until the spectrum plot has laid out (its graph div has a measurable bottom edge) — right
  // after a remount the div can be momentarily absent.
  await poll(async () => (await readOverflow(tauriPage)).spectrumAxisBottom > 0, 15_000, "spectrum plot to lay out");
  const o = await readOverflow(tauriPage);
  console.log(`[layout] scrollHeight=${o.scrollHeight} innerHeight=${o.innerHeight} axisBottom=${o.spectrumAxisBottom}`);
  // No vertical page scroll (allow a couple px for subpixel rounding).
  expect(o.scrollHeight).toBeLessThanOrEqual(o.innerHeight + 2);
  // The spectrum plot (its bottom margin carries the x-axis title) sits inside the viewport — not clipped.
  expect(o.spectrumAxisBottom).toBeGreaterThan(0);
  expect(o.spectrumAxisBottom).toBeLessThanOrEqual(o.innerHeight + 2);
});

// ------------------------------------------------------------ Task 2: MS1 peak m/z + charge labels
test("Task 2 — the MS1 spectrum's prominent peaks carry m/z (and charge) labels", async ({
  tauriPage,
}) => {
  const labels = await readSpectrumAnnotations(tauriPage);
  console.log(`[annotations] ${JSON.stringify(labels)}`);
  // Top-N prominent peaks are labeled.
  expect(labels.length).toBeGreaterThanOrEqual(1);
  // Every label starts with an m/z value (e.g. "356.19" or "356.19 · z1").
  for (const t of labels) expect(t).toMatch(/^\d+\.\d{2}/);
  // At least one peak's charge was inferred from its isotope spacing ("· z2").
  expect(labels.some((t) => /·\s*z\d/.test(t))).toBe(true);
});

// ------------------------------------------------------------------- Task 3: resizable drawer
test("Task 3 — the right-side drawer is resizable by dragging its handle, clamped to [280,760]", async ({
  tauriPage,
}) => {
  // Open a drawer (the walkthrough drawer appears immediately; no backend needed).
  await clickButtonByText(tauriPage, "Walkthrough");
  await tauriPage.waitForSelector("[data-testid=drawer]");
  await tauriPage.waitForSelector("[data-testid=drawer-resize]");

  // Measured width is ~1px over the style width (border box); allow a 2px slack throughout. The
  // starting width isn't asserted — the running dev app persists drawerWidth across test runs — so
  // this proves the drag sets a controlled width and clamps, independent of the start state.
  const initial = await readDrawerWidth(tauriPage);
  console.log(`[drawer] initial width=${initial}`);
  expect(initial).not.toBeNull();

  const innerWidth = await tauriPage.evaluate<number>("window.innerWidth");

  // Drag so the pointer lands 520px from the right edge → width ≈ 520.
  await dragDrawerResizeTo(tauriPage, innerWidth - 520);
  const widened = await readDrawerWidth(tauriPage);
  console.log(`[drawer] after drag→520 width=${widened}`);
  expect(Math.abs((widened as number) - 520)).toBeLessThanOrEqual(2); // dragging controls the width

  // Drag past the right edge → clamps to the max (760).
  await dragDrawerResizeTo(tauriPage, 0);
  const maxed = await readDrawerWidth(tauriPage);
  console.log(`[drawer] after drag→max width=${maxed}`);
  expect(Math.abs((maxed as number) - 760)).toBeLessThanOrEqual(2);

  // Drag past the left/right so target width < min → clamps to the min (280).
  await dragDrawerResizeTo(tauriPage, innerWidth);
  const minned = await readDrawerWidth(tauriPage);
  console.log(`[drawer] after drag→min width=${minned}`);
  expect(Math.abs((minned as number) - 280)).toBeLessThanOrEqual(2);

  await ensureWalkthroughOff(tauriPage);
});

// -------------------------------------------------------- Task 4: persisted spectrum y-zoom
test("Task 4 — stepping scans holds the spectrum y-range fixed; Reset zoom refits", async ({
  tauriPage,
}) => {
  // stepScan and ladder scoring both need the peak index; wait for indexing to finish.
  await poll(async () => !(await isIndexing(tauriPage)), 60_000, "indexing to finish");

  // Zoom in on a peak's isotope comb via the walkthrough seed (the app's real x-zoom path).
  await clickButtonByText(tauriPage, "Walkthrough");
  const base = await readSpectrumBasePeak(tauriPage);
  await seedLadderAt(tauriPage, base.mz, base.intensity);

  // Wait until the spectrum has actually zoomed to a narrow window (the comb around the seed).
  await poll(
    async () => {
      const r = await readSpectrumXRange(tauriPage);
      return !!r && r[1] - r[0] > 0 && r[1] - r[0] < 60;
    },
    30_000,
    "spectrum to x-zoom onto the seed comb",
  );
  const zoomedX = (await readSpectrumXRange(tauriPage))!;
  const zoomedWidth = zoomedX[1] - zoomedX[0];

  const scanBefore = await readSpectrumScan(tauriPage);
  const y0 = await readSpectrumYRange(tauriPage);
  console.log(`[y-zoom] scan=${scanBefore} y0=${JSON.stringify(y0)} xWidth=${zoomedWidth.toFixed(2)}`);
  expect(y0).not.toBeNull();

  // Step to the next MS1 scan — the y-range must stay frozen (envelope height stable).
  await pressArrow(tauriPage, "ArrowRight");
  await poll(async () => (await readSpectrumScan(tauriPage)) !== scanBefore, 15_000, "scan to advance");
  const y1 = await readSpectrumYRange(tauriPage);
  console.log(`[y-zoom] after step y1=${JSON.stringify(y1)}`);
  expect(y1).not.toBeNull();
  // Frozen y is copied verbatim across the step → identical bounds (allow tiny fp slack).
  expect(Math.abs((y1 as number[])[0] - (y0 as number[])[0])).toBeLessThan(1e-3 * ((y0 as number[])[1] || 1) + 1e-6);
  expect(Math.abs((y1 as number[])[1] - (y0 as number[])[1])).toBeLessThan(1e-3 * ((y0 as number[])[1] || 1) + 1e-6);

  // Reset zoom refits: the x-viewport widens back out (persisted zoom cleared).
  await clickButtonByText(tauriPage, "Reset zoom");
  await poll(
    async () => {
      const r = await readSpectrumXRange(tauriPage);
      return !!r && r[1] - r[0] > zoomedWidth * 3;
    },
    15_000,
    "Reset zoom to widen the x-range back out",
  );
  const resetX = (await readSpectrumXRange(tauriPage))!;
  console.log(`[y-zoom] after reset xWidth=${(resetX[1] - resetX[0]).toFixed(2)} (was ${zoomedWidth.toFixed(2)})`);
  expect(resetX[1] - resetX[0]).toBeGreaterThan(zoomedWidth * 3);

  await ensureWalkthroughOff(tauriPage);
});
