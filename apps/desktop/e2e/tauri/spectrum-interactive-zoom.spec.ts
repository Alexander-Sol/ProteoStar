import { expect } from "@playwright/test";
import { test } from "../tauri.fixtures";
import {
  dataPath,
  openDatasetByPath,
  waitForTic,
  readTicApex,
  clickTicAtRt,
  waitForSpectrum,
  readSpectrumXRange,
  readSpectrumScan,
  clickButtonByText,
  clickSpectrumResetZoom,
  isIndexing,
} from "../tauri.helpers";
import type { TauriPage } from "@srsholmes/tauri-playwright";

// Exercises the INTERACTIVE spectrum-zoom path (mouse box-zoom), which Task-4 never covers — it
// zooms via the walkthrough seed (a programmatic reframe that bumps uirevision). A real box-zoom
// does NOT bump uirevision; it fires `plotly_relayout`, and the app captures that range into the
// spectrum viewport (the single source of truth driving the `range` prop). Regression guard for the
// reported bug: the zoom worked once, but was lost on the next new-spectrum selection and then
// stopped working entirely. Closes the interactive-zoom coverage gap noted in
// agent_info/Spectrum-Zoom-Persistence.md.
const YEAST = "SmallCalibratible_Yeast.mzML";

const FIND_SPECTRUM_GD = `
  var gd = Array.prototype.slice.call(document.querySelectorAll('.js-plotly-plot'))
    .find(function(d){
      var t = d.layout && d.layout.xaxis && d.layout.xaxis.title;
      var text = t && (t.text != null ? t.text : t);
      return text === 'm/z';
    });
`;

async function ev<T>(page: TauriPage, body: string): Promise<T> {
  return (await page.evaluate(`(function(){ ${body} })()`)) as T;
}

/** Drive a real interactive relayout on the spectrum: mutate the graph div's x-range exactly as a
 *  mouse box-zoom does, firing `plotly_relayout` → the app's onRelayout. */
async function interactiveZoom(page: TauriPage, frac0: number, frac1: number): Promise<void> {
  const full = (await readSpectrumXRange(page))!;
  const w = full[1] - full[0];
  const lo = full[0] + w * frac0;
  const hi = full[0] + w * frac1;
  await ev(
    page,
    `${FIND_SPECTRUM_GD}
     if (!gd) throw new Error('spectrum gd not found');
     window.Plotly.relayout(gd, { 'xaxis.range[0]': ${lo}, 'xaxis.range[1]': ${hi} });
     return true;`,
  );
}

async function stableXWidth(page: TauriPage): Promise<number> {
  let last: [number, number] | null = null;
  for (let i = 0; i < 12; i++) {
    await new Promise((r) => setTimeout(r, 200));
    last = await readSpectrumXRange(page);
  }
  return last![1] - last![0];
}

async function readTicApexStable(page: TauriPage) {
  let lastErr: unknown;
  for (let i = 0; i < 20; i++) {
    try {
      return await readTicApex(page);
    } catch (e) {
      lastErr = e;
      await new Promise((r) => setTimeout(r, 250));
    }
  }
  throw lastErr;
}

async function openApexSpectrum(page: TauriPage): Promise<void> {
  await openDatasetByPath(page, dataPath(YEAST));
  await waitForTic(page);
  for (let i = 0; i < 8; i++) {
    const apex = await readTicApexStable(page);
    await clickTicAtRt(page, apex.apexRt);
    try {
      await waitForSpectrum(page, 6000);
      return;
    } catch {
      /* TIC remounted / click lost — retry */
    }
  }
  throw new Error("spectrum did not load after retrying the TIC apex click");
}

test("interactive zoom persists across a new-spectrum selection and still works afterward", async ({
  tauriPage,
}) => {
  test.setTimeout(150_000);
  await openApexSpectrum(tauriPage);
  const body = await tauriPage.innerText("body");
  if (body.includes("Walkthrough: on")) await clickButtonByText(tauriPage, "Walkthrough: on");
  for (let i = 0; i < 240 && (await isIndexing(tauriPage)); i++) {
    await new Promise((r) => setTimeout(r, 250));
  }

  // Establish a true full-width baseline (the live dev app persists zoom across test runs).
  await clickSpectrumResetZoom(tauriPage);
  const fullWidth = await stableXWidth(tauriPage);
  expect(fullWidth).toBeGreaterThan(0);

  // (0) Initial interactive zoom holds.
  await interactiveZoom(tauriPage, 0.45, 0.55);
  const w0 = await stableXWidth(tauriPage);
  console.log(`[zoom] initial full=${fullWidth.toFixed(1)} -> held ${w0.toFixed(1)}`);
  expect(w0).toBeLessThan(fullWidth * 0.5);

  // (1) Select a NEW spectrum (TIC click at a different RT). The zoom must persist.
  const apex = await readTicApexStable(tauriPage);
  const scanBefore = await readSpectrumScan(tauriPage);
  await clickTicAtRt(tauriPage, apex.apexRt * 0.6 + 0.01);
  for (let i = 0; i < 40; i++) {
    await new Promise((r) => setTimeout(r, 200));
    if ((await readSpectrumScan(tauriPage)) !== scanBefore) break;
  }
  const w1 = await stableXWidth(tauriPage);
  console.log(`[zoom] after new-spectrum select -> width ${w1.toFixed(1)}`);
  expect(w1).toBeLessThan(fullWidth * 0.5);

  // (2) Interactive zoom AGAIN (to a different, narrower window) must still work.
  await interactiveZoom(tauriPage, 0.2, 0.28);
  const w2 = await stableXWidth(tauriPage);
  console.log(`[zoom] second zoom -> width ${w2.toFixed(1)}`);
  expect(w2).toBeLessThan(fullWidth * 0.35);

  // (3) Double-click autorange reset (the `xaxis.autorange` relayout path) widens back to full.
  await ev(
    tauriPage,
    `${FIND_SPECTRUM_GD}
     if (!gd) throw new Error('spectrum gd not found');
     window.Plotly.relayout(gd, { 'xaxis.autorange': true });
     return true;`,
  );
  const w3 = await stableXWidth(tauriPage);
  console.log(`[zoom] after autorange reset -> width ${w3.toFixed(1)}`);
  expect(w3).toBeGreaterThan(fullWidth * 0.5);
});
