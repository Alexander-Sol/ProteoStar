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
} from "../tauri.helpers";

// End-to-end against the REAL MsViewer webview + Rust backend (requires the app
// running via `bun run e2e:app`). Opens a real mzML, drives the actual TIC/spectrum
// Plotly handlers, and checks a known scientific value.
const YEAST = "SmallCalibratible_Yeast.mzML";

test.describe.configure({ mode: "serial" });
test.beforeEach(async ({ tauriPage }) => {
  test.setTimeout(120_000);
  await openDatasetByPath(tauriPage, dataPath(YEAST));
  await waitForTic(tauriPage);
});

test("opening the file renders its TIC", async ({ tauriPage }) => {
  const apex = await readTicApex(tauriPage);
  expect(apex.points).toBeGreaterThan(0);
  expect(apex.apexIntensity).toBeGreaterThan(0);
  // The app shows the loaded file name in its toolbar once the dataset is open.
  const body = await tauriPage.innerText("body");
  expect(body).toContain(YEAST);
});

test("clicking the TIC apex shows the MS1 spectrum whose base peak is m/z 356.1919", async ({
  tauriPage,
}) => {
  // Click the tallest point of the TIC (deterministic target) exactly as a user would.
  const apex = await readTicApex(tauriPage);
  await clickTicAtRt(tauriPage, apex.apexRt);
  await waitForSpectrum(tauriPage);

  // It should be an MS1 scan (the spectrum panel header reads "Spectrum · MS1 · …").
  const header = await tauriPage.evaluate<string>(
    `(function(){
       var el = Array.prototype.slice.call(document.querySelectorAll('*'))
         .find(function(n){ return /^Spectrum · MS\\d/.test((n.textContent||'').trim()); });
       return el ? (el.textContent||'').trim() : '';
     })()`,
  );
  // textContent runs the title straight into the subtitle ("…MS1Scan 51…"), so
  // guard the level with a negative lookahead rather than a word boundary.
  expect(header).toMatch(/Spectrum · MS1(?!\d)/);

  // The tallest (base) peak of that scan is m/z 356.1919.
  const basePeak = await readSpectrumBasePeak(tauriPage);
  console.log(`[base peak] m/z ${basePeak.mz} (intensity ${basePeak.intensity})`);
  expect(basePeak.mz).toBeCloseTo(356.1919, 3);
});
