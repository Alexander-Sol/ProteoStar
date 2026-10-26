import { expect } from "@playwright/test";
import type { TauriPage } from "@srsholmes/tauri-playwright";
import { test } from "../tauri.fixtures";
import {
  dataPath,
  openDatasetByPath,
  waitForTic,
  readTicApex,
  clickTicAtRt,
  waitForSpectrum,
  clickButtonByText,
  discardSimulation,
  isIndexing,
} from "../tauri.helpers";

// The workbench on REAL top-down data (requires `bun run e2e:app`): one minute of Jurkat rep2
// fract7 (44–45 min: 22 MS1 scans + the 60 MS2 behind its 157 confident proteoform IDs) and the
// SimSpec mzLib fitted to those IDs on the slice's own MS1 scan times. Written by mzLib's
// `AnalysisExample.ExportRep2WorkbenchFixture`.
const SLICE = "02-18-20_jurkat_td_rep2_fract7.rt44-45";
const MS1_SCANS = 22;
const PROTEOFORMS = 131;
const ID_ROWS = 157;

function ev<T = unknown>(page: TauriPage, body: string): Promise<T> {
  return page.evaluate<T>(`(function(){ ${body} })()`);
}

async function waitFor(check: () => Promise<boolean>, what: string, timeoutMs = 30_000): Promise<void> {
  const end = Date.now() + timeoutMs;
  while (Date.now() < end) {
    if (await check()) return;
    await new Promise((r) => setTimeout(r, 150));
  }
  throw new Error(`timed out waiting for ${what}`);
}

function bodyText(page: TauriPage): Promise<string> {
  return ev<string>(page, `return document.body.textContent || '';`);
}

function drawerText(page: TauriPage): Promise<string> {
  return ev<string>(
    page,
    `var d = document.querySelector('[data-testid=drawer]'); return d ? d.textContent || '' : '';`,
  );
}

interface Peaks {
  real: { mz: number; intensity: number }[];
  sim: { mz: number; intensity: number }[];
}

/** The spectrum plot's peak bars: the real scan upright, the simulated scan mirrored below. */
function readPeaks(page: TauriPage): Promise<Peaks> {
  return ev<Peaks>(
    page,
    `var gd = Array.prototype.slice.call(document.querySelectorAll('.js-plotly-plot'))
       .find(function(d){ var t = d.layout && d.layout.xaxis && d.layout.xaxis.title;
                          return (t && (t.text != null ? t.text : t)) === 'm/z'; });
     var out = { real: [], sim: [] };
     if (!gd) return out;
     (gd.data || []).forEach(function(t){
       if (t.type !== 'bar' || t.width !== 0.001 || !t.x) return;
       for (var i = 0; i < t.x.length; i++) {
         if (t.y[i] > 0) out.real.push({ mz: t.x[i], intensity: t.y[i] });
         else if (t.y[i] < 0) out.sim.push({ mz: t.x[i], intensity: -t.y[i] });
       }
     });
     return out;`,
  );
}

/** RTs of the real TIC (the first trace drawn; XIC overlays may follow it) and of the simulated
 *  TIC (the trace mirrored below zero). */
function readTics(page: TauriPage): Promise<{ real: number[]; sim: number[] }> {
  return ev<{ real: number[]; sim: number[] }>(
    page,
    `var gd = Array.prototype.slice.call(document.querySelectorAll('.js-plotly-plot'))
       .find(function(d){ var t = d.layout && d.layout.xaxis && d.layout.xaxis.title;
                          return (t && (t.text != null ? t.text : t)) === 'Retention time (min)'; });
     var out = { real: [], sim: [] };
     if (!gd) return out;
     var lines = (gd.data || []).filter(function(t){
       return t.y && t.y.length && String(t.mode || '').indexOf('lines') >= 0;
     });
     var rts = function(t){ return Array.prototype.map.call(t.x, Number); };
     var mirrored = lines.find(function(t){ return Array.prototype.some.call(t.y, function(y){ return y < 0; }); });
     if (lines.length > 0 && lines[0] !== mirrored) out.real = rts(lines[0]);
     if (mirrored) out.sim = rts(mirrored);
     return out;`,
  );
}

async function setNoiseEnabled(page: TauriPage, on: boolean): Promise<void> {
  await ev(
    page,
    `var lab = Array.prototype.slice.call(document.querySelectorAll('[data-testid=drawer] label'))
       .find(function(l){ return (l.textContent||'').trim() === 'enabled'; });
     if (!lab) throw new Error('noise toggle not found');
     var box = lab.querySelector('input[type=checkbox]');
     if (box.checked !== ${on}) box.click();
     return true;`,
  );
}

/** Of the `n` tallest simulated peaks, the fraction with a real peak within `ppm`. */
function landingFraction(peaks: Peaks, n: number, ppm: number): number {
  const real = peaks.real.map((p) => p.mz).sort((a, b) => a - b);
  const top = [...peaks.sim].sort((a, b) => b.intensity - a.intensity).slice(0, n);
  const hit = (mz: number) => {
    const tol = (mz * ppm) / 1e6;
    let lo = 0;
    let hi = real.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (real[mid] < mz - tol) lo = mid + 1;
      else hi = mid;
    }
    return lo < real.length && real[lo] <= mz + tol;
  };
  return top.filter((p) => hit(p.mz)).length / top.length;
}

test.describe.configure({ mode: "serial" });

test.afterEach(async ({ tauriPage }) => {
  await discardSimulation(tauriPage);
});

test("a SimSpec fitted to real top-down data lands on that data", async ({ tauriPage }) => {
  test.setTimeout(240_000);
  await openDatasetByPath(tauriPage, dataPath(`${SLICE}.mzML`));
  // Opening is fire-and-forget and an earlier spec's dataset may still be on screen, so wait for
  // this file's own name and indexed MS1 count in the toolbar, not just for "a" TIC.
  await waitFor(
    async () => {
      const text = await bodyText(tauriPage);
      return text.includes(`${SLICE}.mzML`) && text.includes(`MS1 scans${MS1_SCANS}`);
    },
    "the slice to open and index",
    120_000,
  );
  await waitForTic(tauriPage);

  // A real scan, then the fitted spec on top of it.
  const apex = await readTicApex(tauriPage);
  await clickTicAtRt(tauriPage, apex.apexRt);
  await waitForSpectrum(tauriPage);
  await ev(
    tauriPage,
    `if (typeof window.__msviewerLoadSimSpecPath !== 'function')
       throw new Error('E2E hook __msviewerLoadSimSpecPath missing — run the app with: bun run e2e:app');
     window.__msviewerLoadSimSpecPath(${JSON.stringify(dataPath(`${SLICE}.simspec.json`))});
     return true;`,
  );
  if (!(await drawerText(tauriPage)).includes("Simulation workbench")) {
    await clickButtonByText(tauriPage, "Simulat");
  }
  await waitFor(async () => (await readPeaks(tauriPage)).sim.length > 0, "the live simulated scan");
  const status = await drawerText(tauriPage);
  expect(status).toContain(`${MS1_SCANS} scans · ${PROTEOFORMS} proteoforms`);

  // The spec was seeded from this slice's scan times, so the simulated scan is the real one's.
  const realRt = (await bodyText(tauriPage)).match(/RT (\d+\.\d\d) min/)?.[1];
  const simRt = status.match(/viewing scan \d+ \(RT (\d+\.\d\d)\)/)?.[1];
  console.log(`[aligned] real RT ${realRt}, simulated RT ${simRt}`);
  expect(simRt).toBe(realRt);

  // Without noise, the simulated scan is the fitted proteoforms alone: their tallest peaks should
  // sit on real peaks of the same scan.
  const noisy = (await readPeaks(tauriPage)).sim.length;
  await setNoiseEnabled(tauriPage, false);
  let peaks = await readPeaks(tauriPage);
  await waitFor(async () => {
    peaks = await readPeaks(tauriPage);
    return peaks.sim.length > 0 && peaks.sim.length < noisy;
  }, "the noise-free scan");
  console.log(`[noise] ${noisy} simulated peaks with noise, ${peaks.sim.length} without`);
  const landed = landingFraction(peaks, 50, 10);
  console.log(
    `[landing] ${peaks.sim.length} simulated vs ${peaks.real.length} real peaks; ` +
      `${(landed * 100).toFixed(0)}% of the 50 tallest simulated peaks have a real peak within 10 ppm`,
  );
  expect(landed).toBeGreaterThan(0.5);

  // Accept: the simulated run's TIC has a point for every real MS1 scan, at the same RTs.
  await clickButtonByText(tauriPage, "Accept");
  await waitFor(async () => (await readTics(tauriPage)).sim.length > 0, "the simulated TIC", 120_000);
  const tics = await readTics(tauriPage);
  expect(tics.sim.length).toBe(MS1_SCANS);
  expect(tics.sim.map((rt) => rt.toFixed(4))).toEqual(tics.real.map((rt) => rt.toFixed(4)));

  // The fixture's IDs load against the slice.
  await ev(
    tauriPage,
    `window.__msviewerLoadPsmsPath(${JSON.stringify(dataPath(`${SLICE}_Proteoforms.psmtsv`))}); return true;`,
  );
  await waitFor(async () => (await bodyText(tauriPage)).includes(`PSMs (${ID_ROWS})`), "the PSMs");
});
