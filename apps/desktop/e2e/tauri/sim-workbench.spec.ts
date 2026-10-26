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
  pressArrow,
  isIndexing,
} from "../tauri.helpers";

// The simulation workbench against the REAL webview + backend (requires `bun run e2e:app`):
// a spec seeded from the open run renders a live simulated scan mirrored under the real one,
// edits re-render it, and Accept draws the simulated run's TIC mirrored under the real TIC.
const YEAST = "SmallCalibratible_Yeast.mzML";

function ev<T = unknown>(page: TauriPage, body: string): Promise<T> {
  return page.evaluate<T>(`(function(){ ${body} })()`);
}

const FIND_GD = (axisTitle: string) => `
  var gd = Array.prototype.slice.call(document.querySelectorAll('.js-plotly-plot'))
    .find(function(d){
      var t = d.layout && d.layout.xaxis && d.layout.xaxis.title;
      var text = t && (t.text != null ? t.text : t);
      return text === ${JSON.stringify(axisTitle)};
    });
`;

interface Mirror {
  /** Peaks (or points) drawn below zero, and the deepest one's height. */
  count: number;
  depth: number;
}

/** The mirrored (negative-y) trace on the plot with this x-axis title. */
function readMirror(page: TauriPage, axisTitle: string): Promise<Mirror> {
  return ev<Mirror>(
    page,
    `${FIND_GD(axisTitle)}
     if (!gd) return { count: 0, depth: 0 };
     var best = { count: 0, depth: 0 };
     (gd.data || []).forEach(function(t){
       if (!t.y || !t.y.length) return;
       var n = 0, d = 0;
       for (var i = 0; i < t.y.length; i++) { if (t.y[i] < 0) { n++; if (-t.y[i] > d) d = -t.y[i]; } }
       if (n > best.count) best = { count: n, depth: d };
     });
     return best;`,
  );
}

/** Set a workbench number field by its label text, as typing would (React-controlled input). */
async function setField(page: TauriPage, label: string, value: number): Promise<void> {
  await ev(
    page,
    `var lab = Array.prototype.slice.call(document.querySelectorAll('[data-testid=drawer] label'))
       .find(function(l){ var s = l.querySelector('span'); return s && (s.textContent||'').indexOf(${JSON.stringify(label)}) === 0; });
     if (!lab) throw new Error('field not found: ' + ${JSON.stringify(label)});
     var input = lab.querySelector('input[type=number]');
     var set = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
     set.call(input, ${JSON.stringify(String(value))});
     input.dispatchEvent(new Event('input', { bubbles: true }));
     return true;`,
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

function drawerText(page: TauriPage): Promise<string> {
  return ev<string>(
    page,
    `var d = document.querySelector('[data-testid=drawer]'); return d ? d.textContent || '' : '';`,
  );
}

async function waitFor(check: () => Promise<boolean>, what: string, timeoutMs = 30_000): Promise<void> {
  const end = Date.now() + timeoutMs;
  while (Date.now() < end) {
    if (await check()) return;
    await new Promise((r) => setTimeout(r, 150));
  }
  throw new Error(`timed out waiting for ${what}`);
}

test.describe.configure({ mode: "serial" });

// The app instance is shared across specs: leave no simulation behind for them to trip over.
test.afterEach(async ({ tauriPage }) => {
  await discardSimulation(tauriPage);
});

test("the workbench simulates the viewed scan live and the run on Accept", async ({ tauriPage }) => {
  test.setTimeout(240_000);
  await openDatasetByPath(tauriPage, dataPath(YEAST));
  await waitForTic(tauriPage);
  await waitFor(async () => !(await isIndexing(tauriPage)), "the index", 120_000);

  const apex = await readTicApex(tauriPage);
  await clickTicAtRt(tauriPage, apex.apexRt);
  await waitForSpectrum(tauriPage);

  // The app may still hold a simulation from an earlier run of this spec (its toolbar button then
  // reads "Simulation"), and the drawer may already be open; "New from run" restarts either way.
  if (!(await drawerText(tauriPage)).includes("Simulation workbench")) {
    await clickButtonByText(tauriPage, "Simulat");
  }
  await waitFor(async () => (await drawerText(tauriPage)).includes("New from run"), "the panel");
  await clickButtonByText(tauriPage, "New from run");

  // A live simulated scan appears, mirrored below the real one.
  await waitFor(async () => (await readMirror(tauriPage, "m/z")).count > 0, "the live scan");
  const noisy = await readMirror(tauriPage, "m/z");
  console.log(`[live] ${noisy.count} simulated peaks, tallest ${noisy.depth.toExponential(2)}`);
  const status = await drawerText(tauriPage);
  console.log(`[status] ${status.match(/\d+ scans[^A-Z]*/)?.[0]}`);
  expect(status).toMatch(/viewing scan \d+/);

  // Dropping the noise leaves only the proteoform's peaks.
  await setNoiseEnabled(tauriPage, false);
  await waitFor(
    async () => {
      const m = await readMirror(tauriPage, "m/z");
      return m.count > 0 && m.count < noisy.count;
    },
    "the noise-free scan",
  );
  const clean = await readMirror(tauriPage, "m/z");
  console.log(`[noise off] ${clean.count} simulated peaks, tallest ${clean.depth.toExponential(2)}`);

  // Doubling the abundance doubles the (noise-free) peaks.
  await setField(tauriPage, "Abundance A", 2e7);
  await waitFor(
    async () => Math.abs((await readMirror(tauriPage, "m/z")).depth / clean.depth - 2) < 0.01,
    "the scan to double",
  );

  // Accept simulates the run and draws its TIC mirrored under the real one.
  await clickButtonByText(tauriPage, "Accept");
  await waitFor(
    async () => (await readMirror(tauriPage, "Retention time (min)")).count > 0,
    "the simulated TIC",
    120_000,
  );
  expect(await drawerText(tauriPage)).toContain("Run up to date");

  // An edit after Accept makes the run stale again.
  await setField(tauriPage, "Abundance A", 3e7);
  await waitFor(async () => (await drawerText(tauriPage)).includes("Accept — simulate run"), "stale");

  // Stepping the real scan moves the simulated one with it.
  const before = (await drawerText(tauriPage)).match(/viewing scan (\d+)/)?.[1];
  await pressArrow(tauriPage, "ArrowRight");
  await waitFor(
    async () => (await drawerText(tauriPage)).match(/viewing scan (\d+)/)?.[1] !== before,
    "the simulated scan to follow",
  );
  const after = (await drawerText(tauriPage)).match(/viewing scan (\d+)/)?.[1];
  console.log(`[step] simulated scan ${before} → ${after}`);
  expect(Number(after)).toBe(Number(before) + 1);
});
