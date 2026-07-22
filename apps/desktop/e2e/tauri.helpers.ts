// Helpers for driving the real MsViewer webview in tauri mode. Everything goes
// through TauriPage.evaluate (arbitrary JS in the webview) because:
//   - the native "Open file…" dialog can't be automated, but the frontend reaches
//     it via __TAURI_INTERNALS__.invoke('plugin:dialog|open'), which we intercept;
//   - the TIC/Spectrum are Plotly plots — we trigger the app's real click handler
//     with gd.emit('plotly_click', …) and read results straight off the graph div.
import { resolve, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import type { TauriPage } from "@srsholmes/tauri-playwright";

const HERE = dirname(fileURLToPath(import.meta.url));

/** Absolute path to a copied test data file under e2e/data. */
export function dataPath(name: string): string {
  return resolve(HERE, "data", name);
}

/** Small helper: run JS in the webview, JSON round-trip the result. */
function ev<T = unknown>(page: TauriPage, body: string): Promise<T> {
  return page.evaluate<T>(`(function(){ ${body} })()`);
}

/**
 * Open a dataset by absolute path via the app's E2E hook (window.__msviewerOpenPath),
 * skipping the native file dialog — which can't be automated. The hook exists only
 * when the app runs with the control server (--features e2e-testing). Fire-and-forget:
 * pair with `waitForTic`.
 */
export async function openDatasetByPath(page: TauriPage, absPath: string): Promise<void> {
  await ev(
    page,
    `if (typeof window.__msviewerOpenPath !== 'function')
       throw new Error('E2E hook window.__msviewerOpenPath missing — run the app with: bun run e2e:app');
     window.__msviewerOpenPath(${JSON.stringify(absPath)});
     return true;`,
  );
}

// Locate a Plotly graph div by its x-axis title ("Retention time (min)" = TIC,
// "m/z" = spectrum). Returned as a JS expression fragment naming `gd`.
const FIND_GD = (axisTitle: string) => `
  var gd = Array.prototype.slice.call(document.querySelectorAll('.js-plotly-plot'))
    .find(function(d){
      var t = d.layout && d.layout.xaxis && d.layout.xaxis.title;
      var text = t && (t.text != null ? t.text : t);
      return text === ${JSON.stringify(axisTitle)};
    });
`;

/** Wait until the TIC plot has rendered at least one point. */
export async function waitForTic(page: TauriPage, timeoutMs = 90_000): Promise<void> {
  await page.waitForFunction(
    `(function(){ ${FIND_GD("Retention time (min)")}
       return !!(gd && gd.data && gd.data[0] && gd.data[0].x && gd.data[0].x.length > 0);
     })()`,
    timeoutMs,
  );
}

export interface TicApex {
  apexRt: number;
  apexIntensity: number;
  points: number;
}

/** Read the TIC's apex (max-intensity) point from the rendered graph. */
export function readTicApex(page: TauriPage): Promise<TicApex> {
  return ev<TicApex>(
    page,
    `${FIND_GD("Retention time (min)")}
     if (!gd) throw new Error('TIC plot not found');
     var line = (gd.data||[]).find(function(t){ return t.x && t.y && String(t.mode||'').indexOf('lines') >= 0; }) || gd.data[0];
     if (!line || !line.x || !line.x.length) throw new Error('TIC has no points');
     var xi = 0, ym = -Infinity;
     for (var i = 0; i < line.y.length; i++) { if (line.y[i] > ym) { ym = line.y[i]; xi = i; } }
     return { apexRt: line.x[xi], apexIntensity: ym, points: line.x.length };`,
  );
}

/**
 * Fire the TIC's real click handler at retention time `rt` — same path a user
 * click takes (Plotly plotly_click → area-click → loads the nearest scan).
 */
export async function clickTicAtRt(page: TauriPage, rt: number): Promise<void> {
  await ev(
    page,
    `${FIND_GD("Retention time (min)")}
     if (!gd) throw new Error('TIC plot not found');
     gd.emit('plotly_click', { points: [{ customdata: { retentionTime: ${rt} } }] });
     return true;`,
  );
}

/** Wait until the spectrum plot has rendered peaks. */
export async function waitForSpectrum(page: TauriPage, timeoutMs = 30_000): Promise<void> {
  await page.waitForFunction(
    `(function(){ ${FIND_GD("m/z")}
       return !!(gd && gd.data && gd.data.some(function(t){ return t.x && t.x.length > 0; }));
     })()`,
    timeoutMs,
  );
}

export interface Peak {
  mz: number;
  intensity: number;
  peaks: number;
}

/** Read the tallest (base) peak of the currently displayed spectrum. */
export function readSpectrumBasePeak(page: TauriPage): Promise<Peak> {
  return ev<Peak>(
    page,
    `${FIND_GD("m/z")}
     if (!gd) throw new Error('spectrum plot not found');
     var tr = (gd.data||[]).find(function(t){ return t.type === 'bar' && t.x && t.x.length; })
           || (gd.data||[]).find(function(t){ return t.x && t.x.length; });
     if (!tr) throw new Error('spectrum has no peaks');
     var xi = 0, ym = -Infinity;
     for (var i = 0; i < tr.y.length; i++) { if (tr.y[i] > ym) { ym = tr.y[i]; xi = i; } }
     return { mz: tr.x[xi], intensity: ym, peaks: tr.x.length };`,
  );
}

// ----------------------------------------------------------------------------
// Task-verification helpers (no-scroll layout, MS1 peak labels, resizable
// drawer, persisted spectrum y-zoom). All read/drive the real webview DOM.

/** Page-overflow + spectrum-visibility metrics (Task 1: no-scroll layout). */
export interface Overflow {
  scrollHeight: number;
  innerHeight: number;
  /** Bottom edge (px from top) of the spectrum plot's x-axis title; -1 if not found. */
  spectrumAxisBottom: number;
}

export function readOverflow(page: TauriPage): Promise<Overflow> {
  // Measure the spectrum graph div's bottom edge (its bottom margin holds the x-axis title). The div
  // is guaranteed present once `waitForSpectrum` returns, unlike the SVG title <text>, which paints a
  // beat later and can read as missing right after a remount.
  return ev<Overflow>(
    page,
    `${FIND_GD("m/z")}
     var axisBottom = gd ? gd.getBoundingClientRect().bottom : -1;
     return {
       scrollHeight: document.documentElement.scrollHeight,
       innerHeight: window.innerHeight,
       spectrumAxisBottom: axisBottom
     };`,
  );
}

/** The text of every annotation on the spectrum plot (Task 2: m/z + charge labels). */
export function readSpectrumAnnotations(page: TauriPage): Promise<string[]> {
  return ev<string[]>(
    page,
    `${FIND_GD("m/z")}
     if (!gd) throw new Error('spectrum plot not found');
     return ((gd.layout && gd.layout.annotations) || []).map(function(a){ return a.text; });`,
  );
}

/** The spectrum plot's current y-axis range (Task 4: persisted y-zoom). */
export function readSpectrumYRange(page: TauriPage): Promise<[number, number] | null> {
  return ev<[number, number] | null>(
    page,
    `${FIND_GD("m/z")}
     if (!gd) throw new Error('spectrum plot not found');
     var r = gd.layout && gd.layout.yaxis && gd.layout.yaxis.range;
     return r ? [r[0], r[1]] : null;`,
  );
}

/** The spectrum plot's current x-axis range (used to detect an applied x-zoom). */
export function readSpectrumXRange(page: TauriPage): Promise<[number, number] | null> {
  return ev<[number, number] | null>(
    page,
    `${FIND_GD("m/z")}
     if (!gd) throw new Error('spectrum plot not found');
     var r = gd.layout && gd.layout.xaxis && gd.layout.xaxis.range;
     return r ? [r[0], r[1]] : null;`,
  );
}

/** The 1-based scan number currently shown in the spectrum panel header, or -1. */
export function readSpectrumScan(page: TauriPage): Promise<number> {
  return ev<number>(
    page,
    `var m = (document.body.textContent || '').match(/Scan\\s+(\\d+)/);
     return m ? parseInt(m[1], 10) : -1;`,
  );
}

/** Click a button whose visible text contains `substr`. Throws if none is found. */
export async function clickButtonByText(page: TauriPage, substr: string): Promise<void> {
  await ev(
    page,
    `var btn = Array.prototype.slice.call(document.querySelectorAll('button'))
       .find(function(b){ return (b.textContent || '').indexOf(${JSON.stringify(substr)}) >= 0; });
     if (!btn) throw new Error('button not found: ' + ${JSON.stringify(substr)});
     btn.click();
     return true;`,
  );
}

/** Click the spectrum panel's own "Reset zoom" button (the TIC panel has one too — target the one
 *  in the same panel <section> as the m/z spectrum plot). Throws if it isn't present. */
export async function clickSpectrumResetZoom(page: TauriPage): Promise<void> {
  await ev(
    page,
    `${FIND_GD("m/z")}
     if (!gd) throw new Error('spectrum plot not found');
     var sec = gd.closest('section');
     if (!sec) throw new Error('spectrum panel <section> not found');
     var btn = Array.prototype.slice.call(sec.querySelectorAll('button'))
       .find(function(b){ return (b.textContent || '').indexOf('Reset zoom') >= 0; });
     if (!btn) throw new Error('spectrum Reset zoom button not found');
     btn.click();
     return true;`,
  );
}

/** Seed the walkthrough charge-ladder by firing a peak-click at m/z `mz` on the spectrum (the same
 *  event a user's click produces). Only has an effect while the walkthrough is on. */
export async function seedLadderAt(page: TauriPage, mz: number, intensity: number): Promise<void> {
  await ev(
    page,
    `${FIND_GD("m/z")}
     if (!gd) throw new Error('spectrum plot not found');
     gd.emit('plotly_click', { points: [{ customdata: { mz: ${mz}, intensity: ${intensity} } }] });
     return true;`,
  );
}

/** Current width (px) of the right-side drawer, or null when no drawer is open. */
export function readDrawerWidth(page: TauriPage): Promise<number | null> {
  return ev<number | null>(
    page,
    `var d = document.querySelector('[data-testid=drawer]');
     return d ? d.getBoundingClientRect().width : null;`,
  );
}

/** Drag the drawer's resize handle so the pointer ends at viewport x = `toClientX`. The drawer is
 *  right-anchored, so the resulting width ≈ innerWidth − toClientX (clamped to [280, 760]). */
export async function dragDrawerResizeTo(page: TauriPage, toClientX: number): Promise<void> {
  await ev(
    page,
    `var h = document.querySelector('[data-testid=drawer-resize]');
     if (!h) throw new Error('drawer resize handle not found');
     var rect = h.getBoundingClientRect();
     var y = rect.top + rect.height / 2;
     h.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: rect.left + 2, clientY: y }));
     window.dispatchEvent(new PointerEvent('pointermove', { bubbles: true, clientX: ${toClientX}, clientY: y }));
     window.dispatchEvent(new PointerEvent('pointerup', { bubbles: true, clientX: ${toClientX}, clientY: y }));
     return true;`,
  );
}

/** Press ArrowRight/ArrowLeft on the window to step to the next/previous MS1 scan. */
export async function pressArrow(page: TauriPage, dir: "ArrowRight" | "ArrowLeft"): Promise<void> {
  await ev(
    page,
    `window.dispatchEvent(new KeyboardEvent('keydown', { key: ${JSON.stringify(dir)}, bubbles: true }));
     return true;`,
  );
}

/** Whether the toolbar still reports indexing in progress ("MS1 scans … indexing…"). */
export function isIndexing(page: TauriPage): Promise<boolean> {
  return ev<boolean>(page, `return (document.body.textContent || '').indexOf('indexing…') >= 0;`);
}
