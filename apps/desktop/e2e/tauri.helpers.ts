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
