import { expect } from "@playwright/test";
import { test } from "../tauri.fixtures";
import { dataPath, openDatasetByPath, waitForTic, waitForSpectrum } from "../tauri.helpers";
import type { TauriPage } from "@srsholmes/tauri-playwright";

// End-to-end against the REAL MsViewer webview + Rust backend (requires `bun run e2e:app`).
// Loads real target + decoy resolved-feature TSVs (produced by detect_features_tsv on the bundled
// yeast mzML — target run and a COMB_MODEL=decoy run) and checks that BOTH overlay layers render:
//   - the TIC rug carries charge-coloured target markers AND muted decoy markers, and
//   - selecting a decoy draws its isotope envelope over the spectrum in the muted decoy colour.
const YEAST = "SmallCalibratible_Yeast.mzML";
const TARGET_TSV = "yeast_features.tsv"; // generated into e2e/data before running
const DECOY_TSV = "yeast_decoys.tsv";

const DECOY_COLOR = "#9098a1"; // App.tsx DECOY_COLOR
const CHARGE_COLORS = [
  "#2f6fb0", "#c0392b", "#27ae60", "#8e44ad",
  "#d35400", "#16a085", "#b7950b", "#c2185b"
];

function ev<T = unknown>(page: TauriPage, body: string): Promise<T> {
  return page.evaluate<T>(`(function(){ ${body} })()`);
}

// The TIC rug is one scattergl markers trace (triangle-up) whose marker.color is a per-point array.
const FIND_TIC = `
  var gd = Array.prototype.slice.call(document.querySelectorAll('.js-plotly-plot'))
    .find(function(d){
      var t = d.layout && d.layout.xaxis && d.layout.xaxis.title;
      var text = t && (t.text != null ? t.text : t);
      return text === 'Retention time (min)';
    });`;
const FIND_SPECTRUM = `
  var gd = Array.prototype.slice.call(document.querySelectorAll('.js-plotly-plot'))
    .find(function(d){
      var t = d.layout && d.layout.xaxis && d.layout.xaxis.title;
      var text = t && (t.text != null ? t.text : t);
      return text === 'm/z';
    });`;

async function loadFeaturesByPath(page: TauriPage, absPath: string): Promise<void> {
  await ev(
    page,
    `if (typeof window.__msviewerLoadFeaturesPath !== 'function')
       throw new Error('E2E hook __msviewerLoadFeaturesPath missing — rebuild with e2e-testing');
     window.__msviewerLoadFeaturesPath(${JSON.stringify(absPath)}); return true;`
  );
}
async function loadDecoysByPath(page: TauriPage, absPath: string): Promise<void> {
  await ev(
    page,
    `if (typeof window.__msviewerLoadDecoysPath !== 'function')
       throw new Error('E2E hook __msviewerLoadDecoysPath missing — rebuild with e2e-testing');
     window.__msviewerLoadDecoysPath(${JSON.stringify(absPath)}); return true;`
  );
}

interface Rug {
  count: number;
  colors: string[];
  decoyIndex: number | null; // a customdata.featureIndex belonging to a decoy marker
}
function readRug(page: TauriPage): Promise<Rug | null> {
  return ev<Rug | null>(
    page,
    `${FIND_TIC}
     if (!gd || !gd.data) return null;
     var rug = gd.data.find(function(t){ return t.marker && t.marker.symbol === 'triangle-up'; });
     if (!rug) return null;
     var colors = Array.isArray(rug.marker.color) ? rug.marker.color : [rug.marker.color];
     var decoyIndex = null;
     var cd = rug.customdata || [];
     for (var i = 0; i < cd.length; i++) {
       if (cd[i] && cd[i].featureIndex >= 10000000) { decoyIndex = cd[i].featureIndex; break; }
     }
     return { count: rug.x ? rug.x.length : 0, colors: colors, decoyIndex: decoyIndex };`
  );
}

// Fire the TIC's real feature-click handler for a given rug featureIndex (selects that feature).
async function clickFeature(page: TauriPage, featureIndex: number): Promise<void> {
  await ev(
    page,
    `${FIND_TIC}
     if (!gd) throw new Error('TIC not found');
     gd.emit('plotly_click', { points: [{ customdata: { featureIndex: ${featureIndex} } }] });
     return true;`
  );
}

// Colours of the vertical envelope lines drawn over the spectrum (layout.shapes, type 'line').
function readEnvelopeColors(page: TauriPage): Promise<string[]> {
  return ev<string[]>(
    page,
    `${FIND_SPECTRUM}
     if (!gd) throw new Error('spectrum not found');
     return ((gd.layout && gd.layout.shapes) || [])
       .filter(function(s){ return s.type === 'line'; })
       .map(function(s){ return s.line && s.line.color; });`
  );
}

test.describe.configure({ mode: "serial" });

test("target + decoy feature overlays both render on the TIC and spectrum", async ({ tauriPage }) => {
  test.setTimeout(120_000);
  await openDatasetByPath(tauriPage, dataPath(YEAST));
  await waitForTic(tauriPage);

  await loadFeaturesByPath(tauriPage, dataPath(TARGET_TSV));
  await loadDecoysByPath(tauriPage, dataPath(DECOY_TSV));

  // Wait until both layers have populated the rug (targets + at least one decoy marker).
  await tauriPage.waitForFunction(
    `(function(){ ${FIND_TIC}
       if (!gd || !gd.data) return false;
       var rug = gd.data.find(function(t){ return t.marker && t.marker.symbol === 'triangle-up'; });
       if (!rug) return false;
       var colors = Array.isArray(rug.marker.color) ? rug.marker.color : [rug.marker.color];
       return colors.indexOf(${JSON.stringify(DECOY_COLOR)}) >= 0;
     })()`,
    60_000
  );

  const rug = await readRug(tauriPage);
  if (!rug) throw new Error("feature rug not found on the TIC");
  console.log(`[rug] ${rug.count} markers; distinct colours: ${[...new Set(rug.colors)].join(", ")}`);

  // Muted decoy markers present…
  expect(rug.colors).toContain(DECOY_COLOR);
  // …alongside charge-coloured target markers.
  expect(rug.colors.some((c) => CHARGE_COLORS.includes(c))).toBe(true);
  expect(rug.decoyIndex).not.toBeNull();

  // Selecting a decoy draws its isotope envelope over the spectrum in the muted decoy colour.
  await clickFeature(tauriPage, rug.decoyIndex as number);
  await waitForSpectrum(tauriPage);
  await tauriPage.waitForFunction(
    `(function(){ ${FIND_SPECTRUM}
       if (!gd) return false;
       return ((gd.layout && gd.layout.shapes) || [])
         .some(function(s){ return s.type === 'line' && s.line && s.line.color === ${JSON.stringify(DECOY_COLOR)}; });
     })()`,
    30_000
  );
  const envColors = await readEnvelopeColors(tauriPage);
  console.log(`[decoy envelope] ${envColors.length} lines; colours: ${[...new Set(envColors)].join(", ")}`);
  expect(envColors.length).toBeGreaterThan(0);
  expect(envColors.every((c) => c === DECOY_COLOR)).toBe(true);
});
