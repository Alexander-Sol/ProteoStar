import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { open as openFileDialog, save as saveFileDialog } from "@tauri-apps/plugin-dialog";
import { getCurrentWebview } from "@tauri-apps/api/webview";

import {
  TicPlot,
  SpectrumPlot,
  createDefaultViewport,
  fitSpectrumYRange,
  hasPersistedY
} from "@msbrowser/plot-adapter";
import type {
  TicPlotTrace,
  SpectrumPlotTrace,
  FeatureMarker,
  PeakAnnotation,
  RtRegion,
  EnvelopeLine,
  SpectrumPeakHighlight,
  PlotViewport
} from "@msbrowser/plot-adapter";
import {
  Panel,
  PanelActionButton,
  PanelHeader,
  StatusBanner,
  ViewerShell,
  MetricReadout
} from "@msbrowser/ui";

import { openDataset } from "./tauri-dataset-provider";
import {
  loadFeatures,
  exportMs1Features,
  runFeatureDetection,
  isotopeGrid,
  featureElutesAt
} from "./features";
import { computePeakLabels } from "./annotate";
import { loadPsms, linkPsms, fetchIsotopeXics, sumXics } from "./psms";
import { scoreSeedLadder } from "./walkthrough";
import type {
  DatasetMetadata,
  DatasetProvider,
  Feature,
  LadderCharge,
  Psm,
  PsmLink,
  ScanSummary,
  SeedLadder,
  Spectrum,
  TicPoint,
  XicPoint
} from "./contract";

const SLOT_COLOR = "#2f6fb0";
const SELECTED_COLOR = "#e8830c";
// Charge-keyed colours for the charge-centric views only: the walkthrough ladder comb + its drawer
// swatches, and the per-isotope XIC traces. Features themselves are coloured by FEATURE_COLORS —
// a feature's colour must not depend on which plot it is drawn in.
const CHARGE_COLORS = [
  "#2f6fb0", "#c0392b", "#27ae60", "#8e44ad",
  "#d35400", "#16a085", "#b7950b", "#c2185b"
];

// THE feature palette — the single source of colour for a feature everywhere it is drawn: the TIC
// rug, the spectrum overlay, and the drawer swatch. Each feature's colour is assigned by its rank in
// GLOBAL m/z order (see `featureColorByIndex`), so a feature keeps the same colour across scans and
// across plots (the top priority) while features near each other in m/z still tend to differ. 8
// distinct hues keep adjacent-in-m/z clashes rare.
const FEATURE_COLORS = [
  "#2f6fb0", "#c0392b", "#27ae60", "#8e44ad",
  "#d35400", "#16a085", "#c2185b", "#b7950b"
];
// Render guard: draw at most this many eluting features' highlights per scan. The overlay re-iterates
// on zoom (only features with a peak in the visible m/z window compete), so far fewer are in play when
// zoomed and this cap can stay modest. NOT a load cap — the whole TSV is in memory.
const SCAN_FEATURE_CAP = 80;
const chargeColor = (z: number): string =>
  CHARGE_COLORS[(Math.max(1, z) - 1) % CHARGE_COLORS.length];

// Decoy features (from a separate decoy-model detector run) are drawn in a single muted colour so
// the real, charge-coloured features stand out; the selected item — target or decoy — is always
// SELECTED_COLOR.
const DECOY_COLOR = "#9098a1";
const DECOY_REGION_FILL = "rgba(144,152,161,0.16)";
// Feature-rug markers carry a numeric featureIndex that comes back on click. Decoy markers are
// offset by this base so the one click handler can tell a target (< base) from a decoy (>= base)
// without threading a discriminator through the shared plot-adapter event type.
const DECOY_INDEX_BASE = 10_000_000;

// Mass-spec dataset extensions accepted by both the "Open file…" dialog and drag-and-drop.
const DATASET_EXTENSIONS = ["raw", "mzML", "mzml", "mzMLb", "mgf"] as const;
// Lowercased for case-insensitive matching of dropped paths.
const DATASET_EXTENSIONS_LC = DATASET_EXTENSIONS.map((e) => e.toLowerCase());
const hasDatasetExtension = (path: string): boolean => {
  const dot = path.lastIndexOf(".");
  if (dot < 0) return false;
  return DATASET_EXTENSIONS_LC.includes(path.slice(dot + 1).toLowerCase());
};

// Features are sorted by intensity on load; these cap what's drawn/listed so a
// pathologically large TSV (top-down noise runs can be 100k+ features) stays
// responsive. The strongest features are always the ones kept.
const RUG_CAP = 8000;
const DRAWER_CAP = 800;
// The feature drawer paginates this many rows per page (so every feature is reachable via the pager,
// and the page can auto-follow a selection made by clicking a peak / rug marker).
const DRAWER_PAGE = 30;

// m/z tolerance for marking a comb tooth "detected" against the displayed scan's peaks. Matches the
// detector's default ppm tolerance so the overlay's notion of a hit agrees with detection.
const DETECT_PPM = 10;

// Right-side drawer sizing (features / walkthrough / future PSM panel). User-resizable via the drag
// handle on the drawer's left edge; the width clamps to [MIN, MAX].
const DEFAULT_DRAWER_WIDTH = 360;
const MIN_DRAWER_WIDTH = 280;
const MAX_DRAWER_WIDTH = 760;
const DRAWER_INSET_GAP = 12;
const clampDrawerWidth = (w: number): number =>
  Math.max(MIN_DRAWER_WIDTH, Math.min(MAX_DRAWER_WIDTH, w));

type LoadState =
  | { status: "idle" }
  | { status: "loading"; message: string }
  | { status: "ready" }
  | { status: "error"; message: string };

export function App() {
  const [load, setLoad] = useState<LoadState>({ status: "idle" });
  // Transient confirmation for actions that produce a file rather than changing the view
  // (currently the ms1.feature export). Shown in the side panel, not over the plot, and
  // cleared on a timer so it can't be mistaken for current state.
  const [notice, setNotice] = useState<string | null>(null);
  const [provider, setProvider] = useState<DatasetProvider | null>(null);
  const [handle, setHandle] = useState<number | null>(null);
  const [detecting, setDetecting] = useState<string | null>(null);
  const [metadata, setMetadata] = useState<DatasetMetadata | null>(null);
  const [ticPoints, setTicPoints] = useState<readonly TicPoint[]>([]);
  const [spectrum, setSpectrum] = useState<Spectrum | null>(null);
  // Per-MS1-scan summaries (RT-ordered; `scanIndex` is the MS1-array index) — the basis for
  // arrow-key scan stepping, which navigates by RT to avoid the file-vs-MS1 index mismatch.
  const [scanSummaries, setScanSummaries] = useState<readonly ScanSummary[]>([]);

  const [features, setFeatures] = useState<readonly Feature[]>([]);
  const [featuresFile, setFeaturesFile] = useState<string | null>(null);
  const [selected, setSelected] = useState<number | null>(null);
  const [drawerOpen, setDrawerOpen] = useState(false);

  // Decoy features live in their own layer (drawn muted, kept out of PSM linking). Loaded from a
  // separate decoy-run resolved TSV. `selectedDecoy` is mutually exclusive with `selected`.
  const [decoys, setDecoys] = useState<readonly Feature[]>([]);
  const [selectedDecoy, setSelectedDecoy] = useState<number | null>(null);
  const [decoyDrawerOpen, setDecoyDrawerOpen] = useState(false);

  // Show every feature eluting at the displayed MS1 scan (cycling colours on the spectrum). On by
  // default: viewing an MS1 spectrum shows all co-eluting features; toggle off to inspect just the
  // selected feature's envelope. Runs over the full loaded feature set (the 800 cap is drawer-only).
  const [showScanFeatures, setShowScanFeatures] = useState(true);

  // PSMs loaded from a MetaMorpheus .psmtsv, linked to detected features client-side.
  // Selecting one extracts its isotope XICs (below) and jumps the spectrum to its RT.
  const [psms, setPsms] = useState<readonly Psm[]>([]);
  const [psmLinks, setPsmLinks] = useState<readonly PsmLink[]>([]);
  const [selectedPsm, setSelectedPsm] = useState<number | null>(null);
  const [psmDrawerOpen, setPsmDrawerOpen] = useState(false);

  // XICs of the selected PSM's most-abundant isotopologues, extracted over the whole run and
  // overlaid on the TIC. `xicMode` toggles between one summed-envelope trace and one trace per
  // isotopologue; `xicIndices` are the isotope indices of `xic`'s traces (for labelling).
  const [xic, setXic] = useState<XicPoint[][]>([]);
  const [xicIndices, setXicIndices] = useState<number[]>([]);
  const [xicMode, setXicMode] = useState<"sum" | "isotopes">("sum");
  const [xicLabel, setXicLabel] = useState<string | null>(null);

  const [ticViewport, setTicViewport] = useState<PlotViewport>(createDefaultViewport());
  const [spectrumViewport, setSpectrumViewport] = useState<PlotViewport>(createDefaultViewport());

  // Plotly `uirevision` for the TIC: bumped only when we intentionally reframe (new file, feature
  // or PSM selection, reset-zoom), so the user's interactive zoom survives ordinary re-renders
  // (e.g. clicking an XIC point to load a spectrum) instead of snapping back to the prop range.
  const [ticUiRev, setTicUiRev] = useState(0);
  const reframeTic = useCallback((vp: PlotViewport) => {
    setTicViewport(vp);
    setTicUiRev((r) => r + 1);
  }, []);

  // Same idea for the spectrum: bumped only on an intentional reframe (feature / PSM selection,
  // walkthrough comb framing, reset-zoom), so the user's zoom stays put while stepping scans with
  // the arrow keys or loading a scan by clicking the TIC/XIC — neither of which reframes.
  const [spectrumUiRev, setSpectrumUiRev] = useState(0);
  // The currently visible m/z window (null = full view). Tracked from reframes AND user zoom/pan
  // (SpectrumPlot's xrange-change event) so the scan overlay can re-iterate on zoom — showing the
  // features in the visible window (including weak ones) rather than only the strongest globally.
  const [spectrumXView, setSpectrumXView] = useState<{ min: number; max: number } | null>(null);
  const reframeSpectrum = useCallback((vp: PlotViewport) => {
    setSpectrumViewport(vp);
    setSpectrumXView(vp.xMin !== null && vp.xMax !== null ? { min: vp.xMin, max: vp.xMax } : null);
    setSpectrumUiRev((r) => r + 1);
  }, []);

  const [ticPinned, setTicPinned] = useState(false);
  const [spectrumPinned, setSpectrumPinned] = useState(false);

  // User-resizable width (px) of the right-side drawer (features / walkthrough / future PSM panel).
  const [drawerWidth, setDrawerWidth] = useState(DEFAULT_DRAWER_WIDTH);

  // True between open (fast TIC shown) and the background peak index landing. While set,
  // spectra / range-XIC / detection aren't available yet (backend returns INDEXING).
  const [indexing, setIndexing] = useState(false);

  // ------------------------------------------------ feature-finding walkthrough
  // When active, clicking a spectrum peak picks a seed and scores its full top-down
  // charge-state ladder at the chosen anchoring charge `zSeed`; the ladder's combs
  // are overlaid on the spectrum and its per-charge scores shown in a right drawer.
  const [walkthrough, setWalkthrough] = useState(false);
  const [ladder, setLadder] = useState<SeedLadder | null>(null);
  const [zSeed, setZSeed] = useState(1);
  // Which charge's comb to isolate in the overlay; null = all retained charges.
  const [focusCharge, setFocusCharge] = useState<number | null>(null);
  const [ladderBusy, setLadderBusy] = useState(false);

  // ------------------------------------------------------------ open raw/mzML
  // Open a dataset by path — the shared core of both the "Open file…" button and
  // the E2E hook below.
  const openPath = useCallback(async (picked: string) => {
    setLoad({ status: "loading", message: "Opening…" });
    setSpectrum(null);
    setScanSummaries([]);
    setSelected(null);
    reframeTic(createDefaultViewport());
    try {
      const { handle: h, provider: p } = await openDataset(picked, (progress) => {
        setLoad({ status: "loading", message: `${progress.phase}…` });
      });
      const meta = await p.getMetadata();
      // Paint the fast-path TIC immediately (no peak decode, no index). For mzML it's the smooth
      // MS1-only trace from per-scan metadata — final, no later swap. For Thermo `.raw` it's the
      // provisional native full TIC (all MS levels), which `markReady` (below) swaps for the smooth
      // MS1 trace once the index lands. Empty only for files exposing neither, which fill in from
      // the index via `markReady`; the panel shows a "building index" note until then.
      const tic = await p.getTicTrace({ maxPoints: 4000 });
      setHandle(h);
      setProvider(p);
      setMetadata(meta);
      setTicPoints(tic);
      // ms1ScanCount is 0 until the background index lands; the effect below flips this off.
      setIndexing(meta.ms1ScanCount === 0);
      setLoad({ status: "ready" });
    } catch (err) {
      setLoad({ status: "error", message: errMessage(err) });
    }
  }, [reframeTic]);

  const handleOpenFile = useCallback(async () => {
    const picked = await openFileDialog({
      multiple: false,
      filters: [{ name: "Mass spec", extensions: [...DATASET_EXTENSIONS] }]
    });
    if (typeof picked !== "string") return;
    await openPath(picked);
  }, [openPath]);

  // Drag-and-drop open. Tauri v2 keeps OS-level drag-drop enabled by default (so HTML5
  // drop events never reach the webview); we subscribe to the native webview drag-drop
  // event, which hands us real filesystem paths — exactly what `openPath` wants. `dragOver`
  // drives the drop overlay while a drag hovers the window.
  const [dragOver, setDragOver] = useState(false);
  useEffect(() => {
    let unlisten: (() => void) | undefined;
    let active = true;
    void getCurrentWebview()
      .onDragDropEvent((event) => {
        const { type } = event.payload;
        if (type === "over") {
          setDragOver(true);
        } else if (type === "leave") {
          setDragOver(false);
        } else if (type === "drop") {
          setDragOver(false);
          // Open the first supported file; ignore folders / unsupported drops.
          const picked = event.payload.paths.find(hasDatasetExtension);
          if (picked) void openPath(picked);
        }
      })
      .then((fn) => {
        if (active) unlisten = fn;
        else fn();
      });
    return () => {
      active = false;
      unlisten?.();
    };
  }, [openPath]);

  // E2E hook: open a dataset by path, skipping the native file dialog (which
  // Playwright can't automate). Exposed only when the app is driven by the
  // tauri-plugin-playwright control server — which sets `window.__PW_ACTIVE__`
  // via its init script — so it is completely inert in normal/production runs.
  useEffect(() => {
    const w = window as typeof window & {
      __PW_ACTIVE__?: boolean;
      __msviewerOpenPath?: (p: string) => void;
    };
    if (!w.__PW_ACTIVE__) return;
    w.__msviewerOpenPath = (p: string) => void openPath(p);
    return () => {
      delete w.__msviewerOpenPath;
    };
  }, [openPath]);

  // Flip `indexing` off once the background peak index is ready (spectra / XIC / detection
  // light up). We poll `getMetadata()` rather than wait on a Tauri event: the backend sets
  // `ms1ScanCount > 0` (under lock) the moment the index lands, so it's the authoritative
  // "done" signal and doesn't depend on event delivery. Also surfaces a stuck index (count
  // never leaves 0 → backend never finished).
  useEffect(() => {
    if (handle === null || provider === null) return;
    let active = true;
    let timer: ReturnType<typeof setTimeout> | undefined;

    const markReady = async () => {
      if (!active) return;
      setIndexing(false);
      try {
        setMetadata(await provider.getMetadata());
        // Re-fetch the TIC. For mzML (per-scan MS1 TIC metadata) this returns the same smooth
        // MS1-only trace already painted on open (no visible change). For Thermo `.raw` the fast
        // path only had the native full TIC (or nothing); this swaps in the smooth MS1-only trace
        // from the freshly-built index.
        setTicPoints(await provider.getTicTrace({ maxPoints: 4000 }));
        // MS1 scan summaries power arrow-key scan stepping (available post-index).
        setScanSummaries(await provider.getScanSummaries());
      } catch {
        /* keep provisional metadata / TIC if the refetch fails */
      }
    };

    const poll = async () => {
      if (!active) return;
      try {
        const meta = await provider.getMetadata();
        if (!active) return;
        if (meta.indexError) {
          // Background index build failed — surface it instead of hanging on "indexing…".
          setIndexing(false);
          setLoad({ status: "error", message: `Indexing failed: ${meta.indexError}` });
          return; // done — stop polling
        }
        if (meta.ms1ScanCount > 0) {
          await markReady();
          return; // done — stop polling
        }
      } catch {
        /* transient (e.g. handle races a close); retry below */
      }
      if (active) timer = setTimeout(() => void poll(), 400);
    };
    void poll();

    return () => {
      active = false;
      if (timer) clearTimeout(timer);
    };
  }, [handle, provider]);

  // A side drawer (features or walkthrough) narrows the plot area. react-plotly's
  // `useResizeHandler` only listens to *window* resize, so the plots don't reflow on their own
  // when the drawer opens/closes — nudge them with a synthetic resize once the layout has painted.
  const drawerVisible =
    walkthrough ||
    (psmDrawerOpen && psms.length > 0) ||
    (drawerOpen && features.length > 0) ||
    (decoyDrawerOpen && decoys.length > 0);
  useEffect(() => {
    const id = window.setTimeout(() => window.dispatchEvent(new Event("resize")), 60);
    return () => window.clearTimeout(id);
  }, [drawerVisible]);

  // Re-link PSMs to features whenever either list changes (a fresh detection run, a
  // reloaded feature TSV, or newly loaded PSMs). Unlinked when no features are loaded yet.
  useEffect(() => {
    if (psms.length === 0) {
      setPsmLinks([]);
      return;
    }
    setPsmLinks(
      features.length > 0
        ? linkPsms(psms, features)
        : psms.map(() => ({ featureIndex: null, massPpmError: null, rtDelta: null }))
    );
  }, [psms, features]);

  // ------------------------------------------------------------ load features
  const loadFeaturesFromPath = useCallback(async (picked: string) => {
    try {
      const feats = await loadFeatures(picked);
      // Strongest first, so the caps below keep the most important features.
      const sorted = [...feats].sort((a, b) => b.summedIntensity - a.summedIntensity);
      setFeatures(sorted);
      setFeaturesFile(picked.split(/[\\/]/).pop() ?? picked);
      setSelected(null);
      setDrawerOpen(true);
    } catch (err) {
      setLoad({ status: "error", message: errMessage(err) });
    }
  }, []);

  const handleLoadFeatures = useCallback(async () => {
    const picked = await openFileDialog({
      multiple: false,
      // ".feature" covers TopFD / FLASHDeconv "*_ms1.feature" tables; the backend sniffs the
      // header to tell them from the detector's resolved TSV, so either can be picked here.
      filters: [{ name: "Feature table", extensions: ["tsv", "txt", "feature"] }]
    });
    if (typeof picked !== "string") return;
    await loadFeaturesFromPath(picked);
  }, [loadFeaturesFromPath]);

  // ------------------------------------------------- export as ms1.feature
  // Hands the loaded features to TopPIC / mzLib / anything else that reads the TopFD
  // interchange format. Writes the FLASHDeconv dialect by default — see `exportMs1Features`
  // for why that beats leaving a TopFD Apex_intensity column blank.
  const handleExportMs1Features = useCallback(async () => {
    if (features.length === 0) return;
    // TopFD's own naming convention: "<run>_ms1.feature" beside the data.
    const stem = (metadata?.fileName ?? "features").replace(/\.[^.]+$/, "");
    const picked = await saveFileDialog({
      defaultPath: `${stem}_ms1.feature`,
      filters: [{ name: "MS1 feature table", extensions: ["feature"] }]
    });
    if (typeof picked !== "string") return;
    try {
      const rows = await exportMs1Features(picked, features, metadata?.fileName);
      const name = picked.split(/[\\/]/).pop() ?? picked;
      // Rows can exceed features when a gapped charge set is split across rows, so report
      // both rather than implying one row per feature.
      setNotice(
        `Wrote ${rows} row${rows === 1 ? "" : "s"} from ${features.length} features → ${name}`
      );
    } catch (err) {
      setLoad({ status: "error", message: errMessage(err) });
    }
  }, [features, metadata]);

  useEffect(() => {
    if (notice === null) return;
    const t = setTimeout(() => setNotice(null), 8000);
    return () => clearTimeout(t);
  }, [notice]);

  // ------------------------------------------------------------- load decoys
  // Decoy features come from a separate detector run (a decoy comb model) and share the resolved
  // TSV format, so they load through the same command — we just keep them in their own layer and
  // draw them muted, for target-vs-decoy comparison on the same raw data.
  const loadDecoysFromPath = useCallback(async (picked: string) => {
    try {
      const feats = await loadFeatures(picked);
      // Strongest first, so the rug/drawer caps keep the most important decoys.
      const sorted = [...feats].sort((a, b) => b.summedIntensity - a.summedIntensity);
      setDecoys(sorted);
      setSelectedDecoy(null);
      setDecoyDrawerOpen(true); // one right-side drawer at a time
      setDrawerOpen(false);
      setPsmDrawerOpen(false);
    } catch (err) {
      setLoad({ status: "error", message: errMessage(err) });
    }
  }, []);

  const handleLoadDecoys = useCallback(async () => {
    const picked = await openFileDialog({
      multiple: false,
      // ".feature" covers TopFD / FLASHDeconv "*_ms1.feature" tables; the backend sniffs the
      // header to tell them from the detector's resolved TSV, so either can be picked here.
      filters: [{ name: "Feature table", extensions: ["tsv", "txt", "feature"] }]
    });
    if (typeof picked !== "string") return;
    await loadDecoysFromPath(picked);
  }, [loadDecoysFromPath]);

  // E2E hook: load target/decoy feature TSVs by path, skipping the native file dialog (which
  // Playwright can't automate). Exposed only under the tauri-plugin-playwright control server
  // (`window.__PW_ACTIVE__`), mirroring `__msviewerOpenPath` — inert in normal/production runs.
  useEffect(() => {
    const w = window as typeof window & {
      __PW_ACTIVE__?: boolean;
      __msviewerLoadFeaturesPath?: (p: string) => void;
      __msviewerLoadDecoysPath?: (p: string) => void;
    };
    if (!w.__PW_ACTIVE__) return;
    w.__msviewerLoadFeaturesPath = (p: string) => void loadFeaturesFromPath(p);
    w.__msviewerLoadDecoysPath = (p: string) => void loadDecoysFromPath(p);
    return () => {
      delete w.__msviewerLoadFeaturesPath;
      delete w.__msviewerLoadDecoysPath;
    };
  }, [loadFeaturesFromPath, loadDecoysFromPath]);

  // ---------------------------------------------------------------- load PSMs
  const handleLoadPsms = useCallback(async () => {
    const picked = await openFileDialog({
      multiple: false,
      filters: [{ name: "PSM TSV", extensions: ["psmtsv", "tsv", "txt"] }]
    });
    if (typeof picked !== "string") return;
    try {
      const loaded = await loadPsms(picked);
      setPsms(loaded);
      setSelectedPsm(null);
      setPsmDrawerOpen(true);
      setDrawerOpen(false); // one right-side drawer at a time
    } catch (err) {
      setLoad({ status: "error", message: errMessage(err) });
    }
  }, []);

  // -------------------------------------------------- run detection in-app
  const handleRunDetection = useCallback(async () => {
    if (handle === null || detecting !== null || indexing) return;
    setDetecting("starting…");
    try {
      const feats = await runFeatureDetection(handle, { maxCharge: 25 }, (p) => {
        const count = p.scansTotal > 0 ? ` ${p.scansDone}/${p.scansTotal}` : "";
        setDetecting(`${p.phase}${count}`);
      });
      const sorted = [...feats].sort((a, b) => b.summedIntensity - a.summedIntensity);
      setFeatures(sorted);
      setFeaturesFile("in-app detection");
      setSelected(null);
      setDrawerOpen(true);
    } catch (err) {
      setLoad({ status: "error", message: errMessage(err) });
    } finally {
      setDetecting(null);
    }
  }, [handle, detecting, indexing]);

  // ------------------------------------------------------- select a feature
  // Shared reframe + apex-spectrum load, used for both a target and a decoy selection.
  const focusFeature = useCallback(
    async (f: Feature) => {
      const pad = Math.max(0.2, (f.rtEnd - f.rtStart) * 0.6);
      reframeTic({ xMin: f.rtStart - pad, xMax: f.rtEnd + pad });

      const z = f.primaryCharge || f.chargeStates[0] || 1;
      const grid = isotopeGrid(f.monoisotopicMass, z, 12);
      reframeSpectrum({ xMin: grid[0] - 1.5, xMax: grid[grid.length - 1] + 1.5 });

      // On-demand read by RT — available immediately (no wait on the peak index).
      if (provider && !spectrumPinned) {
        try {
          setSpectrum(await provider.getSpectrumAtRt(f.rtApex));
        } catch (err) {
          setLoad({ status: "error", message: errMessage(err) });
        }
      }
    },
    [provider, spectrumPinned, reframeTic, reframeSpectrum]
  );

  const selectFeature = useCallback(
    async (index: number) => {
      const f = features[index];
      if (!f) return;
      setSelected(index);
      setSelectedDecoy(null);
      await focusFeature(f);
    },
    [features, focusFeature]
  );

  const selectDecoy = useCallback(
    async (index: number) => {
      const f = decoys[index];
      if (!f) return;
      setSelectedDecoy(index);
      setSelected(null);
      await focusFeature(f);
    },
    [decoys, focusFeature]
  );

  // ------------------------------------------------------------- select a PSM
  // Extract the PSM's isotope XICs over the whole run (from the linked feature's observed
  // mass/charge when linked, else the PSM's theoretical precursor) — summed and overlaid on
  // the TIC — and show the identified MS2 scan in the spectrum pane.
  const selectPsm = useCallback(
    async (index: number) => {
      const psm = psms[index];
      if (!psm) return;
      setSelectedPsm(index);
      setSelected(null); // clear any feature selection (its isotope envelope is MS1-only)
      setSelectedDecoy(null);

      const link = psmLinks[index];
      const feature =
        link && link.featureIndex !== null ? features[link.featureIndex] ?? null : null;

      const mass = feature ? feature.monoisotopicMass : psm.monoisotopicMass;
      const charge = feature ? feature.primaryCharge || psm.precursorCharge : psm.precursorCharge;

      // Full-RT view so the summed XIC is visible over the entire chromatogram.
      reframeTic(createDefaultViewport());
      setXicLabel(
        `${psm.fullSequence} · z${charge} · ${mass.toFixed(2)} Da` +
          (feature ? "" : " · no feature (theoretical m/z)")
      );

      // XIC over the entire RT range (no rtRange) so it spans the whole TIC.
      if (provider && charge > 0 && mass > 0) {
        try {
          const iso = await fetchIsotopeXics(provider, mass, charge, {
            numIsotopes: 3,
            maxPoints: 4000
          });
          setXic(iso.traces);
          setXicIndices(iso.indices);
        } catch (err) {
          setLoad({ status: "error", message: errMessage(err) });
        }
      }
      // Identified MS2 spectrum by its Scan Number; fall back to the MS1 at the PSM's RT when
      // the psmtsv gave no scan number.
      if (provider && !spectrumPinned) {
        reframeSpectrum(createDefaultViewport()); // full m/z range for the fragment spectrum
        try {
          if (psm.ms2ScanNumber > 0) {
            setSpectrum(await provider.getMs2Spectrum(psm.ms2ScanNumber));
          } else if (psm.ms2RetentionTime >= 0) {
            setSpectrum(await provider.getSpectrumAtRt(psm.ms2RetentionTime));
          }
        } catch (err) {
          setLoad({ status: "error", message: errMessage(err) });
        }
      }
    },
    [psms, psmLinks, features, provider, spectrumPinned, reframeTic, reframeSpectrum]
  );

  // Click the TIC background → nearest scan's spectrum (unless the spectrum is pinned).
  const handleAreaClick = useCallback(
    async (rt: number) => {
      if (!provider || spectrumPinned) return;
      // Clicking a new RT is a fresh selection — drop any frozen y-range so the new scan auto-fits
      // (arrow-stepping keeps the frozen range; picking a new scan re-frames it).
      setSpectrumViewport((v) => ({ ...v, yMin: null, yMax: null }));
      try {
        // On-demand read by RT — works during indexing (no peak index needed).
        setSpectrum(await provider.getSpectrumAtRt(rt));
      } catch (err) {
        setLoad({ status: "error", message: errMessage(err) });
      }
    },
    [provider, spectrumPinned]
  );

  // Step `delta` MS1 scans from the currently displayed one, holding the zoom constant. Navigates
  // by retention time via the MS1 scan summaries (RT-ordered) rather than by scanIndex: the
  // displayed spectrum's `scanIndex` is a *file* index (counts MS2 scans) when loaded by RT, so it
  // can't be used directly against the MS1-only scan array. Loads the neighbour by its exact RT.
  const stepScan = useCallback(
    async (delta: number) => {
      if (!provider || indexing || !spectrum || scanSummaries.length === 0) return;
      // When zoomed into an envelope, freeze the y-range on the first step so the envelope height
      // stays stable across scans. Otherwise each scan auto-refits y to its own tallest visible peak
      // and the envelope of interest shrinks as a taller neighbor peak enters view. Only while
      // x-zoomed — at full-spectrum view the natural per-scan autoscale is kept.
      const xZoomed = spectrumViewport.xMin !== null && spectrumViewport.xMax !== null;
      if (xZoomed && !hasPersistedY(spectrumViewport)) {
        const fit = fitSpectrumYRange(spectrum.peaks, spectrumViewport.xMin, spectrumViewport.xMax);
        if (fit) setSpectrumViewport((v) => ({ ...v, yMin: fit[0], yMax: fit[1] }));
      }
      const rt = spectrum.retentionTime;
      let pos = 0;
      let best = Infinity;
      for (let i = 0; i < scanSummaries.length; i++) {
        const d = Math.abs(scanSummaries[i].retentionTime - rt);
        if (d < best) {
          best = d;
          pos = i;
        }
      }
      const next = Math.max(0, Math.min(pos + delta, scanSummaries.length - 1));
      if (next === pos) return;
      try {
        setSpectrum(await provider.getSpectrumAtRt(scanSummaries[next].retentionTime));
      } catch (err) {
        setLoad({ status: "error", message: errMessage(err) });
      }
    },
    [provider, indexing, spectrum, scanSummaries, spectrumViewport]
  );

  // ← / → step to the previous / next MS1 scan, holding the zoom constant.
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key !== "ArrowLeft" && e.key !== "ArrowRight") return;
      const t = e.target as HTMLElement | null;
      if (t && (t.tagName === "INPUT" || t.tagName === "TEXTAREA" || t.isContentEditable)) return;
      if (!spectrum || indexing) return;
      e.preventDefault();
      void stepScan(e.key === "ArrowRight" ? 1 : -1);
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [spectrum, indexing, stepScan]);

  // -------------------------------------------------- walkthrough: score a seed
  const LADDER_MAX_CHARGE = 30;

  // Score the ladder for a seed m/z (at the currently displayed scan's RT) and charge `z`. `frame`
  // controls whether the spectrum re-zooms to the anchored comb — true on a fresh seed pick, false
  // when merely re-scoring the same seed at a new anchoring charge (keep the user's zoom).
  const computeLadder = useCallback(
    async (seedMz: number, z: number, frame: boolean) => {
      if (handle === null || !spectrum) return;
      setLadderBusy(true);
      try {
        const l = await scoreSeedLadder(handle, spectrum.retentionTime, seedMz, z, LADDER_MAX_CHARGE);
        setLadder(l);
        setFocusCharge(z);
        if (frame && !spectrumPinned) {
          const focus = l.charges.find((c) => c.charge === z);
          if (focus && focus.teeth.length > 0) {
            const mzs = [l.seedMz, ...focus.teeth.map((t) => t.expectedMz)];
            const lo = Math.min(...mzs);
            const hi = Math.max(...mzs);
            const pad = Math.max(0.5, (hi - lo) * 0.12);
            // Reframe (bump uiRev) so the comb zoom is re-applied even if the user had manually zoomed —
            // consistent with focusChargeView; a fresh viewport (no persisted y) auto-fits the y-axis.
            reframeSpectrum({ xMin: lo - pad, xMax: hi + pad });
          }
        }
      } catch (err) {
        setLoad({ status: "error", message: errMessage(err) });
      } finally {
        setLadderBusy(false);
      }
    },
    [handle, spectrum, spectrumPinned, reframeSpectrum]
  );

  // Click a spectrum peak: in walkthrough, make it the ladder seed at the current zSeed (frame it);
  // otherwise select the feature that owns the peak — highlighting its row in the list and
  // emphasising it on the plots — WITHOUT reframing, so you stay on the scan you're inspecting.
  const handlePeakClick = useCallback(
    (mz: number) => {
      if (indexing) return;
      if (walkthrough) {
        void computeLadder(mz, zSeed, true);
        return;
      }
      if (!spectrum || spectrum.msLevel !== 1 || features.length === 0) return;
      const idx = findFeatureAtPeak(features, spectrum.retentionTime, mz, DETECT_PPM);
      if (idx === null) return;
      setSelected(idx);
      setSelectedDecoy(null);
      setDrawerOpen(true); // surface the list so the highlighted row is visible
      setDecoyDrawerOpen(false);
      setPsmDrawerOpen(false);
    },
    [walkthrough, indexing, zSeed, computeLadder, spectrum, features]
  );

  // Change the anchoring charge; re-score the current seed and zoom in on the new charge's comb.
  const handleZSeed = useCallback(
    (z: number) => {
      setZSeed(z);
      if (walkthrough && ladder) void computeLadder(ladder.seedMz, z, true);
    },
    [walkthrough, ladder, computeLadder]
  );

  // Click a charge row → isolate its comb overlay and recentre the spectrum on that charge's
  // teeth so you actually see the region it occupies (a higher charge sits at lower m/z).
  const focusChargeView = useCallback(
    (z: number | null) => {
      setFocusCharge(z);
      if (z === null || !ladder || spectrumPinned) return;
      const charge = ladder.charges.find((c) => c.charge === z);
      if (!charge || charge.teeth.length === 0) return;
      const mzs = charge.teeth.map((t) => t.expectedMz);
      const lo = Math.min(...mzs);
      const hi = Math.max(...mzs);
      const pad = Math.max(0.5, (hi - lo) * 0.12);
      reframeSpectrum({ xMin: lo - pad, xMax: hi + pad });
    },
    [ladder, spectrumPinned, reframeSpectrum]
  );

  // ----------------------------------------------------------------- overlays
  // TIC plus, when a PSM is selected, its isotope XICs drawn on the *same absolute* axis (true
  // abundance, not normalised). `xicMode` toggles between one summed-envelope trace and one trace
  // per isotopologue. Every XIC trace opts into `yScale`, so the y-axis zooms to the XIC height and
  // the much taller TIC runs off the top of the view.
  const ticTraces = useMemo<TicPlotTrace[]>(() => {
    const traces: TicPlotTrace[] = [
      {
        slotIndex: 0,
        points: ticPoints,
        selectedScanIndex: spectrum?.scanIndex ?? null,
        color: SLOT_COLOR
      }
    ];
    if (xic.length > 0) {
      if (xicMode === "sum") {
        traces.push({
          slotIndex: 0,
          points: sumXics(xic),
          selectedScanIndex: null,
          color: SELECTED_COLOR,
          yScale: true,
          label: "XIC (Σ isotopes)"
        });
      } else {
        xic.forEach((points, k) => {
          traces.push({
            slotIndex: 0,
            points,
            selectedScanIndex: null,
            color: chargeColor(k + 1),
            yScale: true,
            label: `M+${xicIndices[k] ?? k}`
          });
        });
      }
    }
    return traces;
  }, [ticPoints, spectrum, xic, xicIndices, xicMode]);

  // THE per-feature colour, indexed by position in `features`: cycle the palette in GLOBAL m/z
  // order, so a feature keeps one colour across scans AND across plots (the eluting SET changes
  // scan-to-scan, but a feature's m/z rank doesn't — consistent colour is the priority) while
  // features near each other in m/z still tend to differ. Consumed by the TIC rug, the spectrum
  // overlay, and the drawer swatch; the selected feature overrides to SELECTED_COLOR in all three.
  // Declared before `featureRug` because that memo reads it during render.
  const featureColorByIndex = useMemo<string[]>(() => {
    const order = features.map((f, i) => ({ i, mz: f.detectedMz }));
    order.sort((a, b) => a.mz - b.mz);
    const colors = new Array<string>(features.length);
    order.forEach((o, rank) => {
      colors[o.i] = FEATURE_COLORS[rank % FEATURE_COLORS.length];
    });
    return colors;
  }, [features]);
  const featureColor = useCallback(
    (i: number): string => featureColorByIndex[i] ?? FEATURE_COLORS[0],
    [featureColorByIndex]
  );
  // Decoys are deliberately one flat muted colour everywhere (rug, spectrum, drawer swatch) so the
  // real features' hues stay meaningful — this keeps the drawer honest about that.
  const decoyColor = useCallback((): string => DECOY_COLOR, []);

  const featureRug = useMemo<FeatureMarker[]>(() => {
    const targets = features.slice(0, RUG_CAP).map((f, i) => ({
      featureIndex: i,
      retentionTime: f.rtApex,
      color: i === selected ? SELECTED_COLOR : featureColor(i),
      label: `m/z ${f.detectedMz.toFixed(3)} · z${f.primaryCharge} · ${f.monoisotopicMass.toFixed(1)} Da · RT ${f.rtApex.toFixed(2)}`
    }));
    const decoyMarkers = decoys.slice(0, RUG_CAP).map((f, i) => ({
      featureIndex: DECOY_INDEX_BASE + i,
      retentionTime: f.rtApex,
      color: i === selectedDecoy ? SELECTED_COLOR : DECOY_COLOR,
      label: `decoy · m/z ${f.detectedMz.toFixed(3)} · z${f.primaryCharge} · ${f.monoisotopicMass.toFixed(1)} Da · RT ${f.rtApex.toFixed(2)}`
    }));
    // Decoys first so target markers draw on top of them within the single rug trace.
    return [...decoyMarkers, ...targets];
  }, [features, decoys, selected, selectedDecoy, featureColor]);

  const regions = useMemo<RtRegion[]>(() => {
    if (selected !== null) {
      const f = features[selected];
      if (f) return [{ min: f.rtStart, max: f.rtEnd, color: "rgba(232,131,12,0.14)" }];
    }
    if (selectedDecoy !== null) {
      const f = decoys[selectedDecoy];
      if (f) return [{ min: f.rtStart, max: f.rtEnd, color: DECOY_REGION_FILL }];
    }
    return [];
  }, [features, decoys, selected, selectedDecoy]);

  // Sorted m/z of the currently displayed scan's peaks — the lookup for the "detected" comb marker.
  const spectrumMz = useMemo(
    () => (spectrum ? spectrum.peaks.map((p) => p.mz) : []),
    [spectrum]
  );

  // Walkthrough comb overlay: the focused charge's teeth (or all retained charges if none
  // focused), coloured by charge. Each tooth is marked `detected` when a peak sits at its m/z in
  // the *displayed* scan, so the indicator updates as you arrow through scans. Overrides the
  // feature envelope while walkthrough is active.
  const ladderEnvelope = useMemo<EnvelopeLine[]>(() => {
    if (!walkthrough || !ladder) return [];
    const charges =
      focusCharge !== null
        ? ladder.charges.filter((c) => c.charge === focusCharge)
        : ladder.charges.filter((c) => c.retained);
    return charges.flatMap((c) =>
      c.teeth.map((t) => ({
        mz: t.expectedMz,
        color: chargeColor(c.charge),
        label: t.isotopeIndex === 0 ? `z${c.charge}` : undefined,
        detected: hasPeakNear(spectrumMz, t.expectedMz, DETECT_PPM)
      }))
    );
  }, [walkthrough, ladder, focusCharge, spectrumMz]);

  // Every (target) feature eluting at the displayed MS1 scan's RT — the FULL loaded set, not the
  // 800-row drawer cap. Each entry keeps the feature's global index (`gi`) so the overlay can look
  // up its stable colour. Intensity-sorted (features are), so the strongest survive the cap below.
  const scanEluting = useMemo<readonly { f: Feature; gi: number }[]>(() => {
    if (!showScanFeatures || !spectrum || spectrum.msLevel !== 1) return [];
    const out: { f: Feature; gi: number }[] = [];
    features.forEach((f, gi) => {
      if (featureElutesAt(f, spectrum.retentionTime)) out.push({ f, gi });
    });
    return out;
  }, [showScanFeatures, spectrum, features]);

  // Feature overlays drawn as markers sitting on the MATCHED peaks (peak apex, in the feature's
  // colour), instead of full-height comb lines. Scan mode: every eluting feature, stable colour by
  // m/z (the selected one emphasised in SELECTED_COLOR). Otherwise the single selected feature
  // (per-charge colour) or a selected decoy (muted). Empty during walkthrough (its comb takes over)
  // or when the displayed spectrum isn't MS1.
  const peakHighlights = useMemo<SpectrumPeakHighlight[]>(() => {
    if (walkthrough || !spectrum || spectrum.msLevel !== 1) return [];
    const peaks = spectrum.peaks;
    if (showScanFeatures) {
      // Re-iterate on zoom: when zoomed in, only features with a peak in the visible m/z window are
      // candidates, so weak features in that window get shown (and far fewer compete for the cap).
      const inWindow = spectrumXView
        ? scanEluting.filter(({ f }) =>
            featureHasPeakInWindow(f, spectrumXView.min, spectrumXView.max)
          )
        : scanEluting;
      const items = inWindow.slice(0, SCAN_FEATURE_CAP);
      // Always include the selected feature, even if it fell beyond the cap or the window — otherwise
      // clicking a weak feature in the list wouldn't highlight it.
      if (selected !== null && !items.some((it) => it.gi === selected)) {
        const f = features[selected];
        if (f) items.push({ f, gi: selected });
      }
      return items.flatMap(({ f, gi }) => {
        const color = gi === selected ? SELECTED_COLOR : featureColor(gi);
        return matchedPeakHighlights(f, peaks, () => color, 6);
      });
    }
    if (selected !== null) {
      const f = features[selected];
      // The selection is SELECTED_COLOR here for the same reason it is on the rug and in the drawer:
      // one feature, one colour, in every view it appears in.
      return f ? matchedPeakHighlights(f, peaks, () => SELECTED_COLOR, 12) : [];
    }
    if (selectedDecoy !== null) {
      const f = decoys[selectedDecoy];
      return f ? matchedPeakHighlights(f, peaks, () => DECOY_COLOR, 12) : [];
    }
    return [];
  }, [
    walkthrough,
    spectrum,
    showScanFeatures,
    scanEluting,
    spectrumXView,
    selected,
    selectedDecoy,
    features,
    decoys,
    featureColor
  ]);

  // Full-height comb lines are now ONLY the walkthrough charge-ladder (a predicted-position
  // diagnostic that must show even the missing teeth); feature overlays use `peakHighlights`.
  const spectrumEnvelope = walkthrough ? ladderEnvelope : [];

  const spectrumTraces = useMemo<SpectrumPlotTrace[]>(
    () => (spectrum ? [{ slotIndex: 0, peaks: spectrum.peaks, color: SLOT_COLOR }] : []),
    [spectrum]
  );

  // MS1 peak labels: annotate the most prominent peaks in the current x-window with m/z + inferred
  // charge. Only for MS1 spectra, and suppressed during the walkthrough (its comb overlay already
  // annotates the ladder, so peak labels would just clutter it).
  const spectrumAnnotations = useMemo<PeakAnnotation[]>(() => {
    if (!spectrum || spectrum.msLevel !== 1 || walkthrough) return [];
    return computePeakLabels(spectrum.peaks, {
      xMin: spectrumViewport.xMin,
      xMax: spectrumViewport.xMax,
      maxLabels: 6
    });
  }, [spectrum, spectrumViewport, walkthrough]);

  const selectedFeature =
    selected !== null
      ? features[selected] ?? null
      : selectedDecoy !== null
        ? decoys[selectedDecoy] ?? null
        : null;

  return (
    <>
    {dragOver ? <DropOverlay /> : null}
    <ViewerShell
      title="MsViewer — feature finder"
      subtitle="Real raw/mzML over Tauri IPC, with top-down feature-finding results overlaid."
      rightInset={drawerVisible ? drawerWidth + DRAWER_INSET_GAP : 0}
      toolbar={
        <>
          {metadata ? (
            <>
              <MetricReadout label="File" value={metadata.fileName} />
              <MetricReadout
                label="MS1 scans"
                value={indexing ? "indexing…" : metadata.ms1ScanCount}
              />
              <MetricReadout label="Format" value={metadata.format} />
            </>
          ) : null}
          <PanelActionButton onClick={() => void handleOpenFile()}>Open file…</PanelActionButton>
          <PanelActionButton onClick={() => void handleLoadFeatures()}>
            Load features…
          </PanelActionButton>
          <PanelActionButton onClick={() => void handleLoadDecoys()}>
            Load decoys…
          </PanelActionButton>
          <PanelActionButton onClick={() => void handleLoadPsms()}>Load PSMs…</PanelActionButton>
          {features.length > 0 ? (
            <PanelActionButton onClick={() => void handleExportMs1Features()}>
              Export ms1.feature…
            </PanelActionButton>
          ) : null}
          {notice ? <StatusBanner tone="info">{notice}</StatusBanner> : null}
          {handle !== null ? (
            <PanelActionButton onClick={() => void handleRunDetection()}>
              {indexing
                ? "Indexing…"
                : detecting
                  ? `Detecting: ${detecting}`
                  : "Run feature finding"}
            </PanelActionButton>
          ) : null}
          {handle !== null ? (
            <PanelActionButton
              pressed={walkthrough}
              onClick={() => {
                const next = !walkthrough;
                setWalkthrough(next);
                if (!next) {
                  setLadder(null);
                  setFocusCharge(null);
                }
              }}
            >
              {walkthrough ? "Walkthrough: on" : "Walkthrough"}
            </PanelActionButton>
          ) : null}
          {features.length > 0 ? (
            <PanelActionButton
              pressed={drawerOpen}
              onClick={() =>
                setDrawerOpen((v) => {
                  const next = !v;
                  if (next) {
                    setPsmDrawerOpen(false);
                    setDecoyDrawerOpen(false);
                  }
                  return next;
                })
              }
            >
              Features ({features.length})
            </PanelActionButton>
          ) : null}
          {decoys.length > 0 ? (
            <PanelActionButton
              pressed={decoyDrawerOpen}
              onClick={() =>
                setDecoyDrawerOpen((v) => {
                  const next = !v;
                  if (next) {
                    setDrawerOpen(false);
                    setPsmDrawerOpen(false);
                  }
                  return next;
                })
              }
            >
              Decoys ({decoys.length})
            </PanelActionButton>
          ) : null}
          {psms.length > 0 ? (
            <PanelActionButton
              pressed={psmDrawerOpen}
              onClick={() =>
                setPsmDrawerOpen((v) => {
                  const next = !v;
                  if (next) {
                    setDrawerOpen(false);
                    setDecoyDrawerOpen(false);
                  }
                  return next;
                })
              }
            >
              PSMs ({psms.length})
            </PanelActionButton>
          ) : null}
          {features.length > 0 ? (
            <PanelActionButton
              pressed={!showScanFeatures}
              onClick={() => setShowScanFeatures((v) => !v)}
            >
              Selected only
            </PanelActionButton>
          ) : null}
        </>
      }
    >
      <Panel
        active={ticPinned}
        header={
          <PanelHeader
            title="Total ion chromatogram"
            subtitle={
              xicLabel
                ? `${xicMode === "sum" ? "XIC Σ isotopes" : `XIC isotopologues M+${xicIndices.join(", M+")}`} (absolute abundance): ${xicLabel}`
                : selectedFeature
                  ? `Feature: ${selectedFeature.monoisotopicMass.toFixed(2)} Da · z ${selectedFeature.chargeStates.join(",")} · RT ${selectedFeature.rtStart.toFixed(2)}–${selectedFeature.rtEnd.toFixed(2)}${scoreReadout(selectedFeature)}`
                  : featuresFile
                    ? `${features.length} features from ${featuresFile} — click a marker or a row`
                    : "MS1 TIC · click to load a scan"
            }
            actions={
              <>
                {xic.length > 0 ? (
                  <>
                    <PanelActionButton
                      pressed={xicMode === "sum"}
                      onClick={() => setXicMode("sum")}
                    >
                      XIC: Σ
                    </PanelActionButton>
                    <PanelActionButton
                      pressed={xicMode === "isotopes"}
                      onClick={() => setXicMode("isotopes")}
                    >
                      XIC: isotopes
                    </PanelActionButton>
                  </>
                ) : null}
                {ticPoints.length > 0 ? (
                  <PanelActionButton onClick={() => reframeTic(createDefaultViewport())}>
                    Reset zoom
                  </PanelActionButton>
                ) : null}
                <PanelActionButton
                  pressed={ticPinned}
                  onClick={() => {
                    const next = !ticPinned;
                    setTicPinned(next);
                    if (next) setSpectrumPinned(false);
                  }}
                >
                  {ticPinned ? "Pinned" : "Pin"}
                </PanelActionButton>
              </>
            }
          />
        }
      >
        {load.status === "loading" ? (
          <StatusBanner tone="info">{load.message}</StatusBanner>
        ) : load.status === "error" ? (
          <StatusBanner tone="error">Failed: {load.message}</StatusBanner>
        ) : ticPoints.length === 0 ? (
          indexing ? (
            <StatusBanner tone="info">
              Building peak index — the MS1 chromatogram will appear when it’s ready…
            </StatusBanner>
          ) : (
            <StatusBanner tone="muted">
              No dataset — use “Open file…” to load a .raw or mzML.
            </StatusBanner>
          )
        ) : (
          <TicPlot
            traces={ticTraces}
            viewport={ticViewport}
            featureRug={featureRug}
            regions={regions}
            uirevision={ticUiRev}
            rangeSelectionEnabled={false}
            onEvent={(e) => {
              if (e.type === "area-click") void handleAreaClick(e.retentionTime);
              else if (e.type === "feature-click") {
                if (e.featureIndex >= DECOY_INDEX_BASE)
                  void selectDecoy(e.featureIndex - DECOY_INDEX_BASE);
                else void selectFeature(e.featureIndex);
              }
            }}
          />
        )}
      </Panel>

      <Panel
        active={spectrumPinned}
        header={
          <PanelHeader
            title={spectrum ? `Spectrum · MS${spectrum.msLevel}` : "Spectrum"}
            subtitle={
              spectrum
                ? `Scan ${spectrum.oneBasedScanNumber} · RT ${spectrum.retentionTime.toFixed(2)} min${
                    walkthrough
                      ? " · click a peak to seed the charge-ladder"
                      : showScanFeatures
                        ? ` · ${scanEluting.length} feature${scanEluting.length === 1 ? "" : "s"} eluting${scanEluting.length > SCAN_FEATURE_CAP ? ` (showing strongest ${SCAN_FEATURE_CAP})` : ""}`
                        : selectedFeature
                          ? " · shaded columns mark matched isotope peaks"
                          : ""
                  }`
                : walkthrough
                  ? "Click the chromatogram to load a scan, then click a peak to seed"
                  : "Click the chromatogram or a feature to load a scan"
            }
            actions={
              <>
                {spectrumTraces.length > 0 ? (
                  <PanelActionButton onClick={() => reframeSpectrum(createDefaultViewport())}>
                    Reset zoom
                  </PanelActionButton>
                ) : null}
                <PanelActionButton
                  pressed={spectrumPinned}
                  onClick={() => {
                    const next = !spectrumPinned;
                    setSpectrumPinned(next);
                    if (next) setTicPinned(false);
                  }}
                >
                  {spectrumPinned ? "Pinned" : "Pin"}
                </PanelActionButton>
              </>
            }
          />
        }
      >
        {spectrumTraces.length === 0 ? (
          <StatusBanner tone="muted">No spectrum selected.</StatusBanner>
        ) : (
          <SpectrumPlot
            traces={spectrumTraces}
            viewport={spectrumViewport}
            envelope={spectrumEnvelope}
            highlights={peakHighlights}
            annotations={spectrumAnnotations}
            uirevision={spectrumUiRev}
            rangeSelectionEnabled={false}
            onEvent={(e) => {
              if (e.type === "peak-click") {
                handlePeakClick(e.peak.mz);
              } else if (e.type === "xrange-change") {
                const range = e.range;
                if (range === null) {
                  // Double-click autorange reset: clear the persisted x-zoom and re-fit cleanly
                  // (bumps uiRev so Plotly re-applies the full range).
                  reframeSpectrum(createDefaultViewport());
                } else {
                  // Capture the user's interactive zoom/pan into the viewport — the single source of
                  // truth that drives the `range` prop — WITHOUT bumping uiRev. This is what makes a
                  // manual zoom survive a re-render (new-spectrum select, scan step): the range prop
                  // re-applies it deterministically instead of relying on Plotly's `uirevision`
                  // preserving internal state across data changes (which it does not reliably do).
                  setSpectrumXView(range);
                  setSpectrumViewport((v) =>
                    v.xMin === range.min && v.xMax === range.max
                      ? v
                      : { ...v, xMin: range.min, xMax: range.max }
                  );
                }
              }
            }}
          />
        )}
      </Panel>

      {walkthrough ? (
        <LadderDrawer
          ladder={ladder}
          zSeed={zSeed}
          maxCharge={LADDER_MAX_CHARGE}
          focusCharge={focusCharge}
          busy={ladderBusy}
          width={drawerWidth}
          onResize={setDrawerWidth}
          onZSeed={handleZSeed}
          onFocusCharge={focusChargeView}
          onClose={() => {
            setWalkthrough(false);
            setLadder(null);
            setFocusCharge(null);
          }}
        />
      ) : psmDrawerOpen && psms.length > 0 ? (
        <PsmDrawer
          psms={psms}
          links={psmLinks}
          selected={selectedPsm}
          onSelect={(i) => void selectPsm(i)}
          onClose={() => setPsmDrawerOpen(false)}
        />
      ) : drawerOpen && features.length > 0 ? (
        <FeatureDrawer
          features={features}
          selected={selected}
          width={drawerWidth}
          onResize={setDrawerWidth}
          onSelect={(i) => void selectFeature(i)}
          onClose={() => setDrawerOpen(false)}
          colorOf={featureColor}
        />
      ) : decoyDrawerOpen && decoys.length > 0 ? (
        <FeatureDrawer
          features={decoys}
          selected={selectedDecoy}
          width={drawerWidth}
          onResize={setDrawerWidth}
          onSelect={(i) => void selectDecoy(i)}
          onClose={() => setDecoyDrawerOpen(false)}
          colorOf={decoyColor}
          title="Decoys"
        />
      ) : null}
    </ViewerShell>
    </>
  );
}

// Full-window overlay shown while a file is dragged over the app, cueing the user that a
// drop will open the dataset. Purely visual — the real work happens in the webview
// drag-drop listener; `pointerEvents: none` keeps it from interfering with the OS drop.
function DropOverlay() {
  return (
    <div style={dropOverlayStyle}>
      <div style={dropOverlayCardStyle}>Drop a .raw or mzML file to open it</div>
    </div>
  );
}

// A right-side overlay drawer with a draggable left-edge handle, so the user can resize it. The
// features / walkthrough drawers (and the future PSM panel) render their content inside this frame.
// Dragging updates the width via `onResize`; the parent keeps `ViewerShell`'s `rightInset` in sync so
// the plots reflow to the new width.
function ResizableDrawer({
  width,
  onResize,
  children
}: {
  width: number;
  onResize: (width: number) => void;
  children: React.ReactNode;
}) {
  const startDrag = useCallback(
    (e: React.PointerEvent) => {
      e.preventDefault();
      // Drawer is right-anchored, so width = viewport width − pointer x.
      const onMove = (ev: PointerEvent) =>
        onResize(clampDrawerWidth(window.innerWidth - ev.clientX));
      const onUp = () => {
        window.removeEventListener("pointermove", onMove);
        window.removeEventListener("pointerup", onUp);
      };
      window.addEventListener("pointermove", onMove);
      window.addEventListener("pointerup", onUp);
    },
    [onResize]
  );
  return (
    <div style={{ ...drawerStyle, width }} data-testid="drawer">
      <div
        role="separator"
        aria-orientation="vertical"
        title="Drag to resize"
        onPointerDown={startDrag}
        style={drawerResizeHandleStyle}
        data-testid="drawer-resize"
      />
      {children}
    </div>
  );
}

// The walkthrough drawer: pick the anchoring charge zSeed (1…maxCharge), see the
// candidate mass it anchors, and step the resulting charge-state ladder — each
// charge's comb response, teeth count, and whether it's retained. Click a charge
// row to isolate its comb in the spectrum overlay.
function LadderDrawer({
  ladder,
  zSeed,
  maxCharge,
  focusCharge,
  busy,
  width,
  onResize,
  onZSeed,
  onFocusCharge,
  onClose
}: {
  ladder: SeedLadder | null;
  zSeed: number;
  maxCharge: number;
  focusCharge: number | null;
  busy: boolean;
  width: number;
  onResize: (width: number) => void;
  onZSeed: (z: number) => void;
  onFocusCharge: (z: number | null) => void;
  onClose: () => void;
}) {
  const chargeButtons = Array.from({ length: maxCharge }, (_, i) => i + 1);
  return (
    <ResizableDrawer width={width} onResize={onResize}>
      <div style={drawerHeaderStyle}>
        <strong>Charge-ladder walkthrough</strong>
        <button onClick={onClose} style={drawerCloseStyle} type="button">
          ✕
        </button>
      </div>
      <div style={{ ...drawerBodyStyle, padding: "10px 12px" }}>
        <div style={ladderHintStyle}>
          Anchoring charge z<sub>seed</sub> — the charge the seed is assumed to be. It sets the
          candidate mass; the ladder below then scores that mass at every charge.
        </div>
        <div style={chargeGridStyle}>
          {chargeButtons.map((z) => (
            <button
              key={z}
              type="button"
              onClick={() => onZSeed(z)}
              style={z === zSeed ? zSeedButtonSelStyle : zSeedButtonStyle}
            >
              {z}
            </button>
          ))}
        </div>

        {busy ? (
          <div style={ladderNoteStyle}>Scoring…</div>
        ) : !ladder ? (
          <div style={ladderNoteStyle}>
            Load a scan, then click a spectrum peak to seed the ladder.
          </div>
        ) : (
          <>
            <table style={summaryTableStyle}>
              <tbody>
                <tr>
                  <td style={sumKeyStyle}>Seed</td>
                  <td style={sumValStyle}>
                    m/z {ladder.seedMz.toFixed(4)} · scan {ladder.seedScanIndex} · RT{" "}
                    {ladder.seedRt.toFixed(2)}
                  </td>
                </tr>
                <tr>
                  <td style={sumKeyStyle}>Apex screen</td>
                  <td style={sumValStyle}>
                    {ladder.screenTeeth}/6 teeth ·{" "}
                    <span style={{ color: ladder.screenPassed ? "#1e7e34" : "#b02a37" }}>
                      {ladder.screenPassed ? "passes" : "fails"}
                    </span>
                  </td>
                </tr>
                <tr>
                  <td style={sumKeyStyle}>Anchored mass</td>
                  <td style={sumValStyle}>
                    {ladder.monoMass.toFixed(3)} Da · mono m/z {ladder.monoMz.toFixed(4)} (i*=
                    {ladder.iStar}){ladder.massInRange ? "" : " · out of range"}
                  </td>
                </tr>
                <tr>
                  <td style={sumKeyStyle}>Refined mass</td>
                  <td style={sumValStyle}>
                    {ladder.refinedMonoMass.toFixed(3)} Da
                    {Math.abs(ladder.refinedMonoMass - ladder.monoMass) > 0.5
                      ? ` (Δ ${(ladder.refinedMonoMass - ladder.monoMass).toFixed(2)})`
                      : ""}
                  </td>
                </tr>
                <tr>
                  <td style={sumKeyStyle}>Ladder</td>
                  <td style={sumValStyle}>
                    {ladder.numChargeStates} charge{ladder.numChargeStates === 1 ? "" : "s"} · Σ
                    resp {ladder.totalResponse.toExponential(2)} ·{" "}
                    <span style={{ color: ladder.accepted ? "#1e7e34" : "#b02a37" }}>
                      {ladder.accepted ? "ACCEPT" : "reject"}
                    </span>
                  </td>
                </tr>
                <tr>
                  <td style={sumKeyStyle}>RT window</td>
                  <td style={sumValStyle}>{ladder.windowScanCount} scans</td>
                </tr>
              </tbody>
            </table>

            <div style={ladderTableWrapStyle}>
              <table style={tableStyle}>
                <thead>
                  <tr>
                    <th style={thStyle}>z</th>
                    <th style={thStyle}>mono m/z</th>
                    <th style={thStyleRight}>teeth</th>
                    <th style={thStyleRight}>response</th>
                  </tr>
                </thead>
                <tbody>
                  <tr
                    onClick={() => onFocusCharge(null)}
                    style={focusCharge === null ? rowSelectedStyle : rowStyle}
                  >
                    <td style={tdStyle} colSpan={4}>
                      Show all retained charges
                    </td>
                  </tr>
                  {ladder.charges.map((c) => (
                    <LadderRow
                      key={c.charge}
                      charge={c}
                      isSeed={c.charge === ladder.zSeed}
                      focused={focusCharge === c.charge}
                      onClick={() => onFocusCharge(c.charge)}
                    />
                  ))}
                </tbody>
              </table>
            </div>
          </>
        )}
      </div>
    </ResizableDrawer>
  );
}

function LadderRow({
  charge,
  isSeed,
  focused,
  onClick
}: {
  charge: LadderCharge;
  isSeed: boolean;
  focused: boolean;
  onClick: () => void;
}) {
  const dim = !charge.retained;
  const style: React.CSSProperties = {
    ...(focused ? rowSelectedStyle : rowStyle),
    color: dim ? "#9aa7b6" : "#24364d",
    fontWeight: charge.retained ? 500 : 400
  };
  return (
    <tr onClick={onClick} style={style}>
      <td style={tdStyle}>
        <span style={{ color: chargeColor(charge.charge) }}>■</span> z{charge.charge}
        {isSeed ? " ◄" : ""}
      </td>
      <td style={tdStyle}>{charge.monoMz.toFixed(3)}</td>
      <td style={tdStyleRight}>{charge.numIsotopesObserved}</td>
      <td style={tdStyleRight}>
        {charge.response > 0 ? charge.response.toExponential(1) : "—"}
      </td>
    </tr>
  );
}

function errMessage(err: unknown): string {
  if (err && typeof err === "object" && "message" in err) {
    return String((err as { message: unknown }).message);
  }
  return err instanceof Error ? err.message : String(err);
}

// True if any m/z in the ascending-sorted `sortedMz` lies within `ppm` of `mz`. Binary search for
// the first candidate ≥ (mz − tol), then check it's ≤ (mz + tol).
function hasPeakNear(sortedMz: readonly number[], mz: number, ppm: number): boolean {
  const n = sortedMz.length;
  if (n === 0) return false;
  const tol = (mz * ppm) / 1e6;
  const loTarget = mz - tol;
  let lo = 0;
  let hi = n;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (sortedMz[mid] < loTarget) lo = mid + 1;
    else hi = mid;
  }
  return lo < n && sortedMz[lo] <= mz + tol;
}

// The observed peak within `ppm` of `mz` (lower-bound candidate), or null. Mirrors `hasPeakNear`
// but returns the peak, so an overlay can sit at its true (m/z, intensity). `peaks` must be m/z-sorted.
function nearestPeak(
  peaks: readonly { mz: number; intensity: number }[],
  mz: number,
  ppm: number
): { mz: number; intensity: number } | null {
  const n = peaks.length;
  if (n === 0) return null;
  const tol = (mz * ppm) / 1e6;
  const loTarget = mz - tol;
  let lo = 0;
  let hi = n;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (peaks[mid].mz < loTarget) lo = mid + 1;
    else hi = mid;
  }
  return lo < n && peaks[lo].mz <= mz + tol ? peaks[lo] : null;
}

// The feature that owns a clicked peak in the current scan: among features eluting at `rt`, the one
// whose predicted isotope comb has a tooth closest (in ppm, within tolerance) to `mz`. Returns its
// index in `features`, or null. Lets a peak click select/highlight its feature.
function findFeatureAtPeak(
  features: readonly Feature[],
  rt: number,
  mz: number,
  ppm: number
): number | null {
  let best: number | null = null;
  let bestErr = ppm;
  for (let i = 0; i < features.length; i++) {
    const f = features[i];
    if (!featureElutesAt(f, rt)) continue;
    for (const z of f.chargeStates) {
      for (const tooth of isotopeGrid(f.monoisotopicMass, z, 12)) {
        const err = (Math.abs(tooth - mz) / mz) * 1e6;
        if (err <= bestErr) {
          bestErr = err;
          best = i;
        }
      }
    }
  }
  return best;
}

// True if any predicted isotope tooth of the feature falls within the visible m/z window [lo, hi].
// Used to filter the scan overlay to features visible at the current zoom.
function featureHasPeakInWindow(f: Feature, lo: number, hi: number): boolean {
  for (const z of f.chargeStates) {
    for (const mz of isotopeGrid(f.monoisotopicMass, z, 8)) {
      if (mz >= lo && mz <= hi) return true;
    }
  }
  return false;
}

// Feature-membership highlights: for each predicted isotope tooth that matches an observed peak
// (within DETECT_PPM), a marker at that peak's true (m/z, intensity), coloured by `colorFor(charge)`.
// Predicted-but-absent teeth produce nothing — the overlay marks only real peaks.
function matchedPeakHighlights(
  f: Feature,
  peaks: readonly { mz: number; intensity: number }[],
  colorFor: (z: number) => string,
  count: number
): SpectrumPeakHighlight[] {
  return f.chargeStates.flatMap((z) => {
    const color = colorFor(z);
    return isotopeGrid(f.monoisotopicMass, z, count).flatMap((mz) => {
      const pk = nearestPeak(peaks, mz, DETECT_PPM);
      return pk ? [{ mz: pk.mz, intensity: pk.intensity, color }] : [];
    });
  });
}

// Sortable columns of the feature drawer.
type FeatureSortKey =
  | "mass" | "z" | "mz" | "rt" | "intensity"
  | "decon" | "minDecon" | "maxIso" | "ppm" | "corrAll" | "corr5" | "corr3";

// Numeric sort value for a feature under `key`; NaN for a missing score (sorted last either way).
function featureSortValue(f: Feature, key: FeatureSortKey): number {
  switch (key) {
    case "mass": return f.monoisotopicMass;
    case "z": return f.primaryCharge;
    case "mz": return f.detectedMz;
    case "rt": return f.rtApex;
    case "intensity": return f.summedIntensity;
    case "decon": return f.deconScore ?? NaN;
    case "minDecon": return f.minDeconScore ?? NaN;
    case "maxIso": return f.maxNumIsotopes ?? NaN;
    case "ppm": return f.ppmSpread ?? NaN;
    case "corrAll": return f.isoCorrAll ?? NaN;
    case "corr5": return f.isoCorrTop5 ?? NaN;
    case "corr3": return f.isoCorrTop3 ?? NaN;
  }
}

// A right-side drawer listing resolved features; click a row to select. Sortable by any column and
// paginated; the page auto-follows the current selection.
function FeatureDrawer({
  features,
  selected,
  width,
  onResize,
  onSelect,
  onClose,
  colorOf,
  title = "Features"
}: {
  features: readonly Feature[];
  selected: number | null;
  width: number;
  onResize: (width: number) => void;
  onSelect: (index: number) => void;
  onClose: () => void;
  // The row's swatch colour, by index into `features` — the SAME lookup the rug and the spectrum
  // overlay use, so the drawer reads as a legend for the plots rather than a second scheme.
  colorOf: (index: number) => string;
  title?: string;
}) {
  const [sortKey, setSortKey] = useState<FeatureSortKey>("intensity");
  const [sortDir, setSortDir] = useState<"asc" | "desc">("desc");
  const [filterCharge, setFilterCharge] = useState("");
  const [filterMinIntensity, setFilterMinIntensity] = useState("");
  // Indices into `features`, filtered (charge / min intensity) then sorted by the chosen column
  // (missing scores last). We keep INDICES (not rows) so onSelect / data-row / `selected` stay in the
  // parent's feature-index space.
  const sortedIndices = useMemo(() => {
    const zf = parseInt(filterCharge, 10);
    const minI = parseFloat(filterMinIntensity);
    const hasZ = Number.isFinite(zf);
    const hasI = Number.isFinite(minI);
    const idx: number[] = [];
    for (let i = 0; i < features.length; i++) {
      const f = features[i];
      if (hasZ && !f.chargeStates.includes(zf)) continue;
      if (hasI && f.summedIntensity < minI) continue;
      idx.push(i);
    }
    const dir = sortDir === "asc" ? 1 : -1;
    idx.sort((a, b) => {
      const va = featureSortValue(features[a], sortKey);
      const vb = featureSortValue(features[b], sortKey);
      const na = Number.isNaN(va);
      const nb = Number.isNaN(vb);
      if (na && nb) return 0;
      if (na) return 1;
      if (nb) return -1;
      return va === vb ? 0 : va < vb ? -dir : dir;
    });
    return idx;
  }, [features, sortKey, sortDir, filterCharge, filterMinIntensity]);
  const total = sortedIndices.length;
  const positionOf = useMemo(() => {
    const m = new Map<number, number>();
    sortedIndices.forEach((gi, p) => m.set(gi, p));
    return m;
  }, [sortedIndices]);

  // Rows per page auto-fit the drawer height (measure effect below), in steps of 5. Any feature is
  // still reachable via the pager, which can auto-follow a selection made elsewhere (peak / rug click).
  const [pageSize, setPageSize] = useState(DRAWER_PAGE);
  const pageCount = Math.max(1, Math.ceil(total / pageSize));
  const [page, setPage] = useState(0);
  const clampedPage = Math.min(page, pageCount - 1);
  const start = clampedPage * pageSize;
  const shown = sortedIndices.slice(start, start + pageSize);
  // Score columns show only when the loaded TSV carried them; stable across pages.
  const scored = features.some(hasScores);

  // Reset to the first page on a new feature set or a re-sort. Declared BEFORE the follow effect so
  // that when something is selected, the follow effect's page wins.
  useEffect(() => {
    setPage(0);
  }, [features, sortKey, sortDir, filterCharge, filterMinIntensity]);
  // Jump to the page holding the selection — a peak / rug click (or a re-sort / resize) can move it.
  useEffect(() => {
    if (selected === null) return;
    const p = positionOf.get(selected);
    if (p !== undefined) setPage(Math.floor(p / pageSize));
  }, [selected, positionOf, pageSize]);
  // Snap the highlighted row into view once its page has rendered.
  const bodyRef = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (selected === null) return;
    bodyRef.current
      ?.querySelector(`tr[data-row="${selected}"]`)
      ?.scrollIntoView({ block: "nearest" });
  }, [selected, clampedPage]);
  // Auto-fit rows-per-page to the drawer's height, in steps of 5 (re-measured on resize).
  useEffect(() => {
    const body = bodyRef.current;
    if (!body) return;
    const measure = () => {
      const row = body.querySelector("tbody tr");
      const head = body.querySelector("thead");
      const rowH = row ? row.getBoundingClientRect().height : 0;
      const headH = head ? head.getBoundingClientRect().height : 0;
      const avail = body.clientHeight - headH;
      if (rowH <= 0 || avail <= 0) return;
      const fit = Math.max(5, Math.floor(avail / rowH / 5) * 5);
      setPageSize((prev) => (prev === fit ? prev : fit));
    };
    measure();
    if (typeof ResizeObserver === "undefined") return;
    const ro = new ResizeObserver(measure);
    ro.observe(body);
    return () => ro.disconnect();
  }, []);

  // A clickable, sort-toggling column header (▲/▼ marks the active sort).
  const sortTh = (key: FeatureSortKey, label: string, right = false) => (
    <th
      style={{ ...(right ? thStyleRight : thStyle), cursor: "pointer", userSelect: "none" }}
      onClick={() => {
        if (sortKey === key) setSortDir((d) => (d === "asc" ? "desc" : "asc"));
        else {
          setSortKey(key);
          setSortDir("desc");
        }
      }}
      title={`Sort by ${label}`}
    >
      {label}
      {sortKey === key ? (sortDir === "asc" ? " ▲" : " ▼") : ""}
    </th>
  );

  const first = total === 0 ? 0 : start + 1;
  const last = Math.min(start + pageSize, total);
  const filtered = filterCharge !== "" || filterMinIntensity !== "";
  return (
    <ResizableDrawer width={width} onResize={onResize}>
      <div style={drawerHeaderStyle}>
        <strong>
          {title} (
          {filtered
            ? `${total.toLocaleString()} of ${features.length.toLocaleString()}`
            : total.toLocaleString()}
          )
        </strong>
        <button onClick={onClose} style={drawerCloseStyle} type="button">
          ✕
        </button>
      </div>
      <div style={drawerFilterStyle}>
        <label style={filterLabelStyle}>
          z
          <input
            type="number"
            value={filterCharge}
            onChange={(e) => setFilterCharge(e.target.value)}
            placeholder="any"
            style={filterInputStyle}
          />
        </label>
        <label style={filterLabelStyle}>
          min int.
          <input
            type="number"
            value={filterMinIntensity}
            onChange={(e) => setFilterMinIntensity(e.target.value)}
            placeholder="any"
            style={{ ...filterInputStyle, width: 76 }}
          />
        </label>
        {filtered ? (
          <button
            type="button"
            style={pagerBtnStyle}
            onClick={() => {
              setFilterCharge("");
              setFilterMinIntensity("");
            }}
          >
            clear
          </button>
        ) : null}
      </div>
      {total > pageSize ? (
        <div style={drawerPagerStyle}>
          <button
            style={pagerBtnStyle}
            disabled={clampedPage === 0}
            onClick={() => setPage(0)}
            type="button"
            title="First page"
          >
            ⏮
          </button>
          <button
            style={pagerBtnStyle}
            disabled={clampedPage === 0}
            onClick={() => setPage((p) => Math.max(0, p - 1))}
            type="button"
            title="Previous page"
          >
            ◀
          </button>
          <span style={pagerLabelStyle}>
            {first.toLocaleString()}–{last.toLocaleString()} of {total.toLocaleString()}
          </span>
          <button
            style={pagerBtnStyle}
            disabled={clampedPage >= pageCount - 1}
            onClick={() => setPage((p) => Math.min(pageCount - 1, p + 1))}
            type="button"
            title="Next page"
          >
            ▶
          </button>
          <button
            style={pagerBtnStyle}
            disabled={clampedPage >= pageCount - 1}
            onClick={() => setPage(pageCount - 1)}
            type="button"
            title="Last page"
          >
            ⏭
          </button>
        </div>
      ) : null}
      <div style={drawerBodyStyle} ref={bodyRef}>
        <table style={tableStyle}>
          <thead>
            <tr>
              <th style={{ ...thStyle, width: 18, paddingRight: 0 }} title="Plot colour" />
              {sortTh("mass", "Mono mass")}
              {sortTh("z", "z")}
              {sortTh("mz", "m/z")}
              {sortTh("rt", "RT")}
              {sortTh("intensity", "Intensity", true)}
              {scored ? (
                <>
                  {sortTh("decon", "Decon", true)}
                  {sortTh("minDecon", "MinDec", true)}
                  {sortTh("maxIso", "MaxIso", true)}
                  {sortTh("ppm", "PPM", true)}
                  {sortTh("corrAll", "Corr(all)", true)}
                  {sortTh("corr5", "Corr5", true)}
                  {sortTh("corr3", "Corr3", true)}
                </>
              ) : null}
            </tr>
          </thead>
          <tbody>
            {shown.map((gi) => {
              const f = features[gi];
              return (
                <tr
                  key={gi}
                  data-row={gi}
                  onClick={() => onSelect(gi)}
                  style={gi === selected ? rowSelectedStyle : rowStyle}
                >
                  <td style={{ ...tdStyle, paddingRight: 0 }}>
                    <span
                      style={{
                        color: gi === selected ? SELECTED_COLOR : colorOf(gi),
                        fontSize: 13,
                        lineHeight: 1
                      }}
                    >
                      ■
                    </span>
                  </td>
                  <td style={tdStyle}>{f.monoisotopicMass.toFixed(2)}</td>
                  <td style={tdStyle}>{f.chargeStates.join(",")}</td>
                  <td style={tdStyle}>{f.detectedMz.toFixed(3)}</td>
                  <td style={tdStyle}>{f.rtApex.toFixed(2)}</td>
                  <td style={tdStyleRight}>{f.summedIntensity.toExponential(1)}</td>
                  {scored ? (
                    <>
                      <td style={tdStyleRight}>{fmtScore(f.deconScore)}</td>
                      <td style={tdStyleRight}>{fmtScore(f.minDeconScore)}</td>
                      <td style={tdStyleRight}>{fmtScore(f.maxNumIsotopes, { digits: 0 })}</td>
                      <td style={tdStyleRight}>{fmtScore(f.ppmSpread, { sentinel: 999, digits: 1 })}</td>
                      <td style={tdStyleRight}>{fmtScore(f.isoCorrAll, { sentinel: -2 })}</td>
                      <td style={tdStyleRight}>{fmtScore(f.isoCorrTop5, { sentinel: -2 })}</td>
                      <td style={tdStyleRight}>{fmtScore(f.isoCorrTop3, { sentinel: -2 })}</td>
                    </>
                  ) : null}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </ResizableDrawer>
  );
}

// Format an optional feature score; missing values and the "uncomputable" sentinels render as "—".
function fmtScore(
  v: number | null | undefined,
  opts: { sentinel?: number; digits?: number } = {}
): string {
  if (v == null) return "—";
  if (opts.sentinel != null && v === opts.sentinel) return "—";
  return v.toFixed(opts.digits ?? 2);
}

// True when a feature carries any of the optional resolved-TSV score columns.
function hasScores(f: Feature): boolean {
  return (
    f.deconScore != null ||
    f.minDeconScore != null ||
    f.maxNumIsotopes != null ||
    f.ppmSpread != null ||
    f.isoCorrAll != null ||
    f.isoCorrTop5 != null ||
    f.isoCorrTop3 != null
  );
}

// A compact one-line score summary for the selected-feature readout; "" when the feature is unscored.
function scoreReadout(f: Feature): string {
  if (!hasScores(f)) return "";
  const corr =
    `${fmtScore(f.isoCorrAll, { sentinel: -2 })}/` +
    `${fmtScore(f.isoCorrTop5, { sentinel: -2 })}/` +
    `${fmtScore(f.isoCorrTop3, { sentinel: -2 })}`;
  return (
    ` · Decon ${fmtScore(f.deconScore)} (min ${fmtScore(f.minDeconScore)})` +
    ` · MaxIso ${fmtScore(f.maxNumIsotopes, { digits: 0 })}` +
    ` · PPM ${fmtScore(f.ppmSpread, { sentinel: 999, digits: 1 })}` +
    ` · IsoCorr ${corr}`
  );
}

// Format a q-value compactly: exponential for the very small, fixed otherwise.
function formatQ(q: number): string {
  if (!(q >= 0)) return "—";
  if (q === 0) return "0";
  return q < 0.001 ? q.toExponential(1) : q.toFixed(4);
}

type PsmSortKey = "sequence" | "mass" | "mz" | "charge" | "rt" | "score" | "qValue";

// A right-side drawer listing PSMs with sortable columns and sequence / q-value filters.
// Rows carry their link status to a detected feature; click a row to extract its XIC.
function PsmDrawer({
  psms,
  links,
  selected,
  onSelect,
  onClose
}: {
  psms: readonly Psm[];
  links: readonly PsmLink[];
  selected: number | null;
  onSelect: (index: number) => void;
  onClose: () => void;
}) {
  const [filter, setFilter] = useState("");
  const [maxQ, setMaxQ] = useState("");
  const [minScore, setMinScore] = useState("");
  const [sortKey, setSortKey] = useState<PsmSortKey>("qValue");
  const [sortDir, setSortDir] = useState<"asc" | "desc">("asc");

  // Filtered + sorted view over PSM indices, so onSelect gets the original index.
  const rows = useMemo(() => {
    const q = filter.trim().toLowerCase();
    const qMaxRaw = maxQ.trim();
    const qMax = qMaxRaw === "" ? null : Number(qMaxRaw);
    const scoreMinRaw = minScore.trim();
    const scoreMin = scoreMinRaw === "" ? null : Number(scoreMinRaw);
    const idx = psms
      .map((_, i) => i)
      .filter((i) => {
        const p = psms[i];
        if (q && !p.fullSequence.toLowerCase().includes(q)) return false;
        if (qMax !== null && !Number.isNaN(qMax) && p.qValue > qMax) return false;
        if (scoreMin !== null && !Number.isNaN(scoreMin) && p.score < scoreMin) return false;
        return true;
      });
    const keyVal = (p: Psm): number | string => {
      switch (sortKey) {
        case "sequence":
          return p.fullSequence;
        case "mass":
          return p.monoisotopicMass;
        case "mz":
          return p.precursorMz;
        case "charge":
          return p.precursorCharge;
        case "rt":
          return p.ms2RetentionTime;
        case "score":
          return p.score;
        case "qValue":
          return p.qValue;
      }
    };
    idx.sort((a, b) => {
      const va = keyVal(psms[a]);
      const vb = keyVal(psms[b]);
      const cmp =
        typeof va === "string" ? va.localeCompare(vb as string) : va - (vb as number);
      return sortDir === "asc" ? cmp : -cmp;
    });
    return idx;
  }, [psms, filter, maxQ, minScore, sortKey, sortDir]);

  const shown = rows.slice(0, DRAWER_CAP);
  const toggleSort = (k: PsmSortKey) => {
    if (k === sortKey) {
      setSortDir((d) => (d === "asc" ? "desc" : "asc"));
    } else {
      setSortKey(k);
      // Sequence + q-value read best ascending; the rest default to descending.
      setSortDir(k === "sequence" || k === "qValue" ? "asc" : "desc");
    }
  };
  const arrow = (k: PsmSortKey) => (k === sortKey ? (sortDir === "asc" ? " ▲" : " ▼") : "");

  return (
    <div style={psmDrawerStyle}>
      <div style={drawerHeaderStyle}>
        <strong>
          PSMs{" "}
          {rows.length !== psms.length ? `(${rows.length} of ${psms.length})` : `(${psms.length})`}
        </strong>
        <button onClick={onClose} style={drawerCloseStyle} type="button">
          ✕
        </button>
      </div>
      <div style={psmFilterRowStyle}>
        <input
          type="text"
          placeholder="Filter sequence…"
          value={filter}
          onChange={(e) => setFilter(e.target.value)}
          style={psmInputStyle}
        />
        <input
          type="text"
          inputMode="decimal"
          placeholder="max q"
          value={maxQ}
          onChange={(e) => setMaxQ(e.target.value)}
          style={{ ...psmInputStyle, width: 64, flex: "0 0 auto" }}
        />
        <input
          type="text"
          inputMode="decimal"
          placeholder="min score"
          value={minScore}
          onChange={(e) => setMinScore(e.target.value)}
          style={{ ...psmInputStyle, width: 72, flex: "0 0 auto" }}
        />
      </div>
      <div style={drawerBodyStyle}>
        <table style={tableStyle}>
          <thead>
            <tr>
              <th style={thSortStyle} onClick={() => toggleSort("sequence")}>
                Sequence{arrow("sequence")}
              </th>
              <th style={thSortRightStyle} onClick={() => toggleSort("mass")}>
                Mass{arrow("mass")}
              </th>
              <th style={thSortRightStyle} onClick={() => toggleSort("mz")}>
                m/z{arrow("mz")}
              </th>
              <th style={thSortRightStyle} onClick={() => toggleSort("charge")}>
                z{arrow("charge")}
              </th>
              <th style={thSortRightStyle} onClick={() => toggleSort("rt")}>
                RT{arrow("rt")}
              </th>
              <th style={thSortRightStyle} onClick={() => toggleSort("score")}>
                Score{arrow("score")}
              </th>
              <th style={thSortRightStyle} onClick={() => toggleSort("qValue")}>
                q{arrow("qValue")}
              </th>
              <th style={thStyle}>Feat.</th>
            </tr>
          </thead>
          <tbody>
            {shown.map((i) => {
              const p = psms[i];
              const linked = links[i]?.featureIndex != null;
              return (
                <tr
                  key={i}
                  onClick={() => onSelect(i)}
                  style={i === selected ? rowSelectedStyle : rowStyle}
                >
                  <td style={tdSeqStyle} title={p.fullSequence}>
                    {p.fullSequence}
                  </td>
                  <td style={tdStyleRight}>{p.monoisotopicMass.toFixed(2)}</td>
                  <td style={tdStyleRight}>{p.precursorMz.toFixed(3)}</td>
                  <td style={tdStyleRight}>{p.precursorCharge}</td>
                  <td style={tdStyleRight}>
                    {p.ms2RetentionTime >= 0 ? p.ms2RetentionTime.toFixed(2) : "—"}
                  </td>
                  <td style={tdStyleRight}>{p.score.toFixed(2)}</td>
                  <td style={tdStyleRight}>{formatQ(p.qValue)}</td>
                  <td
                    style={{ ...tdStyle, textAlign: "center", color: linked ? "#1e7e34" : "#b9c2cf" }}
                    title={linked ? "linked to a detected feature" : "no matching feature"}
                  >
                    {linked ? "●" : "—"}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}

const dropOverlayStyle: React.CSSProperties = {
  position: "fixed",
  inset: 0,
  // Above the drawer (zIndex 1000) so the cue is never occluded.
  zIndex: 2000,
  display: "flex",
  alignItems: "center",
  justifyContent: "center",
  backgroundColor: "rgba(47,111,176,0.12)",
  border: "3px dashed #2f6fb0",
  // Never intercept the OS drop; this is a visual cue only.
  pointerEvents: "none"
};
const dropOverlayCardStyle: React.CSSProperties = {
  padding: "16px 28px",
  backgroundColor: "#ffffff",
  border: "1.5px solid #2f6fb0",
  borderRadius: 8,
  boxShadow: "0 6px 24px rgba(0,0,0,0.18)",
  fontSize: 16,
  fontWeight: 600,
  color: "#2f6fb0"
};

const drawerStyle: React.CSSProperties = {
  position: "fixed",
  top: 0,
  right: 0,
  height: "100vh",
  // `width` is supplied by ResizableDrawer (user-resizable); this default is a fallback only.
  width: DEFAULT_DRAWER_WIDTH,
  backgroundColor: "#ffffff",
  borderLeft: "1.5px solid #000",
  boxShadow: "-4px 0 16px rgba(0,0,0,0.12)",
  display: "flex",
  flexDirection: "column",
  zIndex: 1000
};
// Thin draggable strip on the drawer's left edge (straddling the border) for resizing.
const drawerResizeHandleStyle: React.CSSProperties = {
  position: "absolute",
  left: -3,
  top: 0,
  width: 7,
  height: "100%",
  cursor: "col-resize",
  zIndex: 1001,
  touchAction: "none"
};
const drawerHeaderStyle: React.CSSProperties = {
  display: "flex",
  justifyContent: "space-between",
  alignItems: "center",
  padding: "8px 12px",
  borderBottom: "1px solid #d5deeb",
  fontSize: "0.85rem"
};
const drawerCloseStyle: React.CSSProperties = {
  border: "none",
  background: "transparent",
  cursor: "pointer",
  fontSize: "0.9rem"
};
const drawerPagerStyle: React.CSSProperties = {
  display: "flex",
  alignItems: "center",
  justifyContent: "center",
  gap: 6,
  padding: "4px 8px",
  borderBottom: "1px solid #d5deeb",
  fontSize: "0.78rem",
  color: "#3a4a60"
};
const pagerBtnStyle: React.CSSProperties = {
  border: "1px solid #cbd6e5",
  background: "#f4f7fb",
  borderRadius: 4,
  cursor: "pointer",
  padding: "1px 7px",
  lineHeight: 1.4,
  fontSize: "0.8rem"
};
const pagerLabelStyle: React.CSSProperties = { minWidth: 132, textAlign: "center" };
const drawerFilterStyle: React.CSSProperties = {
  display: "flex",
  alignItems: "center",
  gap: 10,
  padding: "5px 10px",
  borderBottom: "1px solid #d5deeb",
  fontSize: "0.76rem",
  color: "#3a4a60"
};
const filterLabelStyle: React.CSSProperties = { display: "flex", alignItems: "center", gap: 4 };
const filterInputStyle: React.CSSProperties = {
  width: 56,
  fontSize: "0.76rem",
  padding: "1px 4px",
  border: "1px solid #cbd6e5",
  borderRadius: 4
};
const drawerBodyStyle: React.CSSProperties = { overflow: "auto", minHeight: 0, flex: 1 };
const tableStyle: React.CSSProperties = {
  width: "100%",
  borderCollapse: "collapse",
  fontSize: "0.72rem"
};
const thStyle: React.CSSProperties = {
  position: "sticky",
  top: 0,
  background: "#eef3fa",
  textAlign: "left",
  padding: "5px 8px",
  borderBottom: "1px solid #cdd8e8"
};
const thStyleRight: React.CSSProperties = { ...thStyle, textAlign: "right" };
const rowStyle: React.CSSProperties = { cursor: "pointer", borderBottom: "1px solid #eef2f7" };
const rowSelectedStyle: React.CSSProperties = { ...rowStyle, background: "#fdecd6" };
const tdStyle: React.CSSProperties = { padding: "4px 8px", whiteSpace: "nowrap" };
const tdStyleRight: React.CSSProperties = { ...tdStyle, textAlign: "right" };

// ---- PSM drawer styles (header + filter row + scrollable table)
const psmDrawerStyle: React.CSSProperties = {
  ...drawerStyle,
  gridTemplateRows: "auto auto 1fr"
};
const psmFilterRowStyle: React.CSSProperties = {
  display: "flex",
  gap: 6,
  padding: "6px 10px",
  borderBottom: "1px solid #e3e9f2"
};
const psmInputStyle: React.CSSProperties = {
  flex: 1,
  minWidth: 0,
  fontSize: "0.72rem",
  padding: "3px 6px",
  border: "1px solid #cdd8e8",
  borderRadius: 4
};
const thSortStyle: React.CSSProperties = { ...thStyle, cursor: "pointer", userSelect: "none" };
const thSortRightStyle: React.CSSProperties = {
  ...thStyleRight,
  cursor: "pointer",
  userSelect: "none"
};
const tdSeqStyle: React.CSSProperties = {
  ...tdStyle,
  maxWidth: 148,
  overflow: "hidden",
  textOverflow: "ellipsis"
};

// ---- walkthrough drawer styles
const ladderHintStyle: React.CSSProperties = {
  fontSize: "0.72rem",
  color: "#5b6b7d",
  lineHeight: 1.4,
  marginBottom: 8
};
const chargeGridStyle: React.CSSProperties = {
  display: "flex",
  flexWrap: "wrap",
  gap: 4,
  marginBottom: 12
};
const zSeedButtonStyle: React.CSSProperties = {
  minWidth: 26,
  padding: "3px 6px",
  fontSize: "0.72rem",
  border: "1px solid #cdd8e8",
  background: "#fff",
  borderRadius: 4,
  cursor: "pointer"
};
const zSeedButtonSelStyle: React.CSSProperties = {
  ...zSeedButtonStyle,
  background: "#2f6fb0",
  color: "#fff",
  borderColor: "#2f6fb0",
  fontWeight: 600
};
const ladderNoteStyle: React.CSSProperties = {
  fontSize: "0.78rem",
  color: "#5b6b7d",
  padding: "12px 0"
};
const summaryTableStyle: React.CSSProperties = {
  width: "100%",
  borderCollapse: "collapse",
  fontSize: "0.73rem",
  marginBottom: 10
};
const sumKeyStyle: React.CSSProperties = {
  padding: "3px 6px 3px 0",
  color: "#5b6b7d",
  verticalAlign: "top",
  whiteSpace: "nowrap",
  width: 92
};
const sumValStyle: React.CSSProperties = { padding: "3px 0", verticalAlign: "top" };
const ladderTableWrapStyle: React.CSSProperties = {
  border: "1px solid #e3e9f2",
  borderRadius: 4,
  overflow: "hidden"
};
