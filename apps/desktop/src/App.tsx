import { useCallback, useEffect, useMemo, useState } from "react";
import { open as openFileDialog } from "@tauri-apps/plugin-dialog";

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
import { loadFeatures, runFeatureDetection, isotopeGrid } from "./features";
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
// Feature-rug + envelope colours keyed by charge state.
const CHARGE_COLORS = [
  "#2f6fb0", "#c0392b", "#27ae60", "#8e44ad",
  "#d35400", "#16a085", "#b7950b", "#c2185b"
];
const chargeColor = (z: number): string =>
  CHARGE_COLORS[(Math.max(1, z) - 1) % CHARGE_COLORS.length];

// Features are sorted by intensity on load; these cap what's drawn/listed so a
// pathologically large TSV (top-down noise runs can be 100k+ features) stays
// responsive. The strongest features are always the ones kept.
const RUG_CAP = 8000;
const DRAWER_CAP = 800;

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
  const reframeSpectrum = useCallback((vp: PlotViewport) => {
    setSpectrumViewport(vp);
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
      // The fast-path TIC is now the smooth MS1-only trace read from per-scan metadata (no peak
      // decode, no index), so paint it immediately — no jagged preview, and no waiting on the
      // full index. Empty only for files without per-scan TIC metadata, which fill in from the
      // index via `markReady` (below); the panel shows a "building index" note until then.
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
      filters: [{ name: "Mass spec", extensions: ["raw", "mzML", "mzml", "mzMLb", "mgf"] }]
    });
    if (typeof picked !== "string") return;
    await openPath(picked);
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
        // Re-fetch the TIC. For files with per-scan TIC metadata this returns the same smooth
        // MS1-only trace already painted on open (no visible change). For files *without* it, the
        // fast-path trace was empty and only now fills in from the freshly-built MS1 index.
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
    (drawerOpen && features.length > 0);
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
  const handleLoadFeatures = useCallback(async () => {
    const picked = await openFileDialog({
      multiple: false,
      filters: [{ name: "Feature TSV", extensions: ["tsv", "txt"] }]
    });
    if (typeof picked !== "string") return;
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
  const selectFeature = useCallback(
    async (index: number) => {
      const f = features[index];
      if (!f) return;
      setSelected(index);

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
    [features, provider, spectrumPinned, reframeTic, reframeSpectrum]
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

  // Click a spectrum peak (walkthrough only) → make it the seed at the current zSeed (frame it).
  const handlePeakClick = useCallback(
    (mz: number) => {
      if (!walkthrough || indexing) return;
      void computeLadder(mz, zSeed, true);
    },
    [walkthrough, indexing, zSeed, computeLadder]
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

  const featureRug = useMemo<FeatureMarker[]>(
    () =>
      features.slice(0, RUG_CAP).map((f, i) => ({
        featureIndex: i,
        retentionTime: f.rtApex,
        color: i === selected ? SELECTED_COLOR : chargeColor(f.primaryCharge),
        label: `m/z ${f.detectedMz.toFixed(3)} · z${f.primaryCharge} · ${f.monoisotopicMass.toFixed(1)} Da · RT ${f.rtApex.toFixed(2)}`
      })),
    [features, selected]
  );

  const regions = useMemo<RtRegion[]>(() => {
    if (selected === null) return [];
    const f = features[selected];
    if (!f) return [];
    return [{ min: f.rtStart, max: f.rtEnd, color: "rgba(232,131,12,0.14)" }];
  }, [features, selected]);

  const envelope = useMemo<EnvelopeLine[]>(() => {
    if (selected === null) return [];
    const f = features[selected];
    if (!f) return [];
    return f.chargeStates.flatMap((z) =>
      isotopeGrid(f.monoisotopicMass, z, 12).map((mz, k) => ({
        mz,
        color: chargeColor(z),
        label: k === 0 ? `z${z} mono` : undefined
      }))
    );
  }, [features, selected]);

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

  // The feature isotope envelope only makes sense over an MS1 scan; a PSM's identified spectrum
  // is MS2 (fragment ions), so suppress the comb there.
  const spectrumEnvelope = walkthrough
    ? ladderEnvelope
    : spectrum?.msLevel === 1
      ? envelope
      : [];

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

  const selectedFeature = selected === null ? null : features[selected] ?? null;

  return (
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
          <PanelActionButton onClick={() => void handleLoadPsms()}>Load PSMs…</PanelActionButton>
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
                  if (next) setPsmDrawerOpen(false);
                  return next;
                })
              }
            >
              Features ({features.length})
            </PanelActionButton>
          ) : null}
          {psms.length > 0 ? (
            <PanelActionButton
              pressed={psmDrawerOpen}
              onClick={() =>
                setPsmDrawerOpen((v) => {
                  const next = !v;
                  if (next) setDrawerOpen(false);
                  return next;
                })
              }
            >
              PSMs ({psms.length})
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
                  ? `Feature: ${selectedFeature.monoisotopicMass.toFixed(2)} Da · z ${selectedFeature.chargeStates.join(",")} · RT ${selectedFeature.rtStart.toFixed(2)}–${selectedFeature.rtEnd.toFixed(2)}`
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
              else if (e.type === "feature-click") void selectFeature(e.featureIndex);
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
                      : selectedFeature
                        ? " · dotted lines = predicted isotope m/z"
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
            annotations={spectrumAnnotations}
            uirevision={spectrumUiRev}
            rangeSelectionEnabled={false}
            onEvent={(e) => {
              if (e.type === "peak-click") handlePeakClick(e.peak.mz);
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
        />
      ) : null}
    </ViewerShell>
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

// A right-side drawer listing resolved features; click a row to select.
function FeatureDrawer({
  features,
  selected,
  width,
  onResize,
  onSelect,
  onClose
}: {
  features: readonly Feature[];
  selected: number | null;
  width: number;
  onResize: (width: number) => void;
  onSelect: (index: number) => void;
  onClose: () => void;
}) {
  const shown = features.slice(0, DRAWER_CAP);
  return (
    <ResizableDrawer width={width} onResize={onResize}>
      <div style={drawerHeaderStyle}>
        <strong>
          Features{" "}
          {features.length > DRAWER_CAP
            ? `(top ${DRAWER_CAP} of ${features.length})`
            : `(${features.length})`}
        </strong>
        <button onClick={onClose} style={drawerCloseStyle} type="button">
          ✕
        </button>
      </div>
      <div style={drawerBodyStyle}>
        <table style={tableStyle}>
          <thead>
            <tr>
              <th style={thStyle}>Mono mass</th>
              <th style={thStyle}>z</th>
              <th style={thStyle}>m/z</th>
              <th style={thStyle}>RT</th>
              <th style={thStyleRight}>Intensity</th>
            </tr>
          </thead>
          <tbody>
            {shown.map((f, i) => (
              <tr
                key={i}
                onClick={() => onSelect(i)}
                style={i === selected ? rowSelectedStyle : rowStyle}
              >
                <td style={tdStyle}>{f.monoisotopicMass.toFixed(2)}</td>
                <td style={tdStyle}>{f.chargeStates.join(",")}</td>
                <td style={tdStyle}>{f.detectedMz.toFixed(3)}</td>
                <td style={tdStyle}>{f.rtApex.toFixed(2)}</td>
                <td style={tdStyleRight}>{f.summedIntensity.toExponential(1)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </ResizableDrawer>
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
  display: "grid",
  gridTemplateRows: "auto 1fr",
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
const drawerBodyStyle: React.CSSProperties = { overflow: "auto", minHeight: 0 };
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
