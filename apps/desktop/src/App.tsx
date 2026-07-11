import { useCallback, useMemo, useState } from "react";
import { open as openFileDialog } from "@tauri-apps/plugin-dialog";

import { TicPlot, SpectrumPlot, createDefaultViewport } from "@msbrowser/plot-adapter";
import type {
  TicPlotTrace,
  SpectrumPlotTrace,
  FeatureMarker,
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
import type {
  DatasetMetadata,
  DatasetProvider,
  Feature,
  Spectrum,
  TicPoint
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

  const [features, setFeatures] = useState<readonly Feature[]>([]);
  const [featuresFile, setFeaturesFile] = useState<string | null>(null);
  const [selected, setSelected] = useState<number | null>(null);
  const [drawerOpen, setDrawerOpen] = useState(false);

  const [ticViewport, setTicViewport] = useState<PlotViewport>(createDefaultViewport());
  const [spectrumViewport, setSpectrumViewport] = useState<PlotViewport>(createDefaultViewport());

  const [ticPinned, setTicPinned] = useState(false);
  const [spectrumPinned, setSpectrumPinned] = useState(false);

  // ------------------------------------------------------------ open raw/mzML
  const handleOpenFile = useCallback(async () => {
    const picked = await openFileDialog({
      multiple: false,
      filters: [{ name: "Mass spec", extensions: ["raw", "mzML", "mzml", "mzMLb", "mgf"] }]
    });
    if (typeof picked !== "string") return;

    setLoad({ status: "loading", message: "Opening…" });
    setSpectrum(null);
    setSelected(null);
    setTicViewport(createDefaultViewport());
    try {
      const { handle: h, provider: p } = await openDataset(picked, (progress) => {
        setLoad({ status: "loading", message: `${progress.phase}…` });
      });
      const meta = await p.getMetadata();
      const tic = await p.getTicTrace({ maxPoints: 4000 });
      setHandle(h);
      setProvider(p);
      setMetadata(meta);
      setTicPoints(tic);
      setLoad({ status: "ready" });
    } catch (err) {
      setLoad({ status: "error", message: errMessage(err) });
    }
  }, []);

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

  // -------------------------------------------------- run detection in-app
  const handleRunDetection = useCallback(async () => {
    if (handle === null || detecting !== null) return;
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
  }, [handle, detecting]);

  // ------------------------------------------------------- select a feature
  const selectFeature = useCallback(
    async (index: number) => {
      const f = features[index];
      if (!f) return;
      setSelected(index);

      const pad = Math.max(0.2, (f.rtEnd - f.rtStart) * 0.6);
      setTicViewport({ xMin: f.rtStart - pad, xMax: f.rtEnd + pad });

      const z = f.primaryCharge || f.chargeStates[0] || 1;
      const grid = isotopeGrid(f.monoisotopicMass, z, 12);
      setSpectrumViewport({ xMin: grid[0] - 1.5, xMax: grid[grid.length - 1] + 1.5 });

      if (provider && !spectrumPinned) {
        try {
          const scan = await provider.getNearestScan(f.rtApex);
          if (scan) setSpectrum(await provider.getSpectrum(scan.scanIndex));
        } catch (err) {
          setLoad({ status: "error", message: errMessage(err) });
        }
      }
    },
    [features, provider, spectrumPinned]
  );

  // Click the TIC background → nearest scan's spectrum (unless the spectrum is pinned).
  const handleAreaClick = useCallback(
    async (rt: number) => {
      if (!provider || spectrumPinned) return;
      try {
        const scan = await provider.getNearestScan(rt);
        if (!scan) return;
        setSpectrum(await provider.getSpectrum(scan.scanIndex));
      } catch (err) {
        setLoad({ status: "error", message: errMessage(err) });
      }
    },
    [provider, spectrumPinned]
  );

  // ----------------------------------------------------------------- overlays
  const ticTraces = useMemo<TicPlotTrace[]>(
    () => [
      {
        slotIndex: 0,
        points: ticPoints,
        selectedScanIndex: spectrum?.scanIndex ?? null,
        color: SLOT_COLOR
      }
    ],
    [ticPoints, spectrum]
  );

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

  const spectrumTraces = useMemo<SpectrumPlotTrace[]>(
    () => (spectrum ? [{ slotIndex: 0, peaks: spectrum.peaks, color: SLOT_COLOR }] : []),
    [spectrum]
  );

  const selectedFeature = selected === null ? null : features[selected] ?? null;

  return (
    <ViewerShell
      title="MsViewer — feature finder"
      subtitle="Real raw/mzML over Tauri IPC, with top-down feature-finding results overlaid."
      toolbar={
        <>
          {metadata ? (
            <>
              <MetricReadout label="File" value={metadata.fileName} />
              <MetricReadout label="MS1 scans" value={metadata.ms1ScanCount} />
              <MetricReadout label="Format" value={metadata.format} />
            </>
          ) : null}
          <PanelActionButton onClick={() => void handleOpenFile()}>Open file…</PanelActionButton>
          <PanelActionButton onClick={() => void handleLoadFeatures()}>
            Load features…
          </PanelActionButton>
          {handle !== null ? (
            <PanelActionButton onClick={() => void handleRunDetection()}>
              {detecting ? `Detecting: ${detecting}` : "Run feature finding"}
            </PanelActionButton>
          ) : null}
          {features.length > 0 ? (
            <PanelActionButton pressed={drawerOpen} onClick={() => setDrawerOpen((v) => !v)}>
              Features ({features.length})
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
              selectedFeature
                ? `Feature: ${selectedFeature.monoisotopicMass.toFixed(2)} Da · z ${selectedFeature.chargeStates.join(",")} · RT ${selectedFeature.rtStart.toFixed(2)}–${selectedFeature.rtEnd.toFixed(2)}`
                : featuresFile
                  ? `${features.length} features from ${featuresFile} — click a marker or a row`
                  : "Summed MS1 intensity · click to load a scan"
            }
            actions={
              <>
                {ticViewport.xMin !== null ? (
                  <PanelActionButton onClick={() => setTicViewport(createDefaultViewport())}>
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
          <StatusBanner tone="muted">No dataset — use “Open file…” to load a .raw or mzML.</StatusBanner>
        ) : (
          <TicPlot
            traces={ticTraces}
            viewport={ticViewport}
            featureRug={featureRug}
            regions={regions}
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
                ? `Scan ${spectrum.oneBasedScanNumber} · RT ${spectrum.retentionTime.toFixed(2)} min${selectedFeature ? " · dotted lines = predicted isotope m/z" : ""}`
                : "Click the chromatogram or a feature to load a scan"
            }
            actions={
              <>
                {spectrumViewport.xMin !== null ? (
                  <PanelActionButton onClick={() => setSpectrumViewport(createDefaultViewport())}>
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
            envelope={envelope}
            rangeSelectionEnabled={false}
            onEvent={() => {}}
          />
        )}
      </Panel>

      {drawerOpen && features.length > 0 ? (
        <FeatureDrawer
          features={features}
          selected={selected}
          onSelect={(i) => void selectFeature(i)}
          onClose={() => setDrawerOpen(false)}
        />
      ) : null}
    </ViewerShell>
  );
}

function errMessage(err: unknown): string {
  if (err && typeof err === "object" && "message" in err) {
    return String((err as { message: unknown }).message);
  }
  return err instanceof Error ? err.message : String(err);
}

// A right-side drawer listing resolved features; click a row to select.
function FeatureDrawer({
  features,
  selected,
  onSelect,
  onClose
}: {
  features: readonly Feature[];
  selected: number | null;
  onSelect: (index: number) => void;
  onClose: () => void;
}) {
  const shown = features.slice(0, DRAWER_CAP);
  return (
    <div style={drawerStyle}>
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
    </div>
  );
}

const drawerStyle: React.CSSProperties = {
  position: "fixed",
  top: 0,
  right: 0,
  height: "100vh",
  width: 360,
  backgroundColor: "#ffffff",
  borderLeft: "1.5px solid #000",
  boxShadow: "-4px 0 16px rgba(0,0,0,0.12)",
  display: "grid",
  gridTemplateRows: "auto 1fr",
  zIndex: 1000
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
