import { useCallback, useEffect, useMemo, useState } from "react";

import { TicPlot, SpectrumPlot, createDefaultViewport } from "@msbrowser/plot-adapter";
import type { TicPlotTrace, SpectrumPlotTrace } from "@msbrowser/plot-adapter";
import {
  Panel,
  PanelActionButton,
  PanelHeader,
  StatusBanner,
  ViewerShell,
  MetricReadout
} from "@msbrowser/ui";

import { openDataset } from "./tauri-dataset-provider";
import type {
  DatasetMetadata,
  DatasetProvider,
  Spectrum,
  TicPoint
} from "./contract";

const SLOT_COLOR = "#2f6fb0";

type LoadState =
  | { status: "idle" }
  | { status: "loading"; message: string }
  | { status: "ready" }
  | { status: "error"; message: string };

export function App() {
  const [load, setLoad] = useState<LoadState>({ status: "idle" });
  const [provider, setProvider] = useState<DatasetProvider | null>(null);
  const [metadata, setMetadata] = useState<DatasetMetadata | null>(null);
  const [ticPoints, setTicPoints] = useState<readonly TicPoint[]>([]);
  const [spectrum, setSpectrum] = useState<Spectrum | null>(null);

  // A pinned panel is frozen: it keeps its current content/view while you explore the
  // other panel. Pin the spectrum to hold a reference scan (MS1 or MS2) while you keep
  // clicking across the chromatogram; the pinned panel is highlighted.
  const [ticPinned, setTicPinned] = useState(false);
  const [spectrumPinned, setSpectrumPinned] = useState(false);

  const openStub = useCallback(async () => {
    setLoad({ status: "loading", message: "Opening stub dataset…" });
    setSpectrum(null);
    try {
      const { provider: p } = await openDataset("stub://demo", (progress) => {
        setLoad({
          status: "loading",
          message: `${progress.phase} (${progress.scansDone}/${progress.scansTotal})`
        });
      });
      const meta = await p.getMetadata();
      const tic = await p.getTicTrace();
      setProvider(p);
      setMetadata(meta);
      setTicPoints(tic);
      setLoad({ status: "ready" });
    } catch (err) {
      const message = err instanceof Error ? err.message : String(err);
      setLoad({ status: "error", message });
    }
  }, []);

  // Auto-open the stub dataset on mount so a real chromatogram renders immediately.
  useEffect(() => {
    void openStub();
  }, [openStub]);

  // Click the TIC to load the nearest scan's spectrum into the bottom panel — unless that
  // panel is pinned, in which case it stays frozen for comparison.
  const handleTicEvent = useCallback(
    async (rt: number) => {
      if (!provider || spectrumPinned) return;
      try {
        const scan = await provider.getNearestScan(rt);
        if (!scan) return;
        setSpectrum(await provider.getSpectrum(scan.scanIndex));
      } catch (err) {
        const message = err instanceof Error ? err.message : String(err);
        setLoad({ status: "error", message });
      }
    },
    [provider, spectrumPinned]
  );

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

  const spectrumTraces = useMemo<SpectrumPlotTrace[]>(
    () =>
      spectrum
        ? [{ slotIndex: 0, peaks: spectrum.peaks, color: SLOT_COLOR }]
        : [],
    [spectrum]
  );

  const viewport = useMemo(() => createDefaultViewport(), []);

  return (
    <ViewerShell
      title="MsViewer (M0 skeleton)"
      subtitle="Tauri + Rust stub — a real Arrow-encoded chromatogram over IPC."
      toolbar={
        <>
          {metadata ? (
            <>
              <MetricReadout label="File" value={metadata.fileName} />
              <MetricReadout label="Scans" value={metadata.scanCount} />
              <MetricReadout label="Format" value={metadata.format} />
            </>
          ) : null}
          <PanelActionButton onClick={() => void openStub()}>
            Open (stub)
          </PanelActionButton>
        </>
      }
    >
      <Panel
        active={ticPinned}
        header={
          <PanelHeader
            title="Total ion chromatogram"
            subtitle={
              ticPinned
                ? "Pinned · view frozen"
                : "Stub data · Rust-produced Arrow IPC → TicPlot"
            }
            actions={
              <PanelActionButton
                pressed={ticPinned}
                onClick={() => setTicPinned((v) => !v)}
              >
                {ticPinned ? "Pinned" : "Pin"}
              </PanelActionButton>
            }
          />
        }
      >
        {load.status === "loading" ? (
          <StatusBanner tone="info">{load.message}</StatusBanner>
        ) : load.status === "error" ? (
          <StatusBanner tone="error">Failed to load: {load.message}</StatusBanner>
        ) : ticPoints.length === 0 ? (
          <StatusBanner tone="muted">No chromatogram loaded.</StatusBanner>
        ) : (
          <TicPlot
            traces={ticTraces}
            viewport={viewport}
            rangeSelectionEnabled={false}
            onEvent={(e) => {
              if (e.type === "area-click") void handleTicEvent(e.retentionTime);
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
              spectrumPinned && spectrum
                ? `Pinned · Scan ${spectrum.oneBasedScanNumber} · RT ${spectrum.retentionTime.toFixed(2)} min`
                : spectrum
                  ? `Scan ${spectrum.oneBasedScanNumber} · RT ${spectrum.retentionTime.toFixed(2)} min · MS${spectrum.msLevel}`
                  : "Click the chromatogram to load a scan"
            }
            actions={
              <PanelActionButton
                pressed={spectrumPinned}
                onClick={() => setSpectrumPinned((v) => !v)}
              >
                {spectrumPinned ? "Pinned" : "Pin"}
              </PanelActionButton>
            }
          />
        }
      >
        {spectrumTraces.length === 0 ? (
          <StatusBanner tone="muted">No spectrum selected.</StatusBanner>
        ) : (
          <SpectrumPlot
            traces={spectrumTraces}
            viewport={viewport}
            rangeSelectionEnabled={false}
            onEvent={() => {}}
          />
        )}
      </Panel>
    </ViewerShell>
  );
}
