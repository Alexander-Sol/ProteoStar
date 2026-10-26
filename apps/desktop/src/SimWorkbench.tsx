// The simulation workbench (Top-Down Simulator Phase 3): a live `SimSpec` edited in a parameter
// panel, a simulated scan that re-renders on every edit, and an explicit Accept that simulates the
// whole run (agent_info/TopDown-Simulator-Rust-Port-Plan.md, interaction model):
//
// - Live = the viewed scan only. Every edit goes to one long-lived backend simulator
//   (`SimSession`), which recomputes only what the edit invalidated and returns the viewed scan.
// - Accept = the full run. The committed run's TIC lags the panel until the next Accept.
//
// The viewed scan follows a target RT (the real spectrum's, or a TIC click) to the nearest
// simulated scan, so a spec seeded from the real run's scan times aligns scan for scan.

import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { CSSProperties, ReactNode } from "react";

import type { DatasetProvider, Spectrum, TicPoint } from "./contract";
import {
  SIM_SCHEMA_VERSION,
  orbitrapNoiseDefaults,
  type SimNoise,
  type SimProteoform,
  type SimSpec
} from "./sim-spec";
import { SimSession, closeDataset } from "./tauri-dataset-provider";

/** The simulated run's colour wherever it is drawn (mirror traces, panel accents). */
export const SIM_COLOR = "#b5482f";

/** A synthetic gradient for a spec with no real run to borrow scan times from: 60 min, 1.2 s. */
const BLANK_RUN_MINUTES = 60;
const BLANK_SCAN_SPACING_MINUTES = 0.02;

/** How many proteoform rows the panel lists; the editor reaches any of them via the list. */
const PROTEOFORM_LIST_CAP = 200;

export interface CommittedSim {
  handle: number;
  provider: DatasetProvider;
  tic: readonly TicPoint[];
  /** The spec this run was simulated from — compared by identity to tell whether the TIC lags. */
  spec: SimSpec;
}

/** A starting spec: the run's own scan times if given, else a blank gradient; one proteoform. */
export function newSimSpec(opts: { scanTimes?: readonly number[]; sourceFile?: string } = {}): SimSpec {
  const scanTimes =
    opts.scanTimes && opts.scanTimes.length > 1
      ? [...opts.scanTimes]
      : Array.from(
          { length: Math.round(BLANK_RUN_MINUTES / BLANK_SCAN_SPACING_MINUTES) + 1 },
          (_, i) => i * BLANK_SCAN_SPACING_MINUTES
        );
  const mid = scanTimes[Math.floor(scanTimes.length / 2)];
  return {
    schemaVersion: SIM_SCHEMA_VERSION,
    generator: { tool: "ProteoStar workbench" },
    sourceRun: opts.sourceFile ? { fileName: opts.sourceFile } : undefined,
    instrument: { polarity: "Positive", analyzer: "Orbitrap", msLevel: 1 },
    chargeRange: { min: 3, max: 50 },
    peakWidthModel: { kind: "Orbitrap", k: 1.2e-7 },
    scanTimes,
    noise: orbitrapNoiseDefaults(0),
    proteoforms: [defaultProteoform(mid, "P1")],
    background: []
  };
}

function defaultProteoform(rt: number, identifier: string): SimProteoform {
  return {
    monoisotopicMass: 15000,
    abundance: 1e7,
    // A scan time carries full float noise (24.512670516967773); the panel shows what it sends.
    rt: { mu: Number(rt.toFixed(2)), sigma: 0.2, tau: 0.1 },
    charge: { kind: "Gaussian", muZ: 15, sigmaZ: 3 },
    identifier
  };
}

/** Index of the scan time nearest `rt` (scan times ascend). */
function nearestScanIndex(scanTimes: readonly number[], rt: number): number {
  let lo = 0;
  let hi = scanTimes.length - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (scanTimes[mid] <= rt) lo = mid;
    else hi = mid;
  }
  return Math.abs(scanTimes[hi] - rt) < Math.abs(scanTimes[lo] - rt) ? hi : lo;
}

function errText(err: unknown): string {
  if (err && typeof err === "object" && "message" in err) {
    return String((err as { message: unknown }).message);
  }
  return String(err);
}

export interface SimWorkbench {
  spec: SimSpec | null;
  /** Replace the spec with an edited copy (the panel's only way to change it). */
  update(edit: (spec: SimSpec) => SimSpec): void;
  /** Start (or restart) the workbench on `spec`. */
  start(spec: SimSpec): Promise<void>;
  /** Drop the spec, the session and the committed run. */
  discard(): void;
  /** Point the viewed scan at the simulated scan nearest this RT. */
  setTargetRt(rt: number): void;
  /** Step the viewed scan; for when no real spectrum drives it. */
  step(delta: number): void;
  viewedIndex: number | null;
  liveSpectrum: Spectrum | null;
  /** Round-trip time of the last live scan, ms. */
  liveMs: number | null;
  error: string | null;
  committed: CommittedSim | null;
  committing: boolean;
  /** The panel has changed since the committed run (its TIC lags). */
  stale: boolean;
  accept(): Promise<void>;
}

export function useSimWorkbench(): SimWorkbench {
  const [spec, setSpec] = useState<SimSpec | null>(null);
  const sessionRef = useRef<SimSession | null>(null);
  // Bumped when a session opens, so the live effect re-runs for it.
  const [sessionVersion, setSessionVersion] = useState(0);
  const [targetRt, setTargetRtState] = useState<number | null>(null);
  const [liveSpectrum, setLiveSpectrum] = useState<Spectrum | null>(null);
  const [liveMs, setLiveMs] = useState<number | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [committed, setCommitted] = useState<CommittedSim | null>(null);
  const committedRef = useRef<CommittedSim | null>(null);
  const [committing, setCommitting] = useState(false);

  const viewedIndex = useMemo(() => {
    if (!spec || spec.scanTimes.length === 0) return null;
    const rt = targetRt ?? spec.scanTimes[Math.floor(spec.scanTimes.length / 2)];
    return nearestScanIndex(spec.scanTimes, rt);
  }, [spec, targetRt]);

  // The live loop: every spec edit or viewed-scan change asks the session for that scan. The
  // session is latest-wins, so superseded requests resolve to null and are skipped.
  useEffect(() => {
    const session = sessionRef.current;
    if (!spec || !session || viewedIndex === null) return;
    const t0 = performance.now();
    session
      .scan(spec, viewedIndex)
      .then((s) => {
        if (!s || sessionRef.current !== session) return;
        setLiveSpectrum(s);
        setLiveMs(performance.now() - t0);
        setError(null);
      })
      .catch((err) => {
        if (sessionRef.current === session) setError(errText(err));
      });
  }, [spec, viewedIndex, sessionVersion]);

  const replaceCommitted = useCallback((next: CommittedSim | null) => {
    const prev = committedRef.current;
    committedRef.current = next;
    setCommitted(next);
    if (prev && prev.handle !== next?.handle) void closeDataset(prev.handle).catch(() => {});
  }, []);

  const start = useCallback(
    async (next: SimSpec) => {
      const old = sessionRef.current;
      sessionRef.current = null;
      if (old) void old.close().catch(() => {});
      replaceCommitted(null);
      setLiveSpectrum(null);
      setError(null);
      setSpec(next);
      try {
        sessionRef.current = await SimSession.open(next);
        setSessionVersion((v) => v + 1);
      } catch (err) {
        setError(errText(err));
      }
    },
    [replaceCommitted]
  );

  const discard = useCallback(() => {
    const old = sessionRef.current;
    sessionRef.current = null;
    if (old) void old.close().catch(() => {});
    replaceCommitted(null);
    setSpec(null);
    setLiveSpectrum(null);
    setLiveMs(null);
    setError(null);
  }, [replaceCommitted]);

  // Close the backend session and committed run with the app.
  useEffect(() => () => discard(), [discard]);

  const update = useCallback((edit: (s: SimSpec) => SimSpec) => {
    setSpec((s) => (s ? edit(s) : s));
  }, []);

  const setTargetRt = useCallback((rt: number) => setTargetRtState(rt), []);

  const step = useCallback(
    (delta: number) => {
      if (!spec || viewedIndex === null) return;
      const next = Math.max(0, Math.min(spec.scanTimes.length - 1, viewedIndex + delta));
      setTargetRtState(spec.scanTimes[next]);
    },
    [spec, viewedIndex]
  );

  const accept = useCallback(async () => {
    const session = sessionRef.current;
    if (!spec || !session || committing) return;
    setCommitting(true);
    try {
      const { handle, provider } = await session.commit(spec);
      const tic = await provider.getTicTrace({ maxPoints: 4000 });
      if (sessionRef.current !== session) {
        // Discarded or restarted mid-commit: this run belongs to nothing any more.
        void closeDataset(handle).catch(() => {});
        return;
      }
      replaceCommitted({ handle, provider, tic, spec });
    } catch (err) {
      setError(errText(err));
    } finally {
      setCommitting(false);
    }
  }, [spec, committing, replaceCommitted]);

  return {
    spec,
    update,
    start,
    discard,
    setTargetRt,
    step,
    viewedIndex,
    liveSpectrum,
    liveMs,
    error,
    committed,
    committing,
    stale: spec !== null && committed?.spec !== spec,
    accept
  };
}

// ------------------------------------------------------------------ the panel

export function SimPanel({
  wb,
  runScanTimes,
  runName,
  onLoad,
  onSave,
  onClose
}: {
  wb: SimWorkbench;
  /** The open real run's MS1 scan times, to seed a spec that aligns scan for scan. */
  runScanTimes: readonly number[] | null;
  runName: string | null;
  onLoad: () => void;
  onSave: () => void;
  onClose: () => void;
}) {
  const { spec } = wb;
  const [selected, setSelected] = useState(0);
  // The noise block last switched off, so switching it back on restores the user's settings.
  const lastNoise = useRef<SimNoise | null>(null);

  const proteoformCount = spec?.proteoforms.length ?? 0;
  const sel = proteoformCount > 0 ? Math.min(selected, proteoformCount - 1) : null;
  const viewedRt =
    spec && wb.viewedIndex !== null ? spec.scanTimes[wb.viewedIndex] : null;

  const setProteoform = (i: number, edit: (p: SimProteoform) => SimProteoform) =>
    wb.update((s) => ({
      ...s,
      proteoforms: s.proteoforms.map((p, j) => (j === i ? edit(p) : p))
    }));

  const setNoise = (edit: (n: SimNoise) => SimNoise) =>
    wb.update((s) => (s.noise ? { ...s, noise: edit(s.noise) } : s));

  const rtRange = spec
    ? { min: spec.scanTimes[0], max: spec.scanTimes[spec.scanTimes.length - 1] }
    : { min: 0, max: 1 };

  return (
    <>
      <div style={headerStyle}>
        <strong>Simulation workbench</strong>
        <button onClick={onClose} style={closeStyle} type="button" aria-label="Close">
          ✕
        </button>
      </div>

      <div style={bodyStyle}>
        <div style={buttonRowStyle}>
          {runScanTimes && runScanTimes.length > 1 ? (
            <SmallButton
              onClick={() =>
                void wb.start(newSimSpec({ scanTimes: runScanTimes, sourceFile: runName ?? undefined }))
              }
              title="A new spec on this run's MS1 scan times, so scans align one to one"
            >
              New from run
            </SmallButton>
          ) : null}
          <SmallButton onClick={() => void wb.start(newSimSpec())}>New blank</SmallButton>
          <SmallButton onClick={onLoad}>Load…</SmallButton>
          {spec ? <SmallButton onClick={onSave}>Save…</SmallButton> : null}
          {spec ? <SmallButton onClick={wb.discard}>Discard</SmallButton> : null}
        </div>

        {!spec ? (
          <div style={noteStyle}>
            Start a spec to simulate. Edits re-render the viewed scan live (mirrored under the real
            spectrum); Accept simulates the whole run and draws its TIC.
          </div>
        ) : (
          <>
            <div style={wb.error ? errorStyle : noteStyle}>
              {wb.error
                ? wb.error
                : `${spec.scanTimes.length} scans · ${proteoformCount} proteoform${proteoformCount === 1 ? "" : "s"}${
                    spec.background && spec.background.length > 0
                      ? ` (+${spec.background.length} background)`
                      : ""
                  }${
                    wb.viewedIndex !== null && viewedRt !== null
                      ? ` · viewing scan ${wb.viewedIndex + 1} (RT ${viewedRt.toFixed(2)})`
                      : ""
                  }${wb.liveMs !== null ? ` · ${wb.liveMs.toFixed(0)} ms` : ""}`}
            </div>

            <Section title="Proteoforms">
              <div style={listStyle}>
                <table style={tableStyle}>
                  <thead>
                    <tr>
                      <th style={thStyle}>id</th>
                      <th style={thRightStyle}>mass</th>
                      <th style={thRightStyle}>RT μ</th>
                      <th style={thRightStyle}>abundance</th>
                    </tr>
                  </thead>
                  <tbody>
                    {spec.proteoforms.slice(0, PROTEOFORM_LIST_CAP).map((p, i) => (
                      <tr
                        key={i}
                        onClick={() => setSelected(i)}
                        style={i === sel ? rowSelectedStyle : rowStyle}
                      >
                        <td style={tdStyle}>{p.identifier ?? `#${i + 1}`}</td>
                        <td style={tdRightStyle}>{p.monoisotopicMass.toFixed(1)}</td>
                        <td style={tdRightStyle}>{p.rt.mu.toFixed(2)}</td>
                        <td style={tdRightStyle}>{p.abundance.toExponential(2)}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
                {proteoformCount > PROTEOFORM_LIST_CAP ? (
                  <div style={noteStyle}>
                    Listing the first {PROTEOFORM_LIST_CAP} of {proteoformCount}.
                  </div>
                ) : null}
              </div>
              <div style={buttonRowStyle}>
                <SmallButton
                  onClick={() => {
                    wb.update((s) => ({
                      ...s,
                      proteoforms: [
                        ...s.proteoforms,
                        defaultProteoform(viewedRt ?? rtRange.min, `P${s.proteoforms.length + 1}`)
                      ]
                    }));
                    setSelected(proteoformCount);
                  }}
                  title="A new proteoform eluting at the viewed scan's RT"
                >
                  + Add at viewed RT
                </SmallButton>
                {sel !== null ? (
                  <>
                    <SmallButton
                      onClick={() => {
                        wb.update((s) => {
                          const copy = { ...s.proteoforms[sel], identifier: `${s.proteoforms[sel].identifier ?? "P"}′` };
                          const proteoforms = [...s.proteoforms];
                          proteoforms.splice(sel + 1, 0, copy);
                          return { ...s, proteoforms };
                        });
                        setSelected(sel + 1);
                      }}
                    >
                      Duplicate
                    </SmallButton>
                    <SmallButton
                      onClick={() =>
                        wb.update((s) => ({
                          ...s,
                          proteoforms: s.proteoforms.filter((_, j) => j !== sel)
                        }))
                      }
                    >
                      Remove
                    </SmallButton>
                  </>
                ) : null}
              </div>

              {sel !== null ? (
                <ProteoformEditor
                  p={spec.proteoforms[sel]}
                  rtRange={rtRange}
                  chargeRange={spec.chargeRange}
                  onChange={(edit) => setProteoform(sel, edit)}
                />
              ) : null}
            </Section>

            <Section title="Run">
              {spec.peakWidthModel.kind === "Orbitrap" ? (
                <NumberField
                  label="Peak width k"
                  hint="σ = k·(m/z)^1.5"
                  value={spec.peakWidthModel.k}
                  slider={{ min: 1e-8, max: 1e-6, log: true }}
                  onChange={(k) => wb.update((s) => ({ ...s, peakWidthModel: { kind: "Orbitrap", k } }))}
                />
              ) : (
                <NumberField
                  label="Peak width σ"
                  value={spec.peakWidthModel.sigma}
                  slider={{ min: 1e-4, max: 0.1, log: true }}
                  onChange={(sigma) =>
                    wb.update((s) => ({ ...s, peakWidthModel: { kind: "Constant", sigma } }))
                  }
                />
              )}
              <div style={pairStyle}>
                <NumberField
                  label="Charge min"
                  value={spec.chargeRange.min}
                  step={1}
                  onChange={(min) =>
                    wb.update((s) => ({ ...s, chargeRange: { ...s.chargeRange, min: Math.round(min) } }))
                  }
                />
                <NumberField
                  label="Charge max"
                  value={spec.chargeRange.max}
                  step={1}
                  onChange={(max) =>
                    wb.update((s) => ({ ...s, chargeRange: { ...s.chargeRange, max: Math.round(max) } }))
                  }
                />
              </div>
              {!spec.noise ? (
                <div style={pairStyle}>
                  <NumberField
                    label="Relative floor"
                    hint="× brightest centroid"
                    value={spec.reduction?.relativeIntensityThreshold ?? 1e-4}
                    onChange={(relativeIntensityThreshold) =>
                      wb.update((s) => ({
                        ...s,
                        reduction: {
                          absoluteIntensityFloor: s.reduction?.absoluteIntensityFloor ?? 0,
                          relativeIntensityThreshold
                        }
                      }))
                    }
                  />
                  <NumberField
                    label="Absolute floor"
                    value={spec.reduction?.absoluteIntensityFloor ?? 0}
                    onChange={(absoluteIntensityFloor) =>
                      wb.update((s) => ({
                        ...s,
                        reduction: {
                          relativeIntensityThreshold: s.reduction?.relativeIntensityThreshold ?? 1e-4,
                          absoluteIntensityFloor
                        }
                      }))
                    }
                  />
                </div>
              ) : null}
            </Section>

            <Section
              title="Noise"
              action={
                <label style={checkLabelStyle}>
                  <input
                    type="checkbox"
                    checked={!!spec.noise}
                    onChange={(e) => {
                      const on = e.target.checked;
                      wb.update((s) => {
                        if (!on && s.noise) lastNoise.current = s.noise;
                        return { ...s, noise: on ? lastNoise.current ?? orbitrapNoiseDefaults(0) : undefined };
                      });
                    }}
                  />
                  enabled
                </label>
              }
            >
              {spec.noise ? (
                <NoiseEditor
                  noise={spec.noise}
                  onChange={setNoise}
                  onReset={() => setNoise((n) => orbitrapNoiseDefaults(n.seed))}
                />
              ) : (
                <div style={noteStyle}>
                  Clean analytical simulation, reduced against the run&apos;s brightest centroid.
                </div>
              )}
            </Section>
          </>
        )}
      </div>

      {spec ? (
        <div style={footerStyle}>
          <button
            type="button"
            style={wb.stale && !wb.committing ? acceptStyle : acceptDisabledStyle}
            disabled={!wb.stale || wb.committing}
            onClick={() => void wb.accept()}
          >
            {wb.committing ? "Simulating run…" : wb.stale ? "Accept — simulate run" : "Run up to date"}
          </button>
          <div style={noteStyle}>
            {wb.committing
              ? "The TIC updates when the run is simulated."
              : wb.committed && wb.stale
                ? "The TIC shows the last accepted run."
                : wb.committed
                  ? "The TIC matches the panel."
                  : "No run simulated yet: only the viewed scan is live."}
          </div>
        </div>
      ) : null}
    </>
  );
}

function ProteoformEditor({
  p,
  rtRange,
  chargeRange,
  onChange
}: {
  p: SimProteoform;
  rtRange: { min: number; max: number };
  chargeRange: { min: number; max: number };
  onChange: (edit: (p: SimProteoform) => SimProteoform) => void;
}) {
  return (
    <div style={editorStyle}>
      <TextField
        label="Identifier"
        value={p.identifier ?? ""}
        onChange={(identifier) => onChange((q) => ({ ...q, identifier: identifier || undefined }))}
      />
      <NumberField
        label="Monoisotopic mass (Da)"
        value={p.monoisotopicMass}
        step={1}
        slider={{ min: 1000, max: 80000 }}
        onChange={(monoisotopicMass) => onChange((q) => ({ ...q, monoisotopicMass }))}
      />
      <NumberField
        label="Abundance A"
        value={p.abundance}
        slider={{ min: 1e3, max: 1e12, log: true }}
        onChange={(abundance) => onChange((q) => ({ ...q, abundance }))}
      />
      <NumberField
        label="RT μ (min)"
        value={p.rt.mu}
        step={0.01}
        slider={{ min: rtRange.min, max: rtRange.max }}
        onChange={(mu) => onChange((q) => ({ ...q, rt: { ...q.rt, mu } }))}
      />
      <div style={pairStyle}>
        <NumberField
          label="RT σ"
          value={p.rt.sigma}
          step={0.01}
          slider={{ min: 0.01, max: 2 }}
          onChange={(sigma) => onChange((q) => ({ ...q, rt: { ...q.rt, sigma } }))}
        />
        <NumberField
          label="RT τ"
          value={p.rt.tau}
          step={0.01}
          slider={{ min: 0, max: 2 }}
          onChange={(tau) => onChange((q) => ({ ...q, rt: { ...q.rt, tau } }))}
        />
      </div>
      <div style={pairStyle}>
        <NumberField
          label="Charge μ_z"
          value={p.charge.muZ}
          step={0.1}
          slider={{ min: chargeRange.min, max: chargeRange.max }}
          onChange={(muZ) => onChange((q) => ({ ...q, charge: { ...q.charge, muZ } }))}
        />
        <NumberField
          label="Charge σ_z"
          value={p.charge.sigmaZ}
          step={0.1}
          slider={{ min: 0.2, max: 15 }}
          onChange={(sigmaZ) => onChange((q) => ({ ...q, charge: { ...q.charge, sigmaZ } }))}
        />
      </div>
    </div>
  );
}

function NoiseEditor({
  noise,
  onChange,
  onReset
}: {
  noise: SimNoise;
  onChange: (edit: (n: SimNoise) => SimNoise) => void;
  onReset: () => void;
}) {
  const floor = (edit: (f: SimNoise["floor"]) => SimNoise["floor"]) =>
    onChange((n) => ({ ...n, floor: edit(n.floor) }));
  const jitter = (edit: (j: SimNoise["jitter"]) => SimNoise["jitter"]) =>
    onChange((n) => ({ ...n, jitter: edit(n.jitter) }));
  const j = noise.jitter;
  return (
    <>
      <div style={pairStyle}>
        <NumberField
          label="Seed"
          value={noise.seed}
          step={1}
          onChange={(seed) => onChange((n) => ({ ...n, seed: Math.round(seed) }))}
        />
        <div style={{ alignSelf: "end" }}>
          <SmallButton onClick={onReset} title="The measured Orbitrap floor and jitter">
            Orbitrap defaults
          </SmallButton>
        </div>
      </div>
      <NumberField
        label="Noise level at reference m/z"
        value={noise.floor.noiseLevelAtReferenceMz}
        slider={{ min: 10, max: 1e6, log: true }}
        onChange={(v) => floor((f) => ({ ...f, noiseLevelAtReferenceMz: v }))}
      />
      <div style={pairStyle}>
        <NumberField
          label="Density scale"
          value={noise.floor.densityScale}
          step={0.05}
          slider={{ min: 0, max: 5 }}
          onChange={(v) => floor((f) => ({ ...f, densityScale: v }))}
        />
        <NumberField
          label="Level dispersion"
          value={noise.floor.levelDispersion}
          step={0.01}
          slider={{ min: 0, max: 1 }}
          onChange={(v) => floor((f) => ({ ...f, levelDispersion: v }))}
        />
      </div>
      <div style={pairStyle}>
        <NumberField
          label="S/N threshold"
          value={noise.floor.signalToNoise.threshold}
          step={0.01}
          slider={{ min: 0.5, max: 10 }}
          onChange={(v) =>
            floor((f) => ({ ...f, signalToNoise: { ...f.signalToNoise, threshold: v } }))
          }
        />
        <NumberField
          label="S/N mean excess"
          value={noise.floor.signalToNoise.meanExcess}
          step={0.01}
          slider={{ min: 0.1, max: 10 }}
          onChange={(v) =>
            floor((f) => ({ ...f, signalToNoise: { ...f.signalToNoise, meanExcess: v } }))
          }
        />
      </div>

      <label style={{ ...checkLabelStyle, marginTop: 8 }}>
        <input
          type="checkbox"
          checked={j.enabled}
          onChange={(e) => jitter((x) => ({ ...x, enabled: e.target.checked }))}
        />
        Peak jitter
      </label>
      {j.enabled ? (
        <>
          <div style={pairStyle}>
            <NumberField
              label="m/z ppm @ S/N 1"
              value={j.mzPpmAtUnitSignalToNoise}
              step={0.1}
              onChange={(v) => jitter((x) => ({ ...x, mzPpmAtUnitSignalToNoise: v }))}
            />
            <NumberField
              label="m/z ppm floor"
              value={j.mzPpmFloor}
              step={0.05}
              onChange={(v) => jitter((x) => ({ ...x, mzPpmFloor: v }))}
            />
          </div>
          <NumberField
            label="Common-mode m/z ppm"
            value={j.commonModeMzPpm}
            step={0.1}
            onChange={(v) => jitter((x) => ({ ...x, commonModeMzPpm: v }))}
          />
          <div style={pairStyle}>
            <NumberField
              label="log-I σ @ S/N 1"
              value={j.logIntensitySigmaAtUnitSignalToNoise}
              step={0.01}
              onChange={(v) => jitter((x) => ({ ...x, logIntensitySigmaAtUnitSignalToNoise: v }))}
            />
            <NumberField
              label="log-I σ floor"
              value={j.logIntensitySigmaFloor}
              step={0.01}
              onChange={(v) => jitter((x) => ({ ...x, logIntensitySigmaFloor: v }))}
            />
          </div>
          <NumberField
            label="Common-mode log-I σ"
            value={j.commonModeLogIntensitySigma}
            step={0.01}
            onChange={(v) => jitter((x) => ({ ...x, commonModeLogIntensitySigma: v }))}
          />
          <label style={checkLabelStyle}>
            <input
              type="checkbox"
              checked={j.dropPeaksFallingBelowFloor}
              onChange={(e) => jitter((x) => ({ ...x, dropPeaksFallingBelowFloor: e.target.checked }))}
            />
            Drop peaks jittered below the floor
          </label>
        </>
      ) : null}
    </>
  );
}

// ------------------------------------------------------------------ inputs

function Section({ title, action, children }: { title: string; action?: ReactNode; children: ReactNode }) {
  return (
    <section style={sectionStyle}>
      <div style={sectionHeaderStyle}>
        <span>{title}</span>
        {action}
      </div>
      {children}
    </section>
  );
}

function SmallButton({
  onClick,
  title,
  children
}: {
  onClick: () => void;
  title?: string;
  children: ReactNode;
}) {
  return (
    <button type="button" onClick={onClick} title={title} style={smallButtonStyle}>
      {children}
    </button>
  );
}

function TextField({
  label,
  value,
  onChange
}: {
  label: string;
  value: string;
  onChange: (v: string) => void;
}) {
  return (
    <label style={fieldStyle}>
      <span style={fieldLabelStyle}>{label}</span>
      <input style={inputStyle} value={value} onChange={(e) => onChange(e.target.value)} />
    </label>
  );
}

/**
 * A number input, optionally with a slider (linear or log). The text is a local draft: a value
 * that does not parse to a finite number is shown but never sent, so a half-typed "1e" does not
 * reach the simulator. Values that parse but the backend rejects (σ ≤ 0, min > max) come back as
 * the panel's error line.
 */
function NumberField({
  label,
  hint,
  value,
  step,
  slider,
  onChange
}: {
  label: string;
  hint?: string;
  value: number;
  step?: number;
  slider?: { min: number; max: number; log?: boolean };
  onChange: (v: number) => void;
}) {
  const [draft, setDraft] = useState(String(value));
  const focused = useRef(false);
  useEffect(() => {
    // Follow outside changes (slider, undo-by-load) unless the user is mid-edit here.
    if (!focused.current) setDraft(String(value));
  }, [value]);

  const SLIDER_STEPS = 1000;
  const toSlider = (v: number): number => {
    if (!slider) return 0;
    const { min, max, log } = slider;
    const f = log
      ? (Math.log(Math.max(v, min)) - Math.log(min)) / (Math.log(max) - Math.log(min))
      : (v - min) / (max - min);
    return Math.round(Math.max(0, Math.min(1, f)) * SLIDER_STEPS);
  };
  const fromSlider = (pos: number): number => {
    const { min, max, log } = slider!;
    const f = pos / SLIDER_STEPS;
    const v = log ? Math.exp(Math.log(min) + f * (Math.log(max) - Math.log(min))) : min + f * (max - min);
    // Three significant figures: slider positions are approximate by nature.
    return Number(v.toPrecision(3));
  };

  return (
    <label style={fieldStyle}>
      <span style={fieldLabelStyle}>
        {label}
        {hint ? <span style={hintStyle}> {hint}</span> : null}
      </span>
      <input
        style={inputStyle}
        type="number"
        step={step ?? "any"}
        value={draft}
        onFocus={() => (focused.current = true)}
        onBlur={() => {
          focused.current = false;
          setDraft(String(value));
        }}
        onChange={(e) => {
          setDraft(e.target.value);
          const v = Number(e.target.value);
          if (e.target.value.trim() !== "" && Number.isFinite(v)) onChange(v);
        }}
      />
      {slider ? (
        <input
          type="range"
          min={0}
          max={SLIDER_STEPS}
          value={toSlider(value)}
          onChange={(e) => onChange(fromSlider(Number(e.target.value)))}
          style={sliderStyle}
        />
      ) : null}
    </label>
  );
}

// ------------------------------------------------------------------ styles

const headerStyle: CSSProperties = {
  display: "flex",
  justifyContent: "space-between",
  alignItems: "center",
  padding: "8px 12px",
  borderBottom: "1px solid #d5deeb",
  fontSize: "0.85rem"
};
const closeStyle: CSSProperties = {
  border: "none",
  background: "transparent",
  cursor: "pointer",
  fontSize: "0.9rem"
};
const bodyStyle: CSSProperties = {
  overflow: "auto",
  minHeight: 0,
  flex: 1,
  padding: "8px 12px",
  fontSize: "0.78rem",
  color: "#24364d"
};
const footerStyle: CSSProperties = {
  borderTop: "1px solid #d5deeb",
  padding: "8px 12px",
  display: "flex",
  flexDirection: "column",
  gap: 4
};
const buttonRowStyle: CSSProperties = { display: "flex", flexWrap: "wrap", gap: 6, margin: "4px 0" };
const smallButtonStyle: CSSProperties = {
  border: "1px solid #cbd6e5",
  background: "#f4f7fb",
  borderRadius: 4,
  cursor: "pointer",
  padding: "2px 8px",
  fontSize: "0.76rem"
};
const acceptStyle: CSSProperties = {
  border: `1px solid ${SIM_COLOR}`,
  background: SIM_COLOR,
  color: "#ffffff",
  borderRadius: 4,
  cursor: "pointer",
  padding: "6px 10px",
  fontWeight: 600,
  fontSize: "0.82rem"
};
const acceptDisabledStyle: CSSProperties = {
  ...acceptStyle,
  background: "#e9edf3",
  borderColor: "#cbd6e5",
  color: "#7a8899",
  cursor: "default"
};
const noteStyle: CSSProperties = { color: "#5a6b80", margin: "4px 0", lineHeight: 1.4 };
const errorStyle: CSSProperties = { ...noteStyle, color: "#b02a37" };
const sectionStyle: CSSProperties = { borderTop: "1px solid #e3e9f2", marginTop: 8, paddingTop: 6 };
const sectionHeaderStyle: CSSProperties = {
  display: "flex",
  justifyContent: "space-between",
  alignItems: "center",
  fontWeight: 600,
  fontSize: "0.8rem",
  marginBottom: 4
};
const editorStyle: CSSProperties = {
  border: "1px solid #e3e9f2",
  borderRadius: 4,
  padding: "4px 8px",
  marginTop: 4,
  background: "#fbfcfe"
};
const pairStyle: CSSProperties = { display: "grid", gridTemplateColumns: "1fr 1fr", gap: 8 };
const fieldStyle: CSSProperties = { display: "flex", flexDirection: "column", gap: 2, margin: "4px 0" };
const fieldLabelStyle: CSSProperties = { color: "#3a4a60" };
const hintStyle: CSSProperties = { color: "#8796a8" };
const inputStyle: CSSProperties = {
  fontSize: "0.78rem",
  padding: "2px 4px",
  border: "1px solid #cbd6e5",
  borderRadius: 4,
  width: "100%",
  boxSizing: "border-box"
};
const sliderStyle: CSSProperties = { width: "100%", accentColor: SIM_COLOR };
const checkLabelStyle: CSSProperties = {
  display: "flex",
  alignItems: "center",
  gap: 4,
  fontWeight: 400,
  fontSize: "0.76rem"
};
const listStyle: CSSProperties = {
  maxHeight: 180,
  overflow: "auto",
  border: "1px solid #e3e9f2",
  borderRadius: 4
};
const tableStyle: CSSProperties = { width: "100%", borderCollapse: "collapse", fontSize: "0.76rem" };
const thStyle: CSSProperties = {
  textAlign: "left",
  padding: "3px 6px",
  position: "sticky",
  top: 0,
  background: "#f4f7fb",
  borderBottom: "1px solid #d5deeb"
};
const thRightStyle: CSSProperties = { ...thStyle, textAlign: "right" };
const rowStyle: CSSProperties = { cursor: "pointer", borderBottom: "1px solid #eef2f7" };
const rowSelectedStyle: CSSProperties = { ...rowStyle, background: "#f7e3dd" };
const tdStyle: CSSProperties = { padding: "2px 6px", whiteSpace: "nowrap" };
const tdRightStyle: CSSProperties = { ...tdStyle, textAlign: "right" };
