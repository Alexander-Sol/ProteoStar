# Top-Down Simulator workbench — handoff

**As of:** 2026-10-20 · **Branch:** `sim-workbench-design` (ProteoStar) · read this first, then the
plan + contract docs below.

## Goal (the vision)

An **interactive simulation workbench** inside ProteoStar: a real MS1 run and a simulated MS1 run
shown **side by side**, a **parameter panel** exposing the simulation knobs, and edits that
re-render the simulated side **live**. Interaction rule (locked with the user):

- **Live = the viewed scan only.** Every parameter edit recomputes just the currently-displayed
  scan's spectrum (fast, index-free). The spectrum panel is the live surface.
- **Accept/Submit = the full run.** TIC/XIC and the whole-run overlay rebuild only on an explicit
  Accept. The TIC is allowed to lag until then.
- **Full noise control** is in v1. Real/sim scan alignment: by index when the sim is seeded from
  the real run, else by RT.

## Read these (design of record)

- `agent_info/TopDown-Simulator-Rust-Port-Plan.md` — phases + locked decisions + interaction model.
- `agent_info/TopDown-SimSpec-Contract-v1.md` — the SimSpec v1 input contract (schema, C# field
  mapping, averagine companion, golden).

## Where things stand

| Phase | State | Where |
|---|---|---|
| 0 — SimSpec contract + emitter + data-free golden | ✅ done | mzLib: `mzLib/TopDownSimulator/Simulation/SimulationSpec.cs`, tests + `SimSpecGoldenGenerator.cs` + `SimSpecFixtures/` in `mzLib/Test/TopDownSimulator/` (mzLib commit `cef4ee75` on branch `TopDownSimClean`) |
| 0.1 — full live-editable noise block | ✅ done | same `SimulationSpec.cs` (`SimNoise`/`SimNoiseFloor`/`SimSignalToNoise`/`SimNoiseJitter`) |
| 1 — `crates/sim-core` (pure-Rust MS1 sim, incl. noise) | ✅ done (another session) | ProteoStar commit `fab92ca`. `simulate_scan` + `simulate_run` → `flashlfq_core::Scan`; 19 tests green (bit-exact signal, statistical noise) |
| 2 — Tauri commands + sim-as-dataset + TS bindings | ✅ done | ProteoStar commits `7dcb659` (backend) + `3827597` (frontend) |
| 2.5 — persistent live session (`Simulator::update`, per-scan max cache, session commands, `SimSession`) | ✅ done, **uncommitted** | `agent_info/TopDown-Simulator-Live-Session-Plan.md` |
| 3 — the workbench UI | ✅ v1 done | `apps/desktop/src/SimWorkbench.tsx` (hook + panel), wired in `App.tsx`; mirror traces in `packages/plot-adapter`; E2E `e2e/tauri/sim-workbench.spec.ts` | `apps/desktop/src/App.tsx` + `packages/plot-adapter` + a new parameter panel |

Branch `sim-workbench-design` history: `e427c65` (docs) → `fab92ca` (sim-core) → `7dcb659`
(backend) → `3827597` (frontend). **Nothing is pushed** — all commits are local.

## The API Phase 3 builds on (already shipped in Phase 2)

Frontend (`apps/desktop/src/tauri-dataset-provider.ts`):
- `openSimulatedDataset(spec: SimSpec) → { handle, provider }` — commit path. `provider` is the same
  `DatasetProvider` a real run has, so a committed sim drops into a viewer slot unchanged.
- `simulateScan(spec: SimSpec, scanIndex, opts?) → Spectrum` — live path, stateless, one scan.
- `SimSpec` type + `orbitrapNoiseDefaults()` in `apps/desktop/src/sim-spec.ts`.

Backend (`apps/desktop/src-tauri/src/commands.rs`): `open_simulated_dataset`, `simulate_scan`;
`OpenDataset.reader` is now `Option` (None = simulated, served from the in-memory index).

## Phase 3 — concrete next steps

The 2-slot compare does **not** exist yet (the old `MsViewer_Architecture.md` described it but it
was never built; only an unused `SlotIndex = 0|1` type survives in `packages/plot-adapter`). So:

1. **Second slot in `App.tsx`.** Today it holds one `provider`/`handle`/`spectrum` as flat state
   (~lines 136–145). Add a second dataset (slot 1 = simulated). Render both via the dormant
   `slotIndex` plumbing in `plot-adapter` (`types.ts` `SlotIndex`, `plots.tsx`) — overlay a second
   trace on the existing `<TicPlot>`/`<SpectrumPlot>`, or a second pane. Nearest working pattern for
   a second data layer: the decoy overlay (`App.tsx` ~153–155, 462–504).
2. **Parameter panel.** Edit a live `SimSpec`: per-proteoform mass/abundance/RT μ,σ,τ/charge
   μ_z,σ_z; run-level width `k` + charge range; full noise (`orbitrapNoiseDefaults()` is the start);
   add/remove proteoforms. Panel state = the live `SimSpec`; load from / save to a `.simspec.json`.
3. **Live loop.** On every edit → `simulateScan(spec, viewedScanIndex)` → repaint the spectrum
   overlay. **Accept** → `openSimulatedDataset(spec)` → swap the committed sim into slot 1 (close the
   previous sim handle) → TIC/XIC/overlay update.

## Phase 3 v1 — what was built (2026-10-20)

- **Layout: a mirror plot.** The real run is drawn upright and the simulation hangs below zero on
  the same axes. The live simulated scan sits under the real spectrum, and the last accepted run's
  TIC under the real TIC. `mirrored` is a trace flag in `plot-adapter`; hover shows the true
  intensity. `fitSpectrumTracesYRange` sizes the y-axis for both halves.
- **Panel** (toolbar "Simulate…", a right drawer):
  - Starting points: New from run (the run's MS1 scan times, so scans align one to one), New
    blank (60 min at 1.2 s), Load… and Save… of `.simspec.json` (backend `read_sim_spec` /
    `write_sim_spec`), and Discard.
  - Proteoforms: a list with add at the viewed RT, duplicate and remove; an editor for mass,
    abundance (log slider), RT μ/σ/τ and charge μ_z/σ_z.
  - Run: `k`, the charge range, and the reduction floors when noise is off.
  - Noise: an on/off toggle that remembers the last settings, plus every floor, S/N and jitter
    knob and Orbitrap defaults.
  - Accept, with an "unaccepted" state while the panel differs from the committed run.
- **Live loop:** `useSimWorkbench` sends every edit through `SimSession.scan` (latest-wins) for
  the scan nearest the target RT. That RT is the real spectrum's, a TIC click, or ←/→ when no
  real scan is shown.
- **Verified:** the E2E spec passes against the real app. With noise the live scan has 17k
  peaks; with noise off it has 499; doubling the abundance doubles the peaks; Accept draws the
  simulated TIC; and the sim follows arrow-key steps. The live round trip is ~9 ms on the small
  yeast run. The full Tauri suite passes 11/12. `feature-overlay` (decoy envelope `shapes`) was
  already stale: it predates this work, since overlays moved to peak highlights in `1134c71`.

**Not in v1 (candidates):** seeding proteoforms from loaded features (abundance calibration is
open), editing background proteoforms in the panel (they are kept and simulated, just not
editable), and a list past 200 proteoforms (the editor only reaches the first 200).

## Gotchas / open items for Phase 3

- **TS `format` union.** The backend reports `format: "simulated"`, but
  `contract.ts` types `DatasetMetadata.format` as `"mzml" | "thermo_raw"`. Extend that union to
  `"simulated"` when the UI starts reading it.
- **Live-path perf: solved by the session (Phase 2.5).** Use `SimSession.scan(spec, idx)` (latest-wins,
  no debounce) and `SimSession.commit(spec)` for Accept. The stateless `simulateScan` still rebuilds
  everything per call.
- **Lane discipline: lifted.** The other (optimizer) session has ended, and this session now owns
  `crates/sim-core` as well.

## Verify the current state

```powershell
cargo test -p sim-core --manifest-path F:\ProteoStar\Cargo.toml     # 19 green
cargo check -p msviewer-tauri --manifest-path F:\ProteoStar\Cargo.toml
cd F:\ProteoStar\apps\desktop; bun run typecheck
```

(Windows / PowerShell box. `cargo`/`bun` on PATH. mzLib checkout at `C:\Users\Alex\source\repos\mzlib`,
branch `TopDownSimClean`, holds the C# `SimulationSpec.cs` + golden generator.)
