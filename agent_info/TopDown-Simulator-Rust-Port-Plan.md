# Top-Down Simulator → Rust port & ProteoStar integration

**Status:** Plan of record · **Owner:** Alex · **Created:** 2026-10-20

## Vision — an interactive simulation workbench

The end goal is not "simulate once and look at it." It is a **workbench**: a real MS1 run and a
simulated MS1 run shown **together**, a panel that exposes the simulation parameters, and edits to
those parameters that re-render the simulated side so you can *see what each knob does*. Port the
mzLib `TopDownSimulator` (C#) MS1 simulation into Rust, run it inside ProteoStar, and drive it from
a live parameter UI compared against real data.

Three pillars:
1. **Side-by-side real + simulated MS1** (spectrum and TIC, same RT/scan).
2. **A parameter panel** exposing the simulation knobs (per-proteoform model params, run-level width
   and charge range, and the full noise controls).
3. **Live feedback while editing.** See the **Interaction model** below for exactly how live.

This doc is the plan; the input-format contract has its own doc:
[`TopDown-SimSpec-Contract-v1.md`](./TopDown-SimSpec-Contract-v1.md).

## Decisions locked

| Area | Decision | Why |
|---|---|---|
| Simulator input | **Seed from a C# parameter file, then edit live in-app.** The C# fitting/extraction pipeline emits a self-describing *simulation spec* (ideally fit from the same real run being compared against); the UI loads it, owns a live editable copy, and `SimSpec` doubles as the save format. Starting from defaults / scratch is also supported. | The fitting/extraction front end (`Extraction/*`, `Fitting/*`) is the most mzLib-coupled part — keep it in C#. But the workbench needs in-app authoring regardless, so the file is a seed, not the only source. |
| First scope | **MS1 only, but with full noise.** | The MS1 forward model is the lightly-coupled core. Noise is in v1 scope because the noise knobs are live-editable (see Interaction model). MS2 reuses mzLib `Omics` fragment chemistry and is deferred. |
| Sim → viewer | **In-memory, two-tier (see Interaction model); optional mzML export later.** | A file round-trip can't be real-time, so the live loop is forward-model → plot, in memory. Export is additive and needs an mzML writer that does not exist yet. |
| Interaction | **Live recompute of the viewed scan only; full dataset rebuilt on explicit Accept/Submit.** | The scan spectrum is where per-scan effects (envelope, charge, width, noise) are visible; recomputing one scan is cheap enough to be continuous. Rebuilding the whole indexed run on every keystroke is not, so it is gated behind a commit. |
| C# sidecar | **Not in Phase 1.** Introduced in the MS2 phase for `Omics` fragmentation; reconsidered for mzML export. | Phase-1 input is an offline file handoff, not a live process. ProteoStar already bundles a .NET 8 runtime (thermorawfilereader), so a sidecar adds no new runtime dependency when we do want one. |

## Interaction model (how "live" works)

Two tiers, split by an explicit commit:

- **Live path — the viewed scan only.** Every parameter edit (slider drag, field change)
  recomputes **just the currently-displayed scan's spectrum**, straight from the forward model +
  centroider (+ that one scan's noise, stable seed), **bypassing the peak index**. One scan's m/z
  span is tiny, so this is continuous/real-time even with noise on. RT-profile edits are visible
  here too, because `g_p(t)` scales the envelope heights at the viewed scan's RT. The TIC/XIC and
  the whole-run overlay do **not** update on the live path.
- **Commit path — Accept/Submit.** On an explicit action, the **full** run is re-simulated (all
  scans), indexed (`PeakIndexingEngine`), and committed as the simulated dataset in its slot — at
  which point TIC/XIC and the full side-by-side overlay update. This is the only time the whole
  dataset and index are rebuilt.

**Rust command surface this implies:**
- `simulate_scan(spec, scanIndex) -> Spectrum` — the live path; no handle, no index, no state.
- `simulate_dataset(spec) -> OpenResult` (handle) — the commit path; builds `Vec<Scan>` +
  `PeakIndexingEngine` and registers/updates a dataset handle, so the sim then behaves like any
  real dataset through the existing `DatasetProvider`.

## Reality check on ProteoStar (verified against code, 2026-10-20)

`agent_info/MsViewer_Architecture.md` describes a 2-slot compare design that is **not built**.
What actually exists:

- **No `viewer-state` package, no `datasetSlots`, no second plot.** The viewer is strictly
  single-dataset (`apps/desktop/src/App.tsx` holds one `provider`/`handle`/`spectrum`). The only
  trace of the 2-slot idea is an **unused** `SlotIndex = 0 | 1` type in `packages/plot-adapter`.
- The clean seam is the **`DatasetProvider`** interface in `apps/desktop/src/contract.ts`
  (`:215-244`). The whole UI depends only on it; the Tauri transport
  (`apps/desktop/src/tauri-dataset-provider.ts`) is one implementation. A new provider drops in
  with no UI rewrite.
- The Rust backend is **hardcoded to open-a-file-from-path**: `open_dataset`
  (`src-tauri/src/commands.rs:316`) → `read_ms1_scans(path)` → `PeakIndexingEngine`. There is **no
  data-source trait**; the per-handle `reader` (`RandomAccessMs1Reader`) has only `open(path)`.
- **No mzML writer** — `mzdata` is compiled read-only (`crates/flashlfq-core/Cargo.toml`).

Consequence: side-by-side is **work to build**, not a feature to plug into. Two layers of work —
(1) make a simulated run display at all, and (2) make the UI show two runs together.

## What the mzLib simulator looks like (port surface)

Project: `mzLib/TopDownSimulator/` (net10.0). Lightly-coupled numeric core to port first:

| Piece | mzLib source | Notes |
|---|---|---|
| EMG elution `g_p(t)` | `Model/EmgProfile.cs` | unit-area; inline Abramowitz–Stegun erf |
| Charge dist `f_p(z)` | `Model/ChargeStateDistribution.cs` | Gaussian density in z — **not** renormalized over `[min,max]` (ported as-is) |
| Peak width `σ_m(mz)` | `Model/IPeakWidthModel.cs` | `Orbitrap(k)`: σ=K·(m/z)^1.5; `Constant(σ)` |
| Isotope envelope `φ` | `Model/IsotopeEnvelopeKernel.cs` | sum of unit-area Gaussians on isotopologue m/z — **uses mzLib `Averagine` tables** (see contract) |
| Forward model | `Model/ForwardModel.cs` | `I(t,mz)=Σ_p A_p·g_p(t)·Σ_z f_p(z)·φ(...)` |
| Centroider | `Simulation/ProfileCentroider.cs` | one centroid per local max; Newton polish on S'(x)=0 |
| Reduction | `Simulation/SimulatedScanReducer.cs` | relative+absolute intensity floor |
| Noise | `Noise/*` | empirical; Poisson density, Box–Muller; SplitMix64 derives a per-scan seed for .NET's legacy seeded `System.Random` (Knuth subtractive) — both ported verbatim |

**Costly edges, deferred:** MS2 (`Ms2/*`, depends on `Omics` fragment chemistry) and the
fitting/extraction front end (`Extraction/*`, `Fitting/*`). Both stay in C# for now.

## Phases

### Phase 0 — Parameter-file contract ✅ *(done 2026-10-20)*
Define and implement the **simulation spec**: a self-describing JSON *sufficient* for Rust to
reproduce an MS1 run deterministically (extends `.groundtruth.tsv`, which omits the scan-time axis
and noise config), plus export the **averagine table** so the isotope math matches.
See [`TopDown-SimSpec-Contract-v1.md`](./TopDown-SimSpec-Contract-v1.md).
**Delivered:** `mzLib/TopDownSimulator/Simulation/SimulationSpec.cs` (emitter/reader);
`mzLib/Test/TopDownSimulator/SimulationSpecTests.cs` (4 round-trip/guard tests, green) +
`SimSpecGoldenGenerator.cs` (`[Explicit]`); the data-free golden in
`mzLib/Test/TopDownSimulator/SimSpecFixtures/` (`golden-v1.simspec.json`, `golden-v1.scans.tsv`,
`averagine.v1.json`). **Phase 0.1 ✅** — the `SimNoise` block is expanded to the full live-editable
set (`SimNoiseFloor` / `SimSignalToNoise` / `SimNoiseJitter`), round-tripping to the real
`NoiseFloorModel` / `SignalToNoiseDistribution` / `PeakJitterModel`, with tests.

### Phase 1 — `crates/sim-core` (pure-Rust MS1 simulator, with noise) ✅ *(done 2026-10-20)*
Port the numeric core: EMG / charge / width / envelope kernel / forward model / centroider, **and
the `Noise/*` subsystem** (noise is v1 because its knobs are live). Input: the Phase-0 spec. Output:
in-memory scans in flashlfq-core's `Scan` type. Expose both granularities the interaction model
needs as library calls: **(a) one scan** (`simulate_scan(spec, scanIndex)`) and **(b) the full run**
(`simulate_run(spec)`).

**Delivered:** `crates/sim-core` — `spec` (SimSpec v1 serde + validation), `averagine` (the
exported table, embedded), `model`, `centroider` (`ProfileCentroider`), `noise` (floor, jitter,
collapse, SplitMix64 + .NET `Random`), and `lib.rs` (`Simulator`, `simulate_scan`, `simulate_run`
= C# `PrepareMs1`: centroid → reduce → noise). `Simulator::centroid_scan`/`centroid_run` stop
after centroiding (= C# `SimulateCentroided`).
**Parity — bit-exact, every peak** (`tests/golden_parity.rs`): `golden-v1` (signal, 44 386 peaks)
plus three C#-generated `PrepareMs1` goldens in `parity/golden/` — `reduced-v1` (8 476),
`noise-jitter-v1` (26 291), `noise-nojitter-v1` (24 290). The live and commit paths agree scan for
scan. Bit-exactness relies on the platform libm (`exp`/`ln`/`cos`/`powf`) matching .NET's; the
assertions use the contract tolerances, so a libm difference elsewhere degrades to ~ulp, not failure.
**Timing** (`examples/simulate_timing.rs`, 1000 proteoforms, z 3–60, full-density noise):
`Simulator::new` ~200 ms (first call parses the 3 MB averagine table), `simulate_scan` ~300 ms,
`simulate_run` ~50 s for 3000 scans (105 M peaks). Centroiding dominates the live path: ~97k
unreduced peaks per scan, most below the floor. **Phase-3 lever:** prune components that cannot
reach the detection floor before centroiding, and keep one `Simulator` alive across edits that
don't touch the proteoform list — trading bit-exactness for speed, so gate it behind the goldens'
tolerance.
**Gotcha for Phase 2:** without noise, the reduction floor is relative to the brightest peak of the
*whole run*, so `simulate_scan` on a fresh `Simulator` centroids every scan once to find it. With
noise on (the realistic case) the floor is the noise model's and the live path touches one scan.

### Phase 2 — Wire both paths through Tauri; sim-as-dataset
`src-tauri` gets the two commands from the Interaction model:
- `simulate_scan(spec, scanIndex) -> Spectrum` — stateless live path; **no index, no handle**.
- `simulate_dataset(spec) -> OpenResult` — the commit path; builds `IndexedData` directly from the
  simulated scans via `PeakIndexingEngine::index_peaks` and registers a handle.

The commit path needs a small source abstraction so a handle's data can be in-memory (today only
`RandomAccessMs1Reader::open(path)` exists; `get_spectrum_at_rt`/`get_ms2_spectrum` read through it).
Result: after Accept, a simulated run behaves like any real dataset through the existing
`DatasetProvider`.

### Phase 3 — The workbench UI *(net-new; does not exist today)*
The payoff phase. Three parts:
1. **Two-slot side-by-side.** Add slot state to `App.tsx` (a second `provider`/`handle`) and render
   real vs simulated via the dormant `slotIndex: 1` plumbing in `plot-adapter` or a second pane. The
   decoy-overlay layer (`App.tsx:153-155, 462-504`) is the nearest working pattern.
2. **Parameter panel.** Edit per-proteoform model params (mass, abundance, RT μ/σ/τ, charge
   μ_z/σ_z), run-level width `k` + charge range, **full noise controls**, and add/remove
   proteoforms. The panel's state is the live `SimSpec`; it loads from / saves to the Phase-0 file.
3. **Live loop.** On every edit, call `simulate_scan` for the viewed scan and repaint its spectrum
   overlay — real-time. An **Accept/Submit** button calls `simulate_dataset` and swaps the committed
   simulated dataset into its slot, updating TIC/XIC and the whole-run overlay.

### Phase 4 — *(optional)* MS2 + mzML export
MS2 via a C# `mzlib-helper` sidecar (batched JSON over stdin/stdout, reusing the bundled .NET
runtime) for `Omics` fragmentation. mzML export via the `mzdata` writer feature or the same sidecar;
decided when we get there.

Phases 1–3 deliver the stated goal. Phase 4 is additive.

## Open risks

1. **Averagine parity.** The envelope kernel's shape comes from mzLib `Averagine`; Rust must use the
   identical table or envelopes drift. Handled by exporting the table in Phase 0.
2. **Reader abstraction (Phase 2).** `get_spectrum_at_rt`/`get_ms2_spectrum` read through the
   per-handle `reader`; an in-memory source must satisfy the same surface. Minimal new trait/enum.
3. ~~**Noise bit-reproducibility.**~~ Resolved in Phase 1: the RNG was ported verbatim and noise
   output is bit-exact against C#.
4. **C# "jitter disabled" ≠ contract.** `SimNoiseJitter.ToModel()` returns null when disabled, and
   `PrepareMs1`/`WriteMzml` replace a null jitter with `PeakJitterModel.Orbitrap` — so in C#,
   `enabled: false` currently means *default* jitter. Rust follows the contract (no jitter). Fix on
   the C# side by making `ToModel()` return `PeakJitterModel.None` when disabled.
