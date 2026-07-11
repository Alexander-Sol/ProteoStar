# ProteoStar

A single Cargo + Bun workspace that holds a pure-Rust mass-spec quantification engine and both
of its frontends:

- **`crates/flashlfq-core`** — a pure-Rust port of the [mzLib](https://github.com/smith-chem-wisc/mzLib)
  **FlashLFQ** label-free quantification core, plus the untargeted feature detector and the MS
  reader/peak-index/XIC pipeline.
- **`crates/flashlfq-py`** — PyO3/maturin Python bindings over the core.
- **`apps/desktop`** — the **MsViewer** Tauri desktop app (TypeScript/Vite UI + a Rust backend
  crate at `apps/desktop/src-tauri`) that reads raw data directly through the core and overlays
  detected features.
- **`packages/*`** — the shared TypeScript packages (`imsp-core`, `viewer-state`, `plot-adapter`,
  `ui`) the desktop UI is built from.

This repo is the merge of the former `flashlfq-rust` (engine + Python bindings) and `MsBrowser`
(the desktop viewer). Both frontends build against the **one** `flashlfq-core` — the desktop
backend depends on it by path, so there is no longer a vendored copy to keep in sync.

## Layout

```
ProteoStar/
├── Cargo.toml                 [workspace] root: crates/* + apps/desktop/src-tauri
├── Cargo.lock                 committed (ships binaries; carries the thermo 0.7.0 pin)
├── crates/
│   ├── flashlfq-core/         pure-Rust core (algorithms, detector, unit + parity tests)
│   │   ├── src/  examples/  tests/
│   │   └── parity/            committed golden .tsv gates + C# golden generator
│   └── flashlfq-py/           PyO3 bindings + Python smoke tests
├── apps/
│   └── desktop/               MsViewer: Vite/React UI …
│       └── src-tauri/         … and its Tauri (Rust) backend — a workspace member
├── packages/                  imsp-core · viewer-state · plot-adapter · ui
├── fixtures/                  small .imsp fixtures for the TS unit tests
├── agent_info/                design docs, architecture notes, task logs
├── package.json               Bun workspaces root (name: proteostar)
└── tsconfig.*.json  vitest.config.ts
```

## Prerequisites

- **Rust** — toolchain 1.94+.
- **.NET 8 runtime** — required to read Thermo `.raw` (the `thermorawfilereader` crate hosts a
  self-contained .NET 8 runtime); mzML stays pure-Rust. Only the .NET **SDK** is needed to
  regenerate the C# goldens.
- **Bun** — for the desktop UI and the TS packages.
- **maturin** + a Python 3.9+ interpreter — only for the Python bindings.

## Build & test

**Engine (Rust):**

```powershell
cargo test -p flashlfq-core        # unit tests + L1–L6 / P2 / MBR parity gates
```

The live parity gates read fixtures from `<MZLIB_DIR>/mzLib/Test/FlashLFQ/TestData/`; the
committed goldens under `crates/flashlfq-core/parity/golden/` let core development and the
golden-diff gates run offline. `MZLIB_DIR` defaults to `F:\mzLib`; override it before running the
live gates or the C# golden generator:

```powershell
$env:MZLIB_DIR = "D:\src\mzLib"
```

**Desktop app (Tauri):**

```powershell
bun install
cd apps/desktop
bun run tauri dev          # or: bun run tauri build
```

A bare `cargo build` / `cargo test` at the root builds the engine + the Tauri backend
(`default-members`); the Python wheel is built separately via maturin.

**Python bindings (maturin):**

```powershell
cd crates/flashlfq-py
maturin develop
python quant_smoke_test.py        # and mbr_/pep_/xic_ smoke tests
```

Pinned stack: pyo3 0.23.5 · numpy 0.23 · arrow 54.x · mzdata 0.65 · thermorawfilereader 0.7.0.

**TypeScript unit tests:**

```powershell
bun run test        # vitest over packages/*/tests and apps/desktop
```

## Where the design lives

- **`agent_info/MsViewer_Architecture.md`** — the desktop viewer architecture (raw data read in
  Rust; the `DatasetProvider` contract).
- **`agent_info/FeatureFinder-Integration.md`** — how the detector is wired into the app.
- **`agent_info/FlashLFQ-Rust-Rewrite-Feasibility.md`** — the original engine porting strategy.
- **`agent_info/ProteoStar-Merge-Plan.md`** — the plan that produced this repo.

## Relationship to mzLib

The engine remains a separate project from mzLib. It **references** an external mzLib checkout
(via `MZLIB_DIR`) only for parity fixtures and golden regeneration; each ported module documents
the exact mzLib source it mirrors, but mzLib is developed independently.
