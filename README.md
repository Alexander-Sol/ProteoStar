# FlashLFQ → Rust

A pure-Rust port of the [mzLib](https://github.com/smith-chem-wisc/mzLib) **FlashLFQ**
label-free quantification core, with PyO3/maturin Python bindings.

This is a **separate project** from mzLib. It does not vendor or build mzLib; instead it
**references an external mzLib checkout** for two purposes only:

1. **Parity-gate test data** — the FlashLFQ golden fixtures under
   `mzLib/Test/FlashLFQ/TestData/` (psmtsv + mzML/raw files).
2. **Golden regeneration** — the C# golden generator (`parity/csharp_golden/`) references
   `mzLib/FlashLFQ/FlashLFQ.csproj` to run the real `FlashLfqEngine` and emit the `.tsv`
   goldens the Rust parity tests diff against.

The generated goldens (`rust/flashlfq-core/parity/golden/*.tsv`) are committed, so **core
development and the parity gates run offline** against them. You only need mzLib present to
*regenerate* goldens or to run the live parity tests that read raw fixtures.

## Locating mzLib: `MZLIB_DIR`

Everything that needs mzLib resolves its location through the **`MZLIB_DIR`** environment
variable, which points at the mzLib repository root (the folder containing the inner
`mzLib/` source tree, i.e. `<MZLIB_DIR>/mzLib/FlashLFQ`, `<MZLIB_DIR>/mzLib/Test`, …).

- **Default when unset:** `F:\mzLib` (the standard checkout on the dev machine).
- **Override:** set `MZLIB_DIR` before running cargo tests, the Python smoke tests, or the
  C# golden generator. Example (PowerShell):

  ```powershell
  $env:MZLIB_DIR = "D:\src\mzLib"
  ```

This single knob is honored by:

| Consumer            | Mechanism                                                            |
|---------------------|---------------------------------------------------------------------|
| Rust tests/examples | `flashlfq_core::mzlib_dir()` / `flashlfq_core::mzlib_test_data(rel)` |
| Python dumpers      | `os.environ.get("MZLIB_DIR", r"F:\mzLib")`                           |
| Python smoke tests  | same                                                                |
| C# golden generator | `$(MZLIB_DIR)` MSBuild property + `Environment.GetEnvironmentVariable` |

## Layout

```
flashlfq-rust/
├── rust/
│   ├── flashlfq-core/         pure-Rust core (algorithms, unit + parity tests)
│   │   ├── src/               the port
│   │   ├── tests/             L1–L6 / P2 / MBR parity gates (read MZLIB_DIR test data)
│   │   ├── examples/
│   │   └── parity/
│   │       ├── golden/        committed golden .tsv files (diffed by the gates)
│   │       ├── dump_*.py      Python golden dumpers (L0–L4)
│   │       └── csharp_golden/ C# real-engine golden generator (L5/L6/MBR)
│   └── flashlfq-py/           PyO3 bindings + Python smoke tests
├── agent_info/                design docs (feasibility study, IMSP design/plan)
├── PLAN.md                    durable task list / progress log for the port
├── ralph-loop.ps1            automation loop that drives PLAN.md task-by-task
└── README.md
```

## Prerequisites

- **Rust** — toolchain 1.94+ (uses `f64::round_ties_even`, stable since 1.77).
- **.NET SDK** — only for regenerating the C# goldens (net8.0). Not needed for normal dev.
- **Python 3.13 venv** — only for the Python bindings and golden dumpers (see below).

## Build & test (Rust)

```powershell
Set-Location F:\flashlfq-rust\rust
# Optional: $env:MZLIB_DIR = "D:\src\mzLib"   # if mzLib is not at F:\mzLib
cargo build
cargo test -p flashlfq-core        # unit tests + L1–L6 / P2 / MBR parity gates
```

The parity gates read fixtures from `<MZLIB_DIR>/mzLib/Test/FlashLFQ/TestData/`; if that
folder is missing they fail with a path error pointing you at `MZLIB_DIR`.

## Python bindings (maturin)

The smoke-test virtualenv (`rust/.venv/`) is **not committed** — recreate it once:

```powershell
py -V:3.13 -m venv F:\flashlfq-rust\rust\.venv
F:\flashlfq-rust\rust\.venv\Scripts\python -m pip install maturin numpy pyarrow scikit-learn matplotlib
```

Then build and exercise the bindings:

```powershell
Set-Location F:\flashlfq-rust\rust\flashlfq-py
F:\flashlfq-rust\rust\.venv\Scripts\Activate.ps1
maturin develop
python quant_smoke_test.py        # (and mbr_/pep_/xic_ smoke tests)
```

Pinned stack: pyo3 0.23.5 · numpy 0.23 · arrow 54.x · mzdata 0.65 (pure-Rust zlib).

## Regenerating goldens

L0–L4 (Python replicas of mzLib chemistry):

```powershell
F:\flashlfq-rust\rust\.venv\Scripts\python F:\flashlfq-rust\rust\flashlfq-core\parity\dump_periodic_table.py
# ...and dump_l1/l2/l3/l4_*.py
```

L5/L6/MBR (real C# `FlashLfqEngine`, needs the .NET SDK + `MZLIB_DIR`):

```powershell
# Optional: $env:MZLIB_DIR = "D:\src\mzLib"
Set-Location F:\flashlfq-rust\rust\flashlfq-core\parity\csharp_golden
dotnet run -c Release
```

The generator writes into this project's `parity/golden/`, reading test data from `MZLIB_DIR`.

## Where the design lives

- **`PLAN.md`** — the phase-by-phase task list and a detailed running log of every ported
  piece (the durable state the `ralph-loop.ps1` automation reads and updates).
- **`agent_info/FlashLFQ-Rust-Rewrite-Feasibility.md`** — the feasibility study and porting
  strategy (isotope-distribution port plan, layered parity harness, Python binding layer).
- **`agent_info/Simulated-IMSP-*.md`** — the simulated-IMSP API/metadata design and plan.

## Relationship to mzLib going forward

The two projects stay separate. This port continues to **reference** mzLib (via `MZLIB_DIR`)
for fixtures and golden regeneration, and each ported module documents the exact mzLib source
(`FlashLfqEngine.cs`, `Chemistry/*.cs`, …) it mirrors, but mzLib is developed independently in
its own repository.
