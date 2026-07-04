# FlashLFQ → Rust Rewrite: Feasibility Analysis & Plan

**Date:** 2026-06-25
**Author:** Claude (Opus 4.8) + Alex
**Status:** Feasibility approved; scope decided; no code written yet.

## Goal

Rewrite the core peak-tracing functionality of mzLib's **FlashLFQ** and **Readers**
projects into **Rust**, exposed to **Python** via bindings. The pipeline:

1. Read identification results (starting with `.psmtsv`).
2. Read raw data files (`.mzML` first; `.raw` later).
3. Detect/trace isotopic-envelope peak traces (XICs) in the data that correspond to
   identified peptides, and integrate them to intensities.

Most logic stays behind the scenes in Rust. Python is for **driving** the pipeline and
**visualizing** detected peak traces.

## Verdict

**Clearly feasible.** The hardest, most tedious part — reading vendor data files — is
already solved by a mature Rust mass-spec ecosystem, and the core peak-tracing algorithm
is conceptually simple. With the scope below, nothing is a research-grade unknown. This is
a favorable first serious Rust project: small, algorithmic, with the gnarly parts (vendor
I/O, ML) delegated to mature crates and to Python.

Main residual risks (each has a concrete mitigation below):
1. Isotope-distribution **numeric parity** → *Isotope-distribution port plan (port-first)*.
2. Discipline on building a **parity harness** early → *Build the parity harness early — layered*.
3. The **Python binding layer** (PyO3/Arrow/GIL), not the algorithm, is where a first Rust
   project actually stalls → *Python binding layer: the real risk* (core/bindings crate split).

## Scope decisions (agreed)

| Area | Decision |
|---|---|
| **Core XIC quant** | In scope, first. psmtsv + mzML → isotope-envelope peak tracing + integration. |
| **MBR (match-between-runs)** | Later. RT alignment + acceptor-peak search in Rust; **ML/PEP model in Python**. |
| **IsoTracker** | Out of scope (for now). |
| **Bayesian protein quant** | Out of scope (for now). |
| **Thermo `.raw`** | Deferred. mzML-only MVP; add `.raw` later via the `mzdata` crate (reader swap only — IndexingEngine is format-agnostic). |
| **Serialization** | Move off custom NetSerializer blob → **Parquet** (or SQLite only if random-access query pattern emerges). Worth doing in mzLib C# too, independently. |

### Why ML-in-Python for MBR

This sidesteps the weakest part of the Rust ecosystem. MBR splits cleanly along a
Rust/Python boundary:

- **Rust half (heavy I/O + search):** RT alignment between runs (`GetRtCalSpline`),
  donor-file selection, acceptor-peak search (`FindAllAcceptorPeaks`) — all reuse the
  IndexingEngine + XIC machinery already built for core quant. Produces a **feature table**
  (one row per candidate transferred peak).
- **Python half (the model):** the PEP step is currently a FastTree gradient-boosted model
  in `Microsoft.ML`. In Python that's a few lines of `lightgbm`/`xgboost`/`scikit-learn`
  on the feature table, returning scores Rust folds back into FDR.

MBR runs **once at the end** on an already-built table, so the Rust→Python→Rust hand-off
is negligible. You get a mature, debuggable, retrainable model and avoid Rust GBM immaturity.

## Architecture

```
              ┌────────────────────────── Python ──────────────────────────┐
              │  driving + visualization + MBR model (lightgbm/sklearn)     │
              └───────────────▲───────────────────────────▲─────────────────┘
                              │ PyO3 / maturin             │ Arrow / NumPy
              ┌───────────────┴───────────────────────────┴─────────────────┐
              │                          Rust                                │
              │  psmtsv parse → isotope dist → binned IndexingEngine →       │
              │  XIC tracing/integration (GetIsotopicEnvelopes + CutPeak) →  │
              │  FDR → Parquet output                                        │
              │  (later) MBR: RT align + acceptor search → feature table     │
              └───────────────▲──────────────────────────────────────────────┘
                              │ mzdata crate
                       .mzML  │  (.raw later)
```

### Python binding surfaces (the public API)

1. **`quant(results_path, raw_paths, params) -> peaks`**
   Returns either an in-memory peak table or a **path to a Parquet/SQLite store**.
   For full runs with many files, return a **Parquet path** (cheap to hand back,
   lazy-loadable in polars/pandas). Keep an option to materialize small result sets in memory.

2. **IndexingEngine calls — the visualization workhorse.**
   Expose: build/load an index from a data file, and query it. The key methods:
   - `GetIndexedPeak(mz, scanIndex, tolerance)` (point query)
   - an **XIC extractor**: "give me (RT, intensity) across scans for this m/z + charge +
     tolerance" → return as **two `f64` NumPy arrays** (via the `numpy` PyO3 crate), the
     natural shape for matplotlib/plotly.
   A detected peak isn't just an intensity number — it's the underlying XIC array plus the
   integration bounds (`CutPeak` gives apex + boundaries), which you'll want to overlay when plotting.

   **Note (verified against the code):** the C# `IndexingEngine<T>` *already* exposes XIC
   extraction — `GetXic`, `GetXicByScanIndex`, and `GetAllXics` — alongside `GetIndexedPeak`
   and `BinarySearchForIndexedPeak`. So this surface is not a new design; it is a faithful
   **port** of existing methods. Add `GetXic`/`GetAllXics` to the port checklist and to the
   parity harness (see below).

3. **MBR** — its own call (or a flag on `quant`), with the Python-side model hook described above.

### Data interchange

Use **Arrow as the interchange spine**: Rust `arrow`/`parquet` ↔ Python `pyarrow`/`polars`
is zero-copy and avoids hand-rolling PyO3 conversions per struct.
- Big result tables → **Parquet** on disk.
- Small interactive things (one peptide's XIC for a plot) → **NumPy** arrays.

## Rust ecosystem leverage (don't start from zero)

- **`mzdata`** — mature reader for mzML, mzMLb, MGF, and **Thermo `.raw`** (bridges the
  Thermo .NET RawFileReader / ThermoRawFileParser). Eliminates the ~4,800-line mzML port
  *and* gives the eventual `.raw` path. Biggest de-risking factor.
- **`sage`** (Michael Lazear's proteomics search engine) — production Rust codebase that
  already does XIC-based label-free quant with the **same** binned-index + isotope-envelope
  approach. Strong reference implementation.
- **`rustyms`** — peptide/proteoform chemistry: chemical formulas from sequence+mods and
  isotopic distribution calculation. Can replace much of the `Chemistry`/`Proteomics` port.
- **`timsrust`** — native Bruker timsTOF `.d` reader, if ever needed.
- **PyO3 + maturin** — standard, mature Python-bindings path. Low risk.

## What the Rust port actually owns (with this scope)

**Port faithfully from mzLib:**
- The binned m/z `IndexingEngine<T>` — `mzLib/MassSpectrometry/PeakIndexing/IndexingEngine.cs`
  (~340 lines). Jagged array bucketed by `m/z × BinsPerDalton`, binary search over scan
  index within each bin. Trivial to port.
- `FlashLfqEngine.GetIsotopicEnvelopes` + `CheckIsotopicEnvelopeCorrelation` + `CutPeak`
  — the XIC tracing / envelope-correlation / peak-integration core. **The algorithmic heart.**
- `.psmtsv` parsing — `Readers/InternalResults/IndividualResultRecords/PsmFromTsv.cs`
  (~144 lines) + `MzLibExtensions.MakeIdentifications` (~66 lines). Trivial with `csv` + `serde`.
- Charge-state aggregation, FDR.

**Borrow from crates (and validate numeric agreement vs mzLib):**
- mzML reading → `mzdata`.
- Isotope distributions + peptide→formula → `rustyms`.
  (mzLib reference: `Chemistry.IsotopicDistribution.GetDistribution`,
  `Proteomics.AminoAcidPolymer.Peptide.GetChemicalFormula`, `PeriodicTable`, `ChemicalFormula`.)

**Defer:**
- MBR Rust-side feature generation until core quant has parity.
- The MBR model lives in Python from the start.

## Key mzLib reference points (for the port)

| Concept | C# location |
|---|---|
| Engine entry / orchestration | `FlashLFQ/FlashLfqEngine.cs` (`Run` at lines 187–343; ~2228 lines total) |
| Theoretical isotope dists | `FlashLfqEngine.CalculateTheoreticalIsotopeDistributions` (~110 lines) |
| XIC tracing core | `FlashLfqEngine.GetIsotopicEnvelopes`, `CheckIsotopicEnvelopeCorrelation`, `CutPeak` |
| Binned index + query | `MassSpectrometry/PeakIndexing/IndexingEngine.cs` (`GetIndexedPeak`, `IndexPeaks`) |
| FlashLFQ index wrapper | `FlashLFQ/PeakIndexingEngine/PeakIndexingEngine.cs` (note: `Serialize/DeserializeIndex` use NetSerializer → replace with Parquet) |
| Identification model | `FlashLFQ/Identification.cs` (BaseSequence, ModifiedSequence, MonoisotopicMass, charge, RT, file, PeakfindingMass) |
| psmtsv → Identification | `FlashLFQ/ResultsReading/MzLibExtensions.cs` (`MakeIdentifications`, `MakeSpectraFileDict`) |
| psmtsv record / reader | `Readers/InternalResults/IndividualResultRecords/PsmFromTsv.cs`, `.../ResultFiles/PsmFromTsvFile.cs` |
| Data file dispatch | `Readers/MsDataFileReader.cs` (`GetDataFile`) → mzML reader in `Readers/MzML/*` |
| MBR (later) | `FlashLfqEngine.GetRtCalSpline`, `FindPeptideDonorFiles`, `FindAllAcceptorPeaks`, `FindIndividualAcceptorPeak`; PEP in `FlashLFQ/MBR/PEP/PepAnalysisEngine.cs` (Microsoft.ML FastTree) |

### External C# dependencies (context for what NOT to port)
- FlashLFQ.csproj: CsvHelper, MathNet.Numerics, **Microsoft.ML + Microsoft.ML.FastTree**
  (MBR/PEP only → Python), NetSerializer (→ Parquet), SharpLearning.Optimization.
- Readers.csproj: CsvHelper, OpenMcdf, System.Data.SQLite, ZstdSharp, **Thermo vendor DLLs**
  (closed-source .NET — the reason `.raw` is deferred to `mzdata`), Bruker native DLLs.

## Serialization change (Parquet)

Two separate things currently use NetSerializer:
- **The peak index** (`PeakIndexingEngine.SerializeIndex`, the `.ind` blob) — columnar
  (m/z, intensity, scan, RT) compresses extremely well → **Parquet** (Rust `arrow` + `parquet`).
- **The results** (peaks per peptide, intensities per file) — natural tabular Parquet output,
  directly readable from Python pandas/polars.

Use **SQLite only** if a *random-access query* pattern emerges ("all peaks for peptide X")
rather than bulk load. FlashLFQ's "load all → process → write all" pattern favors Parquet.
This change is worth making in mzLib's C# too, independent of the Rust rewrite.

## Python binding layer: the real risk, and how to contain it

The core algorithm (`GetIsotopicEnvelopes` → `CutPeak`) is a few hundred lines of arithmetic
and binary search — the part Rust makes *pleasant*. The friction is concentrated at the
Python boundary, and it lands precisely in the two areas a Rust newcomer hasn't internalized:
ownership/lifetimes and the build/ABI story. Specific sharp edges:

- **GIL + ownership interaction.** Modern PyO3 (0.21+) uses `Bound<'py, T>` smart pointers
  carrying a GIL lifetime. The heavy quant wants to run on rayon threads with the GIL
  *released* (`Python::allow_threads`), but nothing inside that closure may touch a Python
  object — the borrow checker enforces this. Get the structure wrong and you either fight a
  wall of lifetime errors or silently serialize everything behind the GIL and lose all
  parallelism.
- **"Zero-copy Arrow" is doing some lifting.** The arrow-rs ↔ pyarrow handoff over the C Data
  Interface is genuinely zero-copy, but (a) you still pay a real transform to *build* Arrow
  columnar arrays from native `Vec<Peak>` structs — only the cross-language handoff is free,
  and (b) arrow-rs and pyarrow version on independent cadences, and **ABI/version skew is a
  classic footgun** that surfaces as opaque runtime errors. Pin a known-compatible pair and
  document it.
- **NumPy copy-vs-move is a silent trap.** `into_pyarray` moves a `Vec` (no copy);
  `to_pyarray` copies. For the interactive XIC surface, you want the move path — and the two
  are one character apart.
- **Errors and panics.** Map `Result` → Python exceptions via `From<MyError> for PyErr`
  (boilerplate). The real trap: a panic unwinding across the FFI boundary is UB unless
  caught, so anything reachable from Python must return `Err`, never `unwrap()` on bad input.
- **Debugging stops at the boundary.** Stack traces don't cross the FFI line cleanly; a
  misuse segfault surfaces as an opaque Python crash. Expect to lean on logging.

**Mitigation — split into two crates.** A pure-Rust **core** crate that knows nothing about
Python (normal types, normal `Result`, all unit tests and the parity harness run here with
zero Python in the loop), and a thin **bindings** crate (PyO3) whose only job is type
translation and error mapping: parse args → call core → convert `Vec<Peak>` to Arrow/NumPy →
map errors. All the hard stuff above is then confined to one small, rarely-changing module,
while the algorithm you care about stays testable in plain Rust. This is the standard
core/bindings layout and it directly defuses the "first serious Rust project" risk: the
learning budget goes to the algorithm, and the FFI weirdness is quarantined.

**Thermo `.raw` has a packaging tail (Phase 2 concern).** `mzdata`'s Thermo path bridges to
.NET (ThermoRawFileParser / RawFileReader), so that route is **not pure Rust** — it needs a
.NET runtime present on the user's machine. Irrelevant for the mzML MVP (pure Rust, clean
static wheel), but it undercuts the "single self-contained wheel" distribution story when
`.raw` lands. Flag it now so it isn't a surprise in Phase 2.

## Isotope-distribution port plan (port-first; `rustyms` deferred)

**Decision: port mzLib's own isotope logic first, even though `rustyms` may replace it
later.** The reason is parity *attribution*. Starting with `rustyms` means a first-pass
disagreement has three simultaneous suspects — your algorithm, your reader, and `rustyms`'s
different abundance table. Porting mzLib's logic collapses that to "did I transcribe it
faithfully," which is mechanically checkable. `rustyms` then becomes a *later* optimization
validated against an already-passing port, not a source of first-pass mystery.

Port the dependency stack bottom-up; each layer is parity-checked before the next is built:

1. **Periodic-table isotope data — the parity anchor, and it's data, not code.**
   `Chemistry/PeriodicTable.cs` loads each element's isotopes as
   `(isotopeNumber, atomicMass, relativeAbundance)` via `AddIsotope`. **Do not re-source
   these from NIST independently** — sub-ppm abundance differences propagate into every
   theoretical envelope and silently break parity. Dump mzLib's loaded table once and embed
   *that* verbatim in Rust. Reproduce `ValidateAbundances` (per-element abundances sum to 1
   within epsilon) as a Rust unit test so a transcription typo fails loudly.
2. **`ChemicalFormula`** (`Chemistry/ChemicalFormula.cs`) — element→count map plus
   `MonoisotopicMass`/`AverageMass`. Parity: exact-integer counts, float on derived masses.
3. **`Peptide.GetChemicalFormula`** — residue→formula table + modification formulas. Same
   discipline as (1): port the table mzLib uses, don't rebuild it.
4. **`IsotopicDistribution.GetDistribution(formula, 0.125, 1e-8)`**
   (`Chemistry/IsotopicDistribution.cs`, ~435 lines) — the Kubinyi-style fine-grained
   polynomial algorithm (`MergeFinePolynomial`, `MultiplyFinePolynomial`,
   `MultiplyFineFinalPolynomial`, `MultipleFinePolynomialRecursiveHelper`, `FactorLn`,
   `CalculateFineGrain`, plus the `Polynomial` struct and `Composition` class). Port
   literally, preserving the three constants (`defaultFineResolution = 0.125`,
   `defaultMinProbability`, `defaultMolecularWeightResolution`) **and the merge/summation
   order** — the polynomial merge is float-addition-order-sensitive, and reordering perturbs
   the low-probability tail.
5. **`CalculateTheoreticalIsotopeDistributions`** (`FlashLfqEngine.cs` 350–460) — averagine
   fill-in (the five `average{C,H,O,N,S}` constants; the `massDiff > 20` threshold; the
   `Math.Round` element counts), normalize-to-most-abundant, the truncation rule
   (`count < NumIsotopesRequired || abundance > 0.1`), then
   `PeakfindingMass = MonoisotopicMass + mostAbundantIsotopeShift`. These thresholds are
   exact behavioral forks — a wrong comparison operator changes which isotopes are retained
   and shifts every downstream peak search.

Port order **is** the test order: layer 4 isn't written until layer 1's table diffs clean.

## Phased plan

- **Phase 0 — spike (days):** Rust crate + maturin/PyO3 skeleton. Read an `.mzML` via
  `mzdata`. **Exercise the whole binding stack, not just an int return** — the toolchain
  surprises live in the data-interchange layer, and scan-count-as-`int` hides them. The
  spike should: (a) return a real `f64` **NumPy** array via `into_pyarray` (the move/no-copy
  path), and (b) round-trip one **Arrow** table Rust→pyarrow over the C Data Interface. If
  both work and the `maturin develop` wheel imports cleanly, the toolchain is validated. See
  *Python binding layer* below for why this is the right thing to de-risk first.
- **Phase 1 — MVP (core ask):** psmtsv parse → isotope dists (`rustyms` or port) → binned
  index (port) → XIC tracing/integration (port `GetIsotopicEnvelopes`/`CutPeak`) → Parquet
  output. **mzML only.** Build the golden-file parity harness here. Expose the three Python
  surfaces (Quant, IndexingEngine/XIC, stub MBR).
- **Phase 2:** Thermo `.raw` via `mzdata`; charge aggregation, FDR hardening.
- **Phase 3:** MBR — Rust RT alignment + acceptor search → feature table → Python ML model → FDR.

## Build the parity harness early — layered, not end-to-end

The thing that consumes time isn't writing Rust — it's verifying numeric parity. A single
end-to-end peptide-by-peptide diff is **not enough**: by the time a final intensity is wrong,
the cause could be in any of six upstream stages, and float-reordering noise is entangled
with real bugs. Build a **layered harness** instead: each stage emits a canonical artifact,
the Rust side emits the parallel artifact, and one diff tool compares stage-by-stage with
per-stage tolerances. The **first** stage that diverges localizes the bug — you never debug
stage 6 while stage 1's table is silently off. (C# emits "golden" artifacts via small dump
hooks in mzLib or a throwaway test harness.)

| # | Stage | Golden artifact (keyed) | Tolerance |
|---|---|---|---|
| L0 | **Periodic-table data** | element → list of `(isotopeNumber, atomicMass, relativeAbundance)` | **exact** (transcribed data; any diff = typo) |
| L1 | **Peptide → ChemicalFormula** | modified sequence → element counts + monoisotopic mass | counts exact; mass rel-1e-9 |
| L2 | **Raw isotope distribution** | formula → `GetDistribution` `(masses[], intensities[])` | rel-1e-6 on intensities; abs on masses |
| L3 | **FlashLFQ theoretical envelope** | modified sequence → final `(massShift, normAbundance)[]` **and** `PeakfindingMass` | massShift abs-1e-6; abundance rel-1e-6; **same array length** (catches truncation-rule bugs) |
| L4 | **Reader peak lists** | per MS1 scan → `(m/z, intensity)[]`, mzLib reader vs `mzdata` | **report-only first** — see note |
| L5 | **Per-(peptide,charge,file) integrated envelope** | `GetIsotopicEnvelopes` → integrated intensity + pearson corr + apex scan | rel-1e-6 |
| L6 | **End-to-end peptide intensity** | peptide × file → intensity | rel-1e-6 |

Make this work:

- **L3 and L5 are the highest-value checkpoints.** L3 isolates the entire isotope port from
  everything downstream — if L0–L3 pass, the chemistry is correct and any later divergence is
  in indexing/tracing, not isotopes. L5 isolates the algorithmic heart
  (`GetIsotopicEnvelopes`/`CutPeak`) from the reader.
- **L4 is a calibration layer, not a pass/fail gate — at first.** `mzdata` and mzLib's mzML
  reader may centroid/filter differently, so expect *some* legitimate per-scan peak
  differences. Run L4 in report-only mode and characterize the magnitude: peaks matching to
  within reader noise → proceed; whole peaks missing or systematically shifted → a real
  upstream problem that would otherwise masquerade as an algorithm bug at L5/L6.
- **Tolerance taxonomy.** Integer-ish quantities (isotope counts, scan indices, charge,
  retained-isotope array length) compare **exact** — drift there is always a bug. Floats
  compare **relative**: `|a−b| / max(|a|,|b|) < 1e-6`, because FlashLFQ quantifies in parallel
  (`AddPeakToConcurrentDict`, concurrent dicts) and a rayon port changes summation order, so
  intensities differ in the last bits. Have the diff tool print the **worst offender per
  stage** (max relative error + which key), not just pass/fail.

**Golden corpus — reuse FlashLFQ's existing test data** (`mzLib/Test/FlashLFQ/TestData/`):
- **Core quant (Phase 1):** `sliced-mzml.mzML` + `AllPSMs.psmtsv` (and the synthetic
  `sliced-*` pair the unit tests build inline) — small, exercises multiple charge states,
  averagine fill-in via unknown mods, and the ±1 off-by-one mass-shift logic in
  `GetIsotopicEnvelopes`. `SmallCalibratibleYeast.mzml` is an additional compact mzML.
- **MBR (Phase 3):** `PSMsForMbrTest.psmtsv` with `f1r1_sliced_mbr.raw` / `f1r2_sliced_mbr.raw`
  and the `20100614_Velos1_TaGe_SA_K562_3/4.mzML` pair.
- Build the L0–L3 corpus in Phase 1 **before** the tracing port — chemistry parity is a
  precondition for L5/L6 meaning anything.

These map onto the phasing cleanly: L0–L3 land with the isotope port early in Phase 1; L4–L6
with the tracing port later in Phase 1. **The full harness passing is the Phase-1 exit
criterion.**
