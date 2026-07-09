# Porting IsoDec deconvolution to Rust — feasibility investigation

Investigation (2026-07-09) into porting **IsoDec** into the Rust FlashLFQ port, motivated by the
top-down monoisotope-assignment bottleneck (our averagine-cosine mono placement leaves a large
strict-recall gap; see `analysis/topdown_bench/RESULTS.md`). IsoDec produced one of the two
ground-truth proteoform sets (`IsoDec_AllProteoforms.psmtsv`), so it is a proven top-down deconvoluter.

## What IsoDec is

A **neural-network isotopic charge-state assignment** method from the Marty lab (UniDec), published as
"A Fast Neural Network for Isotopic Charge State Assignment," *JACS* 2025 (Pavek, Bollis, Grimes,
Shortreed, Smith, Marty). It reports more correct feature assignments than existing tools on complex
top-down spectra, with better speed and accuracy — the gain comes from the NN charge assignment, not
just better scoring/filtering. mzLib/MetaMorpheus call it via FFI to a compiled C library.

## Architecture (the important part — it is small)

**1. Phase encoding (deterministic, no ML).** A peak cluster is turned into a fixed
`(maxz=50, phaseres=8)` = **400-element** "charge phase histogram." For each charge hypothesis
`z = 1..50`, each peak's phase `((m/z · z / 1.00335) mod 1)` is phase-locked to the most intense peak
and its intensity accumulated into one of 8 phase bins. Real isotope spacing `1.00335/z` makes the
true charge's row concentrate intensity in a tight phase band; wrong charges smear. Key constants:
neighborhood `lowmz=-1.5, highmz=+5.5`, alignment `window=16`, `minpeaks=3`.

**2. The network — a single-hidden-layer MLP.** For the default `phaseres=8` model
(`Fast8PhaseNeuralNetwork`):
```
input  400  (flattened 50×8)
Linear 400 → 400  + ReLU
Linear 400 → 50            # 50 charge classes (z = 0..49)
argmax → predicted charge
```
≈ **180k parameters**. The `phaseres=4` variant is ~50k. **No convolutions.** Inference is two
matrix–vector products, a ReLU, and an argmax — reimplementable in pure Rust in tens of lines, no ML
runtime needed.

**3. Post-processing (deterministic).** Given the predicted charge, IsoDec matches the isotope peaks,
scores the envelope, and computes monoisotopic / average / most-abundant masses — the `MPStruct`
(matched-peak struct: `z`, m/z, masses, per-isotope m/z + distribution + match indices [64 each],
score). This is ordinary isotope-envelope bookkeeping, portable directly.

## How it ships in UniDec (`unidec/IsoDec/`)

- **Compiled inference lib:** `isodeclib.{dll,so,lib}` — C, loads model weights and runs inference.
- **Trained weights:** `phase_model_{1,4,8}.pth` (PyTorch) **and `.bin`** (the flat format the C lib
  reads). Small files.
- **C API (via `c_interface.py` ctypes):**
  - `IDSettings DefaultSettings()`
  - `int predict_charge(double* mz, float* intensity, int count, char* modelpath)`
  - `int process_spectrum(double* mz, float* intensity, int count, char* modelpath, MPStruct* out, IDSettings, char* inputtype)` → number of matched peaks
  - `encode(double* mz, float* intensity, int count, float* out, IDConfig, IDSettings)`
  - Structs: `MPStruct` (matched peak), `IDSettings` (~20 thresholds/tolerances), `IDConfig`
    (verbose, phaseres, maxz, dims).
- The C source (`isodeclib.c`) is in the repo and is the deterministic reference spec for both encoding
  and inference — ideal to port from.

## Two port paths

### Path A — FFI to the prebuilt `isodeclib` (fastest to a working, exact result)
Write Rust `extern "C"` bindings for `process_spectrum` + `MPStruct`/`IDSettings`/`IDConfig`, link the
shipped `isodeclib.lib`/`.so`, vendor the `.bin` models, call per candidate cluster.
- **Pros:** byte-for-byte the same result mzLib gets; days not weeks; immediately quantifies the recall
  ceiling IsoDec buys us before committing to a native port.
- **Cons:** a prebuilt binary + model-blob dependency (platform-specific), a C toolchain / `build.rs`
  link step, and it breaks the project's "pure-Rust, parity-gated" discipline. Licensing: UniDec is
  open source (check its exact license before vendoring binaries).

### Path B — native Rust reimplementation (matches project discipline)
Port `isodeclib.c` to Rust: (1) phase encoding, (2) the two-layer MLP forward pass, (3) weight loading
from the `.bin` format (or re-export the `.pth` to a simple layout), (4) the `MPStruct` matching/mass
math.
- **Pros:** no binary/model-runtime dependency; native, portable, embeddable in the detector's
  refine/resolve exactly like `classic_deconvolute`; parity-gateable against the C lib's output and
  against `IsoDec_AllProteoforms`.
- **Cons:** more work than Path A, and must reproduce the C encoding/matching edge cases exactly. The
  network itself is trivial; the effort is the encoding + peak-matching fidelity.
- **Effort:** moderate — the network is ~180k params of two `Linear` layers; the real work is faithful
  encoding + matching, both small deterministic C. Weeks, not months.

## Recommendation

1. **Prototype Path A first** to measure the recall/mono-accuracy IsoDec delivers on our Jurkat +
   golden benchmarks vs the current averagine path. If it does not clearly beat our (improving)
   mono-assignment, stop — no port needed.
2. **If it wins, commit to Path B** as the shippable form, using Path A (and the C source) as the
   parity oracle, following the same golden-diff discipline as the ClassicDeconvolution port.
3. Wire IsoDec as an **alternative refine-stage deconvoluter** (a sibling of `classic_deconvolute`),
   selectable per run — its neural charge assignment plus careful isotope matching directly targets the
   monoisotope off-by-one that dominates our top-down strict-recall gap.

## Relevance to the current bottleneck
Our diagnostic (RESULTS.md) shows the detector *finds* ~99% of proteoforms but strict recall is gated
by monoisotope off-by-one on large envelopes. IsoDec's learned charge assignment + matching is exactly
the component that would replace our fragile averagine-cosine mono placement. This is the strongest
single lever identified, and the port is feasible (small MLP + deterministic C reference).

Sources: UniDec repo (github.com/michaelmarty/UniDec, `unidec/IsoDec/`); "A Fast Neural Network for
Isotopic Charge State Assignment," JACS 2025 (doi:10.1021/jacs.5c03162).
