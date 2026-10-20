# SimSpec v1 — the simulator input contract

**Status:** Draft · **Schema version:** 1 · **Created:** 2026-10-20

The **simulation spec** is the file the C# `TopDownSimulator` emits and the Rust `sim-core` reads.
It is *self-sufficient*: given the spec (and the companion averagine table, §4), Rust can reproduce
an MS1 run deterministically with no access to mzLib or the original raw file.

It supersedes `.groundtruth.tsv` as the Rust-facing input, because that TSV **omits two things Rust
needs** (§3). The TSV stays as the C#-internal round-trip
(`Simulator.WriteGroundTruth`/`ReadGroundTruth`, `Simulation/Simulator.cs:445/489`).

## 1. Format

A single UTF-8 JSON file, extension `.simspec.json`, written with `System.Text.Json`
(`WriteIndented = true`) to match the project's existing JSON idiom
(`Noise/AcquisitionProfile.cs:227`). Doubles use .NET's shortest round-trippable form, which
`serde_json` parses to the identical `f64`. Angles/masses are in mzLib's native units
(neutral monoisotopic mass in Da, m/z in Th, retention time in **minutes**).

## 2. Schema

```jsonc
{
  "schemaVersion": 1,
  "generator": {                 // provenance; informational
    "tool": "TopDownSimulator",
    "mzLibCommit": "733955c0",   // optional
    "writtenUtc": "2026-10-20T00:00:00Z"
  },
  "sourceRun": {                 // provenance; optional, may be null
    "fileName": "Jurkat_rep2.raw",
    "description": "fit from MetaMorpheus proteoform IDs"
  },

  "instrument": {
    "polarity": "Positive",      // mzLib MassSpectrometry.Polarity
    "analyzer": "Orbitrap",      // mzLib MZAnalyzerType
    "msLevel": 1
  },

  "chargeRange": { "min": 3, "max": 60 },   // run-level; z range envelopes are rendered over

  "peakWidthModel": { "kind": "Orbitrap", "k": 1.2e-7 },
  //  OR  { "kind": "Constant", "sigma": 0.01 }

  "scanTimes": [ 30.0000, 30.0085, 30.0170, ... ],  // MS1 RT axis, minutes, strictly ascending,
                                                    // one entry per MS1 scan (reused from real run)

  "reduction": null,            // or { "relativeIntensityThreshold": 1e-4,
                                //       "absoluteIntensityFloor": 0.0 }
                                // null ⇒ these defaults (C# PrepareMs1 always reduces). With noise
                                // non-null the floor is instead S/N-threshold × local noise level.
  "noise": null,                // or the noise block below; null = clean analytical simulation

  "proteoforms": [
    {
      "identifier": "P12345",
      "monoisotopicMass": 10000.123,
      "abundance": 1.0e9,
      "rt":     { "mu": 30.0, "sigma": 0.20, "tau": 0.10 },
      "charge": { "kind": "Gaussian", "muZ": 12.0, "sigmaZ": 2.0 }
    }
  ],

  "background": []              // optional; rendered but excluded from feature truth
}
```

### 2a. Noise block (present when `noise` is non-null)

Noise is **live-editable in v1** (the workbench exposes these as knobs — see the plan's Interaction
model), so the block carries the full parameter set, not just a toggle. The **data-free golden still
uses `noise: null`** so signal-path parity starts on the pure forward model + centroider. Noise
parity is **statistical**: `sim-core` draws the same distributions from its own random streams
(§5), so a seed reproduces a run within one implementation, not across C# and Rust.

Field names below mirror the C# model properties exactly, so a spec round-trips to the real
`NoiseFloorModel` / `SignalToNoiseDistribution` / `PeakJitterModel`. Defaults shown are the measured
Orbitrap values.

```jsonc
"noise": {
  "seed": 0,                                 // floor + jitter reproducible per implementation. C# uses the
                                             // low 32 bits (int noiseSeed); Rust uses all 64
  "floor": {                                 // Noise/NoiseFloorModel.cs
    "noiseLevelAtReferenceMz": 1120.0,       //   amplitude at ReferenceMz (650 Th); per-run knob
    "densityScale": 1.0,                     //   DensityScale (0 = no added floor; ~16900 peaks/scan at 1.0)
    "levelDispersion": 0.25,                 //   lognormal per-peak level scatter (JurkatLevelDispersion)
    "signalToNoise": {                       //   Noise/SignalToNoiseDistribution.cs — S/N = threshold + Exp(meanExcess)
      "threshold": 1.94,                     //     OrbitrapCentroids default
      "meanExcess": 1.165
    }
  },
  "jitter": {                                // Noise/PeakJitterModel.cs; enabled:false ⇒ no jitter applied
                                             // (C# currently disagrees — see the plan's open risk 4)
    "enabled": true,
    "mzPpmAtUnitSignalToNoise": 9.06,        //   per-peak m/z σ = sqrt(a^2/(S/N) + b^2) ppm
    "mzPpmFloor": 0.81,
    "commonModeMzPpm": 1.887,                //   per-scan shared m/z shift σ (calibration drift)
    "logIntensitySigmaAtUnitSignalToNoise": 1.026,
    "logIntensitySigmaFloor": 0.165,         //   irreducible ~17% scan-to-scan scatter
    "commonModeLogIntensitySigma": 0.0,      //   per-scan shared intensity factor σ (unmeasured; default 0)
    "dropPeaksFallingBelowFloor": true
  }
}
```

Implemented in Phase 0.1 (`SimNoise`/`SimNoiseFloor`/`SimSignalToNoise`/`SimNoiseJitter` in
`mzLib/TopDownSimulator/Simulation/SimulationSpec.cs`), with round-trip tests. Not v1: the
`scanConditioned` per-scan floor and its per-bin density/level tables
(`Noise/ScanNoiseConditions.cs`) — a later addition.

## 3. Field → C# source mapping, and the two gaps

| Spec field | C# origin | Note |
|---|---|---|
| `proteoforms[].monoisotopicMass` | `ProteoformModel.MonoisotopicMass` (`Model/ProteoformModel.cs`) | neutral, Da |
| `proteoforms[].abundance` | `ProteoformModel.Abundance` | `A_p` |
| `proteoforms[].rt.{mu,sigma,tau}` | `EmgProfile` (`Model/EmgProfile.cs`) | unit-area EMG |
| `proteoforms[].charge.{muZ,sigmaZ}` | `GaussianChargeDistribution` (`Model/ChargeStateDistribution.cs`) | only Gaussian round-trips today |
| `proteoforms[].identifier` | `ProteoformModel.Identifier` | may be empty |
| `chargeRange.{min,max}` | `minCharge`/`maxCharge` args to `Simulate*` | run-level (the TSV repeats these per row) |
| `peakWidthModel` | `IPeakWidthModel.ToString()` → `Constant(sigma=R)` / `Orbitrap(k=R)` (`Model/IPeakWidthModel.cs:40,82`) | parsed structurally, not from the ToString string |
| `scanTimes` | **`scanTimes` arg to `WriteMzml`** (`Simulation/Simulator.cs:178`) | **GAP 1 — not in `.groundtruth.tsv`.** Reused from the real run. |
| `reduction` | `ScanReductionOptions` (`Simulation/SimulatedScanReducer.cs:12`) | **GAP 2a** — defaults `relative=1e-4` |
| `noise` | `NoiseFloorModel`/`PeakJitterModel` + `noiseSeed` | **GAP 2b** — not in the TSV |

**Gap summary:** `.groundtruth.tsv` has the proteoform vector + charge range + width law, but **not**
the scan-time axis and **not** the noise/reduction config. SimSpec v1 adds exactly those so the file
is self-sufficient.

## 4. The averagine table (companion file — part of the contract)

The isotope-envelope kernel (`Model/IsotopeEnvelopeKernel.cs`) shapes each envelope from the mzLib
`Averagine` static tables (namespace `MassSpectrometry`, file
`mzLib/MassSpectrometry/Deconvolution/AverageResidue/Averagine.cs`): `AllMasses` / `AllIntensities`
(both most-intense-first) / `MostIntenseMasses` / `DiffToMonoisotopic`. For Rust envelopes to match
C# ones, Rust must use the **identical** table. Phase 0 exports it once as `averagine.v1.json`
alongside the golden:

```jsonc
{ "allMasses": [[...],...], "allIntensities": [[...],...],
  "mostIntenseMasses": [...], "diffToMonoisotopic": [...] }
```

Rust embeds this verbatim. The index-selection step Rust must also replicate is
`new Averagine().GetMostIntenseMassIndex(mass)` — it maps a target neutral mass to a table row;
`mostIntenseMasses` is exported so Rust can reproduce that lookup (and its tie-breaking) exactly.
The kernel then anchors the row on the caller's monoisotopic mass by shifting every isotopologue by
`mass − (AllMasses[idx][0] − DiffToMonoisotopic[idx])`, normalizes intensities to sum 1, and
re-sorts ascending by mass.

**Gotcha carried over from the C# side:** averagine's 1e-8 absolute cutoff drops the monoisotopic
peak above ~31 kDa, so `CentroidMzs(z)[0]` is *not* the monoisotopic m/z — compute monoisotopic m/z
from mass directly (via `ToMz`: `mass/z + protonMass`), as the C# kernel does. The golden's
large (>31 kDa) proteoform exists to exercise exactly this path.

## 5. The data-free golden (Phase-1 parity target)

A tiny, self-contained fixture that needs **no** raw data, so the Rust port can be verified offline:

- **Input:** `golden-v1.simspec.json` — 2 proteoforms (one small, ~5 kDa; one large, >31 kDa to
  exercise the monoisotopic-drop gotcha), a short `scanTimes` axis (~20 scans), `Orbitrap` width,
  `noise: null`, `reduction: null`.
- **Expected output:** `golden-v1.scans.tsv` — the centroided result from
  `Simulator.SimulateCentroided` (pure forward model + centroider), one row per peak:
  `scanIndex, oneBasedScanNumber, retentionTime, mz, intensity`.
- **Companion:** `averagine.v1.json`.

Generator lives in `mzLib/Test/TopDownSimulator/` and writes the fixture under
`mzLib/Test/TopDownSimulator/SimSpecFixtures/`. Following the existing parity convention
(ProteoStar reads mzLib fixtures via `MZLIB_DIR`, with committed goldens copied into the repo), a
copy lives in `crates/sim-core/parity/golden/` (the averagine table in `crates/sim-core/data/`,
embedded at build time).

Parity tolerance (Phase 1): m/z match within 1e-6 Th, intensity within 1e-4 relative. **Achieved:**
bit-exact on this golden (Rust's exact centroider path). Three `Simulator.PrepareMs1` goldens sit
alongside it, same models through reduction and noise ± jitter: `reduced-v1` matches within
tolerance (Rust's faster pipeline prunes far-sub-floor Gaussians); `noise-jitter-v1` and
`noise-nojitter-v1` match statistically (Rust uses its own random streams). An earlier,
deliberately slower port matched all three bit for bit.

## 6. Versioning

`schemaVersion` is an integer; Rust refuses an unknown major. Additive fields (new optional keys)
keep the version; a breaking change (renamed/removed field, changed units) bumps it and ships a new
golden (`golden-v2.*`).
