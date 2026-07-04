# Simulated IMSP API And Metadata Design

## Goal

Expose simulation as a reusable backend service instead of a hard-coded test workflow.

The first production-oriented workflow should:

1. Accept an experimental mzML file path.
2. Accept a MetaMorpheus `.psmtsv` results file path.
3. Load top-down proteoform candidates from the results file.
4. Fit/simulate those candidates against the real MS1 data.
5. Write one combined simulated `.imsp` file.
6. Write a sidecar metadata file that links simulated proteoforms to approximate signal regions in the IMSP.
7. Return output paths and summary information to the TauriCS backend/frontend.

## Non-Goals For V1

- Do not write one `.imsp` file per proteoform.
- Do not encode proteoform provenance directly into the IMSP binary format yet.
- Do not support Thermo RAW as the first API target.
- Do not require exact peak-to-proteoform attribution in the first metadata schema.
- Do not move the current simulator into `MassSpectrometry`; simulation should remain in `TopDownSimulator` or a dedicated application/service layer.

## Why Not One IMSP Per Proteoform?

One file per proteoform sounds simple, but it causes practical problems:

- Many duplicate scan tables and m/z/bin directories.
- Thousands of small files for large searches.
- More frontend file coordination.
- Harder overlay rendering.
- No natural way to represent coeluting or overlapping proteoforms.
- More expensive packaging and cleanup.

The better default is one signal file per simulation run plus sidecar metadata.

```text
sample.simulated.imsp
sample.simulated.imsp.metadata.json
```

The `.imsp` file remains a compact render target. The `.json` file carries semantic provenance and UI-friendly search/filter data.

## Current Relevant Code

Existing pieces already support most of the workflow:

- `Readers.PsmFromTsvFile`
  - Location: `mzLib/Readers/InternalResults/ResultFiles/PsmFromTsvFile.cs`
  - Loads MetaMorpheus `.psmtsv` results.

- `Readers.PsmFromTsv`
  - Location: `mzLib/Readers/InternalResults/IndividualResultRecords/PsmFromTsv.cs`
  - Provides fields such as accession, full sequence, precursor charge, monoisotopic mass, retention time, score, q-value, and scan numbers.

- `TopDownSimulator.Extraction.MmResultLoader`
  - Location: `mzLib/TopDownSimulator/Extraction/MmResultLoader.cs`
  - Already converts `PsmFromTsv` records into `MmResultRecord`.

- `TopDownSimulator.Extraction.GroundTruthExtractor`
  - Fits observed real MS1 signal around candidate proteoforms.

- `TopDownSimulator.Fitting.ParameterFitter`
  - Builds fitted `ProteoformModel` objects.

- `TopDownSimulator.Simulation.Simulator`
  - Produces simulated scans from fitted models.

- `MassSpectrometry.ImspExportService`
  - Location: `mzLib/MassSpectrometry/PeakIndexing/ImspExportService.cs`
  - Writes IMSP from `MsDataScan` collections.

## Service Boundary

Use the same two-layer pattern as the mzML IMSP export service:

```text
Readers
  Loads mzML and MetaMorpheus .psmtsv files.

TopDownSimulator
  Owns fitting, simulation, and simulation metadata generation.

MassSpectrometry
  Owns generic IMSP writing from scans.

TauriCS backend
  Owns app endpoints, cancellation, progress reporting, and output path policy.
```

Avoid putting path-based `.psmtsv` + simulation APIs in `MassSpectrometry`. It would pull in `Readers` and simulator dependencies that do not belong in the lightweight IMSP writer layer.

## Proposed Public API

Add a service in `TopDownSimulator`, for example:

```csharp
namespace TopDownSimulator.Export;

public interface ISimulatedImspExportService
{
    SimulatedImspExportResult ExportFromMetaMorpheusPsmTsv(
        string mzmlPath,
        string psmTsvPath,
        string outputDirectory,
        SimulatedImspExportOptions? options = null);
}
```

Result DTO:

```csharp
public sealed record SimulatedImspExportResult(
    string ImspPath,
    string MetadataPath,
    string SummaryPath,
    int InputRecordCount,
    int CandidateRecordCount,
    int FittedProteoformCount,
    int SimulatedScanCount,
    int SimulatedPeakCount,
    IReadOnlyList<SimulatedImspExportWarning> Warnings);
```

Options DTO:

```csharp
public sealed record SimulatedImspExportOptions
{
    public double QValueThreshold { get; init; } = 0.01;
    public bool DeduplicateProteoforms { get; init; } = true;
    public bool CentroidOutput { get; init; } = true;
    public double IntensityThreshold { get; init; } = 10000;
    public int BinsPerDalton { get; init; } = 100;
    public double RtHalfWidthMinutes { get; init; } = 0.25;
    public double PpmTolerance { get; init; } = 20.0;
    public double MzWindowHalfWidth { get; init; } = 0.05;
    public int? MaxProteoforms { get; init; } = null;
    public bool UseGlobalAbundanceRefit { get; init; } = false;
}
```

V1 should keep global abundance refit disabled by default unless performance is well bounded.

## Backend Endpoint Shape

For TauriCS, expose a thin command/API wrapper:

```csharp
public SimulatedImspExportResult SimulateMetaMorpheusToImsp(
    string mzmlPath,
    string psmTsvPath,
    string outputDirectory,
    SimulatedImspExportOptions? options = null)
```

The frontend should receive paths, not large binary payloads:

```json
{
  "imspPath": "D:/run/sample.simulated.imsp",
  "metadataPath": "D:/run/sample.simulated.imsp.metadata.json",
  "summaryPath": "D:/run/sample.simulated.summary.json",
  "inputRecordCount": 12000,
  "candidateRecordCount": 943,
  "fittedProteoformCount": 821,
  "simulatedScanCount": 3500,
  "simulatedPeakCount": 321000,
  "warnings": []
}
```

## Metadata V1: Region-Level Provenance

The current simulator rasterizes all proteoforms into one total intensity grid:

```csharp
RasterizedScanGrid(double[] ScanTimes, double[] MzGrid, double[,] Intensities)
```

That structure does not retain exact per-proteoform contributions. V1 metadata should therefore map each proteoform to an approximate expected region:

- retention time span
- scan index span
- m/z span
- bin span
- charge range
- model parameters
- source `.psmtsv` record references

This supports frontend features like:

- search proteoform by accession/sequence
- highlight approximate region in the IMSP heatmap
- list simulated candidates visible in the current viewport
- show model/fitting summaries

It does not yet answer: "the user clicked this exact peak; which proteoforms contributed to it?"

### Metadata File Example

```json
{
  "format": "mzlib-simulated-imsp-metadata",
  "version": 1,
  "createdUtc": "2026-04-16T00:00:00Z",
  "imspFile": "sample.simulated.imsp",
  "sourceMzmlFile": "sample.mzML",
  "sourceResultsFile": "sample.psmtsv",
  "options": {
    "qValueThreshold": 0.01,
    "deduplicateProteoforms": true,
    "centroidOutput": true,
    "intensityThreshold": 10000,
    "binsPerDalton": 100,
    "rtHalfWidthMinutes": 0.25,
    "ppmTolerance": 20.0
  },
  "summary": {
    "inputRecordCount": 12000,
    "candidateRecordCount": 943,
    "fittedProteoformCount": 821,
    "simulatedScanCount": 3500,
    "simulatedPeakCount": 321000
  },
  "proteoforms": [
    {
      "id": "P12345|PEPTIDEFORM|12",
      "sourceRecordIds": ["P12345|PEPTIDEFORM|12|scan=12345"],
      "accession": "P12345",
      "fullSequence": "PEPTIDEFORM",
      "precursorCharge": 12,
      "monoisotopicMass": 12345.6789,
      "score": 142.2,
      "qValue": 0.004,
      "model": {
        "abundance": 1234567.0,
        "rtApexMinutes": 32.45,
        "rtStartMinutes": 31.95,
        "rtEndMinutes": 32.95,
        "sigmaMz": 0.012,
        "minCharge": 10,
        "maxCharge": 18
      },
      "region": {
        "scanIndexStart": 503,
        "scanIndexEnd": 542,
        "mzStart": 1015.2,
        "mzEnd": 1238.6,
        "binStart": 101520,
        "binEnd": 123860
      }
    }
  ],
  "warnings": []
}
```

## Metadata V2: Exact Contribution Provenance

If the frontend needs click-to-identify behavior, add exact provenance later.

Possible approaches:

1. Per-proteoform sparse contribution ranges.
   - Store compact scan/bin rectangles or sparse ranges per proteoform.
   - Good for "which proteoforms overlap this viewport?"

2. Per-peak contributor index.
   - For each IMSP peak/bin/scan point, store contributing proteoform IDs and contribution fractions.
   - Best for click-to-identify.
   - Potentially very large.

3. On-demand recomputation.
   - Given a clicked m/z/RT region, recompute expected contributors from fitted models.
   - Smaller files, more CPU during interaction.

V1 should not block V2. Use stable proteoform IDs and keep metadata schema versioned.

## Proteoform Identity

Use a stable string ID derived from the search result identity:

```text
{Accession}|{FullSequence}|z{PrecursorCharge}|mass{MonoisotopicMass rounded}
```

If multiple PSMs map to the same proteoform, deduplicate by default and retain:

- best score
- best q-value
- all source record IDs
- representative retention time, likely from best score or apex fitting

The previous dedup key used:

```text
(Accession, FullSequence, PrecursorCharge)
```

That is still a reasonable V1 default.

## Output Files

For input:

```text
sample.mzML
sample.psmtsv
```

Default outputs:

```text
sample.simulated.imsp
sample.simulated.imsp.metadata.json
sample.simulated.summary.json
```

The summary file can duplicate high-level metadata but remain small enough for quick frontend loading before opening the larger metadata file.

## Error Handling

The service should distinguish:

- invalid input path
- unsupported input extension
- mzML read failure
- psmtsv parse failure
- no candidate records after filtering
- no fitted models
- simulation produced zero peaks above threshold
- output write failure

The API result should include warnings for non-fatal skips:

```csharp
public sealed record SimulatedImspExportWarning(
    string Code,
    string Message,
    string? ProteoformId = null);
```

## Progress Reporting

Simulation can be slow. The eventual app service should support progress:

```csharp
public sealed record SimulatedImspExportProgress(
    string Stage,
    int? Completed,
    int? Total,
    string? Message);
```

Stages:

- `loading-mzml`
- `loading-results`
- `filtering-results`
- `indexing-ms1`
- `fitting-proteoforms`
- `simulating-scans`
- `writing-imsp`
- `writing-metadata`

This can be added after the first synchronous service is stable.

## Recommended V1 Decision

Build V1 around one combined IMSP plus region-level JSON metadata. This is the best balance of usability, file size, and implementation risk.

Exact peak-to-proteoform provenance should be a V2 feature after the frontend has proven it needs click-level attribution.
