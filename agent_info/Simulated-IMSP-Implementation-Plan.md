# Simulated IMSP Implementation Plan

This plan turns `AnalysisExample.ExportRep2SliceAndQValueSimulations` into reusable services.

The design target is documented in:

```text
agent_info/Simulated-IMSP-API-and-Metadata-Design.md
```

## Current Problem

The working simulation path lives in a test class and is hard-coded around local Jurkat file paths:

```text
mzLib/Test/FileReadingTests/AnalysisExample.cs
```

That test does several real pieces of production logic:

- load mzML/RAW data
- load MetaMorpheus `.psmtsv` results
- q-value filter
- deduplicate proteoforms
- build an MS1 peak index
- extract ground truth signal
- fit proteoform models
- simulate scans
- optionally centroid simulated scans
- write IMSP

The first refactor should preserve behavior while moving these responsibilities into small services.

## Proposed New Namespaces

```text
TopDownSimulator.Export
TopDownSimulator.Metadata
TopDownSimulator.Results
```

Suggested files:

```text
mzLib/TopDownSimulator/Export/ISimulatedImspExportService.cs
mzLib/TopDownSimulator/Export/SimulatedImspExportService.cs
mzLib/TopDownSimulator/Export/SimulatedImspExportOptions.cs
mzLib/TopDownSimulator/Export/SimulatedImspExportResult.cs
mzLib/TopDownSimulator/Export/SimulatedImspExportWarning.cs

mzLib/TopDownSimulator/Metadata/SimulatedImspMetadata.cs
mzLib/TopDownSimulator/Metadata/SimulatedProteoformMetadata.cs
mzLib/TopDownSimulator/Metadata/SimulatedImspMetadataWriter.cs
mzLib/TopDownSimulator/Metadata/ProteoformRegionEstimator.cs

mzLib/TopDownSimulator/Results/MetaMorpheusProteoformCandidateLoader.cs
mzLib/TopDownSimulator/Results/ProteoformCandidate.cs
mzLib/TopDownSimulator/Results/ProteoformCandidateFilter.cs
```

## Phase 1: Candidate Loading

Goal: turn MetaMorpheus `.psmtsv` into simulation candidates without involving IMSP yet.

Existing code:

```text
TopDownSimulator.Extraction.MmResultLoader
Readers.PsmFromTsvFile
Readers.PsmFromTsv
```

Proposed DTO:

```csharp
public sealed record ProteoformCandidate(
    string Id,
    string FileNameWithoutExtension,
    int PrecursorScanNumber,
    int Ms2ScanNumber,
    int PrecursorCharge,
    double MonoisotopicMass,
    double RetentionTime,
    double Score,
    double? QValue,
    string FullSequence,
    string Accession);
```

Implementation notes:

- Reuse `PsmFromTsvFile`.
- Preserve `MmResultLoader` temporarily or make it delegate to the new candidate loader.
- Keep parser options explicit:

```csharp
new SpectrumMatchParsingParameters
{
    ParseMatchedFragmentIons = false,
}
```

Filtering:

- require positive monoisotopic mass
- require non-negative retention time
- require positive precursor charge
- apply `QValueThreshold` if q-values are present
- support optional file-name filtering if one `.psmtsv` contains multiple files

Tests:

- load a small `.psmtsv` fixture
- verify candidate count
- verify q-value filtering
- verify stable candidate IDs

## Phase 2: Deduplication

Goal: avoid simulating duplicate PSMs for the same proteoform.

Default dedup key:

```text
(Accession, FullSequence, PrecursorCharge)
```

Retention policy:

- keep candidate with best q-value if present
- otherwise keep highest score
- retain all source candidate IDs in metadata

Proposed API:

```csharp
public sealed class ProteoformCandidateFilter
{
    public IReadOnlyList<ProteoformCandidateGroup> FilterAndGroup(
        IReadOnlyList<ProteoformCandidate> candidates,
        SimulatedImspExportOptions options);
}
```

Group DTO:

```csharp
public sealed record ProteoformCandidateGroup(
    string ProteoformId,
    ProteoformCandidate Representative,
    IReadOnlyList<ProteoformCandidate> SourceCandidates);
```

Tests:

- duplicate groups collapse
- highest score wins when q-value missing
- best q-value wins when q-value present
- all source IDs are retained

## Phase 3: Fitting Service

Goal: move `FitProteoforms` out of `AnalysisExample`.

Proposed service:

```csharp
public sealed class ProteoformModelFitService
{
    public ProteoformModelFitResult Fit(
        IReadOnlyList<ProteoformCandidateGroup> candidateGroups,
        GroundTruthExtractor extractor,
        SimulatedImspExportOptions options);
}
```

Result:

```csharp
public sealed record ProteoformModelFitResult(
    IReadOnlyList<FittedProteoform> FittedProteoforms,
    IReadOnlyList<ProteoformGroundTruth> GroundTruths,
    int MinCharge,
    int MaxCharge,
    double SigmaMz,
    IReadOnlyList<SimulatedImspExportWarning> Warnings);
```

Implementation notes:

- This is mostly current `AnalysisExample.FitProteoforms`.
- Keep invalid abundance handling.
- Keep median sigma fallback.
- Keep global abundance refit optional.

Tests:

- small synthetic ground truth fit
- invalid candidates become warnings, not hard crashes
- output min/max charge and sigma are stable

## Phase 4: Simulation Export Service

Goal: orchestrate full mzML + `.psmtsv` to IMSP + metadata.

Service:

```csharp
public sealed class SimulatedImspExportService : ISimulatedImspExportService
{
    public SimulatedImspExportResult ExportFromMetaMorpheusPsmTsv(
        string mzmlPath,
        string psmTsvPath,
        string outputDirectory,
        SimulatedImspExportOptions? options = null);
}
```

Steps:

1. Validate paths.
2. Load mzML with `MsDataFileReader.GetDataFile(mzmlPath)`.
3. Get MS1 scans with `GetMS1Scans()`.
4. Load candidates from `.psmtsv`.
5. Filter/deduplicate candidates.
6. Build peak index with `PeakIndexingEngine.InitializeIndexingEngine`.
7. Build `GroundTruthExtractor`.
8. Fit models.
9. Simulate scans with `Simulator.Simulate`.
10. Centroid simulated scans if requested.
11. Write combined IMSP with `ImspExportService`.
12. Build metadata.
13. Write metadata JSON.
14. Return result.

Dependencies:

- This service can live in `TopDownSimulator` because that project already references `Readers` and `MassSpectrometry`.
- Do not put this service in `Readers`; simulation logic should not be owned by file readers.
- Do not put this service in `MassSpectrometry`; it would make the lightweight IMSP writer depend on search result parsing and simulation.

Tests:

- use synthetic or tiny mzML/test records if available
- otherwise test orchestration pieces independently first
- add one integration test behind an explicit category if it needs larger files

## Phase 5: Metadata V1

Goal: write a compact region-level sidecar file.

Proposed metadata records:

```csharp
public sealed record SimulatedImspMetadata(
    string Format,
    int Version,
    DateTime CreatedUtc,
    string ImspFile,
    string SourceMzmlFile,
    string SourceResultsFile,
    SimulatedImspExportOptions Options,
    SimulatedImspSummary Summary,
    IReadOnlyList<SimulatedProteoformMetadata> Proteoforms,
    IReadOnlyList<SimulatedImspExportWarning> Warnings);
```

Proteoform metadata:

```csharp
public sealed record SimulatedProteoformMetadata(
    string Id,
    IReadOnlyList<string> SourceRecordIds,
    string Accession,
    string FullSequence,
    int PrecursorCharge,
    double MonoisotopicMass,
    double Score,
    double? QValue,
    SimulatedProteoformModelMetadata Model,
    SimulatedProteoformRegion Region);
```

Region estimation:

- Use fitted RT profile to estimate `rtStart`, `rtApex`, `rtEnd`.
- Convert RT boundaries to scan index boundaries using simulated scan times.
- Estimate m/z range from monoisotopic mass, min/max charge, isotope envelope padding, and sigma m/z.
- Convert m/z range to bins using `binsPerDalton`.

V1 does not need to know exact generated peaks.

Tests:

- metadata serializes and deserializes
- schema version is present
- regions contain valid ordered boundaries
- source candidate IDs are preserved

## Phase 6: TauriCS Backend Wrapper

Add an app-level wrapper in the TauriCS project, not mzLib:

```csharp
public sealed class SimulationBackendService
{
    private readonly ISimulatedImspExportService exportService;

    public SimulationBackendService()
        : this(new SimulatedImspExportService())
    {
    }

    public SimulatedImspExportResult SimulateMetaMorpheusToImsp(
        string mzmlPath,
        string psmTsvPath,
        string outputDirectory,
        SimulatedImspExportOptions? options = null)
    {
        return exportService.ExportFromMetaMorpheusPsmTsv(
            mzmlPath,
            psmTsvPath,
            outputDirectory,
            options);
    }
}
```

Frontend should receive paths to files and load/stream them as needed.

## API Evolution Path

V1:

- one combined IMSP
- region-level metadata
- synchronous API
- mzML + MetaMorpheus `.psmtsv`

V2:

- progress reporting
- cancellation token
- exact contribution metadata option
- on-demand contributor lookup
- real + simulated export bundle

V3:

- native app cache/index
- compressed metadata
- viewport-specific metadata loading
- multi-result overlays

## Open Questions

1. Should metadata live in JSON initially, or MessagePack/SQLite for larger datasets?
2. Should region metadata use scan indices, one-based scan numbers, or both?
3. Does the frontend need click-to-identify immediately, or is proteoform search/highlight enough?
4. Should q-value filtering default to `<= 0.01` for all result files?
5. Should simulation output include both centroided and profile modes?
6. Should output directory be caller-specified only, or should service default next to the input mzML?
7. Should the result service write a small summary JSON separately from full metadata?

## Recommended First PR

Keep the first PR narrow:

1. Add DTOs/options/result types.
2. Add MetaMorpheus candidate loader and deduper.
3. Add metadata DTOs and JSON writer.
4. Add tests for loader/deduper/metadata serialization.

Do not implement the full simulation endpoint until these low-level pieces are stable.

## Recommended Second PR

Move the current `AnalysisExample` simulation orchestration into `SimulatedImspExportService`.

Use the existing hard-coded test as a validation harness, but have it call the new service.

## Recommended Third PR

Expose the TauriCS backend command and wire the frontend to:

- start simulation
- receive output paths
- load the `.imsp`
- load sidecar metadata
- highlight proteoform regions
