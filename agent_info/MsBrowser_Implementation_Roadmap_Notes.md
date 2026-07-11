# MsBrowser Implementation Roadmap Notes

## Sample File Sizes

For v1, aim for 4 fixture tiers rather than 3:

1. Small: `100 KB` to `1 MB`
2. Medium: `5 MB` to `25 MB`
3. Large: `50 MB` to `150 MB`
4. Stress: `250 MB+` if feasible

Why these ranges:

- Small lets us validate correctness quickly and makes tests fast.
- Medium is where UI responsiveness and worker boundaries start to matter.
- Large exposes real-world rendering and memory issues.
- Stress helps us avoid designing ourselves into a corner.

If you only want three files, use:

1. Small: around `500 KB`
2. Medium: around `10 MB`
3. Large: around `75 MB`

What matters more than exact size is coverage of:

- low scan count vs high scan count
- sparse vs dense peaks
- narrow RT range vs long run
- simple vs noisy spectra

## Fixture Files

Fixture files are stable sample inputs you keep in the repo or in a known test-data location to drive development and testing.

For this project, fixtures would be:

- real `.imsp` sample files
- tiny synthetic `.imsp` files with known values
- expected outputs for those files, such as:
  - expected header values
  - expected scan count
  - expected TIC points
  - expected spectrum for a chosen scan
  - expected m/z range query results

Why they matter:

- they let us verify the parser automatically
- they prevent regressions when we optimize later
- they give us confidence when integrating the future C# backend

Ideal fixture set:

1. `tiny-known.imsp`
   - manually understandable
   - very small
   - used for exact parser/query tests

2. `small-real.imsp`
   - real-world structure
   - used for UI development

3. `medium-real.imsp`
   - used for performance smoke tests

4. `large-real.imsp`
   - used for profiling and scalability checks

## Versioned Note

A versioned note is a short changelog entry attached to the format spec whenever the format or its interpretation changes.

This matters because binary format misunderstandings become expensive fast, especially once:

- the frontend depends on the format
- the C# generator writes files
- tests and production files already exist

A simple way to do this inside `IMSP_Format.md` is:

```md
## Version Notes

### v1
- Initial format definition
- Scan table entries are 16 bytes
- Canonical extension is `.imsp`

### v1 Clarifications
- 2026-04-11: Corrected section summary that incorrectly listed scan table entries as 12 bytes
- 2026-04-11: Confirmed `OneBasedScanNumber` is uint32
```

This does not need to be elaborate. The goal is to record:

- what changed
- when
- whether it changed the binary format or only clarified docs

Useful categories:

- `format change`
- `documentation correction`
- `compatibility note`

For now, even a small section at the bottom of the spec is enough.

## C# Backend Implications

Your future C# backend should shape the frontend architecture now, but it does not block v1.

What we should do now:

1. Keep `.imsp` parsing/query logic isolated behind interfaces
2. Design app data access around queries, not file internals
3. Make local file loading just one adapter

That way later we can support:

- local browser file
- backend-served file metadata
- backend query endpoints
- desktop local file access via Tauri

Longer term, your frontend should speak in terms of operations like:

- `getScanSummaries`
- `getSpectrumForScan`
- `getPeaksInMzRange`
- `getTicTrace`

Those can be backed by:

- local TypeScript parser now
- C# service later

This is why the UI should depend on a dataset/query interface rather than directly on `DataView`.

## Modern UX Direction

Target:

- modern scientific app, not QualBrowser clone
- cleaner spacing
- lighter chrome
- clearer hierarchy
- compact but readable controls
- desktop-first layout
- light theme only

Good v1 style direction:

- top app bar
- stacked resizable panels
- compact panel headers
- pin/reset controls in the header
- subtle gridlines
- clear hover readouts
- restrained color palette

Preserve the workflow from QualBrowser, not its visual style.

## What You Need To Do Next

Before implementation starts:

1. Scan-table size clarified
2. Canonical `.imsp` extension reflected in the spec
3. A few fixture files prepared
4. If possible, one tiny synthetic file with known expected values
5. A note in the spec indicating current version and clarifications

That is enough to start efficiently.

## Concrete Suggestion For Your Fixture Set

Create these 4 if possible:

1. `tiny-known.imsp`
   - `< 50 KB`
   - handcrafted or generated from minimal known data

2. `small-real.imsp`
   - `300 KB` to `1 MB`

3. `medium-real.imsp`
   - `5 MB` to `15 MB`

4. `large-real.imsp`
   - `50 MB` to `100 MB`

If storage or sharing is inconvenient, keep only the tiny file in-repo and store larger files outside the repo.

## Recommendation

Once the spec is updated and fixtures are prepared, the next step is to produce a concrete build sequence with:

1. repo structure
2. package responsibilities
3. exact first implementation milestones
4. suggested interfaces
5. testing plan tied to the fixtures
