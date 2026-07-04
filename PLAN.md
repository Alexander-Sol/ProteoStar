# FlashLFQ → Rust Port: Execution Plan (ralph loop task list)

This is the durable state for the automated loop in `ralph-loop.ps1`. Each fresh instance
reads this file, does the **single next unchecked task**, records progress here, commits, and
stops. Design rationale lives in `agent_info/FlashLFQ-Rust-Rewrite-Feasibility.md` — read the
relevant section before starting a task.

## Rules for each instance

- Do the **first unchecked `- [ ]` task** in order. Tasks are dependency-ordered; do not skip
  ahead. Port order **is** test order (a chemistry layer must pass parity before the layer
  that depends on it is written).
- A task is **done** only when its acceptance criterion is met (it builds, and its parity/test
  check passes where one is specified).
- When done: change `- [ ]` to `- [x]`, add a one-line result note under the task, update
  **Current focus** below, then commit (do not push).
- If context fills before finishing: write a precise resume note under the task and in
  **Current focus**, commit the partial work, and stop. Do not start a new task.
- If every task is checked, output the single token `ALL_TASKS_COMPLETE` and change nothing.
- Use serena to query the mzLib codebase instead of grep, whenever possible

## Status legend

- `- [ ]` not started · `- [~]` in progress (see note) · `- [x]` done

## Current focus

> _Next up:_ **ALL PLAN TASKS COMPLETE — the next instance should output `ALL_TASKS_COMPLETE` and
> change nothing.** P3.4 was the final task and is now **done and green** (see below). Every checkbox
> from Phase 0 through Phase 3 is `- [x]`. Optional follow-ups, none required by the plan: (1) a faithful
> C# replica / `InternalsVisibleTo` golden for the private+parallel MBR orchestration to upgrade
> P3.2d/P3.3/P3.4 from structural parity to cell-parity; (2) surface `mbr_q_value` (and `mbr_pep`) into
> the Python-returned feature table as trailing nullable columns + wire `apply_mbr_fdr` into the binding.
>
> _P3.4 result (just completed — **Phase 3 / the whole plan complete**):_ **MBR peaks now carry
> q-values.** Ported `FlashLfqEngine.CalculateFdrForMbrPeaks`/`EstimateFdr`/`EstimateDecoyPeptideErrors`/
> `CorrectQValues` (`FlashLfqEngine.cs:1426–1528`) into `src/mbr_search.rs`:
> `calculate_fdr_for_mbr_peaks(&mut Vec<MbrChromatographicPeak>, use_pep)` orders the acceptor's peaks
> (`use_pep=true` → dedup to the best acceptor per donor by `OrderBy(MbrPep).ThenByDescending(MbrScore)`,
> dropping the rest; `use_pep=false` → `OrderByDescending(MbrScore)`, drop nothing), walks the 4-way
> `(decoy_peptide, random_rt)` bucket counts, and assigns each peak the corrected (monotone, 6-dp,
> half-even-rounded) FDR q-value `(1 + decoyPeaks + max(0, decoyPeptides − doubleDecoys)) / total`.
> `apply_mbr_fdr(&mut MbrResult, use_pep)` runs it per acceptor file; `mbr_pep_analysis_succeeded` is the
> `RunPEPAnalysis` `>100`-peaks/`>20`-decoys gate that sets `use_pep`. q-values live on the peaks
> (`MbrChromatographicPeak.mbr_q_value`), a **pure-Rust** step over `mbr_peaks_by_file` exactly as the
> prior note specified; `feature_rows` (the pre-FDR PEP-training table) is untouched. **137 core unit (6
> new) + L1–L6 + P2 + both MBR gates + corpus smoke (now 2 tests, incl.
> `mbr_fdr_assigns_monotone_qvalues_to_every_peak`) all green; full workspace (incl flashlfq-py) builds.**
> **Caveat:** still **no C# golden** for the private+parallel MBR orchestration, so this is faithful
> structural/algorithmic parity (line-by-line port of the formula + ordering), not live-dump cell-parity.
>
> **Caveat carried forward (still applies for any P3.2 parity follow-up):** there is **no C# golden for
> the MBR orchestration yet** — `QuantifyMatchBetweenRunsPeaks` is private + parallel, so exact row
> parity needs a faithful C# replica (or `InternalsVisibleTo`/`Test`-named assembly instrumenting a full
> MBR engine run). The corpus smoke test + P3.2e Python acceptance assert the table is well-formed and
> serializes correctly, but not cell-parity. Fold a parity golden into a dedicated follow-up if wanted.
>
> **Golden-generator note for a P3.2d parity follow-up:** there is **no C# golden for the orchestration
> yet** — `QuantifyMatchBetweenRunsPeaks` is private + parallel, so exact row parity needs a faithful C#
> replica (or `InternalsVisibleTo`/`Test`-named assembly instrumenting a full MBR engine run) the way
> P3.1/P3.2a did for the spline/predicted-RT. The corpus smoke test asserts the table is well-formed
> (targets+decoys, component scores in (0,1], combined in [0,100]) but not cell-parity. Fold this into
> P3.2e or a dedicated follow-up if exact parity is wanted.
>
> **Golden-generator note (still applies):** the MBR methods (`GetRtCalSpline`, `PredictRetentionTime`)
> are private/internal, so the C# generator (`parity/csharp_golden/Program.cs`) **replicates** them
> over the real `results.Peaks` (genuine engine ground truth) rather than invoking them — P3.2c/d
> will likely need the same replica approach (or instrument a full MBR engine run via an
> `InternalsVisibleTo`/`Test`-named assembly). The corpus-as-two-files (K562_3/_4) drives a real
> donor/acceptor pair (245 donor best-peaks, 113 shared anchors) without needing the
> `PSMsForMbrTest.psmtsv` + `f1r*_sliced_mbr.raw` fixtures. The MBR Python switch is already wired
> (`flashlfq_py.quant(..., match_between_runs=True)` → `PyNotImplementedError`, P1.20); Phase 3
> makes it real, then the Python ML model (P3.3) and MBR q-value FDR (P3.4).
>
> _P2.2 result (just completed — **Phase 2 complete**):_ **Decoy / quantify-set semantics made
> faithful to C#, gated by a real-C# decoy parity test.** Split the two conflated C# sets in
> `src/results.rs`: **`engine_quantify_set`** (all modseqs incl decoys = `FlashLfqEngine.cs:94`
> default, the membership gate for `run_error_checking` + `calculate_peptide_results`) and
> **`output_peptide_sequences`** (non-decoy modseqs in the quantify set = `FlashLfqResults.cs:42`
> `PeptideModifiedSequences` keys, the emitted rows). `calculate_peptide_results` now takes both;
> `engine::run_msms` wires them. `default_quantify_set` retained as the documented `FlashLfqResults`
> null-set *fallback* only. **q-value:** MS2 path computes none (FDR = the externally-supplied
> quantify set; engine default applies no q-filter) — MBR q-value FDR is P3.4. **91 core unit** (4 new
> decoy tests) + new **`l6_decoy_flip_matches_csharp_golden`** gate: golden generator flips a
> deterministic third of modseqs to decoy (118/354, sorted-Ordinal `index%3==0`), re-runs the real
> `FlashLfqEngine` → `golden/L6_peptide_intensities_decoy.tsv` (472 rows = 236×2); Rust applies the
> identical flip and **all 472 cells match** + no decoy modseq emitted. All-target L5/L6 goldens
> regenerated byte-identical; all gates green; full workspace (incl flashlfq-py) builds.
>
> _P1.20 result (just completed — **Phase 1 complete**):_ **MBR stub binding green.** Added a
> `match_between_runs=False` kwarg to `flashlfq_py.quant` (`rust/flashlfq-py/src/lib.rs`), the
> permanent home for the MBR switch (mirrors C# `FlashLfqParameters.MatchBetweenRuns`). `True` →
> `PyNotImplementedError` ("…not implemented until Phase 3; pass match_between_runs=False…"), it
> does not silently fall back to MS2-only. Default `False` leaves the P1.18 path untouched (guard
> sits ahead of the `params`/`core_quant` logic). Verified after `maturin develop`. Thin binding;
> no core changes. (Full details in the P1.20 task below.)
>
> _P1.19 result (just completed):_ **IndexingEngine/XIC Python bindings green — plots a real
> peptide's XIC.** Two new PyO3 `#[pyclass]`es in `rust/flashlfq-py/src/lib.rs`: **`PeakIndex`**
> (wraps the core `PeakIndexingEngine`, built once + reused) — `PeakIndex.from_mzml(path)`,
> `num_ms1_scans`, `get_indexed_peak(mz, scan_index, ppm=5.0) -> tuple|None`,
> `get_xic(mz, rt, ppm=20.0, missed_scans_allowed=1, max_peak_half_width=None) -> Xic`; and
> **`Xic`** — the RT-ascending trace as three parallel `f64` NumPy arrays
> (`retention_times`/`intensities`/`mzs`) + integration bounds (`start_rt`/`apex_rt`/`end_rt`,
> `apex_intensity`, `apex_scan_index`; empty → NaN/-1), `__len__`/`__repr__`. Acceptance
> `rust/flashlfq-py/xic_smoke_test.py` builds the index for `…K562_3` (73 MS1 scans), traces an
> 18-peak XIC for `AHQLVMEGYNWC[Carbamidomethyl]HDR` z=3 (apex RT 79.13, 2.06e6), validates the
> arrays/bounds + point-query round-trip + empty-XIC behaviour, and writes `xic_smoke_plot.png`
> → `ALL XIC SMOKE CHECKS PASSED`. (Full details in the P1.19 task below.)
>
> _P1.18 result (previously completed):_ **`quant()` Python binding green — callable end-to-end.** New
> core `engine::quant(results_path, raw_paths) -> Result<Ms2QuantResult, QuantError>` (parse psmtsv →
> resolve spectra files by **bare file name** → delegate to `run_msms`) + `quant_to_parquet`; new
> `QuantError` enum; `psm_tsv::file_name_without_extension` made `pub`. PyO3
> `flashlfq_py.quant(results_path, raw_paths, output_path=None, params=None)` → Parquet path `str`
> (when `output_path` given) or in-memory `pyarrow.RecordBatch` (when `None`); `params` rejected as
> unsupported in Phase 1. End-to-end acceptance `rust/flashlfq-py/quant_smoke_test.py` over the K562
> corpus → **708 rows, detection mix 450/223/33/2**, Parquet == in-memory, error surfaces fire. 87
> core unit + L1–L6 all green. (Full details in the P1.18 task below.)
>
> _P1.17b result (just completed — **Phase 1 complete**):_ **L5 + L6 parity gate GREEN against the
> real C# golden.** New `rust/flashlfq-core/src/engine.rs` ports the MS2 path of `FlashLfqEngine.Run`:
> `run_msms(ids, file_to_mzml)` wires `calculate_theoretical_isotope_distributions` (per-modseq
> envelope cache + engine-global charge range `[min..=max]` + per-id peakfinding mass; `OptionalChemicalFormula`
> is always `None` — confirmed `MakeIdentifications` never sets it) → per-file `quantify_ms2_identified_peptides`
> (per id: `from_identification`, per charge `get_xic(pfm.to_mz(z), ms2_rt, PpmTol(20), missed=1)` →
> 10ppm mass `retain` → `get_isotopic_envelopes`; then `calculate_intensity_for_this_feature(false)`,
> `cut_peak(ms2_rt, false)`, precursor-charge scan-range trim, recompute — empty peaks retained like C#)
> → per-file `run_error_checking` → `calculate_peptide_results`. Constants verbatim (PEAKFINDING_PPM=20,
> PPM=10, ISOTOPE_PPM=5, NUM_ISOTOPES=2, MISSED_SCANS=1, INTEGRATE=false, ID_SPECIFIC_CHARGE=false,
> maxPeakHalfWidth=`i32::MAX`). New gate `tests/l5_l6_parity.rs` drives the live `AllPSMs.psmtsv` + two
> K562 mzMLs, sorts rows with the **same comparator as the C# generator** (L5: file/modseq/apexRT(NaN
> last)/intensity; L6: modseq/file), diffs vs `golden/L5_chromatographic_peaks.tsv` (507) +
> `golden/L6_peptide_intensities.tsv` (708): ints/scans/charges/counts + strings exact, floats rel-1e-6
> (NaN==NaN). **L6 detection-type mix asserted exactly: 450 MSMS / 223 NotDetected / 33
> MSMSIdentifiedButNotQuantified / 2 MSMSAmbiguousPeakfinding.** 85 core unit + L1–L6 all green.
>
> **THE BUG that the gate caught (root cause + fix — important for later readers):** the first run
> diffed at 506/507 L5 + 6/708 L6 cells, all cases of Rust finding *extra/larger* peaks. Traced via a
> C# probe (public `engine.GetXic`/`GetIsotopicEnvelopes`, assembly named `Test` for
> `InternalsVisibleTo`): mzdata returns the **raw profile** mzML including zero-intensity sample
> points, but **mzLib's MzML reader drops peaks with intensity `< 0.01`** ("Remove Zero Intensity
> Peaks", `Readers/MzML/Mzml.cs`). Those retained zeros let the Rust XIC walk **cross gaps mzLib stops
> at** (a kept zero counts as an in-tolerance hit, resetting the missed-scan counter), so Rust traced
> ~33 scans past the real elution into a co-eluting interference cluster and mis-apexed there. **Fix:**
> `read_ms1_scans` now applies the identical `< 0.01` drop (`remove_zero_intensity_peaks`, faithful to
> mzLib's all-zero-scan-kept guard). The L4 reader-calibration test was reconciled to filter its
> raw-decode golden the same way (it now asserts "Rust reader == mzLib reader", not "== raw file").
> **L4's old claim that mzLib does "no zero-trimming" was wrong** — it does, and that filter is
> load-bearing for tracing parity. Probe project + a throwaway debug test were deleted after use.
>
> _Earlier next-up note (P1.17b plan, retained for reference):_ **Rust driver + L5/L6 parity gate.** The C# real-engine golden now exists
> (P1.17a, below) — **the no-.NET-SDK caveat is retired: `dotnet` (SDKs 3.1/6/8/10) IS available in
> this env now.** What remains: **(1) write the Rust driver** = a faithful port of the MS2 path of
> `FlashLfqEngine.Run` wiring the already-ported pieces, and **(2) a parity test** diffing it against
> the two new goldens at rel-1e-6 → **Phase 1 complete.**
>
> **Driver order (port of `FlashLfqEngine.cs:188 Run` MS2 path, default params — read these exact
> methods):** `CalculateTheoreticalIsotopeDistributions` (`FlashLfqEngine.cs:351`) then per-file
> `QuantifyMs2IdentifiedPeptides` (`:472`) then `RunErrorChecking` (already = `results::run_error_checking`)
> then `CalculatePeptideResults` (already = `results::calculate_peptide_results`). Concretely:
> 1. Parse `AllPSMs.psmtsv` → `Vec<Identification>` (`psm_tsv::read_identifications`). The corpus
>    spans **two** files (`20100614_Velos1_TaGe_SA_K562_3` and `_4`); group ids by `file_name`.
> 2. Build per-modseq theoretical envelope = `expected_isotope_peaks(None, base_seq, mono, 2)` keyed
>    by modified sequence (cache; one per distinct modseq). Compute `peakfinding_mass = mono +
>    most_abundant_isotope_shift` per id (the `mostAbundantIsotopeShift` = the massShift where
>    normalized abundance == 1.0). Compute the charge-state range `_chargeStates =
>    [min..=max precursor charge]` across ALL ids (engine-global, not per file).
> 3. For each file: `PeakIndexingEngine::from_mzml(<that file's mzML in TestData>)`. Then for each id
>    in that file build a `ChromatographicPeak::from_identification`; for each charge in `_chargeStates`
>    (IdSpecificChargeState defaults **false** → all charges, not just the id's): `get_xic(peakfinding_mass.to_mz(charge),
>    ms2_rt, PpmTolerance(20), missed=1)` sorted by RT; **RemoveAll where `!PpmTolerance(10).within(peak.M.to_mass(charge),
>    peakfinding_mass)`**; `get_isotopic_envelopes(engine, &xic, &envelope, mono, peakfinding_mass,
>    charge, num_isotopes_required=2, isotope_ppm_tol=PpmTolerance(5))`; append envelopes to the peak.
>    Then `calculate_intensity_for_this_feature(integrate=false)`; `cut_peak(ms2_rt, false)`; if no
>    envelopes → skip; `precursor_xic` = envelopes with charge == id.precursor_charge; if empty → clear
>    + skip; min/max scan over precursor_xic; drop envelopes outside `[min,max]`; recalc intensity.
>    Collect the per-file `Vec<ChromatographicPeak>`.
> 4. `run_error_checking(file_peaks, &quantify_set)` per file (quantify_set = `default_quantify_set`
>    = non-decoy modseqs). Then `calculate_peptide_results` over all files → `PeptideResults` (L6).
> **Watch:** GetXic/envelope APIs already exist — check signatures in `peak_indexing.rs`,
> `isotopic_envelope.rs`, `chromatographic_peak.rs`, `results.rs` before wiring (some take `&engine`).
> The C# `Identification.PeakfindingMass` is set in step 2; the Rust `Identification` may need a
> peakfinding-mass field or a parallel map. The driver likely belongs in a new
> `src/engine.rs` (e.g. `fn run_msms(ids, file_to_mzml_path) -> (Vec<peaks per file>, PeptideResults)`).
>
> **Golden schema (already written by P1.17a):**
> - `parity/golden/L5_chromatographic_peaks.tsv` — 507 rows, header: `file_name base_sequence
>   modified_sequence monoisotopic_mass ms2_retention_time precursor_charge peak_intensity
>   peak_rt_start peak_rt_apex peak_rt_end peak_mz apex_charge apex_scan_index apex_pearson
>   num_charge_states mass_error split_rt num_ids_by_full_seq detection_type decoy`. Rows sorted by
>   (file_name, modified_sequence, peak_rt_apex, peak_intensity) — **sort the Rust peaks the same way**.
>   Compare: ints/scan-index/charge exact; floats rel-1e-6 (intensity, RTs, mass_error, pearson, mz).
> - `parity/golden/L6_peptide_intensities.tsv` — 708 rows (354 modseq × 2 files), header:
>   `modified_sequence file_name intensity retention_time detection_type`. = exactly
>   `parquet_output` columns. Compare intensity rel-1e-6, detection_type string-exact.
> **Detection-type mix in golden (sanity targets):** L5 all 507 MSMS; L6 = 450 MSMS / 223 NotDetected /
> 33 MSMSIdentifiedButNotQuantified / 2 MSMSAmbiguousPeakfinding.
>
> _P1.17a result (just completed):_ **Real C# golden generated — the no-SDK caveat is retired.**
> `.NET SDK is now present` (`dotnet --list-sdks` → 3.1/6.0/8.0/10.0). New **C# console golden
> generator** `rust/flashlfq-core/parity/csharp_golden/` (`Golden.csproj` + `Program.cs`,
> net8.0, references `..\..\..\..\mzLib\FlashLFQ\FlashLFQ.csproj`, **not** in mzLib.sln so it never
> touches the normal build; bin/obj already covered by the root .gitignore). It reads `AllPSMs.psmtsv`
> via the real mzLib path (`FileReader.ReadQuantifiableResultFile` → `MzLibExtensions.MakeIdentifications`
> with `SpectraFileInfo`s for the two K562 mzMLs), runs the **real `FlashLfqEngine`** with default
> `FlashLfqParameters` (MS2 path: no MBR/IsoTracker/normalize; `MaxThreads=1` for full reproducibility),
> and dumps `parity/golden/L5_chromatographic_peaks.tsv` (507 peaks) + `parity/golden/L6_peptide_intensities.tsv`
> (708 cells) with round-trippable (`"R"`) double formatting, deterministically sorted. **594 ids /
> 354 distinct modseqs / 2 files** confirmed. Regenerate:
> `cd rust/flashlfq-core/parity/csharp_golden; dotnet run -c Release`. **This is genuine C# ground
> truth** (not a Python replica like L0–L4) — reconcile L0–L4 against real C# too if ever desired, but
> not required. _Driver + gate (P1.17b) deferred to the next instance to keep this commit clean._
>
> _P1.16 result (just completed):_ **Parquet output — green, loads in Python via pyarrow.** New
> `rust/flashlfq-core/src/parquet_output.rs` writes `results::PeptideResults` to Parquet in **tidy/long
> form** (one row per (modified_sequence, file_name): `modified_sequence`/`file_name`/`intensity`/
> `retention_time`/`detection_type` columns), rows **sorted by (modseq, file)** for reproducibility.
> `write_peptide_results_parquet(results, path)` uses `parquet::arrow::ArrowWriter`. Added
> `DetectionType::as_str()` (verbatim C# enum names) for the detection-type column, and dep
> `parquet = "54"` (lockstep with `arrow = "54"`). **2 new tests (84 core unit + L1–L4 green):** sort/
> shape + an Arrow-reader round-trip. **Python acceptance:** example
> `examples/write_demo_parquet.rs` + reading with `rust/.venv` pyarrow returns exact schema/values.
>
> _P1.15 result (previously completed):_ **Charge aggregation + FDR ported — green, 82 core unit +
> L1–L4.** `ChromatographicPeak` now carries `identifications`/`detection_type`/distinct-seq counts
> (P1.14 envelope fields untouched, CutPeak unaffected); ports `ResolveIdentifications`,
> `MergeFeatureWith` (structural `EnvelopePeakKey` for the C# `IIndexedPeak` reference identity),
> `ApexRetentionTime`, `DecoyPeptide`. New `detection_type.rs` (verbatim enum) + `results.rs`:
> `run_error_checking` = MSMS path of `FlashLfqEngine.RunErrorChecking` (apex grouping + merge/
> overwrite/drop by quantify-set membership, insertion-ordered output); `calculate_peptide_results`
> = MSMS path of `FlashLfqResults.CalculatePeptideResults` → `PeptideResults` (modseq→file→
> `PeptideQuant`), the peptide×file table for L6; `default_quantify_set` = non-decoy modseqs.
> **FDR (Phase 1) = the quantify-set filter**; no per-peak q-value is computed (MBR q-value FDR is
> Phase 3). MBR/IsoTracker branches + `HandleAmbiguityInFractions` intentionally deferred. 10 new
> tests. _Caveat:_ synthetic fixtures (no .NET SDK); peptide×file corpus parity folded into L6 at
> P1.17.
>
> _P1.14 result (previously completed):_ **`CutPeak` ported — green, builds and splits bimodal peaks.**
> New module `rust/flashlfq-core/src/chromatographic_peak.rs`. Introduces `ChromatographicPeak`
> (the `Vec<IsotopicEnvelope>` container) = partial port of `FlashLFQ.ChromatographicPeak`: fields
> `isotopic_envelopes`, `identification_peakfinding_masses` (for the mass-error term),
> `intensity`, `apex: Option<IsotopicEnvelope>`, `split_rt`, `mass_error`, `num_charge_states_observed`.
> `calculate_intensity_for_this_feature(integrate)` ports `CalculateIntensityForThisFeature`
> (apex = first max-intensity envelope per LINQ `MaxBy`; intensity = sum if integrate else apex;
> mass error = smallest-|·| `((ToMass(Apex.M,z) − pfm)/pfm)·1e6` over the id peakfinding masses,
> seeded NaN; distinct charge-state count). `cut_peak(identification_time, integrate)` ports
> `FlashLfqEngine.CutPeak` faithfully: ≥5-envelope guard, apex-charge `timePointsForApexZ` subset,
> the `{+1,−1}` valley walk (running lowest-intensity valley; cut when discrimination factor vs the
> valley **and** vs the point adjacent to the valley both exceed `DISCRIMINATION_FACTOR_TO_CUT_PEAK
> = 0.6`, **or** no envelope in the scan immediately past the valley), then `RemoveAll` of envelopes
> on the far side of the valley RT from the MS2 id (`identificationTime > valleyRT` → drop ≤ valley;
> else drop ≥ valley), intensity recompute, `SplitRT = valleyRT`, recursive re-cut. **Fidelity
> calls:** valley index = the loop index `i` directly (C# `IndexOf` is reference equality, the valley
> *is* `time_points[i]`); apex index found by value equality (the global-max envelope is unique by
> scan index); `RemoveAll` filters **all** charge states by RT even though the valley came from the
> apex-charge subset (matches C#); RT compared as `f32`-widened-to-`f64` (lossless). Made
> `isotopic_envelope::mz_to_mass_f32` `pub` and reused it for the apex mass error. **`ExtractedIonChromatogram.CutPeak`/`FindPeakBoundaries`/`RemovePeaks` were NOT needed** — the
> active FlashLFQ peak-finding path (`FlashLfqEngine.cs:528`) uses the engine's envelope-level
> `CutPeak`, not the XIC-level one; the XIC `CutPeak` is a separate (iso-tracker) path, left for a
> later task if needed. **7 new tests (72 core unit + L1+L2+L3+L4 green):** discrimination constant;
> apex/sum intensity + charge count; empty-resets; no-cut <5 envelopes; no-cut clean unimodal;
> bimodal split keeping the id (left) side; mirror split keeping the right side. _Caveat:_ synthetic-
> fixture ground truth (no .NET SDK); corpus/parity folded into L5 at P1.17.
>
> _P1.13 result (just completed):_ **Isotopic-envelope search + correlation gate — green, builds and
> returns envelopes.** New module `rust/flashlfq-core/src/isotopic_envelope.rs` ports
> `FlashLfqEngine.GetIsotopicEnvelopes` + `CheckIsotopicEnvelopeCorrelation`. `IsotopicEnvelope`
> struct (`indexed_peak`, `charge_state`, `intensity` = summed/charge, `pearson_correlation`) =
> `FlashLFQ.IsotopicEnvelope`. Free fn `get_isotopic_envelopes(engine, xic, isotope_mass_shifts,
> mono, peakfinding_mass, charge, num_isotopes_required, isotope_ppm_tol)` — takes a single
> `&PeakIndexingEngine` (vs C#'s `IndexingEngineDictionary[spectraFile]`) and the
> `&[ExpectedIsotopePeak]` distribution (vs `ModifiedSequenceToIsotopicDistribution[modSeq]`).
> **Ported faithfully:** the `{-1,0,+1}` `massShiftToIsotopePeaks` hypotheses (negative/accurate/
> positive ¹³C off-by-one) as a `[Vec<IsotopePeakRecord>; 3]`; `peakfindingMassIndex =
> round_ties_even(peakfinding − mono)`; the down-then-up peak walk (`direction = -1` from
> `peakfindingMassIndex-1`, then `+1` from `peakfindingMassIndex`), breaking on a missing peak or the
> **`/4`..`*4` intensity-ratio gate** (`theoreticalAbundance * observedPeakIntensity`); the
> `< NumIsotopesRequired` accurate-peak check; correlation gate `pearson>0.7 && corrLeft−corrPadded<0.1
> && corrRight−corrPadded<0.1` with NaN→−1; the "unexpected" padding peak (one ¹³C below each
> hypothesis's lowest expected mass); and the zero-intensity **imputation** + summed-intensity
> envelope. **Conversions:** `peak.M.ToMass` in `f32` then widened (C# `float` overload — `peak.M` is
> f32), `mass.ToMz` all-`f64` (C# `double` overload). **Pearson** = the MathNet `Correlation.Pearson`
> Welford-style online algorithm (running mean/var/cov via differences, not product sums), ported in
> `pearson(&[f64],&[f64])`; empty/constant input → NaN (mapped to −1 by the gate). Constants
> `PROTON_MASS = 1.007276466879`, `C13_MINUS_C12 = 1.00335483810` (mzLib `Chemistry/Constants.cs`).
> Vestigial C# `isotopologuePeaks` list (written, never read) omitted, as with P1.7's
> `highestAbundanceIndex`. **Tests (7 new, 65 core unit + L1+L2+L3+L4 green):** `pearson` identical→1
> / anti→−1 / constant→NaN; mz↔mass round-trip; a clean 3-isotope synthetic envelope (mono most
> abundant, 1.0/0.5/0.2 spaced one ¹³C) traced across a 7-scan elution profile → **1 envelope/scan,
> all pearson>0.7**, apex at the max-multiplier scan with summed intensity `(1+0.5+0.2)*base*mult/z`;
> fewer-isotopes-than-required → empty; same at **charge 2**. _Caveat:_ ground truth is synthetic
> fixtures (no live C# dump; no .NET SDK) — true envelope/correlation parity is folded into the L5
> gate at P1.17 against the corpus.
>
> _P1.12 result (just completed):_ **L4 reader calibration — bit-identical, zero drift.** Golden
> dumper `rust/flashlfq-core/parity/dump_l4_reader_peaks.py` raw-decodes the mzML binary arrays
> (base64 → zlib → 32-bit floats) → `parity/golden/L4_reader_peaks.tsv` (10 MS1 scans, 67,001
> peaks). Report-only test `rust/flashlfq-core/tests/l4_reader_calibration.rs` diffs
> `read_ms1_scans` (mzdata) vs golden: **worst |Δm/z| = 0, worst rel |Δintensity| = 0, worst
> |Δrt| = 0** — exactly bit-for-bit (both decode the same f32 zlib payload and widen f32→f64
> losslessly). File is **profile** data (not centroided), so reader fidelity is total. **No
> reader-level slop needs budgeting into L5/L6 tolerances.** 58 core + L1 + L2 + L3 + L4 green.
>
> _P1.10 result:_ **binned m/z index green.** New module
> `rust/flashlfq-core/src/peak_indexing.rs` ports `IndexingEngine<IndexedMassSpectralPeak>` /
> `PeakIndexingEngine` (the m/z, charge-less flavour): `IndexedMassSpectralPeak` (mz/intensity/RT
> as **f32** per the C# `(float)` casts; `m() == mz`), `ScanInfo`, an input `Scan` struct
> (mz/intensity arrays + scan metadata), and `PeakIndexingEngine` holding the jagged
> `Vec<Option<Vec<peak>>>` index + `Vec<ScanInfo>`. Ported faithfully: `index_peaks` (= C#
> `IndexPeaks`: index length `ceil(maxLastX * 100) + 1`, bin = `(mz*100).round_ties_even()` —
> C# `Math.Round` is **banker's**; peaks appended in scan order so each bin is scan-index-ascending;
> returns `None`/`false` when nothing to index), `get_indexed_peak` (= `GetIndexedPeak`),
> `get_bins_in_range`, `get_best_peak_from_bins`, `get_peak_from_bin`,
> `binary_search_for_indexed_peak` (the binary-search-then-linear-backstep, ported with `isize`),
> and `BINS_PER_DALTON = 100`. Charge handling and the mass-indexing engine were **not** ported
> (m/z peaks only). New module `rust/flashlfq-core/src/tolerance.rs` = `PpmTolerance`
> (`get_minimum_value`/`get_maximum_value`/`within`, faithful to `MzLibUtil/PpmTolerance.cs`).
> mzML loader `read_ms1_scans(path)` + `PeakIndexingEngine::from_mzml(path)` read MS1 scans via
> `mzdata` (`raw_arrays().mzs()/.intensities()`, `start_time()` for RT in minutes, `index()+1`
> for one-based scan number). **Tests:** 10 new peak_indexing unit tests (synthetic 9-scan fixture
> mirroring the C# `IndexingEngineTests` setup: exact-m/z point query, closest-of-adjacent
> resolution, out-of-range → None, absent-scan → None, per-scan presence, bin assignment,
> scan_info) **+ a real-mzML integration test** (`from_mzml("sliced-mzml.mzML")` → 10 MS1 scans,
> pull a real peak out of the index and round-trip-query it back at 5 ppm) + 2 tolerance tests.
> **53 core lib tests + L1 + L2 + L3 integration all green.** _Caveat:_ ground truth for the point
> queries is the synthetic fixture + round-trip (the C# index tests' XIC counts are P1.11's job);
> no live C# dump (no .NET SDK).
>
> _P1.9 result (just completed):_ **psmtsv parse green.** New module
> `rust/flashlfq-core/src/psm_tsv.rs` — `Identification` + `read_identifications(_from_str)` over
> the `csv` crate (`quoting(false)`, tab delim, flexible) faithfully mirroring mzLib
> `SpectrumMatchTsvReader`/`SpectrumMatchFromTsv`/`PsmFromTsv` (plain tab split + own quote/ws
> trim, NOT RFC-4180). Header→column resolution copies `ParseHeader` (no `Accession` column ⇒
> `Peptide Monoisotopic Mass` col 22). Field semantics per the C# ctor (FileName ext-strip,
> `RemoveParentheses`, mono = first `|`-token else `-1.0`, charge `(int)f64`, RT default `-1.0`,
> decoy = DCT contains `'D'`). Corpus test parses the **live** `AllPSMs.psmtsv` → **594 rows / 354
> distinct Full Sequences, all targets** + row-0 spot check. Deps added: `csv = "1.3"`,
> `serde` (derive). **41 core tests + L1+L2+L3 all green.** Section **1b** is underway.
>
> _P1.8 result (just completed):_ **L3 parity gate green — the chemistry gate is complete.** New
> golden dumper `rust/flashlfq-core/parity/dump_l3_theoretical_distributions.py` (faithful Python
> replica of `FlashLfqEngine.CalculateTheoreticalIsotopeDistributions`, layered on the L0 loader +
> L1 residue table + **L2 `get_distribution`** replica, called at `fine_resolution=0.125,
> min_probability=1e-8`) emits `parity/golden/L3_theoretical_distributions.tsv` (D/S/A records:
> descriptor + massShift[] + normAbundance[]) for the **354 distinct Full Sequences** (max 6 peaks).
> **Key fidelity point the port + golden both honor:** the `None`-formula path builds the formula
> from the **base sequence only (no mods)** — N-term H + C-term OH + residues — and tops up with
> averagine when `|mono − baseFormulaMono| > 20` (fired on **44 of 354** peptides, e.g. the
> carbamidomethyl-C set: the +57 mod mass is approximated by averagine, *not* read from the
> `Mods Combined Chemical Formula` column that drives L1/L2). `mono` is the psmtsv
> **`Peptide Monoisotopic Mass`** column (col 22), parsed to `f64` identically on both sides, used
> for re-centering *and* the 20 Da test so averagine fires identically. Rust gate
> `rust/flashlfq-core/tests/l3_parity.rs` drives the distinct Full-Sequence set straight from the
> live psmtsv (asserts set + base + mono match golden), computes
> `expected_isotope_peaks(None, base, mono, 2)` + `peakfinding_mass`, and asserts **array lengths
> exact**, **massShift abs < 1e-6 Da**, **abundance rel < 1e-6**, **peakfinding mass abs < 1e-6 Da**
> — observed **worst massShift abs 9.09e-13 Da, worst abundance rel 5.71e-16, worst peakfinding
> mass abs 4.55e-13 Da** across 1223 total peaks (near bit-identical). **34 core unit tests + L1 +
> L2 + L3 integration all green.** Regenerate golden:
> `& "F:\flashlfq-rust\rust\.venv\Scripts\python.exe" rust/flashlfq-core/parity/dump_l3_theoretical_distributions.py`.
> _Caveat (same as L0–L2):_ ground truth is the Python replica, not a live C# dump (no .NET SDK);
> reconcile against real C# if the SDK ever appears.
>
> _P1.7 result:_ Ported `CalculateTheoreticalIsotopeDistributions` into
> `rust/flashlfq-core/src/theoretical_isotope_distribution.rs` (**34 core unit tests green**, 8
> new; L1 + L2 integration tests still green). Public API: `expected_isotope_peaks(optional_formula:
> Option<&ChemicalFormula>, base_sequence, monoisotopic_mass, num_isotopes_required) ->
> Vec<ExpectedIsotopePeak { mass_shift, normalized_abundance }>`, plus `averagine_mass()`,
> `all_sequence_residues_are_valid()`, `most_abundant_isotope_shift()`, and `peakfinding_mass()`.
> Constants exported verbatim: `AVERAGE_{C,H,O,N,S}` (4.9384/7.7583/1.4773/1.3577/0.0417),
> `FINE_RESOLUTION=0.125`, `MIN_PROBABILITY=1e-8`, `DEFAULT_NUM_ISOTOPES_REQUIRED=2`, 20 Da
> averagine threshold. **Key fidelity decisions:** (1) `Math.Round(x,0)` is banker's rounding →
> used Rust `f64::round_ties_even` (stable since 1.77; toolchain 1.94.1), unit-tested 2.5→2,
> 3.5→4, 0.5→0. (2) The re-center keeps the C# **two-step** `masses[i] += (mono - formulaMono)`
> then `masses[i] -= mono` (not folded to `dist - formulaMono`) so the FP matches bit-for-bit for
> the L3 gate. (3) Distribution uses the **non-default** `(0.125, 1e-8)` params with mwRes at the
> default `1e-12` — matching the C# `GetDistribution(formula, 0.125, 1e-8)` 3-arg overload.
> (4) `averagine_mass()` ≈ 111.1236 Da (computed from L0 average masses — the classic Senko
> averagine residue). Truncation rule (`peaks.len() < num_required || abundance > 0.1`) and the
> normalize-to-max are faithful. The vestigial C# `highestAbundanceIndex` (computed, unused) is
> omitted. _Caveat (same as L0–L2):_ no .NET SDK in env; ground truth will be the Python replica
> (P1.8), reconcile against real C# if the SDK appears.
>
> _P1.6 result (just completed):_ L2 parity gate **green** — `IsotopicDistribution::
> get_distribution` matches the Python Kubinyi replica to **worst mass abs 1.36e-12 Da, worst
> intensity rel 1.86e-15** with **array lengths exact** across all 354 peptides / 257,162 peaks.
> Key decision that made index-aligned comparison work: the golden dumper
> (`parity/dump_l2_isotopic_distributions.py`) iterates elements in **atomic-number order**
> (`sorted(counts)`) to match Rust's BTreeMap, so the FP-order-sensitive convolution lines up
> peak-for-peak. Golden: `parity/golden/L2_isotopic_distributions.tsv` (D/M/I records). Test:
> `tests/l2_parity.rs`. **27 core tests green** (25 unit + L1 + L2).
>
> _P1.5 result:_ Ported `Chemistry/IsotopicDistribution.cs` into
> `rust/flashlfq-core/src/isotopic_distribution.rs` (**25 core unit tests green**, 5 new; L1
> integration test still green). Public API: `IsotopicDistribution { masses, intensities }` with
> `get_distribution(&ChemicalFormula)` and `get_distribution_with(formula, fine_res, min_prob,
> mw_res)`; constants ported exactly (`fine 0.01`, `min_prob 1e-200`, `mw_res 1e-12`; the
> `fineResolution/2` split). Faithful ports of `MultiplyFinePolynomial`,
> `MultipleFinePolynomialRecursiveHelper`, `MultiplyFineFinalPolynomial`, `MergeFinePolynomial`
> (9-pass k-scaled threshold — including the subtle detail that the inner-loop mass comparison
> reads the **running merged centroid** `t_polynomial[i].power`, updated in place, not the
> original), `CalculateFineGrain`, and `FactorLn` (memoised thread-local cumulative log-sum,
> bit-identical summation order to the C# `factorLnArray` cache). `Polynomial`/`Composition` use
> `f64::NAN` sentinels matching the C# structs. Added `ChemicalFormula::elements()`/`isotopes()`
> iterators. **Known limitation (matches C#):** a formula with *only* isotope-specified atoms and
> no plain elements (e.g. `"C{13}5"`) panics in `MultiplyFinePolynomial` — C# throws
> `IndexOutOfRange` on the identical `tPolynomial = fPolynomial[0]`. Real peptide formulas always
> carry backbone elements, so this is never hit; left as-is for parity.
>
> _P1.4 result (just completed):_ The feared "resolve a modified-sequence string to its mod
> formulas" problem is **moot for this corpus** — `AllPSMs.psmtsv` already carries a
> **`Mods Combined Chemical Formula`** column (col 19; distinct values: ``''``, `C2H3NO`, `O`,
> `C2H3NO2`, `C4H6N2O2`) that is the summed mod formula mzLib's `FullChemicalFormula` adds. So
> L1 full formula = `peptide_base_formula(Base Sequence)` + `parse(Mods Combined Chemical
> Formula)`. Golden generator `rust/flashlfq-core/parity/dump_l1_peptide_formulas.py` (Python
> replica reusing the L0 loader) emits `parity/golden/L1_peptide_formulas.tsv` (354 distinct
> peptides; element counts + mono mass + the C# `Peptide Monoisotopic Mass` column). Rust gate is
> the integration test `rust/flashlfq-core/tests/l1_parity.rs`: drives the distinct peptide set
> straight from the live psmtsv, recomputes each formula, and asserts **counts exact**, **mono
> rel-1e-9 vs golden** (worst 2.85e-16), and **abs < 1e-5 vs the C# psmtsv column** (worst
> 4.98e-6 Da = the 5-decimal print-rounding bound, confirming the port matches real C# output).
> **21 core tests green** (20 unit + 1 integration). Regenerate golden if the psmtsv changes:
> `& "F:\flashlfq-rust\rust\.venv\Scripts\python.exe" rust/flashlfq-core/parity/dump_l1_peptide_formulas.py`.
> _Resume notes (env, all validated in P0.5):_ **`.NET SDK is NOT installed`** in this env (`dotnet`
> missing) and **`python`/`py` are not on PATH for the Bash tool** — use the P0.5 venv interpreter
> directly: `& "F:\flashlfq-rust\rust\.venv\Scripts\python.exe" <script>` from PowerShell. Run cargo tests
> with PowerShell: `Set-Location F:\flashlfq-rust\rust; cargo test -p flashlfq-core`. Toolchain
> cargo/rustc 1.94.1. **Python binding
> stack fully works** via `rust/.venv` (Python 3.13.7; maturin 1.14.1, numpy 2.5.0, pyarrow
> 24.0.0). Rebuild bindings any time with `maturin develop` from `rust/flashlfq-py/` with the
> venv active; re-run `rust/flashlfq-py/smoke_test.py` to re-verify the FFI paths. Pinned deps:
> pyo3 0.23.5 (`default-features = false`, bindings list `macros`+`extension-module`+`abi3-py39`),
> numpy 0.23, **arrow 54 → 54.3.1** (last line pinning `pyo3 ^0.23`; 55+ needs pyo3 0.24+),
> **mzdata 0.65** (`default-features = false, features = ["mzml", "miniz_oxide"]` — pure-Rust
> zlib, no C/CMake). `Cargo.lock` and `rust/.venv/` are gitignored. For Phase 1, work happens in
> `flashlfq-core` (pure Rust) — no Python needed until a binding task.

## Conventions

- Repo layout for the port: `rust/` workspace with `flashlfq-core/` (pure Rust, no Python) and
  `flashlfq-py/` (PyO3 bindings only). Parity harness + unit tests live in `flashlfq-core/`.
- Golden test data already exists: `mzLib/Test/FlashLFQ/TestData/` (`sliced-mzml.mzML`,
  `AllPSMs.psmtsv`, `SmallCalibratibleYeast.mzml`; MBR set: `PSMsForMbrTest.psmtsv`,
  `f1r1_sliced_mbr.raw`, `f1r2_sliced_mbr.raw`, `20100614_Velos1_TaGe_SA_K562_3/4.mzML`).
- Parity tolerance: integers/counts/scan-indices/array-lengths **exact**; floats **relative**
  `|a-b| / max(|a|,|b|) < 1e-6`. Diff tool prints worst offender per stage.

## Phase 0 — Toolchain spike

Validate the whole binding stack cheaply before any real logic. Ref: *Python binding layer*
and *Phase 0* in the feasibility doc.

- [x] **P0.1 — Scaffold workspace.** Create `rust/` with `flashlfq-core` (lib) and
  `flashlfq-py` (PyO3 + maturin) crates. _Done when:_ `cargo build` succeeds for both.
  _Result:_ `rust/` workspace builds clean (pyo3 0.23.5, abi3-py39). `flashlfq-core` exposes a
  placeholder `version()` (test passes); `flashlfq-py` exposes `core_version()` via `#[pymodule]`.
  maturin `pyproject.toml` added for P0.5.
- [x] **P0.2 — Read mzML via `mzdata`.** In core, open `sliced-mzml.mzML` and count MS1 scans.
  _Done when:_ a core unit test asserts the MS1 scan count.
  _Result:_ `flashlfq-core::count_ms1_scans` opens the (zlib-compressed) golden mzML via
  `mzdata 0.65` and counts spectra with `ms_level() == 1`. Test `counts_ms1_scans_in_sliced_mzml`
  asserts **10 MS1** (parity number: file has 100 spectra = 10 MS1 + 90 MS2, from the `ms level`
  cvParams; matches the C# reader). `cargo test` green (2 passed); full workspace builds.
  mzdata uses `default-features = false, features = ["mzml", "miniz_oxide"]` to stay pure-Rust.
- [x] **P0.3 — NumPy return.** Expose a function returning a real `f64` array to Python via
  `into_pyarray` (the move/no-copy path, not `to_pyarray`). _Done when:_ Python receives a
  `numpy.ndarray` of the expected length.
  _Result:_ `flashlfq-core::iota_f64(n)` produces an owned `Vec<f64>` ramp `[0.0..(n-1)]`
  (pure Rust, unit-tested: 3 core tests green). `flashlfq-py::iota_array(n)` moves it into a
  `numpy.ndarray` via `IntoPyArray::into_pyarray` (no copy) and is registered in the
  `#[pymodule]`. Added `numpy = "0.23"` to `flashlfq-py` (pairs with pyo3 0.23). Full
  workspace `cargo build` clean; both crates compile. Python-side import deferred to P0.5
  (maturin not yet installed) — Rust binding verified to build here.
- [x] **P0.4 — Arrow round-trip.** Hand one Arrow table Rust→pyarrow over the C Data Interface.
  Pin a known-compatible `arrow-rs` / `pyarrow` version pair and record it here. _Done when:_
  pyarrow reads the table and versions are noted.
  _Result:_ `flashlfq-core::demo_record_batch()` builds a 2-col × 3-row Arrow `RecordBatch`
  `(mz: f64, intensity: f64)` in pure Rust (unit test `demo_record_batch_has_expected_shape_and_values`
  asserts shape, column names, and values — **4 core tests green**). `flashlfq-py::demo_table()`
  hands it to pyarrow via `arrow::pyarrow::IntoPyArrow::into_pyarrow` (zero-copy C Data Interface)
  and is registered in the `#[pymodule]`; full workspace `cargo build` clean. **Version pairing:**
  `arrow = "54"` in `[workspace.dependencies]` (resolves to **54.3.1**) — arrow 54.x is the last
  release line whose `pyarrow` feature pins **`pyo3 = "^0.23"`** (matches our **pyo3 0.23.5**);
  arrow 55+ moves to pyo3 0.24+ and would clash. The handoff rides the Arrow C Data Interface
  (ABI-stable), so any recent **pyarrow (>= 14)** reads the table. End-to-end pyarrow import is
  verified in P0.5 (maturin not yet installed) — Rust binding verified to build here.
- [x] **P0.5 — Wheel smoke test.** `maturin develop` builds and the module imports cleanly;
  a tiny Python script exercises P0.2–P0.4. _Done when:_ the script runs green.
  _Result:_ Built into a **Python 3.13 venv** at `rust/.venv` (`py -V:3.13`; chosen over the
  system 3.14 because pyarrow ships cp313 wheels but no 3.14 yet — and the abi3-py39 wheel
  loads on any ≥3.9 anyway). Installed `maturin 1.14.1`, `numpy 2.5.0`, `pyarrow 24.0.0`.
  `maturin develop` (run from `flashlfq-py/` with the venv on PATH) compiled arrow 54.3.1 +
  pyo3 0.23.5 + numpy 0.23 and produced `flashlfq-0.1.0-cp39-abi3-win_amd64.whl`, installed
  editable. `rust/flashlfq-py/smoke_test.py` exercises all three FFI paths and prints
  **ALL SMOKE CHECKS PASSED**: `core_version()` → `'0.1.0'` (str); `iota_array(5)` →
  `numpy.ndarray` f64 len-5 == `[0..4]`; `demo_table()` → `pyarrow.RecordBatch` 2 cols
  (`mz`, `intensity`, both f64) × 3 rows. **Phase 0 complete — the whole binding stack is
  validated.** `.venv/`, `*.egg-info/`, `__pycache__/` added to `rust/.gitignore`.
  _Repro:_ `py -V:3.13 -m venv rust/.venv`; `rust/.venv/Scripts/python -m pip install maturin
  numpy pyarrow`; from `rust/flashlfq-py/` with the venv active run `maturin develop`; then
  `rust/.venv/Scripts/python rust/flashlfq-py/smoke_test.py`.

## Phase 1 — MVP core quant (mzML only)

The core ask. Isotope chemistry is ported **first** and gated by parity (L0–L3) before the
tracing core (L4–L6). Ref: *Isotope-distribution port plan* and *layered parity harness*.

### 1a. Isotope chemistry port + parity gates L0–L3

- [x] **P1.1 — Periodic-table data (L0).** Dump mzLib's loaded isotope table
  `(element → [(isotopeNumber, atomicMass, relativeAbundance)])` and embed it **verbatim** in
  Rust (do not re-source from NIST). _Done when:_ a Rust `validate_abundances` test (per-element
  sum ≈ 1 within epsilon) passes **and** an L0 diff vs the C# dump is exact.
  _Result:_ **8 core tests green** (4 new). The `.NET SDK is not available` in this env, so the
  C# dumper called for by the plan is substituted by a faithful Python replica of the
  deterministic `Chemistry/PeriodicTable.cs` `static PeriodicTable()` loader
  (`rust/flashlfq-core/parity/dump_periodic_table.py`) — run with the P0.5 venv:
  `rust/.venv/Scripts/python rust/flashlfq-core/parity/dump_periodic_table.py`. It extracts the
  **verbatim** `thePeriodicTable` literal from PeriodicTable.cs into
  `rust/flashlfq-core/data/periodic_table.nist.txt` (single source of truth, embedded by Rust via
  `include_str!`) and emits the golden artifact `parity/golden/L0_periodic_table.tsv`
  (floats as Python `repr()` = shortest round-trip). Rust module `flashlfq_core::periodic_table`
  independently re-implements the same loader over the same embedded data, exposing
  `periodic_table()`, `element_by_symbol/number`, `Element::principal_isotope`, and
  `validate_abundances(epsilon)`. **Loaded shape:** 84 naturally-occurring elements, 288 isotopes
  (only isotopes with a measured abundance are kept; e.g. H folds in deuterium and drops tritium;
  Tc/Pm/transuranics with no abundance are dropped entirely). **L0 gate:** test
  `l0_diff_vs_golden_is_exact` diffs the Rust-loaded table against the golden tsv **exact** (ints
  exact; floats compared with `==` — bit-identical because both sides parse the same decimal
  tokens via correctly-rounded f64 parsing). `validate_abundances(1e-15)` passes and
  `validate_abundances(0.0)` is false, mirroring the C# `ValidatePeriodicTable` test (worst
  per-element sum error is Si at 1.11e-16). Derived avg masses for `[a,b]`-range weights computed
  as `(a+b)/2` identically on both sides. **Caveat for L1:** if the .NET SDK becomes available,
  re-run the real C# `PeriodicTable` loader and diff against this golden to retire the Python-replica
  substitution.
- [x] **P1.2 — `ChemicalFormula`.** Port element→count map + `MonoisotopicMass`/`AverageMass`.
  _Done when:_ core unit tests for a few known formulas pass.
  _Result:_ New module `rust/flashlfq-core/src/chemical_formula.rs` (**14 core tests green**, 6 new).
  Faithful port of `Chemistry/ChemicalFormula.cs`: two stores — `elements` (`atomic_number ->
  count`, BTreeMap) and `isotopes` (`(atomic_number, mass_number) -> count`) — with the C# removal
  thresholds (`== 0` for elements, `<= 0` for isotopes). `monoisotopic_mass()` sums element
  principal-isotope masses + isotope exact masses; `average_mass()` sums element average masses +
  isotope exact masses; plus `atom_count()`, `add_element`/`add_isotope`/`add_formula`. Masses sum
  in deterministic sorted order (L1 tol is rel-1e-9, so order vs C# Dictionary enumeration is
  immaterial). `parse_formula` is a hand-written tokenizer (no `regex` dep) faithful to the C#
  grammar `\s*([A-Z][a-z]*)(?:\{([0-9]+)\})?(-)?([0-9]+)?\s*` repeated over the whole string:
  symbol, optional `{n}` isotope, optional `-` sign, optional count (default 1) — e.g. `"C-2C{13}5"`
  subtracts 2 C then adds 5 ¹³C. Errors via `FormulaParseError` {BadFormat, UnknownElement,
  UnknownIsotope}. Added `Element::isotope_by_mass_number` to `periodic_table` (mirrors the C#
  `Element this[int]` indexer). Tests assert known literature masses (counts exact; H2O mono
  18.0105646863, glucose C6H12O6 180.0633881, ¹³C₆ glucose 186.0835171 — all rel-1e-9) plus
  sign/multi-letter-symbol parsing, formula combination, and malformed/unknown rejection.
- [x] **P1.3 — Peptide → formula.** Port the residue→formula table + modification formulas
  backing `Peptide.GetChemicalFormula`. _Done when:_ formula for several sequences matches C#.
  _Result:_ New module `rust/flashlfq-core/src/peptide.rs` (**20 core tests green**, 6 new).
  Ports the formula side of `Proteomics/AminoAcidPolymer/`: `residue_formula(char)` is the
  **verbatim** 22-entry residue→formula table from `Residue.cs` `ResiduesDictionary` (20 standard
  + `O` pyrrolysine `C12H19N3O2` + `U` selenocysteine `C3H5NOSe`; letters mzLib leaves null —
  B/J/X/Z — return `None`). `peptide_base_formula(seq)` mirrors
  `new Peptide(seq).GetChemicalFormula()`: N-terminus `H` + C-terminus `OH` (= one water) + every
  residue formula, built by summing `ChemicalFormula::parse_formula` results via `add_formula`.
  `peptide_formula_with_mods(seq, &[ChemicalFormula])` mirrors `FullChemicalFormula` (backbone then
  `Add` each mod formula); resolving a *modified-sequence string* to its mod formulas is deferred to
  the P1.4 harness. Added `ChemicalFormula::count_of_element`/`count_of_isotope` accessors so parity
  asserts **exact** element counts. Tests verify exact element counts + known monoisotopic masses
  (rel-1e-7): free Gly `C2H5NO2` 75.03203, Ala `C3H7NO2` 89.04768, canonical **PEPTIDE**
  `C34H53N7O15` 799.35996, sulfur-bearing `MCK`, Se-bearing `U`, PEPTIDE+phospho (`HO3P`, Δ
  +79.96633), and rejection of unknown residues (`B`, `J`). _Caveat (same as L0/L1):_ ground truth
  is literature masses, not a live C# dump (no .NET SDK in this env); reconcile against the real
  `Peptide.GetChemicalFormula` if the SDK becomes available.
- [x] **P1.4 — L1 parity gate.** Harness: modified sequence → element counts + mono mass.
  _Done when:_ counts exact, mass rel-1e-9, over all peptides in `AllPSMs.psmtsv`.
  _Result:_ **21 core tests green** (20 unit + 1 new integration test). The plan's hard part —
  resolving a modified-sequence string to its mod formulas — turned out to be **already done in
  the data**: `AllPSMs.psmtsv` carries a `Mods Combined Chemical Formula` column (col 19) holding
  the summed modification formula mzLib's `PeptideWithSetModifications.FullChemicalFormula` adds
  on top of the backbone. Only 4 distinct non-empty mod formulas appear (`C2H3NO` carbamidomethyl,
  `O` oxidation, `C2H3NO2`, `C4H6N2O2`), all plain (no isotope notation). So the L1 full formula is
  `peptide_base_formula(Base Sequence)` + `parse(Mods Combined Chemical Formula)` — exactly
  `peptide::peptide_formula_with_mods`. New golden generator
  `rust/flashlfq-core/parity/dump_l1_peptide_formulas.py` (Python replica reusing the L0
  `dump_periodic_table` loader as the single element-data source) emits
  `parity/golden/L1_peptide_formulas.tsv`: **354 distinct peptides**, each with exact element
  counts (symbol:count, atomic-number order), the replica mono mass (`repr()` round-trip), and the
  C# `Peptide Monoisotopic Mass` column verbatim. Rust gate
  `rust/flashlfq-core/tests/l1_parity.rs` drives the distinct peptide set **straight from the live
  psmtsv** (so the gate tracks the corpus file, not a frozen copy), recomputes each formula via
  `peptide_formula_with_mods` over `chemical_formula` + L0 `periodic_table`, and asserts three
  things: **element counts exact** (every golden element matched + `atom_count` equals the sum, so
  no extra element can hide); **mono mass rel-1e-9 vs the Python replica** (worst 2.85e-16 —
  bit-identical); **mono mass abs < 1e-5 vs the C# psmtsv column** (worst 4.98e-6 Da = the
  5-decimal print-rounding bound, a genuine tie to real C# output despite no .NET SDK). _Caveat
  (same as L0/L1 parts):_ the per-formula ground truth on element counts is the replica, not a live
  C# `FullChemicalFormula` dump; the mono-mass column keeps the C# anchor. Reconcile against the
  real C# if the SDK becomes available.
- [x] **P1.5 — `IsotopicDistribution.GetDistribution`.** Port the Kubinyi polynomial algorithm
  (`MergeFinePolynomial`, `MultiplyFinePolynomial`, `MultiplyFineFinalPolynomial`,
  `MultipleFinePolynomialRecursiveHelper`, `FactorLn`, `CalculateFineGrain`, `Polynomial`,
  `Composition`), preserving constants (`0.125`, `1e-8`, mol-weight resolution) **and merge
  order**. _Done when:_ it builds and produces `(masses[], intensities[])`.
  _Result:_ New module `rust/flashlfq-core/src/isotopic_distribution.rs` (**25 core unit tests
  green**, 5 new; L1 integration test still passes). `IsotopicDistribution { masses, intensities }`
  with `get_distribution`/`get_distribution_with`. Verbatim port of all five algorithm functions +
  `FactorLn` (memoised thread-local cumulative log-sum, bit-identical summation order to the C#
  `factorLnArray`) and the `Polynomial`/`Composition` structs (using `f64::NAN` sentinels). Default
  constants exact (`fine 0.01`, `min_prob 1e-200`, `mw_res 1e-12`, `fineResolution/2` split). The
  `MergeFinePolynomial` 9-pass merge reads the **running merged centroid** of bin `i` (updated in
  place) for its threshold comparison, matching C#. Added `ChemicalFormula::elements()`/`isotopes()`
  iterators. Tests: water envelope (ascending, base near 18.0106), C100 M/M+1 ~1.0034 apart with
  M+1 dominant, PEPTIDE base within 2.5 Da of mono, ¹³C-label lightest-peak shift, `factor_ln` ==
  direct log-sum. **Known limit (matches C#):** pure-isotope-only formulas (no plain elements, e.g.
  `"C{13}5"`) panic in `MultiplyFinePolynomial`, exactly as C# throws on `tPolynomial =
  fPolynomial[0]`; not hit by real peptide formulas. L2 parity gate is **P1.6**.
- [x] **P1.6 — L2 parity gate.** Harness: formula → `(masses[], intensities[])`. _Done when:_
  intensities rel-1e-6, masses abs-tol, across the corpus.
  _Result:_ **27 core tests green** (25 unit + L1 + new L2 integration test). Golden generator
  `rust/flashlfq-core/parity/dump_l2_isotopic_distributions.py` is a faithful Python replica of
  the full `IsotopicDistribution.cs` Kubinyi algorithm (ports `FactorLn`, `MultiplyFinePolynomial`,
  `MultipleFinePolynomialRecursiveHelper`, `MultiplyFineFinalPolynomial`, `MergeFinePolynomial`,
  `CalculateFineGrain`), reusing the L0 element loader + L1 `full_formula_counts` (backbone + the
  psmtsv `Mods Combined Chemical Formula` column). It emits
  `parity/golden/L2_isotopic_distributions.tsv` (D/M/I records: descriptor + masses + intensities,
  repr() round-trip) for the **354 distinct peptides**; max envelope is 1095 peaks. **The
  element-order caveat is resolved by aligning the replica to *our* order:** the dumper iterates
  elements `for an in sorted(counts)` (atomic-number order) to match Rust's `ChemicalFormula::
  elements()` BTreeMap, making the convolution index-aligned. Rust gate
  `rust/flashlfq-core/tests/l2_parity.rs` drives the distinct set straight from the live psmtsv
  (asserts the set matches golden), recomputes each envelope via `IsotopicDistribution::
  get_distribution`, and asserts: **array lengths exact**, **masses abs < 1e-6 Da**, **intensities
  rel-1e-6** — observed **worst mass abs 1.36e-12 Da, worst intensity rel 1.86e-15** across 257,162
  total peaks, i.e. near-bit-identical. Regenerate:
  `& "F:\flashlfq-rust\rust\.venv\Scripts\python.exe" rust/flashlfq-core/parity/dump_l2_isotopic_distributions.py`.
  _Caveat (same as L0/L1):_ ground truth is the Python replica, not a live C# dump (no .NET SDK);
  reconcile against real C# `IsotopicDistribution.GetDistribution` if the SDK ever appears.
- [x] **P1.7 — `CalculateTheoreticalIsotopeDistributions`.** Port averagine fill-in (the five
  `average{C,H,O,N,S}`, `massDiff > 20`, `Math.Round` counts), normalize-to-most-abundant, the
  truncation rule (`count < NumIsotopesRequired || abundance > 0.1`), and
  `PeakfindingMass = MonoisotopicMass + mostAbundantIsotopeShift`. _Done when:_ it builds.
  _Result:_ New module `rust/flashlfq-core/src/theoretical_isotope_distribution.rs` (**34 core
  unit tests green**, 8 new; L1 + L2 integration still green). `expected_isotope_peaks(optional_formula,
  base_sequence, mono, num_isotopes_required) -> Vec<ExpectedIsotopePeak { mass_shift,
  normalized_abundance }>` is a faithful port of the per-id body: resolve formula
  (explicit → use as-is; valid sequence → bare peptide formula + averagine top-up when
  `|mono - formulaMono| > 20`; invalid sequence → pure averagine sized by `mono`), compute the
  fine distribution at `(0.125, 1e-8)`, re-center on `mono` (two-step `+= (mono-formulaMono)` then
  `-= mono`, kept separate for bit-parity), normalize to the most abundant peak, truncate by
  `len < NumIsotopesRequired || abundance > 0.1`. Plus `averagine_mass()` (≈111.1236 Da from L0
  avg masses), `all_sequence_residues_are_valid()`, `most_abundant_isotope_shift()`,
  `peakfinding_mass()`. Constants verbatim (`AVERAGE_{C,H,O,N,S}`, `FINE_RESOLUTION=0.125`,
  `MIN_PROBABILITY=1e-8`, `DEFAULT_NUM_ISOTOPES_REQUIRED=2`). C# `Math.Round(x,0)` banker's
  rounding → Rust `f64::round_ties_even`. The L3 parity gate is **P1.8**.
- [x] **P1.8 — L3 parity gate (chemistry gate).** Harness: modified sequence → final
  `(massShift, normAbundance)[]` **and** `PeakfindingMass`. _Done when:_ massShift abs-1e-6,
  abundance rel-1e-6, **array lengths exact**, across the corpus. ← all later work depends on
  this passing.
  _Result:_ **Green — chemistry gate complete.** Golden dumper
  `parity/dump_l3_theoretical_distributions.py` (replica of
  `CalculateTheoreticalIsotopeDistributions` over L0+L1+L2 replicas, at `0.125 / 1e-8`) →
  `parity/golden/L3_theoretical_distributions.tsv` (354 distinct Full Sequences, max 6 peaks).
  Formula built from **base sequence only**, averagine top-up when `|mono−baseFormulaMono|>20`
  (fired on 44/354); `mono` is the psmtsv `Peptide Monoisotopic Mass` column. Rust gate
  `tests/l3_parity.rs` drives the live psmtsv, asserts **lengths exact / massShift abs-1e-6 /
  abundance rel-1e-6 / peakfinding mass abs-1e-6** — worst **massShift 9.09e-13 Da, abundance
  rel 5.71e-16, peakfinding 4.55e-13 Da** over 1223 peaks. 34 unit + L1 + L2 + L3 all green.

### 1b. psmtsv + index + tracing + parity gates L4–L6

- [x] **P1.9 — psmtsv parse.** Port `PsmFromTsv` → `Identification` (BaseSequence,
  ModifiedSequence, MonoisotopicMass, charge, RT, file) with `csv` + `serde`. _Done when:_
  parsing `AllPSMs.psmtsv` yields the same record count/fields as C#.
  _Result:_ New module `rust/flashlfq-core/src/psm_tsv.rs` (**41 core tests green**, 7 new; L1+L2+L3
  integration still green). `Identification { file_name, base_sequence, modified_sequence,
  monoisotopic_mass, ms2_retention_time_in_minutes, precursor_charge_state, score, q_value,
  is_decoy }` (`#[derive(Serialize)]` for later output) + `read_identifications(path)` /
  `read_identifications_from_str(text)` / `remove_parentheses` / `distinct_modified_sequences`.
  Faithful port of the slice of mzLib `SpectrumMatchTsvReader`/`SpectrumMatchFromTsv`/`PsmFromTsv`
  + `MzLibExtensions.MakeIdentifications` that FlashLFQ needs. **Tabular model:** drives the `csv`
  crate with `quoting(false)`, `delimiter(b'\t')`, `flexible(true)` to mirror mzLib's plain
  `line.Split('\t')` + its own `Trim('"').Trim()` (NOT RFC-4180 quoting — cells hold `[y1+1:..]`
  text). **Header→column resolution** mirrors `ParseHeader`: an exact `Accession` column selects
  the `Monoisotopic Mass` layout, else (the FlashLFQ psmtsv layout, confirmed: no `Accession`
  column) `Peptide Monoisotopic Mass` (col 22). **Per-field semantics match the C# ctor:** File
  Name → strip accepted spectra extensions + bare file name; Base Sequence → `RemoveParentheses`
  (SILAC `(...)`); Full Sequence verbatim; mono = first `|`-token as f64 else `-1.0`; Precursor
  Charge = `(int)f64` (truncate); Scan Retention Time → minutes default `-1.0`; Score/QValue f64;
  decoy = `Decoy/Contaminant/Target` contains `'D'`. **Corpus parity:** the integration-style unit
  test `parses_all_psms_corpus_record_count_and_distinct_sequences` parses the **live**
  `AllPSMs.psmtsv` → **594 PSM rows / 354 distinct Full Sequences, all targets**, and spot-checks
  row 0 (file `20100614_Velos1_TaGe_SA_K562_3`, base `AHQLVMEGYNWCHDR`, charge 3, mono 1914.82537,
  RT 79.12191, score 23.287). _Divergence from C# (intentional):_ a required-but-missing column or
  unparseable required number is a **hard error** (`PsmTsvError`); mzLib drops the row with a
  warning. For a parity gate we want a hard error over silently fewer records. Added deps
  `csv = "1.3"`, `serde = { version = "1", features = ["derive"] }` to `flashlfq-core/Cargo.toml`.
  _Caveat:_ protein-group resolution (the cross-record `HashSet<ProteinGroup>` dedup) is deferred
  — it is a property of the whole id set, not one row, and only feeds protein quant.
- [x] **P1.10 — `IndexingEngine<T>`.** Port the binned m/z index: `IndexPeaks`,
  `GetIndexedPeak`, `BinarySearchForIndexedPeak`, `BinsPerDalton`. _Done when:_ core tests
  cover point queries.
  _Result:_ New module `rust/flashlfq-core/src/peak_indexing.rs` ports
  `IndexingEngine<IndexedMassSpectralPeak>` / `PeakIndexingEngine` (m/z, charge-less flavour):
  `IndexedMassSpectralPeak` (mz/intensity/RT as **f32** per the C# `(float)` casts; `m()==mz`),
  `ScanInfo`, an input `Scan` struct, and `PeakIndexingEngine` (jagged `Vec<Option<Vec<peak>>>`
  index + `Vec<ScanInfo>`). Faithful ports: `index_peaks` (index length `ceil(maxLastX*100)+1`,
  bin = `(mz*100).round_ties_even()` — C# `Math.Round` banker's; peaks appended scan-ascending;
  `None` when nothing to index), `get_indexed_peak`, `get_bins_in_range`,
  `get_best_peak_from_bins`, `get_peak_from_bin`, `binary_search_for_indexed_peak` (binary search +
  linear backstep via `isize`), `BINS_PER_DALTON=100`. New `rust/flashlfq-core/src/tolerance.rs` =
  `PpmTolerance` (faithful to `MzLibUtil/PpmTolerance.cs`). mzML loader `read_ms1_scans` +
  `PeakIndexingEngine::from_mzml` via `mzdata` (`raw_arrays()`, `start_time()`, `index()+1`).
  Charge + mass-indexing engine **not** ported (m/z peaks only). **53 core lib tests + L1+L2+L3
  integration green** (10 new peak_indexing tests incl. a real `sliced-mzml.mzML` round-trip + 2
  tolerance tests).
- [x] **P1.11 — XIC extraction.** Port `GetXic`, `GetXicByScanIndex`, `GetAllXics`. _Done
  when:_ XIC for a known m/z matches C# scan-by-scan.
  _Result:_ Extended `rust/flashlfq-core/src/peak_indexing.rs`. `get_xic` (finds start scan by RT —
  first scan whose RT is `>= retention_time`, leaving start `-1` if none precede it, exactly as C#)
  delegates to `get_xic_by_scan_index`, a faithful port of `GetXicByScanIndex` (m/z, charge-less
  path): seeds a per-bin pointer array via `binary_search_for_indexed_peak`, grabs the initial peak,
  then walks backward (`-1`) and forward (`+1`) advancing each bin pointer to the first peak of the
  new scan index (the C# `do/while` decrement-then-`++` for `-1`, the `while` increment for `+1`),
  stopping after `missed_scans_allowed` consecutive misses, out-of-range, or RT past
  `max_peak_half_width` from the initial peak. A peak in `matched_peaks` counts as a miss. Result
  sorted by RT ascending. `get_all_xics` ports `GetAllXics` **minus the optional `CutPeak`** (deferred
  to P1.14): flattens the index in bin-then-scan order, stable-sorts by intensity descending (= C#
  `OrderByDescending`), greedily traces each unclaimed peak, keeps XICs of `>= num_peak_threshold`
  peaks, marks their peaks claimed. New partial port of `ExtractedIonChromatogram` (the container +
  `SetXicInfo`/`AverageM` summary: `apex_peak`, `start_rt`/`end_rt`, `start_scan_index`/
  `end_scan_index`, intensity-weighted `averaged_mass_or_mz`); `CutPeak`/`FindPeakBoundaries`/
  `RemovePeaks` deferred to P1.14. **Peak identity:** C# threads a `Dictionary<IIndexedPeak,…>` keyed
  by object reference; our peaks are `Copy`, so dedup uses a structural `PeakKey = (scan_index,
  mz.to_bits(), intensity.to_bits())` in a `HashSet`. **Tests (5 new, 58 core + L1+L2+L3 green):**
  `get_xic_by_scan_index_traces_all_scans` (= C# `TestGetXicWithStartIndex`, start 0 → **9**),
  `get_xic_by_retention_time_traces_all_scans` (= `TestGetXicWithRetentionTime`, RT 1.0 → **9**),
  `get_xic_returns_empty_for_absent_mz` (= `TestMissingXic`), `get_xic_stops_after_missed_scans`
  (= `TestXicStops`, 10-scan fixture with a charge-shift gap at scans 2/3, start 7, 1 miss allowed
  → **6**, scans [4,5,6,7,8,9]), `get_all_xics_traces_every_mz_trace` (synthetic 5-trace × 9-scan
  fixture at 20 ppm → **5 XICs × 9 peaks = 45**, each apex at scan 4). _Caveat:_ ground truth is the
  C# `IndexingEngineTests` expected counts replicated as Rust fixtures (no live C# dump; no .NET SDK).
  `GetAllXics` parity is structural (no C# count fixture exists for it); the `CutPeak` trimming it
  can apply lands with P1.14.
- [x] **P1.12 — L4 reader calibration (report-only).** Diff `mzdata` vs mzLib per-scan
  `(m/z, intensity)` peak lists for `sliced-mzml.mzML`. _Done when:_ the magnitude of any
  reader delta is characterized and recorded here (not a pass/fail gate yet).
  _Result:_ **Bit-identical — zero reader drift.** Golden dumper
  `rust/flashlfq-core/parity/dump_l4_reader_peaks.py` is an independent raw decode of the mzML
  binary data arrays (base64 → zlib-decompress → reinterpret as 32-bit little-endian floats per
  the `MS:1000521`/`MS:1000574` cvParams), emitting `parity/golden/L4_reader_peaks.tsv`
  (`S` scan-summary + `P` per-peak records, `repr()` round-trip) for the **10 MS1 scans**. This is
  exactly what mzLib's mzML reader does (no re-centroiding / sorting / zero-trimming), so the
  golden values are the file's ground-truth peaks = what any faithful reader (mzLib included)
  yields. Report-only test `rust/flashlfq-core/tests/l4_reader_calibration.rs` reads MS1 scans via
  `read_ms1_scans` (mzdata) and diffs against golden, printing per-scan + summary magnitudes.
  **Observed across all 67,001 MS1 peaks / 10 scans: worst |Δm/z| = 0, worst rel |Δintensity| = 0,
  worst |Δintensity| = 0, worst |Δrt| = 0 — exactly bit-for-bit.** Both sides decode the same
  32-bit-float zlib payload and widen f32→f64 identically (an exact, lossless conversion). Note the
  file is **profile** data (8156/6990/… peaks per scan, including zero-intensity points), not
  centroided — reader fidelity is total regardless. **Implication for L5/L6:** no reader-level slop
  needs to be budgeted into the tracing parity tolerances for `sliced-mzml.mzML`; any L5/L6
  divergence will be algorithmic, not from the reader. The test asserts peak-count parity exactly
  (a structural invariant) + loose float sanity bounds (`<1e-2` Da m/z, `<1e-3` rel int) so a
  catastrophic reader regression still fails; the recorded magnitude is the deliverable. **58 core
  unit tests + L1 + L2 + L3 + L4 all green.** Regenerate golden:
  `& "F:\flashlfq-rust\rust\.venv\Scripts\python.exe" rust/flashlfq-core/parity/dump_l4_reader_peaks.py`.
  _Caveat:_ characterized for `sliced-mzml.mzML` only (the Phase-1 mzML); `.raw` reader fidelity is
  a Phase-2 concern (P2.1).
- [x] **P1.13 — `GetIsotopicEnvelopes` + `CheckIsotopicEnvelopeCorrelation`.** Port the
  envelope search incl. the ±1 mass-shift off-by-one logic and the intensity ratio gate
  (`/4`..`*4`). _Done when:_ it builds and returns envelopes.
  _Result:_ **Green — builds and returns envelopes.** New module
  `rust/flashlfq-core/src/isotopic_envelope.rs`: `IsotopicEnvelope` struct + free fn
  `get_isotopic_envelopes(engine, xic, isotope_mass_shifts, mono, peakfinding_mass, charge,
  num_isotopes_required, isotope_ppm_tol)` faithfully porting `FlashLfqEngine.GetIsotopicEnvelopes`
  and `check_isotopic_envelope_correlation` porting `CheckIsotopicEnvelopeCorrelation`. The
  `{-1,0,+1}` off-by-one hypotheses are a `[Vec<IsotopePeakRecord>; 3]`; `peakfindingMassIndex =
  round_ties_even(peakfinding−mono)`; down-then-up peak walk breaking on missing peak or the
  `/4`..`*4` intensity-ratio gate; correlation gate `pearson>0.7 && corrLeft−corrPadded<0.1 &&
  corrRight−corrPadded<0.1` (NaN→−1) with the "unexpected" padding peak; zero-intensity imputation +
  summed/charge intensity. `peak.M.ToMass` in f32-then-widened (C# float overload), `mass.ToMz`
  all-f64; `pearson(&[f64],&[f64])` = MathNet `Correlation.Pearson` Welford online algorithm
  (empty/constant→NaN). Constants `PROTON_MASS`, `C13_MINUS_C12` (mzLib `Chemistry/Constants.cs`).
  Vestigial C# `isotopologuePeaks` omitted. **7 new tests (65 core unit + L1+L2+L3+L4 green)**:
  pearson identical/anti/constant, mz↔mass round-trip, clean 3-isotope synthetic envelope across a
  7-scan profile (1 env/scan, all pearson>0.7, correct apex + summed intensity), fewer-isotopes→
  empty, charge 2. _Caveat:_ synthetic-fixture ground truth (no .NET SDK); corpus parity folded into
  L5 at P1.17.
- [x] **P1.14 — `CutPeak`.** Port apex + boundary integration. _Done when:_ it builds.
  _Result:_ New module `rust/flashlfq-core/src/chromatographic_peak.rs`. Introduces
  `ChromatographicPeak` (the `Vec<IsotopicEnvelope>` container) + `calculate_intensity_for_this_feature`
  (apex/sum intensity, mass error, distinct-charge count) + `cut_peak` (the recursive valley
  split, `DISCRIMINATION_FACTOR_TO_CUT_PEAK = 0.6`, ≥5-envelope guard, far-side `RemoveAll`,
  `SplitRT`). Made `isotopic_envelope::mz_to_mass_f32` `pub` for the apex mass error. The
  engine-level `FlashLfqEngine.CutPeak` is the active FlashLFQ path; `ExtractedIonChromatogram.CutPeak`
  (iso-tracker path) was not needed. **7 new tests, 72 core unit + L1+L2+L3+L4 green.**
- [x] **P1.15 — Charge aggregation + FDR.** Port charge-state aggregation and FDR.
  _Result:_ **Green — builds, 82 core unit + L1+L2+L3+L4 all pass.** Extended
  `rust/flashlfq-core/src/chromatographic_peak.rs` so a `ChromatographicPeak` now carries its
  `identifications: Vec<Identification>` + `detection_type: DetectionType` +
  `num_identifications_by_base_seq/full_seq` (alongside the existing P1.14 envelope/intensity
  fields, kept untouched so CutPeak is unaffected). New constructors `from_identification` /
  `from_identifications` (run `resolve_identifications` like the C# ctor); the bare
  `new(peakfinding_masses)` is retained for the CutPeak fixtures (empty ids). Ported
  `ChromatographicPeak.ResolveIdentifications` (distinct base/full-seq counts),
  `MergeFeatureWith` (union ids by value-equality + peakfinding masses in lockstep, add
  non-duplicate envelopes via a structural `EnvelopePeakKey` = scan-index + m/z-bits +
  intensity-bits standing in for C# `IIndexedPeak` reference equality, recompute intensity,
  IsoTracker→Ambiguous promotion), `ApexRetentionTime`, `DecoyPeptide`. New module
  `rust/flashlfq-core/src/detection_type.rs` = verbatim `FlashLFQ.DetectionType` enum. New module
  `rust/flashlfq-core/src/results.rs`: `run_error_checking` ports the **MSMS path** of
  `FlashLfqEngine.RunErrorChecking` (group peaks by apex peak; same-apex MSMS peaks both in the
  quantify set → `merge_feature_with`; only-try in set → overwrite stored in place preserving
  iteration order; neither → drop; apex-less peaks appended first, matching C# `errorCheckedPeaks`
  then grouped-values order via an insertion-ordered key-vec + HashMap standing in for C# Dictionary
  order). `calculate_peptide_results` ports the MSMS path of
  `FlashLfqResults.CalculatePeptideResults`: per file, seed every quantify-set sequence to
  `NotDetected`/0, group non-decoy single-full-seq peaks by modseq, pick max-intensity peak (best =
  first attaining it) → `MSMS` if >0 else `MSMSIdentifiedButNotQuantified`; ambiguous (>1 full-seq)
  peaks zero a sequence + mark `MSMSAmbiguousPeakfinding` when `peak.Intensity/(already+peak)>0.3`
  (the `quantifyAmbiguousPeptides=false` branch). `PeptideResults` = `modseq → file → PeptideQuant
  {intensity, retention_time, detection_type}` — the (peptide × file) table L6 will check.
  `default_quantify_set` = non-decoy modified sequences (= C# `_peptideModifiedSequencesToQuantify`
  fallback). **FDR for Phase 1 = the quantify-set filter** (decoys excluded; only quantify-set
  sequences reported); MSMS peaks carry the search q-value, there is no per-peak q-value computed
  in Phase 1 — the MBR q-value FDR (`CalculateFdrForMbrPeaks`/`CorrectQValues`) is Phase 3.
  **Intentional Phase-1 omissions (noted in the module doc):** the MBR/IsoTracker branches of both
  routines (unreachable when every peak is MSMS) and `HandleAmbiguityInFractions` (fraction-only).
  **10 new tests:** resolve-counts, decoy read, merge unions ids+dedups envelopes+integrates,
  error-checking merges same-apex / keeps distinct-apex / overwrites unquantified-stored, peptide
  results pick-most-intense / undetected-stays-zero / ambiguous-zeroed, default-set-excludes-decoys.
  _Caveat:_ synthetic-fixture ground truth (no .NET SDK); corpus/parity for the aggregated
  peptide×file table is folded into the L6 gate at P1.17.
- [x] **P1.16 — Parquet output.** Write the results table (peptide × file → intensity) to
  Parquet via `arrow`/`parquet`. _Done when:_ output loads in Python (polars/pandas).
  _Result:_ **Green — loads in Python via pyarrow.** New module
  `rust/flashlfq-core/src/parquet_output.rs` serializes `results::PeptideResults` to Parquet in
  **tidy / long form**: one row per `(modified_sequence, file_name)` cell with columns
  `modified_sequence:Utf8`, `file_name:Utf8`, `intensity:Float64`, `retention_time:Float64`,
  `detection_type:Utf8` (all non-nullable). `peptide_results_schema()` / `peptide_results_record_batch()`
  build an Arrow `RecordBatch` (rows **sorted by (modseq, file)** for reproducible output —
  `HashMap` order is otherwise nondeterministic); `write_peptide_results_parquet(results, path)`
  writes it via `parquet::arrow::ArrowWriter` (default props) and returns the row count. Added
  `DetectionType::as_str()` = the **verbatim C# enum member names** (incl. underscores:
  `IsoTrack_MSMS`/`IsoTrack_MBR`/`IsoTrack_Ambiguous`) used for the `detection_type` column. Long
  form chosen over a wide peptide×file matrix because it's the simplest faithful serialization of
  the nested `modseq→file→quant` map and pivots back to wide in polars/pandas. Added dep
  `parquet = "54"` to `flashlfq-core/Cargo.toml` (pinned to the **same 54.x line as `arrow`** —
  the two crates version in lockstep over shared `arrow-*` sub-crates; keeps us off pyo3 0.24).
  **2 new tests (84 core unit + L1+L2+L3+L4 green):** `record_batch_is_sorted_and_complete`
  (sort order + shape), `parquet_round_trips_through_the_arrow_reader` (write to a temp file,
  read back via `ParquetRecordBatchReaderBuilder`, assert every cell intensity/RT/detection
  survived). **Python acceptance evidence:** example
  `rust/flashlfq-core/examples/write_demo_parquet.rs` writes a 3-row demo file; reading it with
  the `rust/.venv` **pyarrow** (`pyarrow.parquet.read_table`) returns the exact schema + values
  (`MSMS`/`MSMSAmbiguousPeakfinding`/`NotDetected`, intensities `1e6`/`2.5e5`/`0.0`). Repro:
  `cargo run -p flashlfq-core --example write_demo_parquet -- <out.parquet>` then
  `& "F:\flashlfq-rust\rust\.venv\Scripts\python.exe" -c "import pyarrow.parquet as pq; print(pq.read_table(r'<out.parquet>').to_pydict())"`.
  _Note:_ the `Serialize` derive on `Identification`/`PeptideQuant` is unrelated — Parquet goes
  through Arrow, not serde.
- [x] **P1.17 — L5 + L6 parity gate (Phase 1 exit).** L5: per-(peptide,charge,file) integrated
  intensity + pearson corr + apex scan. L6: end-to-end peptide × file intensity. _Done when:_
  both rel-1e-6 across the corpus. ← **Phase 1 complete when this passes.**
  _Result (P1.17b):_ **GREEN — Phase 1 complete.** Rust driver `src/engine.rs` (`run_msms`) ports the
  MS2 path of `FlashLfqEngine.Run`; gate `tests/l5_l6_parity.rs` diffs against the real C# golden
  (507 L5 peaks + 708 L6 cells) at rel-1e-6, detection-type mix exact (450/223/33/2). Root-cause bug
  found + fixed: mzLib's reader drops `<0.01`-intensity peaks; the Rust reader now matches (see
  **Current focus → P1.17b result** for the full diagnosis). 85 core unit + L1–L6 all green.
  _Progress (P1.17a done):_ **Real C# golden generated** (the no-.NET-SDK caveat is retired —
  `dotnet` is now installed). New C# console generator `parity/csharp_golden/` runs the real
  `FlashLfqEngine` (MS2 path, default params, `MaxThreads=1`) over the corpus (`AllPSMs.psmtsv` +
  the two K562 mzMLs) and wrote `parity/golden/L5_chromatographic_peaks.tsv` (507 peaks) +
  `parity/golden/L6_peptide_intensities.tsv` (708 cells). Schema + driver port plan + sanity
  targets are in **Current focus → P1.17b**. _Remaining (P1.17b):_ port the MS2 path of
  `FlashLfqEngine.Run` into a Rust driver and add the parity test diffing it against the goldens.
  Regenerate golden: `cd rust/flashlfq-core/parity/csharp_golden; dotnet run -c Release`.

### 1c. Python surfaces

- [x] **P1.18 — `quant()` binding.** `quant(results_path, raw_paths, params)` → Parquet path
  (or in-memory table for small sets). _Done when:_ callable from Python end-to-end.
  _Result:_ **Green — callable from Python end-to-end over the real corpus.** New high-level core
  entry `engine::quant(results_path, raw_paths) -> Result<Ms2QuantResult, QuantError>`: parses the
  psmtsv (`read_identifications`), resolves each id's spectra file to a supplied path by **bare file
  name** (the same extension+directory strip the psmtsv `File Name` column applies — `raw_paths`
  order is immaterial, extra paths ignored, unmatched files error), then delegates to the existing
  `run_msms` (no per-file-loop duplication). New `QuantError` enum (`PsmTsv`/`Io`/
  `MissingSpectraFile{file_name,supplied}`/`Parquet`) with `Display` + `From` for clean error text.
  Helper `engine::quant_to_parquet(results_path, raw_paths, output_path)` wraps `quant` +
  `write_peptide_results_parquet`. Made `psm_tsv::file_name_without_extension` `pub` so the resolver
  keys raw paths exactly as the parser keys ids. **PyO3 binding** `flashlfq_py.quant(results_path,
  raw_paths, output_path=None, params=None)`: with `output_path` → writes the L6 Parquet (tidy long
  form) and returns the path `str`; with `None` → returns the results as an in-memory
  `pyarrow.RecordBatch` (via `peptide_results_record_batch` + `into_pyarrow`). `params` is accepted
  for a stable signature but **must be None/empty** (Phase 1 = default MS2 path only) else
  `ValueError`; a bad psmtsv / unmatched spectra file / write failure also map to `ValueError`.
  **2 new core unit tests** (`resolve_file_to_mzml` bare-name matching incl. mixed dirs/exts +
  dedup + extra-path-ignored; missing-file error) → **87 core unit + L1–L6 all green**. **End-to-end
  acceptance** `rust/flashlfq-py/quant_smoke_test.py` runs `maturin develop`'s module over
  `AllPSMs.psmtsv` + the two K562 mzMLs: **708 rows**, detection mix **450/223/33/2** (matches L6
  golden), Parquet-path == in-memory-RecordBatch identical, and both error surfaces fire → prints
  `ALL QUANT SMOKE CHECKS PASSED`. Repro: from `rust/flashlfq-py/` with `rust/.venv` active run
  `maturin develop`, then `& "F:\flashlfq-rust\rust\.venv\Scripts\python.exe" rust/flashlfq-py/quant_smoke_test.py`.
- [x] **P1.19 — IndexingEngine/XIC bindings.** Expose index build/load, `GetIndexedPeak`, and
  an XIC extractor returning two `f64` NumPy arrays (RT, intensity) + integration bounds.
  _Done when:_ a Python script plots one peptide's XIC.
  _Result:_ **Green — a Python script plots a real peptide's XIC end-to-end.** Two new PyO3
  `#[pyclass]`es in `rust/flashlfq-py/src/lib.rs`: **`PeakIndex`** wraps the core
  `PeakIndexingEngine` (built once, reused — it is not cheap to rebuild), exposing
  `PeakIndex.from_mzml(path)` (staticmethod; `ValueError` on unreadable file / no MS1 peaks),
  `num_ms1_scans` getter, `get_indexed_peak(mz, scan_index, ppm=5.0) -> (mz, intensity, rt,
  scan_index) | None`, and `get_xic(mz, retention_time, ppm=20.0, missed_scans_allowed=1,
  max_peak_half_width=None) -> Xic` (defaults mirror the engine constants; `None` half-width →
  `i32::MAX` minutes as in C#). **`Xic`** carries the RT-ascending trace as three parallel `f64`
  NumPy arrays (`retention_times`/`intensities`/`mzs`, returned via `PyArray1::from_slice` so the
  getters can be called repeatedly) plus the **integration bounds** `start_rt`/`apex_rt`/`end_rt` +
  `apex_intensity`/`apex_scan_index` (computed in the binding from the core peak list; empty trace →
  NaN bounds, `apex_scan_index = -1`), with `__len__`/`__repr__`. Both registered in the
  `#[pymodule]`. The bindings stay thin (translation only); all tracing logic stays in core.
  **Acceptance** `rust/flashlfq-py/xic_smoke_test.py`: reads the first usable id from the live
  `AllPSMs.psmtsv` (`AHQLVMEGYNWC[Carbamidomethyl]HDR`, z=3, file `…K562_3`), builds the index
  (**73 MS1 scans**), computes monoisotopic m/z 639.28240, traces an **18-peak XIC** (apex RT
  79.1327, intensity 2.06e6, scan 2), asserts the arrays are f64 ndarrays of equal length,
  RT-ascending, bounds == array min/max/argmax, the apex peak round-trips through
  `get_indexed_peak`, an absent m/z → empty XIC with NaN bounds, and **writes `xic_smoke_plot.png`**
  (a clean elution peak with apex + integration-bound shading) → prints `ALL XIC SMOKE CHECKS
  PASSED`. Repro: from `rust/flashlfq-py/` with `rust/.venv` active run `maturin develop`, then
  `& "F:\flashlfq-rust\rust\.venv\Scripts\python.exe" rust/flashlfq-py/xic_smoke_test.py` (installs
  matplotlib for the PNG; arrays are plot-ready without it). The PNG is gitignored
  (`rust/.gitignore`).
- [x] **P1.20 — Stub MBR binding.** Placeholder call/flag so the API shape is fixed. _Done
  when:_ it exists and raises a clear "not implemented" until Phase 3.
  _Result:_ **Green — the MBR switch exists and refuses clearly.** Added a
  `match_between_runs=False` kwarg to the existing `flashlfq_py.quant(...)` binding (PyO3
  signature `(results_path, raw_paths, output_path=None, params=None, match_between_runs=false)`
  in `rust/flashlfq-py/src/lib.rs`). This is the **permanent home** for the MBR switch, mirroring
  C#'s `FlashLfqParameters.MatchBetweenRuns` flag (MBR is a mode of the quant run, not a separate
  entrypoint), so the API shape is fixed before Phase 3 fills it in. When `True` the binding
  returns early with `PyNotImplementedError` (Python `NotImplementedError`): *"match_between_runs
  is not implemented until Phase 3; pass match_between_runs=False (the default) for MS2-only
  quantification"* — it does **not** silently fall back to a plain MS2 quant. The default
  (`False`) leaves the existing P1.18 path untouched (the MBR guard sits ahead of the
  `params` check and `core_quant` call). Verified after `maturin develop`: `quant(...,
  match_between_runs=True)` raises `NotImplementedError`; `quant(...)` with a default flag still
  reaches the real path (raising `ValueError` on a bad file, not `NotImplementedError`). Thin
  binding only — no core changes. **Phase 1 §1c (and all of Phase 1) is now complete.**

## Phase 2 — Thermo `.raw` + hardening

- [x] **P2.1 — `.raw` via `mzdata`.** Swap the reader for Thermo `.raw`; the IndexingEngine is
  format-agnostic. Document the **.NET runtime requirement** (mzdata bridges to
  ThermoRawFileParser — not a pure-Rust wheel). _Done when:_ a `.raw` from TestData quants and
  matches the mzML result within tolerance.
  _Result:_ **Green — `.raw` quant is bit-identical to the matching mzML.** The reader boundary
  was the only change: `flashlfq-core/src/peak_indexing.rs` now opens spectra via mzdata's
  **`MZReader::open_path`** (format-sniffing) instead of `MzMLReader`, so `read_ms1_scans` reads
  **mzML and Thermo `.raw`** transparently (the `< 0.01` zero-intensity filter and everything
  downstream are unchanged). Added the `thermo` feature to the `mzdata` dep in
  `flashlfq-core/Cargo.toml`. **.NET runtime requirement documented** in both the Cargo.toml
  comment and the `read_ms1_scans` doc-comment: `.raw` rides the `thermorawfilereader` crate, which
  hosts a **self-contained .NET 8 runtime** (verified present: `dotnet --list-runtimes` shows
  `Microsoft.NETCore.App 8.0.x`), so a build/run touching `.raw` is **not a pure-Rust wheel** and
  needs a .NET 8 runtime on the machine — **mzML stays pure-Rust**. The `thermo` feature's
  `nethost-download` fetched the .NET host at build time (network was available; build took ~44 s
  the first time, pulling `thermorawfilereader`/`netcorehost`/`nethost-sys`/`dotnetrawfilereader-sys`).
  `PeakIndexingEngine::from_spectra_file` is the new canonical entry; `from_mzml` is kept as a thin
  alias (now also reads `.raw`) so the P1.x callers + the Python `PeakIndex.from_mzml` binding keep
  working. `engine::run_msms` now calls `from_spectra_file`; **`quant()` already supports `.raw`
  with no binding change** (it resolves files by bare name regardless of extension, then reads via
  `run_msms`). **Parity gate** `tests/p2_raw_vs_mzml_parity.rs`: the head-to-head fixtures
  `sliced-raw.raw` / `sliced-mzml.mzML` are the *same* spectra in two containers (mzLib's own
  `TestFlashLFQ.TestFlashLfq` runs `EGFQVADGPLYR` against both and asserts the rounded intensities
  equal, `TestFlashLFQ.cs:93-95`). The test feeds the same 4-PSM set to `run_msms` keyed once to the
  `.raw` and once to the `.mzML` → **both peptide intensities = 3367919.3125, relative difference
  0e0** (bit-identical, not just within tolerance). **87 core unit + L1–L6 + the new P2 gate all
  green; full workspace (incl. flashlfq-py) builds.** _Note for later:_ the `thermo` feature is now
  always-on for `flashlfq-core`, so the whole workspace build now requires network on a clean build
  (nethost download) and the Python wheel now carries the .NET-bridge native deps. If a pure-Rust /
  offline mzML-only build is ever wanted, gate `thermo` behind an opt-in cargo feature.
- [x] **P2.2 — FDR hardening.** Edge cases, decoy handling, q-value correctness.
  _Result:_ **Green — decoy / quantify-set semantics now faithful to C#, gated by a real-C# decoy
  parity test.** The Phase-1 port conflated two distinct C# sets (legal only because the corpus is
  all targets); P2.2 splits them in `src/results.rs`:
  - **`engine_quantify_set(ids)`** = every distinct modseq **incl decoys** = the C# `FlashLfqEngine`
    ctor default `allIdentifications.Select(modseq)` (`FlashLfqEngine.cs:94`). This is the membership
    gate for `run_error_checking` merge/overwrite and the `calculate_peptide_results` filters. Decoys
    belong here so they drive error-checking exactly as C# (the MSMS-vs-MSMS branch checks only
    quantify-set membership, not `DecoyPeptide`).
  - **`output_peptide_sequences(ids, quantify_set)`** = `ids.Where(!is_decoy &&
    quantify_set.contains(modseq)).map(modseq)` = the C# `FlashLfqResults` ctor
    `PeptideModifiedSequences` keys (`FlashLfqResults.cs:42`) — the rows that get emitted. Decoy-only
    modseqs never appear even though they sit in the quantify set.
  - `default_quantify_set` (non-decoy modseqs) kept but redocumented as the *`FlashLfqResults` ctor
    fallback* (`FlashLfqResults.cs:33`, null-set path only) — not used by the engine path anymore.
  `calculate_peptide_results` now takes both sets (`quantify_set` gate + `peptide_sequences` output
  rows; NotDetected cells seeded from the output set). `engine::run_msms` builds both and wires them
  through. **q-value:** the MS2 path computes no per-peak q-value (confirmed against C#: FDR control
  is entirely the externally-supplied quantify set; the engine default applies *no* q-filter — it
  quantifies every modseq). The MBR q-value FDR (`CalculateFdrForMbrPeaks`/`CorrectQValues`) stays
  Phase 3 (P3.4). **Tests:** 4 new decoy unit tests in `results.rs` (engine set incl decoys; output
  set excl decoys/unquantified; decoy-only modseq produces no row; decoy in quantify set overwrites
  an unquantified target in error checking) → **91 core unit** (was 87). **Real-C# decoy parity:** the
  golden generator (`parity/csharp_golden/Program.cs`) now also flips a deterministic third of the
  modseqs to decoy (distinct modseqs sorted Ordinal, `index % 3 == 0` → decoy; 118 of 354), re-runs
  the **real `FlashLfqEngine`**, and dumps `golden/L6_peptide_intensities_decoy.tsv` (**472 rows** =
  236 non-decoy modseqs × 2 files). New gate `tests/l5_l6_parity.rs::l6_decoy_flip_matches_csharp_golden`
  applies the identical flip (Rust `str` Ord == Ordinal for these ASCII seqs), runs `run_msms`, asserts
  no decoy modseq appears as a row, and diffs every surviving cell vs the golden — **all 472 match**
  (intensity/RT rel-1e-6, detection-type exact). L5/L6 all-target goldens regenerated **byte-identical**
  (deterministic), so the existing gates stay green. Regenerate:
  `cd rust/flashlfq-core/parity/csharp_golden; dotnet run -c Release`.

## Phase 3 — MBR (Rust search + Python model)

- [x] **P3.1 — RT alignment.** Port `GetRtCalSpline`. _Done when:_ alignment spline matches C#.
  _Result:_ **Green — the donor→acceptor RT calibration spline matches C# over real engine peaks
  (113/113 anchor points exact).** New module `rust/flashlfq-core/src/mbr.rs` ports
  `FlashLfqEngine.GetRtCalSpline` (`FlashLfqEngine.cs:559`) + `ChooseBestPeak` (`:670`) +
  `RetentionTimeCalibDataPoint`. `get_rt_cal_spline(donor_peaks, acceptor_peaks,
  donor_q_value_threshold, criterion) -> RtCalSpline { calibration_curve, anchor_rt_diffs,
  donor_best_peaks_ordered_by_mass }`: from each run keep anchor candidates
  (`num_identifications_by_full_seq==1 && DetectionType==MSMS && !envelopes.is_empty() &&
  min(id.q_value) < threshold`), group by first modseq (first-seen order preserved like LINQ
  `GroupBy`), `choose_best_peak` per group, pair acceptor↔donor on shared modseq into
  `RetentionTimeCalibDataPoint`s (`rt_diff = acceptorApexRT − donorApexRT`), collect
  `donorRT−acceptorRT` anchor diffs (both apex RTs > 0), order donor best-peaks by first id
  peakfinding mass, and sort the curve by donor apex RT. `choose_best_peak` ports the default
  `DonorCriterion::Score` (MaxBy max-id-score, first on ties; fall through to Intensity when the
  chosen peak's first id score ≤ 0 — verified all 594 K562 scores > 0, so Score always wins) +
  `Intensity`; `Neighbors` is `unimplemented!` (default path is Score). Constants verbatim:
  `DONOR_Q_VALUE_THRESHOLD=0.01` (`FlashLfqParameters.cs:36`), `NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR=3`,
  `MBR_ALIGNMENT_WINDOW=2.5`. **The q-filter is load-bearing:** 18 of 594 K562 ids have QValue ≥ 0.01
  and are excluded. **Golden:** the C# generator (`parity/csharp_golden/Program.cs`) gained
  `WriteMbrSpline` — a faithful replica of the *private* `GetRtCalSpline`/`ChooseBestPeak` (Score)
  over the **real** engine `results.Peaks` (genuine C# peaks; the method is private and needs an
  `MbrScorer`, a P3.3 type, so it is replicated not invoked) for donor=K562_3, acceptor=K562_4 →
  `golden/MBR_rt_cal_spline.tsv` (113 points: donor_modseq, donor_apex_rt, acceptor_apex_rt,
  rt_diff, donor_peakfinding_mass). **Gate** `tests/mbr_parity.rs::rt_cal_spline_matches_csharp_golden`
  runs `run_msms` over both K562 files, builds the spline, sorts both sides by (donor_apex_rt,
  modseq) (tie-robust — peptides can share a scan RT), and diffs every point: modseq exact, RTs/
  rt_diff/mass rel-1e-6 → **all 113 match**. **95 core unit** (4 new mbr tests) + L1–L6 + P2 +
  the new MBR gate all green; L5/L6/decoy goldens regenerated byte-identical. Regenerate:
  `cd rust/flashlfq-core/parity/csharp_golden; dotnet run -c Release`.
- [~] **P3.2 — Donor selection + acceptor search.** Port `FindPeptideDonorFiles`,
  `FindAllAcceptorPeaks`, `FindIndividualAcceptorPeak` → emit a **feature table** (one row per
  candidate transferred peak) to Parquet. **Decomposed into sub-tasks (do in order):**
  - [x] **P3.2a — `RtInfo` + `PredictRetentionTime` (local RT alignment) + stats helpers.**
    _Result:_ **Green — predicted acceptor RTs match C# 245/245 over the real K562 donor/acceptor
    spline.** Added to `src/mbr.rs`: **`RtInfo`** (`predicted_rt`/`width` + `rt_start_hypothesis`/
    `rt_end_hypothesis` = `MBR/RtInfo.cs`); **`predict_retention_time(curve, donor_peak,
    max_mbr_rt_window, number_of_anchor_peptides) -> RtInfo`** = faithful port of
    `FlashLfqEngine.PredictRetentionTime` (`FlashLfqEngine.cs:716`): a `.NET Array.BinarySearch`
    replica (`binary_search_donor_rt`, exact `lo/hi` midpoint algorithm so the chosen index on
    duplicate donor RTs matches C#; compares donor apex RT per `RetentionTimeCalibDataPoint.CompareTo`),
    the forward-then-backward anchor gather (≤`NUMBER_OF_ANCHOR_PEPTIDES_FOR_MBR=3` per side, break
    when a neighbour's donor-RT diff > 0.5 min, **the element at the insertion index is excluded
    from both directions** — faithful C# quirk), then `donorRt`/`donorRt−medianDiff` with width
    `0.25` (0/1 anchors) or `min(6·stddev, MaxMbrRtWindow)`. **Stats helpers** `median` (MathNet
    order-statistic, even = avg of two central) and `standard_deviation` (sample SD via
    `StreamingStatistics.Variance`'s running `j`/`t` accumulation — bit-faithful order for a
    `List<double>`). New consts `MAX_MBR_RT_WINDOW=1.0`, `MBR_PPM_TOLERANCE=10.0`
    (`FlashLfqParameters.cs:31-32`). **The fraction gate is omitted on purpose** (core peak has no
    `SpectraFileInfo.Fraction`; that early-return belongs to the P3.2d orchestration — both K562
    files are unfractionated anyway). **Golden:** the C# generator (`parity/csharp_golden/Program.cs`)
    gained `PredictRtReplica` (faithful replica of the *internal* `PredictRetentionTime` over real
    `RetentionTimeCalibDataPoint[]` objects — internal, so replicated not invoked) + `BestBySeq`
    helper + `WriteMbrPredictedRt`: predicts every donor best-peak (ordered by peakfinding mass)
    for donor=K562_3/acceptor=K562_4 → `golden/MBR_predicted_rt.tsv` (245 rows: donor_modseq,
    predicted_rt, width). **Gate** `tests/mbr_parity.rs::predicted_rt_matches_csharp_golden` runs
    `run_msms` over both files, builds the spline, predicts each donor best-peak, sorts by modseq
    (unique key) and diffs predicted_rt + width rel-1e-6 → **all 245 match**. **102 core unit** (7
    new mbr tests) + L1–L6 + P2 + both MBR gates green; all other goldens regenerated byte-identical.
    Regenerate: `cd rust/flashlfq-core/parity/csharp_golden; dotnet run -c Release`.
  - [x] **P3.2b — `MbrChromatographicPeak` + peak feature fields.** _Result:_ **Green — the MBR peak
    type + base accessors are in place; 105 core unit (3 new) + L1–L6 + P2 + both MBR gates all
    pass.** Added `scan_count()` (= `IsotopicEnvelopes.Count`) and `isotopic_pearson_correlation()`
    (= `Apex?.PearsonCorrelation ?? -1`) accessors to `ChromatographicPeak` (`src/chromatographic_peak.rs`)
    — the two feature getters `MbrScorer` reads. New module `src/mbr_chromatographic_peak.rs` ports
    `MBR/MbrChromatographicPeak.cs`: since Rust has no inheritance, **`MbrChromatographicPeak`
    *composes*** a base `ChromatographicPeak` (`pub peak`, `detection_type = DetectionType::MBR`)
    plus the MBR-exclusive fields: `predicted_retention_time` (C# `init`), `random_rt` (C# read-only
    `RandomRt`), `rt_prediction_error`, `mbr_q_value`, `mbr_pep: Option<f64>`, `mbr_score`, and the
    five component scores `ppm_score`/`intensity_score`/`rt_score`/`scan_count_score`/
    `isotopic_distribution_score`, plus `charge_list: Vec<i32>` (C# base `ChargeList`, set by the
    P3.2d acceptor search). Constructor `new(id, spectra_file_name, peakfinding_mass,
    predicted_retention_time, random_rt)` mirrors the C# ctor (builds the base via
    `from_identifications([id], [pfm], MBR)`; all scores/q-value default `0.0`, `mbr_pep` `None`,
    `charge_list` empty). **`spectra_file_name: String`** stands in for the C# `SpectraFileInfo`
    reference (the base peak carries no file info in this port) — the acceptor-file key the P3.2d
    loop/scorer group by. `is_decoy_peak()` returns `random_rt` (C# doc: a randomized RT "implies
    this peak is a decoy peak identified by the MBR algorithm"), distinct from `peak.decoy_peptide()`
    (decoy *peptide* id). Delegating accessors `intensity()`/`apex_retention_time()`/`mass_error()`/
    `scan_count()`/`isotopic_pearson_correlation()` forward to the base for the scorer. Module
    registered in `lib.rs`. 3 new unit tests (defaults + MBR detection type; `random_rt`→decoy peak;
    empty-base-peak accessors). _Read:_ `MBR/MbrChromatographicPeak.cs`, `ChromatographicPeak.cs`
    (ScanCount/IsotopicPearsonCorrelation getters), `MBR/MbrScorer.cs` (confirmed which fields the
    scorer reads/writes).
  - [x] **P3.2c — `MbrScorer` (statistical distributions + `ScoreMbr`).** _Result:_ **Green — the MBR
    scorer + factory are ported; 120 core unit (15 new) + L1–L6 + P2 + both MBR gates all pass; full
    workspace (incl flashlfq-py) builds.** Two new modules:
    - **`src/special_functions.rs`** — the MathNet special functions the distributions reduce to.
      `gamma_ln` (Lanczos, constants verbatim from `SpecialFunctions.GammaLn`),
      `gamma_lower_regularized`/`gamma_upper_regularized` (Cephes series + continued-fraction port of
      `GammaLowerRegularized`/`GammaUpperRegularized`, incl. `big`/`bigInv` rescaling, ε=1e-15), and
      `erf`/`erfc` expressed via the **exact identities** `erf(x)=P(½,x²)` / `erfc(x)=Q(½,x²)` (the
      `Q` route avoids large-x cancellation) rather than MathNet's Boost rational polynomial — both
      accurate to ≈1e-13, far inside the rel-1e-6 budget, and the scorer's Normal CDF only ever
      evaluates `erfc` of a non-negative arg. Value types **`Normal`** (`cumulative_distribution` =
      `0.5·erfc((mean−x)/(stddev·√2))`, `is_valid_parameter_set` = `stddev≥0 && !mean.is_nan()`) and
      **`Gamma`** (rate-parameterized; CDF = `P(shape, x·rate)`). A C#-faithful `log2(x)=ln(x)/ln(2)`
      (NOT fused `f64::log2`). 6 tests (known erf/erfc/gammaln/Normal-CDF/Gamma-CDF values + validity).
    - **`src/mbr_scorer.rs`** — `MbrScorer` + `build_mbr_scorer` (= `MbrScorerFactory.BuildMbrScorer`).
      Fits the four acceptor-file distributions in `initialize_scorer` exactly as C#: ppm→`Normal`
      (median center, IQR/1.36 or stddev spread, raw stats stashed for `get_ppm_error_tolerance`),
      log2-intensity→`Normal`, scan-count→`Normal` (**LINQ `Average` mean, not incremental `Mean`**),
      `1−isoCorr`→`Gamma` (method-of-moments α=mean²/var, β=mean/var via MathNet incremental
      `Mean`/`Variance`). `calculate_fold_change_between_files` (incl. the verbatim "replace with the
      *less* intense peak" quirk, ≥100-pair gate) and `add_rt_pred_error_distribution` add the
      per-donor-file dists. **`ScoreMbr`** writes the five component scores + `rt_prediction_error`
      and returns `100·(∏ scores)^0.2`; `CalculateScore(Normal)`=`2·CDF(mean−|mean−v|)` floored at
      `_minScore=3e-7`, `CalculateScore(Gamma)`=`1−CDF(v)`. **The `sigma=1` RT-dist bug
      (`MbrScorer.cs:288`) is replicated** — `sigma` is computed then discarded; dist built
      `Normal(medianRtError, 1)`. 8 tests (init gate, factory build/clamp/None, score-in-range,
      min-score on NaN/null/negative, score=1 at the mean, **sigma=1-bug assertion**, default RT dist).
    - **Stats helpers added to `src/mbr.rs`** (`pub(crate)`, MathNet-faithful, siblings of the existing
      `median`/`standard_deviation`): `variance` (StreamingStatistics, = stddev²), `mean` (incremental
      `Mean`), `average` (LINQ sum/count — kept distinct from `mean`), `quantile` (type-8
      `QuantileInplace`) + `interquartile_range`.

    **Divergence from C# (documented):** the base `ChromatographicPeak` has no `SpectraFileInfo`, so
    the per-donor-file dicts (logFC, RT-pred-error) are keyed by donor **file name** (`&str`);
    `calculate_fold_change_between_files`/`add_rt_pred_error_distribution`/`score_mbr`/`is_valid_for_donor`
    take a `donor_file: &str` (mirrors P3.2b's `spectra_file_name` stand-in). The apex→peak dict uses
    `EnvelopePeakKey` (structural apex identity) as the `DistinctBy(Apex.IndexedPeak)` key. **Read:**
    `MBR/MbrScorer.cs` + `MBR/MbrScorerFactory.cs` in full.
  - [x] **P3.2d — Acceptor search + MBR orchestration → feature table.** _Result:_ **Green — the full
    MBR transfer runs end-to-end over the real K562 corpus and emits a 156-row feature table (129
    target / 27 decoy peaks); 132 core unit (6 new mbr_search) + L1–L6 + P2 + both MBR gates + the new
    corpus smoke test all pass; full workspace (incl flashlfq-py) builds.** New module
    `src/mbr_search.rs` ports the remaining MBR pieces of `FlashLfqEngine`:
    `FindPeptideDonorFiles` → [`find_peptide_donor_files`] (pick one donor peak per quantifiable
    peptide, key by source file; `ChooseBestPeak` Score+Intensity-fallback via `choose_best_index`,
    tracking the file the peak came from); `GetRandomPeak` → [`get_random_peak`] (the 5–11-H mass
    window, ×10 widening to 1e5, verbatim `(int)(1e5·(pfm%1)·(rt%1)) % count` pseudo-random draw —
    H mass from `periodic_table().element_by_symbol("H").principal_isotope().atomic_mass`);
    `FindIndividualAcceptorPeak` → `find_individual_acceptor_peak` (seed from the **least-intense**
    envelope — the OrderBy-ascending+First C# quirk — GetXic→GetIsotopicEnvelopes→CalculateIntensity→
    CutPeak, claimed-peak removal incl. the seed, `ApexToAcceptorFilePeakDict` skip, ScoreMbr);
    `FindAllAcceptorPeaks` → `find_all_acceptor_peaks` (scan-snip to the predicted-RT window, charge
    set = distinct donor precursor charges + apex charge, per-charge GetIndexedPeak XIC, the
    seed-and-remove while loop, best-by-MbrScore + ChargeList); and the
    `QuantifyMatchBetweenRunsPeaks` driver → `quantify_mbr_for_acceptor` + `finalize_acceptor_peaks`
    + top-level `run_mbr`. The driver: build the scorer from the acceptor's quantifiable MS/MS peaks
    (`build_mbr_scorer`, then **mbrTol overwritten with the flat `MbrPpmTolerance=10`** per
    `FlashLfqEngine.cs:895`), per donor file build the spline + `add_rt_pred_error_distribution` +
    `is_valid_for_donor` gate, per donor best-peak `predict_retention_time` → target search + decoy
    (`get_random_peak`→predict→search) + the window-widening retry loop, then dedup (best-per-apex),
    `msmsImsPeaks` conflict skip, GroupBy(RandomRt) best-result selection, and different-charge merge.
    **Feature table** = `FeatureRow` (donor mod/base seq, acceptor file, predicted/apex RT, intensity,
    the five component scores, combined `mbr_score`, mass error, scan count, isotopic correlation,
    `rt_prediction_error`, `random_rt`/decoy-peptide flags) — one row per surviving transferred peak.
    **Documented divergences (all in the module header):** single-condition/unfractionated path (the
    fraction gate + cross-condition fold-change branch + `RequireMsmsIdInCondition` never fire under
    default params — matches the golden generator's setup, and the core peak has no `SpectraFileInfo`);
    the `Parallel.ForEach` is run **sequentially** (matches `MaxThreads=1`); the C# decoy-call bug
    (`FlashLfqEngine.cs:986/:1004` passes the *target* `rtInfo` width with the decoy's randomRt centre)
    is replicated; apex-less peaks are dropped (a C# `Dictionary` can't key a null apex either); the
    C# double-add (lines 1092+1110) is collapsed to one add; the final per-acceptor MBR `RunErrorChecking`
    pass is **not** re-run (the apex-claim guards already remove the MS/MS conflicts that pass handles).
    Engine wiring: `Ms2QuantResult` now carries `engines_by_file` (the `IndexingEngineDictionary` —
    `run_msms` keeps each `PeakIndexingEngine` alive) + `peptide_sequences_to_quantify`. New corpus
    smoke test `tests/mbr_search_corpus.rs` runs `run_msms`→`run_mbr` over the two K562 mzMLs and
    asserts the table is well-formed (targets+decoys present, every component score in (0,1], combined
    in [0,100], one row per stored MBR peak, all `DetectionType::MBR`). 6 new unit tests
    (donor-file selection by score + file, q/quantify exclusion, random-peak mass/seq window + none,
    apexless-drop, MBR detection type). **No C# golden yet:** the orchestration is private + parallel,
    so exact row parity needs a faithful C# replica (or `InternalsVisibleTo` instrumentation) — folded
    into **P3.2e** or a follow-up; this task validates the port runs end-to-end on real data.
  - [x] **P3.2e — Parquet feature table + Python read.** _Result:_ **Green — the MBR feature table
    serializes to Parquet and round-trips through pyarrow, both from Rust and through the Python
    binding.** Extended `src/parquet_output.rs` with `feature_table_schema()` (18 non-nullable
    columns in `FeatureRow` field order: 3 Utf8, 12 Float64, 1 UInt64 `scan_count`, 2 Boolean
    `random_rt`/`decoy_peptide`), `feature_table_record_batch(&[FeatureRow])` (sorts rows by
    `(acceptor_file, donor_modified_sequence, random_rt, predicted_rt, apex_rt, intensity)` —
    `total_cmp` on floats — for reproducible output), and `write_feature_table_parquet`. New core
    helper `engine::quant_mbr(results_path, raw_paths) -> Result<MbrResult, QuantError>` (= `quant`
    then `run_mbr`, reusing the kept-alive per-file engines). **Binding (P1.20 stub made real):**
    `flashlfq_py.quant(..., match_between_runs=True)` now runs MBR and returns the **feature table**
    (Parquet path `str` with `output_path`, else in-memory `pyarrow.RecordBatch`) instead of raising
    `NotImplementedError`; the MS2 (`match_between_runs=False`) path is untouched and still returns
    the peptide×file table. New example `examples/write_mbr_feature_table.rs` runs the real corpus →
    **156 feature rows (129 target / 27 decoy)** to Parquet. Python acceptance
    `rust/flashlfq-py/mbr_smoke_test.py` (run with `rust/.venv`) drives `quant(match_between_runs=True)`
    both ways, asserts the 18-col schema + 156/129/27 + component scores in (0,1] + combined in [0,100]
    + Parquet==in-memory + MS2 path still returns 708 peptide rows → **ALL MBR SMOKE CHECKS PASSED**.
    2 new parquet_output unit tests (sorted batch; Arrow round-trip with bool/uint/float cells); 130
    core unit + L1–L6 + P2 + both MBR gates + corpus smoke all green; full workspace (incl flashlfq-py
    via `maturin develop`) builds. **No cell-parity golden** (the orchestration still has no C# golden —
    see the P3.2d note); this validates serialization + the Python read path end-to-end.
- [x] **P3.3 — Python ML model.** Train/score the PEP model in Python (`lightgbm`/`sklearn`) on
  the feature table, replacing `Microsoft.ML` FastTree. _Done when:_ scores returned to Rust.
  _Result:_ **Green — the Python PEP model trains on the MBR feature table, scores every peak, and
  the scores round-trip back into the Rust peaks; the corpus smoke test passes end-to-end.**
  New Python module `rust/flashlfq-py/pep_model.py` ports `FlashLFQ.PEP.PepAnalysisEngine.
  ComputePEPValuesForAllPeaks` onto **scikit-learn** (`HistGradientBoostingClassifier`, the
  FastTree-shaped histogram GBDT): `compute_pep(table) -> list[float]` takes the feature table (a
  `pyarrow.RecordBatch`/`Table`) and returns one PEP (`1 - P(target)`) per row **in table order**.
  Faithful structural port: the 10 features (`PpmErrorScore`/`IntensityScore`/`RtScore`/
  `ScanCountScore`/`IsotopicDistributionScore` + `PpmErrorRaw=|MassError|`/`IntensityRaw=log2(Int)`/
  `RtPredictionErrorRaw=|·|`/`ScanCountRaw`/`IsotopicPearsonCorrelation` — = `ChromatographicPeakData`
  "standard" set), donor grouping (`_build_donor_groups`/`_order_donor_groups` = the
  `OrderByDescending(MbrScore).GroupBy(donor).` + target-count/decoy-count/best-score ordering),
  the PIP-score cutoff (`floor(N·0.25)`), the 3-fold round-robin split **+ the target/decoy
  equalization swap** (`_get_donor_group_indices`/`_equalize_donor_group_indices`/`_group_swap`,
  incl. the C# `>=`-in-surplus vs `>`-in-candidate asymmetry and the `continue`-skips-min-index
  control flow), and the **10 training passes** (iteration 0 selects positives by MBR score,
  iterations 1-9 by the previous round's PEP; predict-on-held-out-fold → `pep = 1 - prob`). FastTree
  hyper-params mapped: 100 trees→`max_iter`, 20 leaves→`max_leaf_nodes`, min-10/leaf→
  `min_samples_leaf`, LR 0.2, `UnbalancedSets`→`class_weight="balanced"`, seed 42. **Scores returned
  to Rust:** new core `mbr_search::apply_mbr_pep(&mut MbrResult, &[f64])` writes the peps onto
  `MbrChromatographicPeak.mbr_pep` (the C# `Compute_PEP_For_All_Peaks` sink), mapping table-order
  peps back onto peaks via the shared `parquet_output::feature_table_sort_order` (extracted from
  `feature_table_record_batch` so the table Python receives and the write-back agree exactly). New
  parquet helpers `feature_table_with_pep_schema`/`feature_table_record_batch_with_pep`/
  `write_feature_table_with_pep_parquet` (append a trailing **nullable** `mbr_pep` Float64 col → 19
  cols). **Binding:** `flashlfq_py.quant(..., match_between_runs=True, pep_model=<callable>)` now does
  the genuine Rust→Python→Rust round trip in one call — build the 18-col table, hand it to
  `pep_model`, take back the peps, `apply_mbr_pep`, return the **19-col** table (or 19-col Parquet);
  with `pep_model=None` the path is byte-unchanged (18 cols, P3.2e safe). Acceptance
  `rust/flashlfq-py/pep_smoke_test.py` over the real K562 corpus → **156 PEPs, all in [0,1], decoys
  separate from targets (mean PEP target 0.589 vs decoy 0.936, ROC-AUC 0.716), Rust `mbr_pep` ==
  a direct `compute_pep` recompute (deterministic), Parquet==in-memory, MS2 path untouched →
  ALL PEP SMOKE CHECKS PASSED**. 3 new Rust unit tests (apply_mbr_pep maps table-order/ rejects
  length mismatch; with-pep batch schema+null). **131 core unit** + L1-L6 + P2 + both MBR gates +
  corpus smoke all green; full workspace (incl flashlfq-py via `maturin develop --release`) builds.
  **Env note:** `scikit-learn 1.9.0` installed into `rust/.venv`
  (`& "F:\flashlfq-rust\rust\.venv\Scripts\python.exe" -m pip install scikit-learn`); `lightgbm` not used.
  **Divergences (documented in the module header):** donor grouping keys on
  `donor_modified_sequence` (the table's available proxy for the C# donor `Identification` ref); the
  C# `swappedDonors`-never-populated latent bug is replicated; exact FastTree cell-parity is not
  attainable across ML frameworks and there is **no C# golden** for the MBR orchestration that feeds
  this — the smoke test asserts the model learns + the round trip is faithful, not cell-parity.
- [x] **P3.4 — Fold scores into FDR.** Combine model scores with `CalculateFdrForMbrPeaks`.
  _Done when:_ MBR peaks carry q-values.
  _Result:_ **Green — every transferred MBR peak now carries an `mbr_q_value`; 137 core unit (6 new) +
  L1–L6 + P2 + both MBR gates + corpus smoke (now 2 tests, incl. the new FDR assertion) all pass; full
  workspace (incl flashlfq-py) builds.** Ported `FlashLfqEngine.CalculateFdrForMbrPeaks` /
  `EstimateFdr` / `EstimateDecoyPeptideErrors` / `CorrectQValues` (`FlashLfqEngine.cs:1426–1528`) into
  `src/mbr_search.rs`:
  - `estimate_fdr(doubleDecoys, decoyPeptides, decoyPeaks, total) = (1 + decoyPeaks + max(0,
    decoyPeptides − doubleDecoys)) / total` (the double-decoy correction is `usize::saturating_sub`);
    `correct_q_values` = the bottom-up running-min so q-values are monotone non-decreasing down the
    score-ordered list; `round6` = `Math.Round(_, 6)` via `round_ties_even` (half-to-even, the C#
    default `MidpointRounding`).
  - **`calculate_fdr_for_mbr_peaks(&mut Vec<MbrChromatographicPeak>, use_pep)`** — per acceptor file.
    `use_pep=true`: dedup to the single best acceptor **per donor** (group by donor modified sequence —
    the available proxy for the C# `Identifications.First()` ref; one donor per peptide in this port, so
    it *is* the donor id), keeping the `OrderBy(MbrPep).ThenByDescending(MbrScore)` first, then re-sort
    the survivors the same way; the dropped peaks are removed (mirrors C#
    `_results.Peaks[acceptorFile] = filteredMbr.Concat(nonMbr)`). Target + decoy (`random_rt`)
    hypotheses of the same donor compete in this group ("acceptor can be target or decoy!").
    `use_pep=false`: just `OrderByDescending(MbrScore)`, drop nothing ("err on the safe side").
    Then walk the ordered list accumulating the `(decoy_peptide, random_rt)` 4-way bucket counts and
    assign each peak its corrected FDR q-value.
  - **`apply_mbr_fdr(&mut MbrResult, use_pep)`** runs it over every acceptor file (sorted-name order;
    C# loops `CalculateFdrForMbrPeaks(spectraFile, pepSuccesful)` at `:291`), and
    **`mbr_pep_analysis_succeeded(&MbrResult)`** = the `RunPEPAnalysis` gate (`>100` peaks **and** `>20`
    `random_rt` decoys) that decides `use_pep`. `MbrResult.feature_rows` is the *pre-FDR* table (already
    consumed by the P3.3 PEP model) and is left unchanged; q-values live on the peaks
    (`MbrChromatographicPeak.mbr_q_value`), which the PLAN note specified as a **pure-Rust** step over
    `mbr_peaks_by_file`.
  **Tests:** 6 new unit tests (`estimate_fdr` formula incl. double-decoy correction; `correct_q_values`
  monotone + empty; no-pep score-ordering + hand-computed q-walk; use_pep per-donor dedup keeps the
  better-pep target over the higher-score decoy; empty + double-decoy bucket; the `>`-strict success
  thresholds) + a new corpus integration test (`mbr_fdr_assigns_monotone_qvalues_to_every_peak`):
  `run_msms`→`run_mbr`→`apply_mbr_fdr(false)` over the two K562 mzMLs asserts the corpus clears the PEP
  gate, every peak's q-value ∈ [0,1], peaks are score-descending, and q-values are non-decreasing.
  **Caveat (carried from P3.2/P3.3):** there is still **no C# golden** for the MBR orchestration (private
  + parallel), so this is structural/algorithmic parity (the formula + ordering are a faithful line-by-
  line port), not cell-parity against a live C# dump. **Not wired into the Python binding's returned
  table** (no `mbr_q_value` column added to the 19-col feature table) — the q-values are core-side on the
  peaks, matching the PLAN's "pure Rust over `MbrResult.mbr_peaks_by_file`" scope; surfacing them to
  Python would be a small follow-up (append a nullable `mbr_q_value` Float64 col, like P3.3's `mbr_pep`).
