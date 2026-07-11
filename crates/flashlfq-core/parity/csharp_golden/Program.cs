using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using Readers;

// =============================================================================
// L5/L6 golden generator (PLAN.md P1.17 -- Phase 1 exit gate).
//
// Runs the real FlashLfqEngine over the Phase-1 corpus and writes two TSVs that
// the Rust parity gate diffs against:
//
//   L5_chromatographic_peaks.tsv : one row per ChromatographicPeak (post error
//       checking) -- per-(peptide, file) integrated intensity, apex scan index,
//       apex pearson correlation, apex charge, mass error, etc.
//   L6_peptide_intensities.tsv   : one row per (modified sequence, file) cell --
//       the end-to-end peptide x file intensity table (CalculatePeptideResults).
//
// Engine config = FlashLfqParameters defaults (the MS2-only path):
//   Normalize=false, PpmTolerance=10, IsotopePpmTolerance=5, Integrate=false,
//   NumIsotopesRequired=2, IdSpecificChargeState=false, QuantifyAmbiguousPeptides=false,
//   MatchBetweenRuns=false, IsoTracker=false.
// MaxThreads is forced to 1 so the run is fully reproducible (the numbers are
// per-id-independent and thus deterministic regardless, but single-threaded
// removes any doubt). Output rows are sorted deterministically on this side and
// must be sorted identically on the Rust side before diffing.
// =============================================================================

static class Program
{
    static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // Round-trippable double formatting so the Rust side parses the exact f64.
    static string D(double v) => v.ToString("R", Inv);

    static int Main()
    {
        // The port and mzLib are now separate projects, so there are two roots:
        //   * mzLib TestData comes from the EXTERNAL mzLib checkout ($MZLIB_DIR, default F:\mzLib).
        //   * the golden output dir lives inside THIS port (parity/golden), located by walking
        //     up from the compiled exe to flashlfq-core.
        string mzlibDir = MzLibDir();
        string testData = Path.Combine(mzlibDir, "mzLib", "Test", "FlashLFQ", "TestData");
        string goldenDir = Path.Combine(FindFlashLfqCore(), "parity", "golden");
        Directory.CreateDirectory(goldenDir);

        string psmPath = Path.Combine(testData, "AllPSMs.psmtsv");
        var mzml3 = new SpectraFileInfo(Path.Combine(testData, "20100614_Velos1_TaGe_SA_K562_3.mzML"), "a", 0, 0, 0);
        var mzml4 = new SpectraFileInfo(Path.Combine(testData, "20100614_Velos1_TaGe_SA_K562_4.mzML"), "a", 1, 0, 0);
        var spectraFiles = new List<SpectraFileInfo> { mzml3, mzml4 };

        Console.WriteLine($"Reading identifications from {psmPath}");
        IQuantifiableResultFile quant = FileReader.ReadQuantifiableResultFile(psmPath);
        List<Identification> ids = MzLibExtensions.MakeIdentifications(quant, spectraFiles);
        Console.WriteLine($"  {ids.Count} identifications, " +
                          $"{ids.Select(i => i.ModifiedSequence).Distinct().Count()} distinct modified sequences, " +
                          $"{ids.Select(i => i.FileInfo).Distinct().Count()} files");

        var p = new FlashLfqParameters
        {
            MaxThreads = 1,
            Silent = true,
        };
        var engine = new FlashLfqEngine(p, ids);
        Console.WriteLine("Running FlashLfqEngine (MS2 path, default params, single-threaded)...");
        FlashLfqResults results = engine.Run();

        WriteL5(results, Path.Combine(goldenDir, "L5_chromatographic_peaks.tsv"));
        WriteL6(results, Path.Combine(goldenDir, "L6_peptide_intensities.tsv"));

        // ---------------------------------------------------------------------
        // MBR RT calibration spline (PLAN.md P3.1 -- GetRtCalSpline).
        //
        // Build the donor=K562_3, acceptor=K562_4 RT calibration spline from the
        // *real* engine peaks (results.Peaks, genuine C# ground truth). The
        // spline-assembly logic here is a faithful replica of the private
        // FlashLfqEngine.GetRtCalSpline + ChooseBestPeak (Score criterion, the
        // FlashLfqParameters default) over those public peak lists -- the method
        // itself is private and takes an MbrScorer (a P3.3 dependency), so it is
        // replicated rather than invoked. The Rust port (mbr::get_rt_cal_spline)
        // reproduces the identical logic and the gate diffs the two splines.
        // ---------------------------------------------------------------------
        WriteMbrSpline(results, mzml3, mzml4, Path.Combine(goldenDir, "MBR_rt_cal_spline.tsv"));

        // ---------------------------------------------------------------------
        // MBR predicted retention times (PLAN.md P3.2 -- PredictRetentionTime).
        //
        // For each donor best-peak (donor=K562_3), predict where it would elute
        // in the acceptor run (K562_4) via the local-alignment strategy of
        // FlashLfqEngine.PredictRetentionTime (FlashLfqEngine.cs:716). That method
        // is internal and not callable from this assembly, so it is replicated
        // here over real RetentionTimeCalibDataPoint[] objects (built from the
        // same real engine peaks as WriteMbrSpline). The fraction gate is omitted
        // because both K562 files are unfractionated (single fraction per
        // condition/biorep), exactly as the Rust port documents.
        // ---------------------------------------------------------------------
        WriteMbrPredictedRt(results, mzml3, mzml4, Path.Combine(goldenDir, "MBR_predicted_rt.tsv"));

        // ---------------------------------------------------------------------
        // Decoy parity variant (PLAN.md P2.2 -- FDR / decoy hardening).
        //
        // The corpus is all targets, so the decoy code paths are otherwise
        // never exercised against C# ground truth. We mark a deterministic
        // subset of modified sequences as decoys -- distinct modseqs sorted
        // Ordinal, every 3rd one (index % 3 == 0) -- rebuild the Identification
        // list with that decoy flag, and re-run the *real* engine. The Rust
        // parity test reproduces the identical flip rule (Rust String Ord ==
        // Ordinal for these ASCII sequences) and diffs the L6 table.
        //
        // Expected effect: decoy modseqs vanish from the output peptide table
        // (FlashLfqResults excludes !IsDecoy ids), while the quantify set still
        // contains them (FlashLfqEngine default includes every modseq), so the
        // error-checking merges -- and therefore the surviving target rows --
        // match the all-target run except where a now-decoy id changes a peak's
        // first-id / ambiguity classification.
        // ---------------------------------------------------------------------
        var distinctModseqs = ids.Select(i => i.ModifiedSequence).Distinct()
            .OrderBy(s => s, StringComparer.Ordinal).ToList();
        var decoySet = new HashSet<string>();
        for (int i = 0; i < distinctModseqs.Count; i++)
        {
            if (i % 3 == 0) decoySet.Add(distinctModseqs[i]);
        }

        var decoyIds = ids.Select(id => new Identification(
            id.FileInfo, id.BaseSequence, id.ModifiedSequence, id.MonoisotopicMass,
            id.Ms2RetentionTimeInMinutes, id.PrecursorChargeState, id.ProteinGroups.ToList(),
            id.OptionalChemicalFormula, id.UseForProteinQuant, id.PsmScore, id.QValue,
            decoy: decoySet.Contains(id.ModifiedSequence))).ToList();

        var pDecoy = new FlashLfqParameters { MaxThreads = 1, Silent = true };
        var decoyEngine = new FlashLfqEngine(pDecoy, decoyIds);
        Console.WriteLine($"Running decoy-flip variant ({decoySet.Count} of {distinctModseqs.Count} modseqs marked decoy)...");
        FlashLfqResults decoyResults = decoyEngine.Run();
        WriteL6(decoyResults, Path.Combine(goldenDir, "L6_peptide_intensities_decoy.tsv"));

        Console.WriteLine("Golden generation complete.");
        return 0;
    }

    // Root of the external mzLib checkout, from the MZLIB_DIR environment variable
    // (default F:\mzLib). This port references mzLib but does not vendor it -- see README.
    static string MzLibDir()
    {
        string v = Environment.GetEnvironmentVariable("MZLIB_DIR");
        if (string.IsNullOrEmpty(v)) v = @"F:\mzLib";
        if (!Directory.Exists(Path.Combine(v, "mzLib", "Test", "FlashLFQ", "TestData")))
        {
            throw new DirectoryNotFoundException(
                $"MZLIB_DIR='{v}' does not contain mzLib\\Test\\FlashLFQ\\TestData. " +
                "Set the MZLIB_DIR environment variable to your mzLib checkout.");
        }
        return v;
    }

    // Locate this port's flashlfq-core directory (which holds parity/golden) by walking up
    // from the compiled exe (parity/csharp_golden/bin/<cfg>/net8.0/).
    static string FindFlashLfqCore()
    {
        string d = AppContext.BaseDirectory;
        for (int i = 0; i < 12 && d != null; i++)
        {
            if (Directory.Exists(Path.Combine(d, "parity", "golden")) &&
                Directory.Exists(Path.Combine(d, "src")))
            {
                return d;
            }
            d = Path.GetDirectoryName(d);
        }
        throw new DirectoryNotFoundException("Could not locate flashlfq-core (expected parity/golden + src).");
    }

    static void WriteL5(FlashLfqResults results, string path)
    {
        var rows = new List<string[]>();
        foreach (var kvp in results.Peaks)
        {
            string fileName = kvp.Key.FilenameWithoutExtension;
            foreach (ChromatographicPeak peak in kvp.Value)
            {
                var id = peak.Identifications.First();
                FlashLFQ.IsotopicEnvelope apex = peak.Apex;

                string baseSeqs = string.Join("|", peak.Identifications.Select(x => x.BaseSequence).Distinct());
                string modSeqs = string.Join("|", peak.Identifications.Select(x => x.ModifiedSequence).Distinct());

                double rtStart = apex != null ? peak.IsotopicEnvelopes.Min(e => e.IndexedPeak.RetentionTime) : double.NaN;
                double rtApex = apex != null ? apex.IndexedPeak.RetentionTime : double.NaN;
                double rtEnd = apex != null ? peak.IsotopicEnvelopes.Max(e => e.IndexedPeak.RetentionTime) : double.NaN;
                double peakMz = apex != null ? apex.IndexedPeak.M : double.NaN;
                int apexCharge = apex != null ? apex.ChargeState : -1;
                int apexScanIndex = apex != null ? apex.IndexedPeak.ZeroBasedScanIndex : -1;
                double apexPearson = apex != null ? apex.PearsonCorrelation : double.NaN;

                rows.Add(new[]
                {
                    fileName,
                    baseSeqs,
                    modSeqs,
                    D(id.MonoisotopicMass),
                    D(id.Ms2RetentionTimeInMinutes),
                    id.PrecursorChargeState.ToString(Inv),
                    D(peak.Intensity),
                    D(rtStart),
                    D(rtApex),
                    D(rtEnd),
                    D(peakMz),
                    apexCharge.ToString(Inv),
                    apexScanIndex.ToString(Inv),
                    D(apexPearson),
                    peak.NumChargeStatesObserved.ToString(Inv),
                    D(peak.MassError),
                    D(peak.SplitRT),
                    peak.NumIdentificationsByFullSeq.ToString(Inv),
                    peak.DetectionType.ToString(),
                    peak.DecoyPeptide.ToString(Inv),
                });
            }
        }

        // Deterministic order: file, modified sequence, apex RT, intensity.
        rows = rows
            .OrderBy(r => r[0], StringComparer.Ordinal)
            .ThenBy(r => r[2], StringComparer.Ordinal)
            .ThenBy(r => ParseSortable(r[8]))
            .ThenBy(r => ParseSortable(r[6]))
            .ToList();

        string[] header =
        {
            "file_name", "base_sequence", "modified_sequence", "monoisotopic_mass",
            "ms2_retention_time", "precursor_charge", "peak_intensity", "peak_rt_start",
            "peak_rt_apex", "peak_rt_end", "peak_mz", "apex_charge", "apex_scan_index",
            "apex_pearson", "num_charge_states", "mass_error", "split_rt",
            "num_ids_by_full_seq", "detection_type", "decoy",
        };

        var sb = new StringBuilder();
        sb.Append(string.Join("\t", header)).Append('\n');
        foreach (var r in rows) sb.Append(string.Join("\t", r)).Append('\n');
        File.WriteAllText(path, sb.ToString());
        Console.WriteLine($"Wrote {rows.Count} peak rows -> {path}");
    }

    static void WriteL6(FlashLfqResults results, string path)
    {
        var files = results.SpectraFiles.OrderBy(f => f.FilenameWithoutExtension, StringComparer.Ordinal).ToList();
        var rows = new List<string[]>();

        foreach (var kvp in results.PeptideModifiedSequences)
        {
            string modSeq = kvp.Key;
            Peptide pep = kvp.Value;
            foreach (var file in files)
            {
                rows.Add(new[]
                {
                    modSeq,
                    file.FilenameWithoutExtension,
                    D(pep.GetIntensity(file)),
                    D(pep.GetRetentionTime(file)),
                    pep.GetDetectionType(file).ToString(),
                });
            }
        }

        rows = rows
            .OrderBy(r => r[0], StringComparer.Ordinal)
            .ThenBy(r => r[1], StringComparer.Ordinal)
            .ToList();

        string[] header = { "modified_sequence", "file_name", "intensity", "retention_time", "detection_type" };
        var sb = new StringBuilder();
        sb.Append(string.Join("\t", header)).Append('\n');
        foreach (var r in rows) sb.Append(string.Join("\t", r)).Append('\n');
        File.WriteAllText(path, sb.ToString());
        Console.WriteLine($"Wrote {rows.Count} peptide x file rows -> {path}");
    }

    // Default DonorQValueThreshold (FlashLfqParameters.cs:36).
    const double DonorQValueThreshold = 0.01;

    // Faithful replica of FlashLfqEngine.ChooseBestPeak for the default DonorCriterion.Score:
    // the peak whose highest-scoring identification is greatest (MaxBy = first on ties); if that
    // peak's first id has a positive score it stands, else fall through to most-intense.
    static ChromatographicPeak ChooseBestPeakScore(List<ChromatographicPeak> peaks)
    {
        ChromatographicPeak best = peaks.MaxBy(p => p.Identifications.Max(id => id.PsmScore));
        if (best.Identifications.First().PsmScore > 0)
            return best;
        return peaks.MaxBy(p => p.Intensity);
    }

    // Replica of FlashLfqEngine.GetRtCalSpline (FlashLfqEngine.cs:559) over the real engine peaks.
    static void WriteMbrSpline(FlashLfqResults results, SpectraFileInfo donor, SpectraFileInfo acceptor, string path)
    {
        Func<SpectraFileInfo, Dictionary<string, ChromatographicPeak>> bestBySeq = file =>
        {
            var bySeq = results.Peaks[file]
                .Where(peak => peak.NumIdentificationsByFullSeq == 1
                    && peak.DetectionType == DetectionType.MSMS
                    && peak.IsotopicEnvelopes.Any()
                    && peak.Identifications.Min(id => id.QValue) < DonorQValueThreshold)
                .GroupBy(peak => peak.Identifications.First().ModifiedSequence)
                .ToDictionary(g => g.Key, g => g.ToList());
            var best = new Dictionary<string, ChromatographicPeak>();
            foreach (var kvp in bySeq)
            {
                if (!kvp.Value.Any()) continue;
                var bp = ChooseBestPeakScore(kvp.Value);
                if (bp == null) continue;
                best[kvp.Key] = bp;
            }
            return best;
        };

        var donorBest = bestBySeq(donor);
        var acceptorBest = bestBySeq(acceptor);

        var curve = new List<(string seq, double donorRt, double acceptorRt, double rtDiff, double donorMass)>();
        foreach (var kvp in acceptorBest)
        {
            if (donorBest.TryGetValue(kvp.Key, out ChromatographicPeak donorPeak))
            {
                double donorRt = donorPeak.Apex.IndexedPeak.RetentionTime;
                double acceptorRt = kvp.Value.Apex.IndexedPeak.RetentionTime;
                double rtDiff = acceptorRt - donorRt; // RetentionTimeCalibDataPoint.RtDiff
                double donorMass = donorPeak.Identifications.First().PeakfindingMass;
                curve.Add((kvp.Key, donorRt, acceptorRt, rtDiff, donorMass));
            }
        }

        // OrderBy donor apex RT (the C# return sort).
        curve = curve.OrderBy(c => c.donorRt).ToList();

        string[] header = { "donor_modified_sequence", "donor_apex_rt", "acceptor_apex_rt", "rt_diff", "donor_peakfinding_mass" };
        var sb = new StringBuilder();
        sb.Append(string.Join("\t", header)).Append('\n');
        foreach (var c in curve)
        {
            sb.Append(string.Join("\t", new[] { c.seq, D(c.donorRt), D(c.acceptorRt), D(c.rtDiff), D(c.donorMass) })).Append('\n');
        }
        File.WriteAllText(path, sb.ToString());
        Console.WriteLine($"Wrote {curve.Count} RT-cal-spline points -> {path}");
    }

    // Builds the per-sequence best MSMS peak dictionary for one file (the GetRtCalSpline filter +
    // ChooseBestPeak Score grouping). Shared by the spline + predicted-RT golden writers.
    static Dictionary<string, ChromatographicPeak> BestBySeq(FlashLfqResults results, SpectraFileInfo file)
    {
        var bySeq = results.Peaks[file]
            .Where(peak => peak.NumIdentificationsByFullSeq == 1
                && peak.DetectionType == DetectionType.MSMS
                && peak.IsotopicEnvelopes.Any()
                && peak.Identifications.Min(id => id.QValue) < DonorQValueThreshold)
            .GroupBy(peak => peak.Identifications.First().ModifiedSequence)
            .ToDictionary(g => g.Key, g => g.ToList());
        var best = new Dictionary<string, ChromatographicPeak>();
        foreach (var kvp in bySeq)
        {
            if (!kvp.Value.Any()) continue;
            var bp = ChooseBestPeakScore(kvp.Value);
            if (bp == null) continue;
            best[kvp.Key] = bp;
        }
        return best;
    }

    // Number of anchor peptides used per side for local RT alignment (FlashLfqEngine.cs:35).
    const int NumberOfAnchorPeptidesForMbr = 3;
    // FlashLfqParameters.MaxMbrRtWindow default (FlashLfqParameters.cs:32).
    const double MaxMbrRtWindow = 1.0;

    // Faithful replica of FlashLfqEngine.PredictRetentionTime (FlashLfqEngine.cs:716), minus the
    // fraction gate (both K562 files are unfractionated). Operates over real
    // RetentionTimeCalibDataPoint[] objects so Array.BinarySearch uses the real CompareTo.
    static RtInfo PredictRtReplica(RetentionTimeCalibDataPoint[] curve, ChromatographicPeak donorPeak)
    {
        double donorRt = donorPeak.Apex.IndexedPeak.RetentionTime;

        var testPoint = new RetentionTimeCalibDataPoint(donorPeak, null);
        int index = Array.BinarySearch(curve, testPoint);
        if (index < 0) index = ~index;
        if (index >= curve.Length && index >= 1) index = curve.Length - 1;

        var nearby = new List<RetentionTimeCalibDataPoint>();

        int forward = 0;
        for (int r = index + 1; r < curve.Length; r++)
        {
            double rtDiff = curve[r].DonorFilePeak.Apex.IndexedPeak.RetentionTime - donorRt;
            if (curve[r].AcceptorFilePeak != null && curve[r].AcceptorFilePeak.ApexRetentionTime > 0)
            {
                if (Math.Abs(rtDiff) > 0.5) break;
                nearby.Add(curve[r]);
                forward++;
                if (forward >= NumberOfAnchorPeptidesForMbr) break;
            }
        }

        int backward = 0;
        for (int r = index - 1; r >= 0; r--)
        {
            double rtDiff = curve[r].DonorFilePeak.Apex.IndexedPeak.RetentionTime - donorRt;
            if (curve[r].AcceptorFilePeak != null && curve[r].AcceptorFilePeak.ApexRetentionTime > 0)
            {
                if (Math.Abs(rtDiff) > 0.5) break;
                nearby.Add(curve[r]);
                backward++;
                if (backward >= NumberOfAnchorPeptidesForMbr) break;
            }
        }

        if (!nearby.Any())
            return new RtInfo(donorRt, 0.25);

        var rtDiffs = nearby
            .Select(p => p.DonorFilePeak.ApexRetentionTime - p.AcceptorFilePeak.ApexRetentionTime)
            .ToList();

        double medianRtDiff = rtDiffs.Median();
        if (rtDiffs.Count == 1)
            return new RtInfo(donorRt - medianRtDiff, 0.25);

        double rtRange = rtDiffs.StandardDeviation() * 6;
        rtRange = Math.Min(rtRange, MaxMbrRtWindow);
        return new RtInfo(donorRt - medianRtDiff, rtRange);
    }

    static void WriteMbrPredictedRt(FlashLfqResults results, SpectraFileInfo donor, SpectraFileInfo acceptor, string path)
    {
        var donorBest = BestBySeq(results, donor);
        var acceptorBest = BestBySeq(results, acceptor);

        // Build the same RT calibration curve as GetRtCalSpline (paired shared sequences, sorted by donor apex RT).
        var curveList = new List<RetentionTimeCalibDataPoint>();
        foreach (var kvp in acceptorBest)
        {
            if (donorBest.TryGetValue(kvp.Key, out ChromatographicPeak donorPeak))
            {
                curveList.Add(new RetentionTimeCalibDataPoint(donorPeak, kvp.Value));
            }
        }
        var curve = curveList.OrderBy(p => p.DonorFilePeak.Apex.IndexedPeak.RetentionTime).ToArray();

        // Predict the acceptor RT for every donor best-peak (ordered by peakfinding mass, the
        // donorPeaksMassOrdered set), one row per donor modified sequence.
        var rows = new List<(string seq, double predictedRt, double width)>();
        foreach (var donorPeak in donorBest.Values.OrderBy(p => p.Identifications.First().PeakfindingMass))
        {
            RtInfo info = PredictRtReplica(curve, donorPeak);
            rows.Add((donorPeak.Identifications.First().ModifiedSequence, info.PredictedRt, info.Width));
        }

        rows = rows.OrderBy(r => r.seq, StringComparer.Ordinal).ToList();

        string[] header = { "donor_modified_sequence", "predicted_rt", "width" };
        var sb = new StringBuilder();
        sb.Append(string.Join("\t", header)).Append('\n');
        foreach (var r in rows)
        {
            sb.Append(string.Join("\t", new[] { r.seq, D(r.predictedRt), D(r.width) })).Append('\n');
        }
        File.WriteAllText(path, sb.ToString());
        Console.WriteLine($"Wrote {rows.Count} predicted-RT rows -> {path}");
    }

    static double ParseSortable(string s)
    {
        if (double.TryParse(s, NumberStyles.Float, Inv, out double v)) return double.IsNaN(v) ? double.MaxValue : v;
        return double.MaxValue;
    }
}
