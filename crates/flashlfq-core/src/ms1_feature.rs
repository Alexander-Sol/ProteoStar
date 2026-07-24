//! Reading and writing TopFD / FLASHDeconv **`_ms1.feature`** files — the tab-delimited
//! MS1 feature table that top-down tools exchange (TopFD and FLASHDeconv write it; TopPIC
//! and mzLib read it).
//!
//! Unlike [`crate::feature_export`]'s `.msalign` (a *per-scan deconvoluted spectrum*
//! format), `_ms1.feature` is genuinely feature-shaped: one row per proteoform feature,
//! carrying a neutral mass, an RT range, a charge-state range, and intensities. That makes
//! it the natural interchange format for ProteoStar's resolved features, and the one the
//! viewer can load from other tools' output.
//!
//! # The three dialects
//!
//! There is no single spec; three column vocabularies are in the wild, and this module
//! reads all of them (see [`Ms1FeatureDialect`]). Column aliasing follows mzLib's
//! `Readers/ExternalResults/IndividualResultRecords/Ms1Feature.cs`, so anything mzLib
//! accepts is accepted here.
//!
//! ```text
//! TopFD v1.6.2   Sample_ID ID Mass Intensity Time_begin Time_end Apex_time
//!                Apex_intensity Minimum_charge_state Maximum_charge_state
//!                Minimum_fraction_id Maximum_fraction_id
//!
//! FLASHDeconv    same as above, but Time_apex instead of Apex_time and no
//!                Apex_intensity column at all
//!
//! TopFD v1.7.0   File_name Fraction_ID Feature_ID Mass Intensity Min_time Max_time
//!                Min_scan Max_scan Min_charge Max_charge Apex_time Apex_scan
//!                Apex_intensity Rep_charge Rep_average_mz Envelope_num EC_score
//! ```
//!
//! # Retention-time units — the one real trap
//!
//! The dialects disagree, and nothing in the file declares which is in use:
//!
//! - **TopFD v1.6.2 and FLASHDeconv write seconds.** In the reference files a feature spans
//!   `Time_begin 2375.98 → Time_end 2398.21`; read as seconds that is a 22-second elution
//!   peak (right for top-down LC), read as minutes it would be 22 minutes (absurd for one
//!   feature).
//! - **TopFD v1.7.0 writes minutes.** There, `Min_time 44.458 → Max_time 45.280` spans
//!   `Min_scan 1539 → Max_scan 1593`; 54 scans across 0.82 minutes is ~1.1 scans/s, which is
//!   right. As seconds it would be 54 scans in 0.82 s.
//!
//! ProteoStar works in **minutes** throughout, so [`read_ms1_feature`] normalises to minutes
//! on the way in and records what it found in [`Ms1FeatureSet::rt_unit`];
//! [`write_ms1_feature`] converts back to the target dialect's native unit on the way out.
//! Note mzLib does *not* normalise — it reads the raw doubles — so a file round-tripped
//! through here is not byte-identical to one round-tripped through mzLib unless the dialect
//! is preserved (which [`Ms1FeatureSet::write`] does).
//!
//! Unit selection is by dialect, with a magnitude guard for mislabelled files: a
//! "minutes" dialect whose largest RT exceeds [`MAX_PLAUSIBLE_RUN_MINUTES`] is re-read as
//! seconds, and a "seconds" dialect whose largest RT falls below
//! [`MIN_PLAUSIBLE_RUN_SECONDS`] is re-read as minutes. Those windows are far outside any
//! real LC run, so the guard cannot fire on well-formed input.

use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::feature_refinement::ResolvedFeature;

/// An RT larger than this in a minutes-dialect file means the file is really in seconds.
/// A 12-hour LC gradient does not exist; 720 s (12 min) runs are routine.
pub const MAX_PLAUSIBLE_RUN_MINUTES: f64 = 720.0;

/// An RT smaller than this in a seconds-dialect file means the file is really in minutes.
/// A 12-second acquisition does not exist; 12-minute runs are routine.
pub const MIN_PLAUSIBLE_RUN_SECONDS: f64 = 12.0;

/// Which column vocabulary a file uses. Detected from the header on read; chosen by the
/// caller on write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ms1FeatureDialect {
    /// TopFD ≤ 1.6.x: `Time_begin`/`Time_end`/`Apex_time` + `Apex_intensity`, RT in seconds.
    TopFdV1_6,
    /// TopFD ≥ 1.7.0: the wide schema (`File_name`, `Feature_ID`, scan bounds, `EC_score`),
    /// RT in minutes.
    TopFdV1_7,
    /// FLASHDeconv (OpenMS): the v1.6 column set minus `Apex_intensity`, RT in seconds.
    FlashDeconv,
}

impl Ms1FeatureDialect {
    /// The retention-time unit this dialect natively writes.
    pub fn native_rt_unit(self) -> RtUnit {
        match self {
            Ms1FeatureDialect::TopFdV1_7 => RtUnit::Minutes,
            Ms1FeatureDialect::TopFdV1_6 | Ms1FeatureDialect::FlashDeconv => RtUnit::Seconds,
        }
    }
}

impl fmt::Display for Ms1FeatureDialect {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Ms1FeatureDialect::TopFdV1_6 => "TopFD v1.6",
            Ms1FeatureDialect::TopFdV1_7 => "TopFD v1.7",
            Ms1FeatureDialect::FlashDeconv => "FLASHDeconv",
        };
        f.write_str(s)
    }
}

/// The retention-time unit found in (or written to) a file. Records always hold minutes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RtUnit {
    Minutes,
    Seconds,
}

impl RtUnit {
    /// Converts a value in this unit to minutes.
    fn to_minutes(self, v: f64) -> f64 {
        match self {
            RtUnit::Minutes => v,
            RtUnit::Seconds => v / 60.0,
        }
    }

    /// Converts a value in minutes to this unit.
    fn from_minutes(self, v: f64) -> f64 {
        match self {
            RtUnit::Minutes => v,
            RtUnit::Seconds => v * 60.0,
        }
    }
}

/// One `_ms1.feature` row: a proteoform feature as a neutral mass with RT and charge extents.
///
/// **All retention times are minutes**, regardless of the source dialect's unit — the reader
/// normalises and the writer denormalises. Fields absent from the source dialect are `None`
/// (they are omitted, not defaulted, so a round trip does not fabricate values).
#[derive(Debug, Clone, PartialEq)]
pub struct Ms1FeatureRecord {
    /// Source spectra file (TopFD v1.7 `File_name`).
    pub file_name: Option<String>,
    /// Sample index (`Sample_ID`); absent in TopFD v1.7.
    pub sample_id: Option<i32>,
    /// Fraction index. Written to both `Minimum_fraction_id`/`Maximum_fraction_id` in the
    /// narrow dialects and to the single `Fraction_ID` column in TopFD v1.7.
    pub fraction_id: Option<i32>,
    /// Row identifier (`ID` / `Feature_ID`).
    pub id: i32,
    /// Neutral monoisotopic mass in daltons.
    pub mass: f64,
    /// Summed feature intensity across the charge states this row covers.
    pub intensity: f64,
    /// Start of the elution range, **minutes**.
    pub rt_begin: f64,
    /// End of the elution range, **minutes**.
    pub rt_end: f64,
    /// Apex retention time, **minutes**.
    pub rt_apex: f64,
    /// Intensity at the apex. `None` for FLASHDeconv, which has no such column.
    pub apex_intensity: Option<f64>,
    /// Lowest charge state covered by this row (inclusive).
    pub charge_min: i32,
    /// Highest charge state covered by this row (inclusive).
    pub charge_max: i32,
    /// First scan of the elution range (TopFD v1.7 only).
    pub scan_min: Option<i32>,
    /// Last scan of the elution range (TopFD v1.7 only).
    pub scan_max: Option<i32>,
    /// Apex scan number (TopFD v1.7 only).
    pub apex_scan: Option<i32>,
    /// The most intense / representative charge state (TopFD v1.7 only).
    pub rep_charge: Option<i32>,
    /// Average m/z of the representative charge's envelope (TopFD v1.7 only).
    pub rep_average_mz: Option<f64>,
    /// Number of isotopic envelopes contributing to the feature (TopFD v1.7 only).
    pub envelope_num: Option<i32>,
    /// TopFD's envelope-collection score in [0, 1] (TopFD v1.7 only). Deliberately *not*
    /// mapped onto ProteoStar's deconvolution score — it is a different quantity on a
    /// different scale.
    pub ec_score: Option<f64>,
}

impl Ms1FeatureRecord {
    /// The charge states this row covers, ascending. The format stores only the bounds, so a
    /// row is by definition a contiguous run — which is why [`resolved_to_ms1_feature_records`]
    /// splits gapped charge sets across rows rather than emitting a span it cannot represent.
    pub fn charge_states(&self) -> Vec<i32> {
        (self.charge_min.min(self.charge_max)..=self.charge_max.max(self.charge_min)).collect()
    }
}

/// A parsed `_ms1.feature` file: the rows plus the provenance needed to write it back out
/// unchanged.
#[derive(Debug, Clone, PartialEq)]
pub struct Ms1FeatureSet {
    /// The column vocabulary detected from the header.
    pub dialect: Ms1FeatureDialect,
    /// The retention-time unit the *source file* used. Records are always in minutes.
    pub rt_unit: RtUnit,
    /// The feature rows, in file order.
    pub records: Vec<Ms1FeatureRecord>,
}

impl Ms1FeatureSet {
    /// Writes this set back out in the dialect it was read as (and therefore its RT unit),
    /// so a read/write cycle is faithful.
    pub fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        write_ms1_feature(writer, &self.records, self.dialect)
    }
}

/// Why an `_ms1.feature` file could not be read.
#[derive(Debug)]
pub enum Ms1FeatureError {
    /// The file could not be opened or read.
    Io(io::Error),
    /// The file was read but is not a well-formed `_ms1.feature` table.
    Parse(String),
}

impl fmt::Display for Ms1FeatureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Ms1FeatureError::Io(e) => write!(f, "{e}"),
            Ms1FeatureError::Parse(m) => f.write_str(m),
        }
    }
}

impl std::error::Error for Ms1FeatureError {}

impl From<io::Error> for Ms1FeatureError {
    fn from(e: io::Error) -> Self {
        Ms1FeatureError::Io(e)
    }
}

// ------------------------------------------------------------------------ header

/// Header column names, in the order the narrow (v1.6 / FLASHDeconv) dialects use them.
/// The first entry of each alias list is what the writer emits.
const A_SAMPLE_ID: &[&str] = &["Sample_ID"];
const A_ID: &[&str] = &["ID", "Feature_ID"];
const A_MASS: &[&str] = &["Mass"];
const A_INTENSITY: &[&str] = &["Intensity"];
const A_RT_BEGIN: &[&str] = &["Time_begin", "Min_time"];
const A_RT_END: &[&str] = &["Time_end", "Max_time"];
const A_RT_APEX: &[&str] = &["Apex_time", "Time_apex"];
const A_APEX_INTENSITY: &[&str] = &["Apex_intensity", "Intensity_Apex"];
const A_CHARGE_MIN: &[&str] = &["Minimum_charge_state", "Min_charge"];
const A_CHARGE_MAX: &[&str] = &["Maximum_charge_state", "Max_charge"];
/// A record holds one fraction id. The narrow dialects split it across a min/max pair that is
/// always the same value in practice, so the reader takes the min column and the writer emits
/// it to both.
const A_FRACTION: &[&str] = &["Minimum_fraction_id", "Fraction_ID"];
const A_FILE_NAME: &[&str] = &["File_name"];
const A_SCAN_MIN: &[&str] = &["Min_scan"];
const A_SCAN_MAX: &[&str] = &["Max_scan"];
const A_APEX_SCAN: &[&str] = &["Apex_scan"];
const A_REP_CHARGE: &[&str] = &["Rep_charge"];
const A_REP_AVG_MZ: &[&str] = &["Rep_average_mz"];
const A_ENVELOPE_NUM: &[&str] = &["Envelope_num"];
const A_EC_SCORE: &[&str] = &["EC_score"];

/// Header name → column index. Trailing empty names (TopFD v1.7 ends its header with a tab)
/// are dropped so they cannot shadow a real column.
struct Header<'a>(HashMap<&'a str, usize>);

impl<'a> Header<'a> {
    fn parse(line: &'a str) -> Self {
        Header(
            line.split('\t')
                .enumerate()
                .filter_map(|(i, c)| {
                    let c = c.trim();
                    (!c.is_empty()).then_some((c, i))
                })
                .collect(),
        )
    }

    /// Index of the first alias present, if any.
    fn find(&self, aliases: &[&str]) -> Option<usize> {
        aliases.iter().find_map(|a| self.0.get(a).copied())
    }

    fn has(&self, aliases: &[&str]) -> bool {
        self.find(aliases).is_some()
    }

    fn require(&self, aliases: &[&str]) -> Result<usize, Ms1FeatureError> {
        self.find(aliases).ok_or_else(|| {
            Ms1FeatureError::Parse(format!("missing column '{}'", aliases[0]))
        })
    }
}

/// True if `header` looks like an `_ms1.feature` table. Used to tell this format apart from
/// ProteoStar's own resolved-feature TSV without trusting the file extension.
///
/// The test is the format's irreducible core: a neutral `Mass`, a charge range, and an
/// elution start — under any dialect's spelling. ProteoStar's resolved TSV has none of these
/// column names (it uses `Monoisotopic Mass`, `Charge States`, `RT Start`), so the two never
/// collide.
pub fn is_ms1_feature_header(header: &str) -> bool {
    let h = Header::parse(header);
    h.has(A_MASS) && h.has(A_CHARGE_MIN) && h.has(A_CHARGE_MAX) && h.has(A_RT_BEGIN)
}

/// Picks the dialect from the columns present.
fn detect_dialect(h: &Header<'_>) -> Ms1FeatureDialect {
    // The wide schema is unmistakable: only v1.7 carries scan bounds and a source file name.
    if h.has(A_FILE_NAME) || h.has(A_SCAN_MIN) || h.has(A_EC_SCORE) {
        Ms1FeatureDialect::TopFdV1_7
    } else if h.has(A_APEX_INTENSITY) {
        // mzLib uses this same discriminator: an apex-intensity column means TopFD.
        Ms1FeatureDialect::TopFdV1_6
    } else {
        Ms1FeatureDialect::FlashDeconv
    }
}

/// Resolves the file's RT unit: the dialect's native unit, overridden when the observed
/// magnitudes make that reading impossible. See the module docs for the thresholds.
fn infer_rt_unit(dialect: Ms1FeatureDialect, max_rt: f64) -> RtUnit {
    match dialect.native_rt_unit() {
        RtUnit::Minutes if max_rt > MAX_PLAUSIBLE_RUN_MINUTES => RtUnit::Seconds,
        RtUnit::Seconds if max_rt > 0.0 && max_rt < MIN_PLAUSIBLE_RUN_SECONDS => RtUnit::Minutes,
        native => native,
    }
}

// ------------------------------------------------------------------------- read

/// Parses an `_ms1.feature` document, normalising retention times to minutes.
///
/// Blank lines are skipped. Rows with fewer columns than the header, or whose required
/// numeric cells do not parse, are rejected with a [`Ms1FeatureError::Parse`] naming the line
/// — a silently-zeroed mass or charge would surface much later as a mysteriously missing
/// feature.
pub fn read_ms1_feature(text: &str) -> Result<Ms1FeatureSet, Ms1FeatureError> {
    let mut lines = text.lines().enumerate().skip_while(|(_, l)| l.trim().is_empty());
    let (_, header_line) = lines
        .next()
        .ok_or_else(|| Ms1FeatureError::Parse("empty ms1.feature file".into()))?;
    let h = Header::parse(header_line);
    let dialect = detect_dialect(&h);

    let c_mass = h.require(A_MASS)?;
    let c_intensity = h.require(A_INTENSITY)?;
    let c_begin = h.require(A_RT_BEGIN)?;
    let c_end = h.require(A_RT_END)?;
    let c_apex = h.require(A_RT_APEX)?;
    let c_zmin = h.require(A_CHARGE_MIN)?;
    let c_zmax = h.require(A_CHARGE_MAX)?;

    let c_id = h.find(A_ID);
    let c_sample = h.find(A_SAMPLE_ID);
    let c_fraction = h.find(A_FRACTION);
    let c_apex_int = h.find(A_APEX_INTENSITY);
    let c_file = h.find(A_FILE_NAME);
    let c_scan_min = h.find(A_SCAN_MIN);
    let c_scan_max = h.find(A_SCAN_MAX);
    let c_apex_scan = h.find(A_APEX_SCAN);
    let c_rep_z = h.find(A_REP_CHARGE);
    let c_rep_mz = h.find(A_REP_AVG_MZ);
    let c_env_num = h.find(A_ENVELOPE_NUM);
    let c_ec = h.find(A_EC_SCORE);

    // Pass 1: parse every row with raw (un-normalised) times, tracking the largest so the
    // unit can be inferred from the file as a whole rather than row by row.
    let mut raw: Vec<Ms1FeatureRecord> = Vec::new();
    let mut max_rt = 0.0_f64;
    for (lineno, line) in lines {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        // 1-based line number for the message, matching what an editor shows.
        let at = lineno + 1;

        let req_f = |i: usize, name: &str| -> Result<f64, Ms1FeatureError> {
            f.get(i)
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .and_then(|s| s.parse::<f64>().ok())
                .ok_or_else(|| {
                    Ms1FeatureError::Parse(format!("line {at}: bad or missing '{name}'"))
                })
        };
        let req_i = |i: usize, name: &str| -> Result<i32, Ms1FeatureError> {
            f.get(i)
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .and_then(|s| s.parse::<i32>().ok())
                .ok_or_else(|| {
                    Ms1FeatureError::Parse(format!("line {at}: bad or missing '{name}'"))
                })
        };
        let opt_f = |i: Option<usize>| -> Option<f64> {
            i.and_then(|i| f.get(i))
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .and_then(|s| s.parse::<f64>().ok())
        };
        let opt_i = |i: Option<usize>| -> Option<i32> {
            i.and_then(|i| f.get(i))
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .and_then(|s| s.parse::<i32>().ok())
        };

        let rt_begin = req_f(c_begin, "Time_begin")?;
        let rt_end = req_f(c_end, "Time_end")?;
        let rt_apex = req_f(c_apex, "Apex_time")?;
        max_rt = max_rt.max(rt_begin).max(rt_end).max(rt_apex);

        raw.push(Ms1FeatureRecord {
            file_name: c_file
                .and_then(|i| f.get(i))
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(str::to_owned),
            sample_id: opt_i(c_sample),
            fraction_id: opt_i(c_fraction),
            id: opt_i(c_id).unwrap_or(raw.len() as i32),
            mass: req_f(c_mass, "Mass")?,
            intensity: req_f(c_intensity, "Intensity")?,
            rt_begin,
            rt_end,
            rt_apex,
            apex_intensity: opt_f(c_apex_int),
            charge_min: req_i(c_zmin, "Minimum_charge_state")?,
            charge_max: req_i(c_zmax, "Maximum_charge_state")?,
            scan_min: opt_i(c_scan_min),
            scan_max: opt_i(c_scan_max),
            apex_scan: opt_i(c_apex_scan),
            rep_charge: opt_i(c_rep_z),
            rep_average_mz: opt_f(c_rep_mz),
            envelope_num: opt_i(c_env_num),
            ec_score: opt_f(c_ec),
        });
    }

    // Pass 2: now that the file's RT range is known, convert to minutes.
    let rt_unit = infer_rt_unit(dialect, max_rt);
    for r in &mut raw {
        r.rt_begin = rt_unit.to_minutes(r.rt_begin);
        r.rt_end = rt_unit.to_minutes(r.rt_end);
        r.rt_apex = rt_unit.to_minutes(r.rt_apex);
    }

    Ok(Ms1FeatureSet { dialect, rt_unit, records: raw })
}

/// Reads an `_ms1.feature` file from disk. See [`read_ms1_feature`].
pub fn read_ms1_feature_file<P: AsRef<Path>>(path: P) -> Result<Ms1FeatureSet, Ms1FeatureError> {
    read_ms1_feature(&std::fs::read_to_string(path)?)
}

// ------------------------------------------------------------------------ write

/// Writes `records` as an `_ms1.feature` document in the given dialect.
///
/// Retention times are converted from the records' minutes to the dialect's native unit
/// (seconds for TopFD v1.6 / FLASHDeconv, minutes for TopFD v1.7). Columns the dialect does
/// not define are omitted; a column it defines but the record lacks is written empty, which
/// is what mzLib's `double?`/`int?` fields expect.
pub fn write_ms1_feature<W: Write>(
    writer: &mut W,
    records: &[Ms1FeatureRecord],
    dialect: Ms1FeatureDialect,
) -> io::Result<()> {
    let unit = dialect.native_rt_unit();
    // Optional numeric cell → its text, or empty when absent.
    fn cell<T: fmt::Display>(v: Option<T>) -> String {
        v.map(|x| x.to_string()).unwrap_or_default()
    }

    match dialect {
        Ms1FeatureDialect::TopFdV1_7 => {
            writeln!(
                writer,
                "File_name\tFraction_ID\tFeature_ID\tMass\tIntensity\tMin_time\tMax_time\t\
                 Min_scan\tMax_scan\tMin_charge\tMax_charge\tApex_time\tApex_scan\t\
                 Apex_intensity\tRep_charge\tRep_average_mz\tEnvelope_num\tEC_score"
            )?;
            for r in records {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    r.file_name.as_deref().unwrap_or(""),
                    cell(r.fraction_id),
                    r.id,
                    r.mass,
                    r.intensity,
                    unit.from_minutes(r.rt_begin),
                    unit.from_minutes(r.rt_end),
                    cell(r.scan_min),
                    cell(r.scan_max),
                    r.charge_min,
                    r.charge_max,
                    unit.from_minutes(r.rt_apex),
                    cell(r.apex_scan),
                    cell(r.apex_intensity),
                    cell(r.rep_charge),
                    cell(r.rep_average_mz),
                    cell(r.envelope_num),
                    cell(r.ec_score),
                )?;
            }
        }
        Ms1FeatureDialect::TopFdV1_6 | Ms1FeatureDialect::FlashDeconv => {
            // The two narrow dialects differ only by the Apex_intensity column.
            let apex_col = dialect == Ms1FeatureDialect::TopFdV1_6;
            writeln!(
                writer,
                "Sample_ID\tID\tMass\tIntensity\tTime_begin\tTime_end\tApex_time\t{}\
                 Minimum_charge_state\tMaximum_charge_state\tMinimum_fraction_id\t\
                 Maximum_fraction_id",
                if apex_col { "Apex_intensity\t" } else { "" }
            )?;
            for r in records {
                write!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
                    cell(r.sample_id.or(Some(0))),
                    r.id,
                    r.mass,
                    r.intensity,
                    unit.from_minutes(r.rt_begin),
                    unit.from_minutes(r.rt_end),
                    unit.from_minutes(r.rt_apex),
                )?;
                if apex_col {
                    write!(writer, "{}\t", cell(r.apex_intensity))?;
                }
                let frac = r.fraction_id.unwrap_or(0);
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    r.charge_min, r.charge_max, frac, frac
                )?;
            }
        }
    }
    Ok(())
}

/// Convenience wrapper: writes an `_ms1.feature` file at `path`. See [`write_ms1_feature`].
pub fn write_ms1_feature_file<P: AsRef<Path>>(
    path: P,
    records: &[Ms1FeatureRecord],
    dialect: Ms1FeatureDialect,
) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_ms1_feature(&mut writer, records, dialect)?;
    writer.flush()
}

// -------------------------------------------------------------------- conversion

/// Maps the pipeline's resolved features onto `_ms1.feature` rows.
///
/// **A gapped charge set becomes several rows.** The format stores only `charge_min` and
/// `charge_max`, so a feature seen at charges {10, 12, 15} written as a single 10–15 row
/// would be read back as 10, 11, 12, 13, 14, 15 — fabricating three charges that were never
/// observed. Splitting into maximal contiguous runs (10–10, 12–12, 15–15) carries the exact
/// set. A feature whose charges are already contiguous yields exactly one row. This mirrors
/// mzLib's `ToMs1Features`.
///
/// Where mzLib repeats the whole feature's summed intensity on every split row, each row here
/// carries only the intensity of the members in *its* run, so summing the rows reproduces the
/// feature's total instead of multiplying it.
///
/// `file_name` is recorded on each row (the spectra file's basename); it is only emitted by
/// the TopFD v1.7 writer.
pub fn resolved_to_ms1_feature_records(
    resolved: &[ResolvedFeature],
    file_name: &str,
) -> Vec<Ms1FeatureRecord> {
    let mut out: Vec<Ms1FeatureRecord> = Vec::new();
    for feature in resolved {
        for (lo, hi) in contiguous_charge_runs(&feature.charge_states) {
            // Members belonging to this run — the ones whose intensity and apex the row reports.
            let members: Vec<_> = feature
                .members
                .iter()
                .filter(|m| m.refined_charge >= lo && m.refined_charge <= hi)
                .collect();
            // A run with no members can only arise if charge_states and members disagree;
            // skip rather than emit a zero-intensity row.
            let Some(dominant) = members.iter().max_by(|a, b| {
                a.detected
                    .summed_intensity
                    .total_cmp(&b.detected.summed_intensity)
            }) else {
                continue;
            };
            let intensity: f64 = members.iter().map(|m| m.detected.summed_intensity).sum();
            // Apex intensity = the dominant member's signal at its apex scan, i.e. the height
            // of the strongest charge state's envelope at the top of its elution. Matches
            // mzLib's "apex envelope of the dominant trace".
            let apex_intensity: f64 = dominant
                .detected
                .peaks
                .iter()
                .filter(|p| p.zero_based_scan_index == dominant.detected.apex_scan_index)
                .map(|p| p.intensity as f64)
                .sum();

            out.push(Ms1FeatureRecord {
                file_name: Some(file_name.to_owned()),
                sample_id: Some(0),
                fraction_id: Some(0),
                id: out.len() as i32,
                mass: feature.monoisotopic_mass,
                intensity,
                // The row's RT extent is the run's, not the whole feature's, so a split row
                // describes only the charges it actually covers.
                rt_begin: members
                    .iter()
                    .map(|m| m.detected.start_rt)
                    .fold(f64::INFINITY, f64::min),
                rt_end: members
                    .iter()
                    .map(|m| m.detected.end_rt)
                    .fold(f64::NEG_INFINITY, f64::max),
                rt_apex: dominant.detected.apex_rt,
                apex_intensity: Some(apex_intensity),
                charge_min: lo,
                charge_max: hi,
                // Detector scan indices are zero-based; the format's scan numbers are one-based.
                scan_min: None,
                scan_max: None,
                apex_scan: Some(dominant.detected.apex_scan_index.max(0) + 1),
                rep_charge: Some(dominant.refined_charge),
                rep_average_mz: Some(dominant.detected.mono_mz),
                envelope_num: Some(members.len() as i32),
                // EC_score is TopFD's own envelope-collection metric; we have no equivalent
                // on the same scale, so leave it unset rather than write a look-alike.
                ec_score: None,
            });
        }
    }
    out
}

/// Distinct charges, ascending, grouped into maximal contiguous `(lo, hi)` runs.
/// `[10, 12, 13, 15]` → `[(10,10), (12,13), (15,15)]`. Empty input yields no runs.
///
/// Public because every producer of `_ms1.feature` rows needs it: the format cannot express
/// a gapped charge set in one row, so a gapped set must be split or it will be read back with
/// the gaps filled in.
pub fn contiguous_charge_runs(charges: &[i32]) -> Vec<(i32, i32)> {
    let mut sorted: Vec<i32> = charges.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let mut runs = Vec::new();
    let mut iter = sorted.into_iter();
    let Some(first) = iter.next() else {
        return runs;
    };
    let (mut start, mut prev) = (first, first);
    for z in iter {
        if z == prev + 1 {
            prev = z;
        } else {
            runs.push((start, prev));
            start = z;
            prev = z;
        }
    }
    runs.push((start, prev));
    runs
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The header TopFD v1.6.2 actually writes (from a reference file).
    const V16_HEADER: &str = "Sample_ID\tID\tMass\tIntensity\tTime_begin\tTime_end\tApex_time\tApex_intensity\tMinimum_charge_state\tMaximum_charge_state\tMinimum_fraction_id\tMaximum_fraction_id";

    fn v16_file() -> String {
        format!(
            "{V16_HEADER}\n\
             0\t0\t10835.85272090354\t10849947123.04\t2372.27\t2401.92\t2390.8\t912795138.8\t7\t17\t0\t0\n\
             0\t1\t11299.38479102015\t21715373919.44001\t2387.1\t2974.11\t2416.74\t395601212.66\t7\t18\t0\t0\n"
        )
    }

    /// FLASHDeconv (OpenMS 3.0.0): no Apex_intensity, and it spells the apex column Time_apex.
    fn flashdeconv_file() -> String {
        "Sample_ID\tID\tMass\tIntensity\tTime_begin\tTime_end\tTime_apex\tMinimum_charge_state\tMaximum_charge_state\tMinimum_fraction_id\tMaximum_fraction_id\n\
         0\t1\t10835.9\t1.21E+10\t2375.98\t2398.21\t2390.8\t7\t18\t0\t0\n\
         0\t2\t13997.9\t7.99E+09\t2546.47\t2561.3\t2561.3\t9\t23\t0\t0\n"
            .to_string()
    }

    /// TopFD v1.7.0: the wide schema, minutes, and a trailing tab on every line.
    fn v17_file() -> String {
        "File_name\tFraction_ID\tFeature_ID\tMass\tIntensity\tMin_time\tMax_time\tMin_scan\tMax_scan\tMin_charge\tMax_charge\tApex_time\tApex_scan\tApex_intensity\tRep_charge\tRep_average_mz\tEnvelope_num\tEC_score\t\n\
         E:/Projects/yeast.mzML\t0\t0\t6983.005647734753\t29168462662.22477\t44.45833333333334\t45.27966666666667\t1539\t1593\t6\t14\t44.67983333333333\t1557\t2503126332.811494\t10\t699.7515683472939\t129\t0.9979024819216775\t\n"
            .to_string()
    }

    #[test]
    fn reads_topfd_v16_and_converts_seconds_to_minutes() {
        let set = read_ms1_feature(&v16_file()).unwrap();
        assert_eq!(set.dialect, Ms1FeatureDialect::TopFdV1_6);
        assert_eq!(set.rt_unit, RtUnit::Seconds);
        assert_eq!(set.records.len(), 2);

        let r = &set.records[0];
        assert!((r.mass - 10835.85272090354).abs() < 1e-9);
        assert!((r.intensity - 10849947123.04).abs() < 1e-3);
        // 2372.27 s -> minutes.
        assert!((r.rt_begin - 2372.27 / 60.0).abs() < 1e-9);
        assert!((r.rt_apex - 2390.8 / 60.0).abs() < 1e-9);
        assert_eq!(r.apex_intensity, Some(912795138.8));
        assert_eq!((r.charge_min, r.charge_max), (7, 17));
        assert_eq!(r.charge_states().len(), 11);
        // Narrow dialect: no wide-schema fields invented.
        assert_eq!(r.scan_min, None);
        assert_eq!(r.ec_score, None);
        assert_eq!(r.file_name, None);
    }

    #[test]
    fn reads_flashdeconv_without_apex_intensity() {
        let set = read_ms1_feature(&flashdeconv_file()).unwrap();
        assert_eq!(set.dialect, Ms1FeatureDialect::FlashDeconv);
        assert_eq!(set.rt_unit, RtUnit::Seconds);
        // The Time_apex spelling resolves through the alias list.
        assert!((set.records[0].rt_apex - 2390.8 / 60.0).abs() < 1e-9);
        // Scientific notation parses.
        assert!((set.records[0].intensity - 1.21e10).abs() < 1.0);
        // No apex-intensity column anywhere in this dialect.
        assert!(set.records.iter().all(|r| r.apex_intensity.is_none()));
    }

    #[test]
    fn reads_topfd_v17_wide_schema_in_minutes() {
        let set = read_ms1_feature(&v17_file()).unwrap();
        assert_eq!(set.dialect, Ms1FeatureDialect::TopFdV1_7);
        // Already minutes — no conversion applied.
        assert_eq!(set.rt_unit, RtUnit::Minutes);
        let r = &set.records[0];
        assert!((r.rt_begin - 44.45833333333334).abs() < 1e-12);
        assert_eq!(r.file_name.as_deref(), Some("E:/Projects/yeast.mzML"));
        assert_eq!(r.id, 0); // read from Feature_ID
        assert_eq!(r.scan_min, Some(1539));
        assert_eq!(r.apex_scan, Some(1557));
        assert_eq!(r.rep_charge, Some(10));
        assert_eq!(r.envelope_num, Some(129));
        assert!((r.ec_score.unwrap() - 0.9979024819216775).abs() < 1e-12);
        // The trailing tab must not become a phantom column or shift any index.
        assert_eq!((r.charge_min, r.charge_max), (6, 14));
    }

    #[test]
    fn rt_unit_guard_overrides_a_mislabelled_dialect() {
        // A v1.7-shaped header (minutes by default) carrying values that can only be seconds.
        let text = "File_name\tFeature_ID\tMass\tIntensity\tMin_time\tMax_time\tApex_time\tMin_charge\tMax_charge\n\
                    a.mzML\t0\t1000.0\t500.0\t2400.0\t2430.0\t2415.0\t5\t9\n";
        let set = read_ms1_feature(text).unwrap();
        assert_eq!(set.rt_unit, RtUnit::Seconds, "2400 'minutes' is a 40-hour run");
        assert!((set.records[0].rt_apex - 2415.0 / 60.0).abs() < 1e-9);

        // ...and the converse: a seconds-dialect file whose values can only be minutes.
        let text = "Sample_ID\tID\tMass\tIntensity\tTime_begin\tTime_end\tApex_time\tApex_intensity\tMinimum_charge_state\tMaximum_charge_state\n\
                    0\t0\t1000.0\t500.0\t1.5\t2.5\t2.0\t100.0\t5\t9\n";
        let set = read_ms1_feature(text).unwrap();
        assert_eq!(set.rt_unit, RtUnit::Minutes, "a 2.5-second run is impossible");
        assert!((set.records[0].rt_apex - 2.0).abs() < 1e-12);
    }

    #[test]
    fn round_trips_through_each_dialect() {
        for (name, text) in [
            ("v1.6", v16_file()),
            ("flashdeconv", flashdeconv_file()),
            ("v1.7", v17_file()),
        ] {
            let set = read_ms1_feature(&text).unwrap();
            let mut buf: Vec<u8> = Vec::new();
            set.write(&mut buf).unwrap();
            let reparsed = read_ms1_feature(&String::from_utf8(buf).unwrap()).unwrap();

            assert_eq!(reparsed.dialect, set.dialect, "{name}: dialect not preserved");
            assert_eq!(reparsed.rt_unit, set.rt_unit, "{name}: rt unit not preserved");
            assert_eq!(reparsed.records.len(), set.records.len(), "{name}");
            for (a, b) in set.records.iter().zip(&reparsed.records) {
                assert!((a.mass - b.mass).abs() < 1e-9, "{name}: mass");
                assert!((a.rt_begin - b.rt_begin).abs() < 1e-9, "{name}: rt_begin");
                assert!((a.rt_end - b.rt_end).abs() < 1e-9, "{name}: rt_end");
                assert!((a.rt_apex - b.rt_apex).abs() < 1e-9, "{name}: rt_apex");
                assert_eq!(a.charge_min, b.charge_min, "{name}: charge_min");
                assert_eq!(a.charge_max, b.charge_max, "{name}: charge_max");
                assert_eq!(a.apex_intensity, b.apex_intensity, "{name}: apex intensity");
            }
        }
    }

    #[test]
    fn header_sniffer_accepts_every_dialect_and_rejects_the_resolved_tsv() {
        for text in [v16_file(), flashdeconv_file(), v17_file()] {
            assert!(is_ms1_feature_header(text.lines().next().unwrap()));
        }
        // ProteoStar's own resolved-feature TSV must not be mistaken for this format.
        let resolved = "Detected m/z (primary)\tRT Start\tRT Apex\tRT End\tCharge States\tPer-Charge Detected m/z\tPrimary Charge\tMonoisotopic Mass\tMono m/z (primary)\tSummed Intensity\tCross-Charge Support\tNum Members";
        assert!(!is_ms1_feature_header(resolved));
    }

    #[test]
    fn rejects_malformed_input_instead_of_zeroing_it() {
        // No Mass column at all.
        let e = read_ms1_feature("Foo\tBar\n1\t2\n").unwrap_err();
        assert!(matches!(e, Ms1FeatureError::Parse(ref m) if m.contains("Mass")));

        // Header is fine, but a row's mass is not a number: must fail loudly, and name the line.
        let text = format!("{V16_HEADER}\n0\t0\tnot_a_mass\t1.0\t10\t20\t15\t1.0\t2\t3\t0\t0\n");
        let e = read_ms1_feature(&text).unwrap_err();
        assert!(matches!(e, Ms1FeatureError::Parse(ref m) if m.contains("line 2") && m.contains("Mass")));

        assert!(read_ms1_feature("").is_err());
    }

    #[test]
    fn header_only_file_reads_as_zero_records() {
        let set = read_ms1_feature(&format!("{V16_HEADER}\n")).unwrap();
        assert_eq!(set.dialect, Ms1FeatureDialect::TopFdV1_6);
        assert!(set.records.is_empty());
    }

    #[test]
    fn charge_runs_split_on_gaps_only() {
        assert_eq!(contiguous_charge_runs(&[]), vec![]);
        assert_eq!(contiguous_charge_runs(&[5]), vec![(5, 5)]);
        assert_eq!(contiguous_charge_runs(&[2, 3, 4]), vec![(2, 4)]);
        assert_eq!(
            contiguous_charge_runs(&[10, 12, 13, 15]),
            vec![(10, 10), (12, 13), (15, 15)]
        );
        // Unsorted and duplicated input normalises to the same runs.
        assert_eq!(contiguous_charge_runs(&[13, 10, 12, 15, 13]), vec![
            (10, 10),
            (12, 13),
            (15, 15)
        ]);
    }

    #[test]
    fn gapped_charges_never_fabricate_a_charge_on_round_trip() {
        // The failure this guards against: {10, 12} written as one 10-12 row reads back as
        // {10, 11, 12}. Split rows carry the exact set.
        let records = vec![
            Ms1FeatureRecord {
                file_name: None,
                sample_id: Some(0),
                fraction_id: Some(0),
                id: 0,
                mass: 5000.0,
                intensity: 100.0,
                rt_begin: 10.0,
                rt_end: 11.0,
                rt_apex: 10.5,
                apex_intensity: Some(50.0),
                charge_min: 10,
                charge_max: 10,
                scan_min: None,
                scan_max: None,
                apex_scan: None,
                rep_charge: None,
                rep_average_mz: None,
                envelope_num: None,
                ec_score: None,
            },
            Ms1FeatureRecord {
                id: 1,
                charge_min: 12,
                charge_max: 12,
                ..Default::default()
            },
        ];
        let mut buf: Vec<u8> = Vec::new();
        write_ms1_feature(&mut buf, &records, Ms1FeatureDialect::TopFdV1_6).unwrap();
        let back = read_ms1_feature(&String::from_utf8(buf).unwrap()).unwrap();
        let charges: Vec<i32> = back.records.iter().flat_map(|r| r.charge_states()).collect();
        assert_eq!(charges, vec![10, 12], "charge 11 must not appear");
    }

    #[test]
    fn writer_emits_empty_cells_for_absent_optional_values() {
        let records = vec![Ms1FeatureRecord { apex_intensity: None, ..Default::default() }];
        let mut buf: Vec<u8> = Vec::new();
        write_ms1_feature(&mut buf, &records, Ms1FeatureDialect::TopFdV1_6).unwrap();
        let out = String::from_utf8(buf).unwrap();
        let header: Vec<&str> = out.lines().next().unwrap().split('\t').collect();
        let row: Vec<&str> = out.lines().nth(1).unwrap().split('\t').collect();
        // The column is present but blank — what mzLib's `double?` reads as null — and the
        // row still has exactly as many fields as the header.
        assert_eq!(row.len(), header.len());
        let i = header.iter().position(|c| *c == "Apex_intensity").unwrap();
        assert_eq!(row[i], "");
        // And it reads back as absent rather than zero.
        assert_eq!(read_ms1_feature(&out).unwrap().records[0].apex_intensity, None);
    }

    #[test]
    fn flashdeconv_writer_omits_the_apex_intensity_column() {
        let records = vec![Ms1FeatureRecord { apex_intensity: Some(1.0), ..Default::default() }];
        let mut buf: Vec<u8> = Vec::new();
        write_ms1_feature(&mut buf, &records, Ms1FeatureDialect::FlashDeconv).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(!out.contains("Apex_intensity"));
        // Losing the column is what makes it FLASHDeconv on re-read.
        assert_eq!(read_ms1_feature(&out).unwrap().dialect, Ms1FeatureDialect::FlashDeconv);
    }

    #[test]
    fn narrow_writer_emits_seconds() {
        let records = vec![Ms1FeatureRecord {
            rt_begin: 10.0,
            rt_end: 12.0,
            rt_apex: 11.0,
            ..Default::default()
        }];
        let mut buf: Vec<u8> = Vec::new();
        write_ms1_feature(&mut buf, &records, Ms1FeatureDialect::TopFdV1_6).unwrap();
        let out = String::from_utf8(buf).unwrap();
        // 11 minutes -> 660 seconds on the wire...
        assert!(out.contains("660"), "apex should be written in seconds: {out}");
        // ...and back to 11 minutes in the record.
        assert!((read_ms1_feature(&out).unwrap().records[0].rt_apex - 11.0).abs() < 1e-9);
    }
}

#[cfg(test)]
impl Default for Ms1FeatureRecord {
    /// A minimal valid record, for tests that only care about a field or two.
    fn default() -> Self {
        Ms1FeatureRecord {
            file_name: None,
            sample_id: Some(0),
            fraction_id: Some(0),
            id: 0,
            mass: 1000.0,
            intensity: 1.0,
            rt_begin: 10.0,
            rt_end: 11.0,
            rt_apex: 10.5,
            apex_intensity: None,
            charge_min: 1,
            charge_max: 1,
            scan_min: None,
            scan_max: None,
            apex_scan: None,
            rep_charge: None,
            rep_average_mz: None,
            envelope_num: None,
            ec_score: None,
        }
    }
}
