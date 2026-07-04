//! psmtsv parsing (PLAN.md P1.9): MetaMorpheus `*.psmtsv` -> [`Identification`].
//!
//! This is the Rust port of the slice of mzLib's `SpectrumMatchTsvReader` / `PsmFromTsv` /
//! `SpectrumMatchFromTsv` (`Readers/InternalResults/...`) plus `FlashLFQ`'s
//! `MzLibExtensions.MakeIdentifications` that FlashLFQ actually needs to seed quantification:
//! turning each PSM row into an [`Identification`] carrying base sequence, full (modified)
//! sequence, monoisotopic mass, precursor charge, MS2 retention time, and source file, along
//! with the score / q-value / decoy flag the later FDR + protein-quant stages read.
//!
//! ## Faithfulness to the C# reader
//!
//! * **Tabular model.** mzLib splits each line on a single `'\t'` (`SpectrumMatchTsvReader.Split`)
//!   and trims `"` then whitespace off every cell — it does *not* use RFC-4180 CSV quoting (a
//!   cell like `[y1+1:175.1]` would break real CSV quote parsing). We mirror this by driving the
//!   `csv` crate with `quoting(false)`, `delimiter(b'\t')`, and `flexible(true)`, then trimming
//!   `"`/whitespace ourselves in [`cell`].
//!
//! * **Header -> column resolution.** mzLib's `ParseHeader` maps logical fields to column
//!   indices by header-text match, choosing between two layouts: when an exact `Accession`
//!   column exists it reads `Monoisotopic Mass`; otherwise (the FlashLFQ psmtsv layout) it reads
//!   `Peptide Monoisotopic Mass`. [`HeaderMap::resolve`] reproduces that branch so the same files
//!   the C# reader accepts resolve the same column.
//!
//! * **Per-field semantics**, matching `SpectrumMatchFromTsv`'s constructor:
//!   - `FileName` -> [`Identification::file_name`]: known spectra-file extensions stripped, then
//!     the bare file name (`FileNameWithoutExtension`). For this corpus the values carry no
//!     extension, so it is a pass-through; the strip logic is kept for fidelity.
//!   - `Base Sequence` -> base sequence with SILAC-style `(...)` parentheses removed
//!     (`RemoveParentheses`).
//!   - `Full Sequence` -> modified sequence verbatim.
//!   - monoisotopic mass: first `|`-split token parsed as `f64`, else `-1.0` (mirrors the C#
//!     `MonoisotopicMass` getter's `TryParse(... .Split('|')[0]) ? value : -1`).
//!   - `Precursor Charge` parsed as `f64` then truncated to `i32` (C# `(int)GetRequiredValue<double>`).
//!   - `Scan Retention Time` -> minutes, defaulting to `-1.0` when absent/blank.
//!   - `Score`, `QValue` parsed as `f64`; decoy flag = `Decoy/Contaminant/Target` contains `'D'`.
//!
//! Required-but-missing columns or unparseable required numbers make the whole parse fail
//! (`PsmTsvError`); mzLib instead drops the offending row with a warning, but for a parity gate
//! we want a hard error rather than silently fewer records.

use std::collections::HashMap;
use std::fmt;
use std::path::Path;

use serde::Serialize;

/// Spectra-file extensions mzLib strips off `File Name` before taking the bare file name.
/// Verbatim from `SpectrumMatchFromTsvHeader.AcceptedSpectraFormats`.
const ACCEPTED_SPECTRA_FORMATS: &[&str] =
    &[".raw", ".mzML", ".mgf", ".d", "ms1.msalign", "ms2.msalign"];

/// A single FlashLFQ identification, ported from `FlashLFQ.Identification`.
///
/// Only the fields FlashLFQ reads from a psmtsv are populated here; protein-group resolution
/// (the cross-record `HashSet<ProteinGroup>` dedup in `MakeIdentifications`) is left to a later
/// task, as it is a property of the whole id set rather than one row.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Identification {
    /// Source spectra file name, extension stripped (`FileInfo` key in C#).
    pub file_name: String,
    /// Unmodified amino-acid sequence (`BaseSequence`), SILAC `(...)` labels removed.
    pub base_sequence: String,
    /// Full / modified sequence (`ModifiedSequence`), verbatim from the `Full Sequence` column.
    pub modified_sequence: String,
    /// Monoisotopic mass in daltons (`MonoisotopicMass`); `-1.0` when the cell is unparseable.
    pub monoisotopic_mass: f64,
    /// MS2 retention time in minutes (`Ms2RetentionTimeInMinutes`); `-1.0` when absent/blank.
    pub ms2_retention_time_in_minutes: f64,
    /// Precursor charge state (`PrecursorChargeState`).
    pub precursor_charge_state: i32,
    /// PSM score (`PsmScore`).
    pub score: f64,
    /// PSM q-value (`QValue`).
    pub q_value: f64,
    /// Decoy flag (`IsDecoy`): the `Decoy/Contaminant/Target` cell contains `'D'`.
    pub is_decoy: bool,
}

/// Errors raised while parsing a psmtsv into [`Identification`]s.
#[derive(Debug)]
pub enum PsmTsvError {
    /// The file could not be read or the `csv` reader failed on a record.
    Io(String),
    /// The header was missing entirely (empty file).
    MissingHeader,
    /// A required logical column was not present under any of its accepted header names.
    MissingColumn(&'static str),
    /// A required numeric cell could not be parsed.
    BadNumber {
        column: &'static str,
        value: String,
        row: usize,
    },
    /// A data row had fewer columns than a resolved index demanded.
    ShortRow { row: usize, needed: usize, got: usize },
}

impl fmt::Display for PsmTsvError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PsmTsvError::Io(e) => write!(f, "psmtsv read error: {e}"),
            PsmTsvError::MissingHeader => write!(f, "psmtsv has no header line"),
            PsmTsvError::MissingColumn(c) => write!(f, "psmtsv missing required column {c:?}"),
            PsmTsvError::BadNumber { column, value, row } => {
                write!(f, "psmtsv row {row}: cannot parse {column:?} value {value:?} as a number")
            }
            PsmTsvError::ShortRow { row, needed, got } => write!(
                f,
                "psmtsv row {row}: needs at least {needed} columns but has {got}"
            ),
        }
    }
}

impl std::error::Error for PsmTsvError {}

/// Maps logical fields to resolved column indices, mirroring mzLib's `ParseHeader`.
struct HeaderMap {
    file_name: usize,
    scan_retention_time: Option<usize>,
    precursor_charge: usize,
    base_sequence: usize,
    full_sequence: usize,
    monoisotopic_mass: usize,
    score: usize,
    q_value: usize,
    decoy_contam_target: usize,
}

impl HeaderMap {
    /// Resolves required/optional columns from the header cells.
    ///
    /// The monoisotopic-mass column follows the C# branch: an exact `Accession` column selects
    /// the `Monoisotopic Mass` layout; otherwise (the FlashLFQ psmtsv layout) `Peptide
    /// Monoisotopic Mass` is used.
    fn resolve(header: &[String]) -> Result<HeaderMap, PsmTsvError> {
        let index_of = |name: &str| header.iter().position(|h| h == name);
        let required = |name: &'static str| index_of(name).ok_or(PsmTsvError::MissingColumn(name));

        let mono_col = if index_of("Accession").is_some() {
            required("Monoisotopic Mass")?
        } else {
            required("Peptide Monoisotopic Mass")?
        };

        Ok(HeaderMap {
            file_name: required("File Name")?,
            scan_retention_time: index_of("Scan Retention Time"),
            precursor_charge: required("Precursor Charge")?,
            base_sequence: required("Base Sequence")?,
            full_sequence: required("Full Sequence")?,
            monoisotopic_mass: mono_col,
            score: required("Score")?,
            q_value: required("QValue")?,
            decoy_contam_target: required("Decoy/Contaminant/Target")?,
        })
    }
}

/// Trims `"` quote chars then surrounding whitespace off a cell, matching mzLib's
/// `p.Trim('"')` followed by `.Trim()`.
fn cell(fields: &csv::StringRecord, idx: usize) -> &str {
    fields.get(idx).unwrap_or("").trim_matches('"').trim()
}

/// Strips a known spectra-file extension and returns the bare file name
/// (`FileNameWithoutExtension` in `SpectrumMatchFromTsv`).
///
/// Public so a caller can key spectra-file paths by the *same* bare name the psmtsv `File Name`
/// column resolves to (e.g. matching `--raw foo.mzML` to ids whose `file_name` is `foo`).
pub fn file_name_without_extension(raw: &str) -> String {
    if !raw.contains('.') {
        return raw.to_string();
    }
    let mut s = raw.to_string();
    for ext in ACCEPTED_SPECTRA_FORMATS {
        // Case-insensitive replace of the extension token (C# uses InvariantCultureIgnoreCase).
        s = replace_ignore_case(&s, ext);
    }
    // C# then takes Path.GetFileName; take the last path component on either separator.
    s.rsplit(['/', '\\']).next().unwrap_or(&s).to_string()
}

/// Case-insensitive removal of every occurrence of `needle` from `haystack`.
fn replace_ignore_case(haystack: &str, needle: &str) -> String {
    if needle.is_empty() {
        return haystack.to_string();
    }
    let hay_lower = haystack.to_ascii_lowercase();
    let needle_lower = needle.to_ascii_lowercase();
    let mut out = String::with_capacity(haystack.len());
    let mut i = 0;
    while i < haystack.len() {
        if hay_lower[i..].starts_with(&needle_lower) {
            i += needle.len();
        } else {
            // advance one full char to stay on a UTF-8 boundary
            let ch = haystack[i..].chars().next().unwrap();
            out.push(ch);
            i += ch.len_utf8();
        }
    }
    out
}

/// Removes SILAC-style `(...)` parenthetical labels from a base sequence, keeping residues
/// outside parentheses (`SpectrumMatchFromTsv.RemoveParentheses`).
pub fn remove_parentheses(base_sequence: &str) -> String {
    if !base_sequence.contains('(') {
        return base_sequence.to_string();
    }
    let mut out = String::with_capacity(base_sequence.len());
    let mut within = false;
    for c in base_sequence.chars() {
        match c {
            ')' => within = false,
            '(' => within = true,
            _ if !within => out.push(c),
            _ => {}
        }
    }
    out
}

/// Parses a `Precursor Charge` cell the C# way: `(int)double.Parse(...)` (truncate toward zero).
fn parse_charge(value: &str, row: usize) -> Result<i32, PsmTsvError> {
    let f: f64 = value.parse().map_err(|_| PsmTsvError::BadNumber {
        column: "Precursor Charge",
        value: value.to_string(),
        row,
    })?;
    Ok(f as i32)
}

/// Parses the monoisotopic-mass cell: first `|`-split token as `f64`, else `-1.0`
/// (`MonoisotopicMass` getter).
fn parse_monoisotopic_mass(value: &str) -> f64 {
    let first = value.split('|').next().unwrap_or("").trim();
    first.parse::<f64>().unwrap_or(-1.0)
}

/// Parses a psmtsv file at `path` into a list of [`Identification`]s, preserving row order.
pub fn read_identifications<P: AsRef<Path>>(path: P) -> Result<Vec<Identification>, PsmTsvError> {
    let text = std::fs::read_to_string(path).map_err(|e| PsmTsvError::Io(e.to_string()))?;
    read_identifications_from_str(&text)
}

/// Parses psmtsv `text` into a list of [`Identification`]s, preserving row order.
pub fn read_identifications_from_str(text: &str) -> Result<Vec<Identification>, PsmTsvError> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .quoting(false) // mzLib splits on raw tabs and trims quotes itself; no CSV quoting.
        .flexible(true) // rows may have trailing-empty / ragged column counts.
        .has_headers(true)
        .from_reader(text.as_bytes());

    let header: Vec<String> = reader
        .headers()
        .map_err(|e| PsmTsvError::Io(e.to_string()))?
        .iter()
        .map(|h| h.to_string())
        .collect();
    if header.is_empty() {
        return Err(PsmTsvError::MissingHeader);
    }
    let map = HeaderMap::resolve(&header)?;

    // Largest index any required column needs, to give a precise short-row error.
    let needed = [
        map.file_name,
        map.precursor_charge,
        map.base_sequence,
        map.full_sequence,
        map.monoisotopic_mass,
        map.score,
        map.q_value,
        map.decoy_contam_target,
    ]
    .into_iter()
    .max()
    .unwrap_or(0);

    let mut out = Vec::new();
    for (i, record) in reader.records().enumerate() {
        let row = i + 2; // 1-based, +1 for the header line, matching mzLib's warning numbering
        let record = record.map_err(|e| PsmTsvError::Io(e.to_string()))?;
        if record.iter().all(|c| c.trim().is_empty()) {
            continue; // skip blank lines
        }
        if record.len() <= needed {
            return Err(PsmTsvError::ShortRow {
                row,
                needed: needed + 1,
                got: record.len(),
            });
        }

        let score = cell(&record, map.score)
            .parse::<f64>()
            .map_err(|_| PsmTsvError::BadNumber {
                column: "Score",
                value: cell(&record, map.score).to_string(),
                row,
            })?;
        let q_value = cell(&record, map.q_value)
            .parse::<f64>()
            .map_err(|_| PsmTsvError::BadNumber {
                column: "QValue",
                value: cell(&record, map.q_value).to_string(),
                row,
            })?;
        let ms2_rt = match map.scan_retention_time {
            Some(idx) => {
                let v = cell(&record, idx);
                if v.is_empty() {
                    -1.0
                } else {
                    v.parse::<f64>().unwrap_or(-1.0)
                }
            }
            None => -1.0,
        };

        out.push(Identification {
            file_name: file_name_without_extension(cell(&record, map.file_name)),
            base_sequence: remove_parentheses(cell(&record, map.base_sequence)),
            modified_sequence: cell(&record, map.full_sequence).to_string(),
            monoisotopic_mass: parse_monoisotopic_mass(cell(&record, map.monoisotopic_mass)),
            ms2_retention_time_in_minutes: ms2_rt,
            precursor_charge_state: parse_charge(cell(&record, map.precursor_charge), row)?,
            score,
            q_value,
            is_decoy: cell(&record, map.decoy_contam_target).contains('D'),
        });
    }

    Ok(out)
}

/// Convenience: counts distinct full (modified) sequences across a list of identifications.
pub fn distinct_modified_sequences(ids: &[Identification]) -> usize {
    let set: HashMap<&str, ()> = ids
        .iter()
        .map(|id| (id.modified_sequence.as_str(), ()))
        .collect();
    set.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn all_psms_path() -> PathBuf {
        crate::mzlib_test_data("AllPSMs.psmtsv")
    }

    #[test]
    fn remove_parentheses_strips_silac_labels() {
        assert_eq!(remove_parentheses("PEPTIDE"), "PEPTIDE");
        assert_eq!(remove_parentheses("PEP(label)TIDE"), "PEPTIDE");
        assert_eq!(remove_parentheses("(a)PEP(bc)"), "PEP");
    }

    #[test]
    fn file_name_strips_known_extensions() {
        assert_eq!(file_name_without_extension("run1"), "run1");
        assert_eq!(file_name_without_extension("run1.mzML"), "run1");
        assert_eq!(file_name_without_extension("run1.RAW"), "run1");
        assert_eq!(file_name_without_extension("C:\\data\\run1.mzML"), "run1");
    }

    #[test]
    fn monoisotopic_mass_takes_first_pipe_token_or_minus_one() {
        assert_eq!(parse_monoisotopic_mass("1914.82537"), 1914.82537);
        assert_eq!(parse_monoisotopic_mass("1914.82537|2000.0"), 1914.82537);
        assert_eq!(parse_monoisotopic_mass(""), -1.0);
        assert_eq!(parse_monoisotopic_mass("not-a-number"), -1.0);
    }

    #[test]
    fn charge_truncates_double_to_int() {
        assert_eq!(parse_charge("3.00000", 1).unwrap(), 3);
        assert_eq!(parse_charge("2", 1).unwrap(), 2);
        assert!(parse_charge("abc", 1).is_err());
    }

    #[test]
    fn small_inline_psmtsv_parses() {
        // Minimal psmtsv with the FlashLFQ layout (no `Accession`, has `Peptide Monoisotopic Mass`).
        let text = "File Name\tScan Retention Time\tPrecursor Charge\tBase Sequence\tFull Sequence\tPeptide Monoisotopic Mass\tScore\tDecoy/Contaminant/Target\tQValue\n\
                    run1.mzML\t79.12\t3.00000\tPEPTIDE\tPEP[mod]TIDE\t799.35996\t23.287\tT\t0.01\n\
                    run2\t81.25\t2.00000\tAA(silac)CC\tAACC\t1000.0|2000.0\t10.0\tD\t0.0\n";
        let ids = read_identifications_from_str(text).expect("parses");
        assert_eq!(ids.len(), 2);

        let a = &ids[0];
        assert_eq!(a.file_name, "run1"); // extension stripped
        assert_eq!(a.base_sequence, "PEPTIDE");
        assert_eq!(a.modified_sequence, "PEP[mod]TIDE");
        assert_eq!(a.monoisotopic_mass, 799.35996);
        assert_eq!(a.ms2_retention_time_in_minutes, 79.12);
        assert_eq!(a.precursor_charge_state, 3);
        assert_eq!(a.score, 23.287);
        assert_eq!(a.q_value, 0.01);
        assert!(!a.is_decoy);

        let b = &ids[1];
        assert_eq!(b.file_name, "run2");
        assert_eq!(b.base_sequence, "AACC"); // (silac) removed
        assert_eq!(b.monoisotopic_mass, 1000.0); // first |-token
        assert_eq!(b.precursor_charge_state, 2);
        assert!(b.is_decoy); // DCT contains 'D'
    }

    #[test]
    fn missing_required_column_errors() {
        let text = "File Name\tBase Sequence\n run1\tPEPTIDE\n";
        assert!(matches!(
            read_identifications_from_str(text),
            Err(PsmTsvError::MissingColumn(_))
        ));
    }

    #[test]
    fn parses_all_psms_corpus_record_count_and_distinct_sequences() {
        // Parity numbers established in PLAN.md: AllPSMs.psmtsv has 594 data rows and
        // 354 distinct Full Sequences; this corpus is all targets.
        let ids = read_identifications(all_psms_path()).expect("AllPSMs.psmtsv parses");
        assert_eq!(ids.len(), 594, "expected 594 PSM rows");
        assert_eq!(
            distinct_modified_sequences(&ids),
            354,
            "expected 354 distinct Full Sequences"
        );
        assert!(ids.iter().all(|id| !id.is_decoy), "corpus is all targets");

        // Spot-check the first row against the known psmtsv values.
        let first = &ids[0];
        assert_eq!(first.file_name, "20100614_Velos1_TaGe_SA_K562_3");
        assert_eq!(first.base_sequence, "AHQLVMEGYNWCHDR");
        assert_eq!(
            first.modified_sequence,
            "AHQLVMEGYNWC[Common Fixed:Carbamidomethyl on C]HDR"
        );
        assert_eq!(first.precursor_charge_state, 3);
        assert!((first.monoisotopic_mass - 1914.82537).abs() < 1e-9);
        assert!((first.ms2_retention_time_in_minutes - 79.12191).abs() < 1e-9);
        assert!((first.score - 23.287).abs() < 1e-9);
    }
}
