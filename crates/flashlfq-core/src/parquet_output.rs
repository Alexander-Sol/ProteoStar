//! Parquet output for the peptide-level results table (PLAN.md P1.16).
//!
//! Serializes [`PeptideResults`] — the (peptide × file) quantification table produced by
//! [`crate::results::calculate_peptide_results`] — to a Parquet file via Arrow.
//!
//! ## Layout: tidy / long form
//!
//! One row per `(modified_sequence, file_name)` cell, with columns:
//!
//! | column              | type    | source                                 |
//! |---------------------|---------|----------------------------------------|
//! | `modified_sequence` | Utf8    | the peptide's full (modified) sequence |
//! | `file_name`         | Utf8    | the spectra file the cell belongs to   |
//! | `intensity`         | Float64 | [`PeptideQuant::intensity`]            |
//! | `retention_time`    | Float64 | [`PeptideQuant::retention_time`]       |
//! | `detection_type`    | Utf8    | [`DetectionType::as_str`]              |
//!
//! Long form is the simplest faithful serialization of the nested `modseq -> file -> quant`
//! map and round-trips cleanly into polars/pandas (`pivot` to recover the wide peptide×file
//! matrix). Rows are emitted in a deterministic order (sorted by modified sequence, then file
//! name) so the output is reproducible across runs regardless of `HashMap` iteration order.
//!
//! [`PeptideQuant::intensity`]: crate::results::PeptideQuant::intensity
//! [`PeptideQuant::retention_time`]: crate::results::PeptideQuant::retention_time
//! [`DetectionType::as_str`]: crate::detection_type::DetectionType::as_str

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{BooleanArray, Float64Array, StringArray, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::errors::ParquetError;

use crate::mbr_search::FeatureRow;
use crate::results::PeptideResults;

/// The Arrow schema of the tidy peptide-results table. All columns are non-nullable.
pub fn peptide_results_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("modified_sequence", DataType::Utf8, false),
        Field::new("file_name", DataType::Utf8, false),
        Field::new("intensity", DataType::Float64, false),
        Field::new("retention_time", DataType::Float64, false),
        Field::new("detection_type", DataType::Utf8, false),
    ]))
}

/// Builds a tidy (long-form) Arrow [`RecordBatch`] from [`PeptideResults`]: one row per
/// `(modified_sequence, file_name)` cell. Rows are sorted by modified sequence then file name
/// for deterministic output.
pub fn peptide_results_record_batch(results: &PeptideResults) -> RecordBatch {
    // Flatten the nested map into a sorted row list so the output is reproducible.
    let mut rows: Vec<(&str, &str)> = Vec::new();
    for (modseq, per_file) in &results.quant {
        for file in per_file.keys() {
            rows.push((modseq.as_str(), file.as_str()));
        }
    }
    rows.sort_unstable();

    let mut modseqs: Vec<&str> = Vec::with_capacity(rows.len());
    let mut files: Vec<&str> = Vec::with_capacity(rows.len());
    let mut intensities: Vec<f64> = Vec::with_capacity(rows.len());
    let mut retention_times: Vec<f64> = Vec::with_capacity(rows.len());
    let mut detection_types: Vec<&str> = Vec::with_capacity(rows.len());

    for (modseq, file) in rows {
        // The cell is guaranteed present: rows were enumerated straight from `results.quant`.
        let cell = results
            .get(modseq, file)
            .expect("row enumerated from results.quant must resolve");
        modseqs.push(modseq);
        files.push(file);
        intensities.push(cell.intensity);
        retention_times.push(cell.retention_time);
        detection_types.push(cell.detection_type.as_str());
    }

    RecordBatch::try_new(
        peptide_results_schema(),
        vec![
            Arc::new(StringArray::from(modseqs)),
            Arc::new(StringArray::from(files)),
            Arc::new(Float64Array::from(intensities)),
            Arc::new(Float64Array::from(retention_times)),
            Arc::new(StringArray::from(detection_types)),
        ],
    )
    .expect("schema and columns are constructed in lockstep")
}

/// Writes [`PeptideResults`] to a Parquet file at `path` (tidy long form). Overwrites any
/// existing file. Returns the number of rows written (= number of (peptide, file) cells).
pub fn write_peptide_results_parquet<P: AsRef<Path>>(
    results: &PeptideResults,
    path: P,
) -> Result<usize, ParquetError> {
    let batch = peptide_results_record_batch(results);
    let num_rows = batch.num_rows();
    let file = File::create(path).map_err(ParquetError::from)?;
    let mut writer = ArrowWriter::try_new(file, batch.schema(), None)?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(num_rows)
}

// ---------------------------------------------------------------------------------------------
// MBR feature table (PLAN.md P3.2e)
// ---------------------------------------------------------------------------------------------

/// The Arrow schema of the tidy MBR feature table — one row per transferred candidate peak
/// ([`FeatureRow`]). This is the training table P3.3's Python PEP model reads: the donor
/// peptide, the acceptor file, the predicted/apex retention times, the five component scores,
/// the combined `mbr_score`, the apex statistics, and the `random_rt`/`decoy_peptide` labels.
/// All columns are non-nullable.
///
/// Column order matches [`FeatureRow`]'s field order so the Parquet schema reads like the struct.
pub fn feature_table_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("donor_modified_sequence", DataType::Utf8, false),
        Field::new("donor_base_sequence", DataType::Utf8, false),
        Field::new("acceptor_file", DataType::Utf8, false),
        Field::new("predicted_retention_time", DataType::Float64, false),
        Field::new("apex_retention_time", DataType::Float64, false),
        Field::new("intensity", DataType::Float64, false),
        Field::new("ppm_score", DataType::Float64, false),
        Field::new("intensity_score", DataType::Float64, false),
        Field::new("rt_score", DataType::Float64, false),
        Field::new("scan_count_score", DataType::Float64, false),
        Field::new("isotopic_distribution_score", DataType::Float64, false),
        Field::new("mbr_score", DataType::Float64, false),
        Field::new("mass_error", DataType::Float64, false),
        Field::new("scan_count", DataType::UInt64, false),
        Field::new("isotopic_pearson_correlation", DataType::Float64, false),
        Field::new("rt_prediction_error", DataType::Float64, false),
        Field::new("random_rt", DataType::Boolean, false),
        Field::new("decoy_peptide", DataType::Boolean, false),
    ]))
}

/// The permutation that puts `rows` into the canonical feature-table order: sorted by
/// `(acceptor_file, donor_modified_sequence, random_rt, predicted_retention_time,
/// apex_retention_time, intensity)`. Floats are compared with `total_cmp` so `NaN`s order
/// deterministically. Both [`feature_table_record_batch`] and the P3.3 PEP write-back
/// ([`crate::mbr_search::apply_mbr_pep`]) sort with this so a pep computed from the emitted table
/// maps unambiguously back onto the peak the row came from.
pub fn feature_table_sort_order(rows: &[FeatureRow]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| {
        let ra = &rows[a];
        let rb = &rows[b];
        ra.acceptor_file
            .cmp(&rb.acceptor_file)
            .then_with(|| ra.donor_modified_sequence.cmp(&rb.donor_modified_sequence))
            .then_with(|| ra.random_rt.cmp(&rb.random_rt))
            .then_with(|| {
                ra.predicted_retention_time
                    .total_cmp(&rb.predicted_retention_time)
            })
            .then_with(|| ra.apex_retention_time.total_cmp(&rb.apex_retention_time))
            .then_with(|| ra.intensity.total_cmp(&rb.intensity))
    });
    order
}

/// Builds a tidy Arrow [`RecordBatch`] from a slice of [`FeatureRow`]s (one row per row).
///
/// Rows are sorted by [`feature_table_sort_order`] for reproducible output regardless of the order
/// [`run_mbr`] happened to emit them in.
///
/// [`run_mbr`]: crate::mbr_search::run_mbr
pub fn feature_table_record_batch(rows: &[FeatureRow]) -> RecordBatch {
    let order = feature_table_sort_order(rows);
    let sorted: Vec<&FeatureRow> = order.iter().map(|&i| &rows[i]).collect();

    let donor_modified_sequence: Vec<&str> =
        sorted.iter().map(|r| r.donor_modified_sequence.as_str()).collect();
    let donor_base_sequence: Vec<&str> =
        sorted.iter().map(|r| r.donor_base_sequence.as_str()).collect();
    let acceptor_file: Vec<&str> = sorted.iter().map(|r| r.acceptor_file.as_str()).collect();
    let predicted_retention_time: Vec<f64> =
        sorted.iter().map(|r| r.predicted_retention_time).collect();
    let apex_retention_time: Vec<f64> = sorted.iter().map(|r| r.apex_retention_time).collect();
    let intensity: Vec<f64> = sorted.iter().map(|r| r.intensity).collect();
    let ppm_score: Vec<f64> = sorted.iter().map(|r| r.ppm_score).collect();
    let intensity_score: Vec<f64> = sorted.iter().map(|r| r.intensity_score).collect();
    let rt_score: Vec<f64> = sorted.iter().map(|r| r.rt_score).collect();
    let scan_count_score: Vec<f64> = sorted.iter().map(|r| r.scan_count_score).collect();
    let isotopic_distribution_score: Vec<f64> =
        sorted.iter().map(|r| r.isotopic_distribution_score).collect();
    let mbr_score: Vec<f64> = sorted.iter().map(|r| r.mbr_score).collect();
    let mass_error: Vec<f64> = sorted.iter().map(|r| r.mass_error).collect();
    let scan_count: Vec<u64> = sorted.iter().map(|r| r.scan_count as u64).collect();
    let isotopic_pearson_correlation: Vec<f64> =
        sorted.iter().map(|r| r.isotopic_pearson_correlation).collect();
    let rt_prediction_error: Vec<f64> = sorted.iter().map(|r| r.rt_prediction_error).collect();
    let random_rt: Vec<bool> = sorted.iter().map(|r| r.random_rt).collect();
    let decoy_peptide: Vec<bool> = sorted.iter().map(|r| r.decoy_peptide).collect();

    RecordBatch::try_new(
        feature_table_schema(),
        vec![
            Arc::new(StringArray::from(donor_modified_sequence)),
            Arc::new(StringArray::from(donor_base_sequence)),
            Arc::new(StringArray::from(acceptor_file)),
            Arc::new(Float64Array::from(predicted_retention_time)),
            Arc::new(Float64Array::from(apex_retention_time)),
            Arc::new(Float64Array::from(intensity)),
            Arc::new(Float64Array::from(ppm_score)),
            Arc::new(Float64Array::from(intensity_score)),
            Arc::new(Float64Array::from(rt_score)),
            Arc::new(Float64Array::from(scan_count_score)),
            Arc::new(Float64Array::from(isotopic_distribution_score)),
            Arc::new(Float64Array::from(mbr_score)),
            Arc::new(Float64Array::from(mass_error)),
            Arc::new(UInt64Array::from(scan_count)),
            Arc::new(Float64Array::from(isotopic_pearson_correlation)),
            Arc::new(Float64Array::from(rt_prediction_error)),
            Arc::new(BooleanArray::from(random_rt)),
            Arc::new(BooleanArray::from(decoy_peptide)),
        ],
    )
    .expect("schema and columns are constructed in lockstep")
}

/// Writes a slice of [`FeatureRow`]s — the MBR feature table from [`run_mbr`] — to a Parquet file
/// at `path` (tidy long form, one row per transferred peak). Overwrites any existing file.
/// Returns the number of rows written.
///
/// [`run_mbr`]: crate::mbr_search::run_mbr
pub fn write_feature_table_parquet<P: AsRef<Path>>(
    rows: &[FeatureRow],
    path: P,
) -> Result<usize, ParquetError> {
    let batch = feature_table_record_batch(rows);
    let num_rows = batch.num_rows();
    let file = File::create(path).map_err(ParquetError::from)?;
    let mut writer = ArrowWriter::try_new(file, batch.schema(), None)?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(num_rows)
}

/// The feature-table schema with a trailing nullable `mbr_pep` Float64 column appended — the
/// shape returned by the binding once the P3.3 Python PEP model has scored the peaks. `mbr_pep`
/// is nullable because an unscored peak carries `None`.
pub fn feature_table_with_pep_schema() -> Arc<Schema> {
    let mut fields: Vec<Field> = feature_table_schema().fields().iter().map(|f| f.as_ref().clone()).collect();
    fields.push(Field::new("mbr_pep", DataType::Float64, true));
    Arc::new(Schema::new(fields))
}

/// Builds the feature-table [`RecordBatch`] with an appended `mbr_pep` column. `peps` must be in
/// **canonical table order** ([`feature_table_sort_order`]) — i.e. parallel to the rows of
/// [`feature_table_record_batch`] — and the same order the binding hands the table to the Python
/// PEP model. A `None` entry serializes as a Parquet null. Panics if `peps.len() != rows.len()`.
pub fn feature_table_record_batch_with_pep(rows: &[FeatureRow], peps: &[Option<f64>]) -> RecordBatch {
    assert_eq!(
        peps.len(),
        rows.len(),
        "peps must be parallel to the feature rows"
    );
    let base = feature_table_record_batch(rows);
    let mut columns: Vec<Arc<dyn arrow::array::Array>> = base.columns().to_vec();
    columns.push(Arc::new(arrow::array::Float64Array::from(peps.to_vec())));
    RecordBatch::try_new(feature_table_with_pep_schema(), columns)
        .expect("base columns + pep column match the with-pep schema")
}

/// Writes the pep-augmented feature table to Parquet at `path`. `peps` is in canonical table order
/// (see [`feature_table_record_batch_with_pep`]). Returns the number of rows written.
pub fn write_feature_table_with_pep_parquet<P: AsRef<Path>>(
    rows: &[FeatureRow],
    peps: &[Option<f64>],
    path: P,
) -> Result<usize, ParquetError> {
    let batch = feature_table_record_batch_with_pep(rows, peps);
    let num_rows = batch.num_rows();
    let file = File::create(path).map_err(ParquetError::from)?;
    let mut writer = ArrowWriter::try_new(file, batch.schema(), None)?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(num_rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detection_type::DetectionType;
    use crate::results::{PeptideQuant, PeptideResults};
    use arrow::array::{Array, BooleanArray, Float64Array, StringArray, UInt64Array};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use std::collections::HashMap;

    fn sample_results() -> PeptideResults {
        // Two peptides across two files; a mix of detected / not-detected cells.
        let mut quant: HashMap<String, HashMap<String, PeptideQuant>> = HashMap::new();

        let mut pep_a: HashMap<String, PeptideQuant> = HashMap::new();
        pep_a.insert(
            "fileA".to_string(),
            PeptideQuant {
                intensity: 1.0e6,
                retention_time: 12.5,
                detection_type: DetectionType::MSMS,
            },
        );
        pep_a.insert(
            "fileB".to_string(),
            PeptideQuant {
                intensity: 0.0,
                retention_time: 0.0,
                detection_type: DetectionType::NotDetected,
            },
        );
        quant.insert("PEP[ox]TIDE".to_string(), pep_a);

        let mut pep_b: HashMap<String, PeptideQuant> = HashMap::new();
        pep_b.insert(
            "fileA".to_string(),
            PeptideQuant {
                intensity: 2.5e5,
                retention_time: 30.1,
                detection_type: DetectionType::MSMSAmbiguousPeakfinding,
            },
        );
        pep_b.insert(
            "fileB".to_string(),
            PeptideQuant {
                intensity: 7.0e5,
                retention_time: 31.0,
                detection_type: DetectionType::MSMS,
            },
        );
        quant.insert("PEPTIDEK".to_string(), pep_b);

        PeptideResults { quant }
    }

    #[test]
    fn record_batch_is_sorted_and_complete() {
        let results = sample_results();
        let batch = peptide_results_record_batch(&results);

        assert_eq!(batch.num_columns(), 5);
        assert_eq!(batch.num_rows(), 4); // 2 peptides * 2 files

        let modseq = batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let file = batch
            .column(1)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();

        // sorted by (modified_sequence, file_name)
        let pairs: Vec<(&str, &str)> = (0..batch.num_rows())
            .map(|i| (modseq.value(i), file.value(i)))
            .collect();
        // Byte order: 'T' (0x54) < '[' (0x5B), so "PEPTIDEK" sorts before "PEP[ox]TIDE".
        assert_eq!(
            pairs,
            vec![
                ("PEPTIDEK", "fileA"),
                ("PEPTIDEK", "fileB"),
                ("PEP[ox]TIDE", "fileA"),
                ("PEP[ox]TIDE", "fileB"),
            ]
        );
    }

    #[test]
    fn parquet_round_trips_through_the_arrow_reader() {
        let results = sample_results();

        // Write to a unique temp path (Date/random are unavailable in this env, so derive a
        // stable-but-unique name from the process id + a fixed tag).
        let mut path = std::env::temp_dir();
        path.push(format!("flashlfq_p116_{}.parquet", std::process::id()));

        let written = write_peptide_results_parquet(&results, &path).expect("write parquet");
        assert_eq!(written, 4);

        // Read it back with the Parquet Arrow reader and verify every cell survived.
        let file = File::open(&path).expect("open parquet");
        let reader = ParquetRecordBatchReaderBuilder::try_new(file)
            .expect("parquet reader builder")
            .build()
            .expect("parquet reader");
        let batches: Vec<RecordBatch> = reader.map(|b| b.expect("batch")).collect();

        // Reassemble the (modseq, file) -> (intensity, rt, detection) map from the read batches.
        let mut got: HashMap<(String, String), (f64, f64, String)> = HashMap::new();
        for batch in &batches {
            let modseq = batch.column(0).as_any().downcast_ref::<StringArray>().unwrap();
            let file = batch.column(1).as_any().downcast_ref::<StringArray>().unwrap();
            let intensity = batch.column(2).as_any().downcast_ref::<Float64Array>().unwrap();
            let rt = batch.column(3).as_any().downcast_ref::<Float64Array>().unwrap();
            let det = batch.column(4).as_any().downcast_ref::<StringArray>().unwrap();
            for i in 0..batch.num_rows() {
                got.insert(
                    (modseq.value(i).to_string(), file.value(i).to_string()),
                    (intensity.value(i), rt.value(i), det.value(i).to_string()),
                );
            }
        }

        assert_eq!(got.len(), 4);
        let a = &got[&("PEP[ox]TIDE".to_string(), "fileA".to_string())];
        assert!((a.0 - 1.0e6).abs() < 1e-9);
        assert!((a.1 - 12.5).abs() < 1e-9);
        assert_eq!(a.2, "MSMS");

        let b = &got[&("PEPTIDEK".to_string(), "fileA".to_string())];
        assert!((b.0 - 2.5e5).abs() < 1e-9);
        assert_eq!(b.2, "MSMSAmbiguousPeakfinding");

        let nd = &got[&("PEP[ox]TIDE".to_string(), "fileB".to_string())];
        assert_eq!(nd.0, 0.0);
        assert_eq!(nd.2, "NotDetected");

        std::fs::remove_file(&path).ok();
    }

    fn feature_row(
        modseq: &str,
        acceptor: &str,
        predicted_rt: f64,
        intensity: f64,
        random_rt: bool,
    ) -> FeatureRow {
        FeatureRow {
            donor_modified_sequence: modseq.to_string(),
            donor_base_sequence: modseq.chars().filter(|c| c.is_ascii_uppercase()).collect(),
            acceptor_file: acceptor.to_string(),
            predicted_retention_time: predicted_rt,
            apex_retention_time: predicted_rt + 0.05,
            intensity,
            ppm_score: 0.8,
            intensity_score: 0.7,
            rt_score: 0.9,
            scan_count_score: 0.6,
            isotopic_distribution_score: 0.95,
            mbr_score: 75.0,
            mass_error: 1.5,
            scan_count: 12,
            isotopic_pearson_correlation: 0.99,
            rt_prediction_error: 0.05,
            random_rt,
            decoy_peptide: false,
        }
    }

    #[test]
    fn feature_batch_is_sorted_and_complete() {
        // Deliberately emit out of sort order: fileB before fileA, decoy before target.
        let rows = vec![
            feature_row("PEPTIDEK", "fileB", 30.0, 5.0e5, false),
            feature_row("PEPTIDEK", "fileA", 12.0, 1.0e6, true),
            feature_row("PEPTIDEK", "fileA", 12.0, 1.0e6, false),
        ];
        let batch = feature_table_record_batch(&rows);

        assert_eq!(batch.num_columns(), 18);
        assert_eq!(batch.num_rows(), 3);
        assert_eq!(batch.schema(), feature_table_schema());

        let acceptor = batch.column(2).as_any().downcast_ref::<StringArray>().unwrap();
        let random = batch.column(16).as_any().downcast_ref::<BooleanArray>().unwrap();

        // Sorted by (acceptor_file, modseq, random_rt, ...): fileA target, fileA decoy, fileB target.
        let got: Vec<(&str, bool)> = (0..3).map(|i| (acceptor.value(i), random.value(i))).collect();
        assert_eq!(
            got,
            vec![("fileA", false), ("fileA", true), ("fileB", false)]
        );
    }

    #[test]
    fn feature_table_with_pep_appends_aligned_pep_column() {
        // Out-of-sort emission order; peps are in canonical (sorted) table order.
        let rows = vec![
            feature_row("PEPTIDEK", "fileB", 30.0, 5.0e5, false), // sorted pos 2
            feature_row("PEPTIDEK", "fileA", 12.0, 1.0e6, true),  // sorted pos 1
            feature_row("PEPTIDEK", "fileA", 12.0, 1.0e6, false), // sorted pos 0
        ];
        let peps = vec![Some(0.1), Some(0.2), None];
        let batch = feature_table_record_batch_with_pep(&rows, &peps);

        assert_eq!(batch.num_columns(), 19);
        assert_eq!(batch.schema(), feature_table_with_pep_schema());
        assert_eq!(batch.schema().field(18).name(), "mbr_pep");

        let acceptor = batch.column(2).as_any().downcast_ref::<StringArray>().unwrap();
        let random = batch.column(16).as_any().downcast_ref::<BooleanArray>().unwrap();
        let pep = batch.column(18).as_any().downcast_ref::<Float64Array>().unwrap();
        // Row order: (fileA target, fileA decoy, fileB target) carries peps (0.1, 0.2, null).
        assert_eq!(acceptor.value(0), "fileA");
        assert!(!random.value(0));
        assert!((pep.value(0) - 0.1).abs() < 1e-12);
        assert!((pep.value(1) - 0.2).abs() < 1e-12);
        assert!(pep.is_null(2));
    }

    #[test]
    fn feature_table_round_trips_through_the_arrow_reader() {
        let rows = vec![
            feature_row("PEPTIDEK", "fileA", 12.0, 1.0e6, false),
            feature_row("ELVISK", "fileB", 25.5, 3.3e5, true),
        ];

        let mut path = std::env::temp_dir();
        path.push(format!("flashlfq_p32e_{}.parquet", std::process::id()));

        let written = write_feature_table_parquet(&rows, &path).expect("write parquet");
        assert_eq!(written, 2);

        let file = File::open(&path).expect("open parquet");
        let reader = ParquetRecordBatchReaderBuilder::try_new(file)
            .expect("parquet reader builder")
            .build()
            .expect("parquet reader");
        let batches: Vec<RecordBatch> = reader.map(|b| b.expect("batch")).collect();

        let total: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total, 2);

        // Spot-check one cell of each type round-tripped. Sorted by acceptor_file, so the
        // fileA target (PEPTIDEK) sorts before the fileB decoy (ELVISK).
        let b = &batches[0];
        assert_eq!(b.schema(), feature_table_schema());
        let modseq = b.column(0).as_any().downcast_ref::<StringArray>().unwrap();
        let mbr = b.column(11).as_any().downcast_ref::<Float64Array>().unwrap();
        let scan = b.column(13).as_any().downcast_ref::<UInt64Array>().unwrap();
        let random = b.column(16).as_any().downcast_ref::<BooleanArray>().unwrap();
        assert_eq!(modseq.value(0), "PEPTIDEK");
        assert!((mbr.value(0) - 75.0).abs() < 1e-9);
        assert_eq!(scan.value(0), 12);
        assert!(!random.value(0));
        // The fileB decoy is the second row.
        assert_eq!(modseq.value(1), "ELVISK");
        assert!(random.value(1));

        std::fs::remove_file(&path).ok();
    }
}
