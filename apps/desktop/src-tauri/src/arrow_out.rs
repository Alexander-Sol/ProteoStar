//! Arrow IPC builders — one per response schema in the IPC contract §3.
//! Each returns a `Vec<u8>` Arrow IPC stream, ready for `Response::new(bytes)`.

use std::collections::HashMap;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, Float32Array, Float64Array, Int32Array, UInt32Array, UInt8Array,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::error::ArrowError;
use arrow::ipc::writer::StreamWriter;
use arrow::record_batch::RecordBatch;

/// Serialize a `RecordBatch` to an Arrow IPC stream (schema `custom_metadata`
/// is carried on the batch's schema).
fn to_ipc(batch: &RecordBatch) -> Result<Vec<u8>, ArrowError> {
    let mut buf = Vec::new();
    {
        let mut writer = StreamWriter::try_new(&mut buf, &batch.schema())?;
        writer.write(batch)?;
        writer.finish()?;
    }
    Ok(buf)
}

/// `trace` schema — shared by get_tic_trace and get_range_xic.
pub fn build_trace(
    retention_time: Vec<f32>,
    intensity: Vec<f32>,
    scan_index: Vec<u32>,
) -> Result<Vec<u8>, ArrowError> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("retentionTime", DataType::Float32, false),
        Field::new("intensity", DataType::Float32, false),
        Field::new("scanIndex", DataType::UInt32, false),
    ]));
    let columns: Vec<ArrayRef> = vec![
        Arc::new(Float32Array::from(retention_time)),
        Arc::new(Float32Array::from(intensity)),
        Arc::new(UInt32Array::from(scan_index)),
    ];
    let batch = RecordBatch::try_new(schema, columns)?;
    to_ipc(&batch)
}

/// `scan_summaries` schema — all MS levels; precursor columns nullable.
pub fn build_scan_summaries(
    one_based_scan_number: Vec<u32>,
    retention_time: Vec<f32>,
    tic: Vec<f32>,
    ms_level: Vec<u8>,
    precursor_mz: Vec<Option<f64>>,
    precursor_scan_index: Vec<Option<i32>>,
    isolation_low: Vec<Option<f64>>,
    isolation_high: Vec<Option<f64>>,
) -> Result<Vec<u8>, ArrowError> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("oneBasedScanNumber", DataType::UInt32, false),
        Field::new("retentionTime", DataType::Float32, false),
        Field::new("tic", DataType::Float32, false),
        Field::new("msLevel", DataType::UInt8, false),
        Field::new("precursorMz", DataType::Float64, true),
        Field::new("precursorScanIndex", DataType::Int32, true),
        Field::new("isolationLow", DataType::Float64, true),
        Field::new("isolationHigh", DataType::Float64, true),
    ]));
    let columns: Vec<ArrayRef> = vec![
        Arc::new(UInt32Array::from(one_based_scan_number)),
        Arc::new(Float32Array::from(retention_time)),
        Arc::new(Float32Array::from(tic)),
        Arc::new(UInt8Array::from(ms_level)),
        Arc::new(Float64Array::from(precursor_mz)),
        Arc::new(Int32Array::from(precursor_scan_index)),
        Arc::new(Float64Array::from(isolation_low)),
        Arc::new(Float64Array::from(isolation_high)),
    ];
    let batch = RecordBatch::try_new(schema, columns)?;
    to_ipc(&batch)
}

/// `spectrum` schema — per-peak columns; scalar metadata rides in `custom_metadata`.
pub fn build_spectrum(
    mz: Vec<f64>,
    intensity: Vec<f32>,
    metadata: HashMap<String, String>,
) -> Result<Vec<u8>, ArrowError> {
    let schema = Arc::new(
        Schema::new(vec![
            Field::new("mz", DataType::Float64, false),
            Field::new("intensity", DataType::Float32, false),
        ])
        .with_metadata(metadata),
    );
    let columns: Vec<ArrayRef> = vec![
        Arc::new(Float64Array::from(mz)),
        Arc::new(Float32Array::from(intensity)),
    ];
    let batch = RecordBatch::try_new(schema, columns)?;
    to_ipc(&batch)
}
