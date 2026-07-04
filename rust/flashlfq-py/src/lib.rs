//! `flashlfq-py` — thin PyO3 bindings over `flashlfq-core`.
//!
//! This crate's only job is type translation and error mapping at the Python boundary
//! (args → core call → NumPy/Arrow conversion → error mapping). All the FFI sharp edges
//! (GIL/ownership, Arrow ABI skew, NumPy copy-vs-move, panic-across-FFI) are quarantined
//! here so the algorithm in `flashlfq-core` stays plain, testable Rust.
//!
//! At P0.1 it exposes a single placeholder function to prove the binding stack links.

use std::path::PathBuf;

use arrow::pyarrow::IntoPyArrow;
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyString;

use flashlfq_core::engine::{quant as core_quant, quant_mbr as core_quant_mbr};
use flashlfq_core::mbr_search::apply_mbr_pep;
use flashlfq_core::parquet_output::{
    feature_table_record_batch, feature_table_record_batch_with_pep, peptide_results_record_batch,
    write_feature_table_parquet, write_feature_table_with_pep_parquet, write_peptide_results_parquet,
};
use flashlfq_core::peak_indexing::{IndexedMassSpectralPeak, PeakIndexingEngine};
use flashlfq_core::tolerance::PpmTolerance;

/// Returns the `flashlfq-core` version string. Smoke check that the bindings reach core.
#[pyfunction]
fn core_version() -> String {
    flashlfq_core::version().to_string()
}

/// Returns `[0.0, 1.0, ..., (n-1)]` to Python as a 1-D `numpy.ndarray` of `f64`.
///
/// P0.3 binding spike: the core hands us an owned `Vec<f64>`, which `into_pyarray`
/// *moves* into NumPy without copying the buffer (contrast `to_pyarray`, which clones).
/// Proves the NumPy return path end-to-end before any real algorithm output flows through it.
#[pyfunction]
fn iota_array(py: Python<'_>, n: usize) -> Bound<'_, PyArray1<f64>> {
    flashlfq_core::iota_f64(n).into_pyarray(py)
}

/// Returns a small `(mz, intensity)` Arrow table to Python as a `pyarrow.RecordBatch`.
///
/// P0.4 binding spike: the core builds a native `arrow` `RecordBatch`, and `into_pyarrow`
/// hands it to pyarrow over the Arrow C Data Interface — a genuinely zero-copy cross-language
/// transfer (only building the columns costs anything). Proves the Arrow handoff end-to-end
/// before real quant output (peak lists, the results matrix) flows through it. The arrow 54 ↔
/// pyo3 0.23 version pairing is the ABI-skew mitigation called out in the feasibility doc.
#[pyfunction]
fn demo_table(py: Python<'_>) -> PyResult<PyObject> {
    flashlfq_core::demo_record_batch().into_pyarrow(py)
}

/// Runs FlashLFQ MS2 quantification end-to-end from Python (PLAN.md P1.18).
///
/// * `results_path` — a MetaMorpheus `*.psmtsv` of identifications.
/// * `raw_paths` — the spectra files to quantify (mzML in Phase 1). Each identification is
///   matched to a spectra file by *bare* file name (directory + extension stripped), so the
///   order is immaterial and extra paths are ignored.
/// * `output_path` — when given, the peptide × file results table is written to a Parquet file
///   there (tidy long form: `modified_sequence`, `file_name`, `intensity`, `retention_time`,
///   `detection_type`) and the function returns that path as a `str`. When `None`, nothing is
///   written and the function returns the results as an in-memory `pyarrow.RecordBatch`.
/// * `params` — reserved for the parameter surface; Phase 1 only supports the default MS2 path,
///   so this is accepted (for a stable signature) but must be `None`/empty for now.
/// * `match_between_runs` — turns on match-between-runs quantification (the C#
///   `FlashLfqParameters.MatchBetweenRuns` flag). When `True`, the donor/acceptor MBR transfer
///   runs (PLAN.md Phase 3) and the function returns the **feature table** instead of the
///   peptide × file table: one row per transferred candidate peak (donor sequence, acceptor file,
///   predicted/apex RT, the five component scores + combined `mbr_score`, apex statistics, and the
///   `random_rt`/`decoy_peptide` labels). With `output_path` it is written to Parquet and the path
///   is returned; with `None` it is returned as an in-memory `pyarrow.RecordBatch`.
///
/// * `pep_model` — only used when `match_between_runs=True` (PLAN.md P3.3). When given, it must be
///   a callable that takes the MBR feature table (a `pyarrow.RecordBatch`) and returns one
///   posterior error probability (PEP) per row, **in the table's row order**, as a sequence of
///   floats (e.g. a Python `list`). The scores are written back onto the Rust MBR peaks
///   (`MbrChromatographicPeak.mbr_pep`) and the returned/serialized feature table gains a trailing
///   `mbr_pep` column. When `None`, the table is returned without PEP scoring (18 columns).
///
/// Returns either the output path (`str`) or a `pyarrow.RecordBatch` (the peptide table for the
/// MS2 path, the feature table for MBR). Raises `ValueError` on a bad psmtsv, an unmatched spectra
/// file, or a write failure.
#[pyfunction]
#[pyo3(signature = (results_path, raw_paths, output_path=None, params=None, match_between_runs=false, pep_model=None))]
fn quant(
    py: Python<'_>,
    results_path: PathBuf,
    raw_paths: Vec<PathBuf>,
    output_path: Option<PathBuf>,
    params: Option<Bound<'_, PyAny>>,
    match_between_runs: bool,
    pep_model: Option<Bound<'_, PyAny>>,
) -> PyResult<PyObject> {
    // Phase 1 exposes no tunables yet; reject a non-empty params payload rather than silently
    // ignoring it, so callers don't think an unsupported option took effect.
    if let Some(p) = &params {
        let non_empty = match p.len() {
            Ok(n) => n > 0,
            Err(_) => !p.is_none(),
        };
        if non_empty {
            return Err(PyValueError::new_err(
                "params are not supported in Phase 1 (default MS2 path only); pass None",
            ));
        }
    }

    if match_between_runs {
        // Phase 3: run the donor/acceptor transfer and return the MBR feature table.
        let mut mbr = core_quant_mbr(&results_path, &raw_paths)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        // P3.3: if a Python PEP model was supplied, hand it the feature table, take back one PEP
        // per row (in table order), and write the scores onto the Rust peaks. The returned table
        // then carries the trailing `mbr_pep` column.
        if let Some(model) = pep_model {
            let table = feature_table_record_batch(&mbr.feature_rows).into_pyarrow(py)?;
            let returned = model.call1((table,)).map_err(|e| {
                PyValueError::new_err(format!("pep_model callable raised: {e}"))
            })?;
            let peps: Vec<f64> = returned.extract().map_err(|e| {
                PyValueError::new_err(format!(
                    "pep_model must return a sequence of floats (one PEP per feature-table row): {e}"
                ))
            })?;
            apply_mbr_pep(&mut mbr, &peps)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            let peps_opt: Vec<Option<f64>> = peps.iter().map(|&p| Some(p)).collect();
            return match output_path {
                Some(path) => {
                    write_feature_table_with_pep_parquet(&mbr.feature_rows, &peps_opt, &path)
                        .map_err(|e| {
                            PyValueError::new_err(format!("failed to write Parquet output: {e}"))
                        })?;
                    Ok(PyString::new(py, &path.to_string_lossy()).into_any().unbind())
                }
                None => {
                    let batch = feature_table_record_batch_with_pep(&mbr.feature_rows, &peps_opt);
                    batch.into_pyarrow(py)
                }
            };
        }

        return match output_path {
            Some(path) => {
                write_feature_table_parquet(&mbr.feature_rows, &path).map_err(|e| {
                    PyValueError::new_err(format!("failed to write Parquet output: {e}"))
                })?;
                Ok(PyString::new(py, &path.to_string_lossy()).into_any().unbind())
            }
            None => {
                let batch = feature_table_record_batch(&mbr.feature_rows);
                batch.into_pyarrow(py)
            }
        };
    }

    let result = core_quant(&results_path, &raw_paths).map_err(|e| PyValueError::new_err(e.to_string()))?;

    match output_path {
        Some(path) => {
            write_peptide_results_parquet(&result.peptide_results, &path)
                .map_err(|e| PyValueError::new_err(format!("failed to write Parquet output: {e}")))?;
            Ok(PyString::new(py, &path.to_string_lossy()).into_any().unbind())
        }
        None => {
            let batch = peptide_results_record_batch(&result.peptide_results);
            batch.into_pyarrow(py)
        }
    }
}

/// A built m/z peak index over one mzML file's MS1 scans (PLAN.md P1.19).
///
/// Wraps the core `PeakIndexingEngine` so a Python script can build the index once and
/// reuse it for interactive point queries (`get_indexed_peak`) and XIC extraction
/// (`get_xic`) — e.g. to pull and plot one peptide's chromatographic trace. The engine is
/// not cheap to rebuild, so construct a `PeakIndex` once with [`PeakIndex::from_mzml`] and
/// hold onto it.
#[pyclass]
struct PeakIndex {
    engine: PeakIndexingEngine,
}

#[pymethods]
impl PeakIndex {
    /// Builds the index from the MS1 scans of an mzML file (the C# `PeakIndexingEngine`
    /// initialization). Raises `ValueError` if the file cannot be read or holds no MS1 peaks.
    #[staticmethod]
    fn from_mzml(path: PathBuf) -> PyResult<Self> {
        match PeakIndexingEngine::from_mzml(path.clone()) {
            Ok(Some(engine)) => Ok(PeakIndex { engine }),
            Ok(None) => Err(PyValueError::new_err(format!(
                "no indexable MS1 peaks in {}",
                path.display()
            ))),
            Err(e) => Err(PyValueError::new_err(format!(
                "failed to read {}: {e}",
                path.display()
            ))),
        }
    }

    /// Number of MS1 scans indexed.
    #[getter]
    fn num_ms1_scans(&self) -> usize {
        self.engine.scan_info().len()
    }

    /// Finds the peak closest to `mz` in the scan at zero-based index `scan_index`, within
    /// `ppm` tolerance. Returns `(mz, intensity, retention_time, scan_index)` or `None`.
    ///
    /// Faithful to the core `get_indexed_peak`. `ppm` defaults to 5 (the engine's isotope
    /// tolerance); pass a wider value for a coarser lookup.
    #[pyo3(signature = (mz, scan_index, ppm=5.0))]
    fn get_indexed_peak(&self, mz: f64, scan_index: i32, ppm: f64) -> Option<(f64, f64, f64, i32)> {
        self.engine
            .get_indexed_peak(mz, scan_index, &PpmTolerance::new(ppm))
            .map(|p| {
                (
                    p.mz as f64,
                    p.intensity as f64,
                    p.retention_time as f64,
                    p.zero_based_scan_index,
                )
            })
    }

    /// Traces an XIC for `mz` across retention time, beginning near `retention_time`.
    ///
    /// Faithful to the core `get_xic` (m/z, charge-less path): the walk stops after
    /// `missed_scans_allowed` consecutive misses or once a scan's RT is more than
    /// `max_peak_half_width` from the initial peak. `ppm` defaults to 20 (the engine's
    /// peak-finding tolerance), `missed_scans_allowed` to 1, and `max_peak_half_width` to
    /// unbounded (`int.MaxValue` minutes, as in the C# default). Returns an [`Xic`] holding
    /// the RT-ascending trace (RT/intensity/m-z arrays) plus its integration bounds.
    #[pyo3(signature = (mz, retention_time, ppm=20.0, missed_scans_allowed=1, max_peak_half_width=None))]
    fn get_xic(
        &self,
        mz: f64,
        retention_time: f64,
        ppm: f64,
        missed_scans_allowed: i32,
        max_peak_half_width: Option<f64>,
    ) -> Xic {
        let half_width = max_peak_half_width.unwrap_or(i32::MAX as f64);
        let peaks = self.engine.get_xic(
            mz,
            retention_time,
            &PpmTolerance::new(ppm),
            missed_scans_allowed,
            half_width,
            None,
        );
        Xic::from_peaks(&peaks)
    }
}

/// An extracted-ion chromatogram returned to Python (PLAN.md P1.19): the traced peaks as
/// parallel `f64` NumPy arrays (`retention_times`, `intensities`, `mzs`, RT-ascending) plus
/// the integration bounds (`start_rt`/`apex_rt`/`end_rt` and the apex intensity + scan index).
///
/// An empty trace (no peak found) has zero-length arrays; its bound scalars are `NaN`
/// (`apex_scan_index` is `-1`).
#[pyclass]
struct Xic {
    retention_times: Vec<f64>,
    intensities: Vec<f64>,
    mzs: Vec<f64>,
    start_rt: f64,
    apex_rt: f64,
    end_rt: f64,
    apex_intensity: f64,
    apex_scan_index: i32,
}

impl Xic {
    /// Builds the Python-facing XIC from a core peak trace (RT-ascending). Computes the
    /// integration bounds: start/end RT are the min/max over the trace, the apex is the
    /// max-intensity peak (ties → first, matching the core's `total_cmp` max).
    fn from_peaks(peaks: &[IndexedMassSpectralPeak]) -> Xic {
        let retention_times: Vec<f64> = peaks.iter().map(|p| p.retention_time as f64).collect();
        let intensities: Vec<f64> = peaks.iter().map(|p| p.intensity as f64).collect();
        let mzs: Vec<f64> = peaks.iter().map(|p| p.mz as f64).collect();

        let (start_rt, apex_rt, end_rt, apex_intensity, apex_scan_index) = if peaks.is_empty() {
            (f64::NAN, f64::NAN, f64::NAN, f64::NAN, -1)
        } else {
            let apex = peaks
                .iter()
                .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
                .expect("non-empty trace has an apex");
            let start = retention_times.iter().copied().fold(f64::INFINITY, f64::min);
            let end = retention_times
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            (
                start,
                apex.retention_time as f64,
                end,
                apex.intensity as f64,
                apex.zero_based_scan_index,
            )
        };

        Xic {
            retention_times,
            intensities,
            mzs,
            start_rt,
            apex_rt,
            end_rt,
            apex_intensity,
            apex_scan_index,
        }
    }
}

#[pymethods]
impl Xic {
    /// Retention times of the traced peaks (minutes), RT-ascending, as a `numpy.ndarray`.
    #[getter]
    fn retention_times<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, &self.retention_times)
    }

    /// Peak intensities, parallel to `retention_times`, as a `numpy.ndarray`.
    #[getter]
    fn intensities<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, &self.intensities)
    }

    /// Peak m/z values, parallel to `retention_times`, as a `numpy.ndarray`.
    #[getter]
    fn mzs<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, &self.mzs)
    }

    /// Lower integration bound: the earliest retention time in the trace (`NaN` if empty).
    #[getter]
    fn start_rt(&self) -> f64 {
        self.start_rt
    }

    /// Apex retention time: the RT of the most intense peak (`NaN` if empty).
    #[getter]
    fn apex_rt(&self) -> f64 {
        self.apex_rt
    }

    /// Upper integration bound: the latest retention time in the trace (`NaN` if empty).
    #[getter]
    fn end_rt(&self) -> f64 {
        self.end_rt
    }

    /// Intensity of the apex (most intense) peak (`NaN` if empty).
    #[getter]
    fn apex_intensity(&self) -> f64 {
        self.apex_intensity
    }

    /// Zero-based scan index of the apex peak (`-1` if empty).
    #[getter]
    fn apex_scan_index(&self) -> i32 {
        self.apex_scan_index
    }

    /// Number of peaks in the trace.
    fn __len__(&self) -> usize {
        self.retention_times.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "Xic(n={}, start_rt={:.4}, apex_rt={:.4}, end_rt={:.4}, apex_intensity={:.1})",
            self.retention_times.len(),
            self.start_rt,
            self.apex_rt,
            self.end_rt,
            self.apex_intensity,
        )
    }
}

/// The `flashlfq_py` Python module.
#[pymodule]
fn flashlfq_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(core_version, m)?)?;
    m.add_function(wrap_pyfunction!(iota_array, m)?)?;
    m.add_function(wrap_pyfunction!(demo_table, m)?)?;
    m.add_function(wrap_pyfunction!(quant, m)?)?;
    m.add_class::<PeakIndex>()?;
    m.add_class::<Xic>()?;
    Ok(())
}
