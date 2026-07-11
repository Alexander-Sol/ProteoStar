"""End-to-end acceptance for the P1.18 `quant()` binding.

Drives the full FlashLFQ MS2 path from Python over the Phase-1 corpus (`AllPSMs.psmtsv` + the
two K562 mzMLs) two ways:

  1. `quant(..., output_path=...)`  -> writes a Parquet file, returns the path (str).
  2. `quant(..., output_path=None)` -> returns an in-memory `pyarrow.RecordBatch`.

Both must agree with the L6 golden: 708 (peptide x file) rows with the detection-type mix
450 MSMS / 223 NotDetected / 33 MSMSIdentifiedButNotQuantified / 2 MSMSAmbiguousPeakfinding.

Run with the rust/.venv interpreter:
  & "F:\\flashlfq-rust\\rust\\.venv\\Scripts\\python.exe" rust/flashlfq-py/quant_smoke_test.py
"""

import collections
import os
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq

import flashlfq_py as fl

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.normpath(
    os.path.join(os.environ.get("MZLIB_DIR", r"F:\mzLib"), "mzLib", "Test", "FlashLFQ", "TestData")
)

PSMTSV = os.path.join(TEST_DATA, "AllPSMs.psmtsv")
RAW = [
    os.path.join(TEST_DATA, "20100614_Velos1_TaGe_SA_K562_3.mzML"),
    os.path.join(TEST_DATA, "20100614_Velos1_TaGe_SA_K562_4.mzML"),
]

EXPECTED_ROWS = 708
EXPECTED_DETECTION_MIX = {
    "MSMS": 450,
    "NotDetected": 223,
    "MSMSIdentifiedButNotQuantified": 33,
    "MSMSAmbiguousPeakfinding": 2,
}


def detection_mix(detection_types):
    return dict(collections.Counter(detection_types))


def check_table(table):
    assert table.num_rows == EXPECTED_ROWS, f"expected {EXPECTED_ROWS} rows, got {table.num_rows}"
    cols = set(table.column_names)
    expected_cols = {
        "modified_sequence",
        "file_name",
        "intensity",
        "retention_time",
        "detection_type",
    }
    assert cols == expected_cols, f"unexpected columns: {cols}"
    mix = detection_mix(table.column("detection_type").to_pylist())
    assert mix == EXPECTED_DETECTION_MIX, f"detection mix mismatch: {mix}"


def main():
    # ---- Path 1: write to Parquet, returns the output path -----------------------------------
    out_path = os.path.join(tempfile.gettempdir(), "flashlfq_quant_smoke.parquet")
    returned = fl.quant(PSMTSV, RAW, out_path)
    assert isinstance(returned, str), f"expected str path, got {type(returned)}"
    assert returned == out_path
    assert os.path.exists(out_path), "Parquet file was not written"
    parquet_table = pq.read_table(out_path)
    check_table(parquet_table)
    print(f"[1] Parquet path OK: {out_path} ({parquet_table.num_rows} rows)")

    # ---- Path 2: in-memory RecordBatch (output_path=None) ------------------------------------
    batch = fl.quant(PSMTSV, RAW, None)
    assert isinstance(batch, pa.RecordBatch), f"expected RecordBatch, got {type(batch)}"
    in_mem = pa.Table.from_batches([batch])
    check_table(in_mem)
    print(f"[2] In-memory RecordBatch OK ({batch.num_rows} rows)")

    # ---- The two paths must produce identical tables -----------------------------------------
    assert parquet_table.equals(in_mem), "Parquet and in-memory tables disagree"
    print("[3] Parquet == in-memory: identical")

    # ---- Error surface: an unmatched spectra file raises ValueError --------------------------
    try:
        fl.quant(PSMTSV, [RAW[0]], None)  # file _4 ids have no spectra path
        raise AssertionError("expected ValueError for a missing spectra file")
    except ValueError as e:
        assert "no spectra file supplied" in str(e), str(e)
    print("[4] Missing-spectra-file raises ValueError")

    # ---- Error surface: unsupported params rejected ------------------------------------------
    try:
        fl.quant(PSMTSV, RAW, None, {"ppm": 5})
        raise AssertionError("expected ValueError for unsupported params")
    except ValueError as e:
        assert "params are not supported" in str(e), str(e)
    print("[5] Unsupported params raise ValueError")

    os.remove(out_path)
    print("ALL QUANT SMOKE CHECKS PASSED")


if __name__ == "__main__":
    main()
