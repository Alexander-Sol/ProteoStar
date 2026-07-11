"""End-to-end acceptance for the P3.2e MBR feature-table binding.

Drives the full FlashLFQ MS2 quant + match-between-runs transfer from Python over the Phase-1
corpus (`AllPSMs.psmtsv` + the two K562 mzMLs) via `quant(..., match_between_runs=True)`, which
returns the **feature table** (one row per transferred candidate peak, targets + `random_rt`
decoys), both ways:

  1. `quant(..., match_between_runs=True, output_path=...)`  -> writes a Parquet file, returns the path.
  2. `quant(..., match_between_runs=True, output_path=None)` -> returns an in-memory `pyarrow.RecordBatch`.

Both must agree and carry the P3.2d feature table: 156 rows (129 target / 27 decoy) with the
18-column schema the P3.3 Python PEP model trains on.

Run with the rust/.venv interpreter:
  & "F:\\flashlfq-rust\\rust\\.venv\\Scripts\\python.exe" rust/flashlfq-py/mbr_smoke_test.py
"""

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

EXPECTED_ROWS = 156
EXPECTED_TARGETS = 129
EXPECTED_DECOYS = 27

EXPECTED_COLS = [
    "donor_modified_sequence",
    "donor_base_sequence",
    "acceptor_file",
    "predicted_retention_time",
    "apex_retention_time",
    "intensity",
    "ppm_score",
    "intensity_score",
    "rt_score",
    "scan_count_score",
    "isotopic_distribution_score",
    "mbr_score",
    "mass_error",
    "scan_count",
    "isotopic_pearson_correlation",
    "rt_prediction_error",
    "random_rt",
    "decoy_peptide",
]

COMPONENT_SCORES = [
    "ppm_score",
    "intensity_score",
    "rt_score",
    "scan_count_score",
    "isotopic_distribution_score",
]


def check_table(table):
    assert table.num_rows == EXPECTED_ROWS, f"expected {EXPECTED_ROWS} rows, got {table.num_rows}"
    assert table.column_names == EXPECTED_COLS, f"unexpected columns: {table.column_names}"

    random_rt = table.column("random_rt").to_pylist()
    decoys = sum(1 for r in random_rt if r)
    targets = sum(1 for r in random_rt if not r)
    assert targets == EXPECTED_TARGETS, f"expected {EXPECTED_TARGETS} targets, got {targets}"
    assert decoys == EXPECTED_DECOYS, f"expected {EXPECTED_DECOYS} decoys, got {decoys}"

    # Acceptors are the two corpus files.
    acceptors = set(table.column("acceptor_file").to_pylist())
    assert acceptors == {
        "20100614_Velos1_TaGe_SA_K562_3",
        "20100614_Velos1_TaGe_SA_K562_4",
    }, f"unexpected acceptors: {acceptors}"

    # Every component score in (0, 1], combined in [0, 100].
    for col in COMPONENT_SCORES:
        vals = table.column(col).to_pylist()
        assert all(0.0 < v <= 1.0 + 1e-9 for v in vals), f"{col} out of (0, 1]"
    mbr = table.column("mbr_score").to_pylist()
    assert all(0.0 <= v <= 100.0 + 1e-9 for v in mbr), "mbr_score out of [0, 100]"

    # Every transferred peak names a donor; scan_count is a non-negative integer.
    assert all(s for s in table.column("donor_modified_sequence").to_pylist()), "empty donor seq"
    assert all(c >= 0 for c in table.column("scan_count").to_pylist()), "negative scan_count"


def main():
    # ---- Path 1: write to Parquet, returns the output path -----------------------------------
    out_path = os.path.join(tempfile.gettempdir(), "flashlfq_mbr_smoke.parquet")
    returned = fl.quant(PSMTSV, RAW, out_path, None, True)
    assert isinstance(returned, str), f"expected str path, got {type(returned)}"
    assert returned == out_path
    assert os.path.exists(out_path), "Parquet file was not written"
    parquet_table = pq.read_table(out_path)
    check_table(parquet_table)
    print(f"[1] Parquet path OK: {out_path} ({parquet_table.num_rows} rows)")

    # ---- Path 2: in-memory RecordBatch (output_path=None) ------------------------------------
    batch = fl.quant(PSMTSV, RAW, None, None, True)
    assert isinstance(batch, pa.RecordBatch), f"expected RecordBatch, got {type(batch)}"
    in_mem = pa.Table.from_batches([batch])
    check_table(in_mem)
    print(f"[2] In-memory RecordBatch OK ({batch.num_rows} rows)")

    # ---- The two paths must produce identical tables -----------------------------------------
    assert parquet_table.equals(in_mem), "Parquet and in-memory tables disagree"
    print("[3] Parquet == in-memory: identical")

    # ---- The MS2 (match_between_runs=False) path is unchanged: peptide table, not feature ----
    ms2 = fl.quant(PSMTSV, RAW, None)
    assert "modified_sequence" in ms2.schema.names, "MS2 path should still return the peptide table"
    assert "mbr_score" not in ms2.schema.names, "MS2 path must not return the feature table"
    print(f"[4] MS2 path unaffected ({ms2.num_rows} peptide x file rows)")

    os.remove(out_path)
    print("ALL MBR SMOKE CHECKS PASSED")


if __name__ == "__main__":
    main()
