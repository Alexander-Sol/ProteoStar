"""End-to-end acceptance for the P3.3 MBR PEP model + Rust write-back.

Drives the full FlashLFQ MS2 quant + match-between-runs transfer over the Phase-1 corpus
(`AllPSMs.psmtsv` + the two K562 mzMLs), this time passing the Python PEP model
(`pep_model.compute_pep`) into `quant(..., match_between_runs=True, pep_model=...)`. The binding:

  1. builds the 18-column MBR feature table,
  2. hands it to `compute_pep`, which trains the cross-validated gradient-boosted-tree model
     (replacing the C# Microsoft.ML FastTree) and returns one PEP per row,
  3. writes those PEPs back onto the Rust MBR peaks (`apply_mbr_pep`), and
  4. returns the feature table with a trailing `mbr_pep` column (19 columns).

Checks:
  * the returned table has 19 columns ending in `mbr_pep`, every PEP in [0, 1];
  * the model separates classes — mean PEP(decoy) > mean PEP(target), ROC-AUC well above 0.5;
  * the `mbr_pep` Rust returned == `compute_pep` recomputed on the 18-col table (deterministic
    round trip: scores really came back from Python and back out of Rust unchanged);
  * the Parquet and in-memory paths agree;
  * the no-model MBR path still returns the 18-column table (P3.2e unaffected).

Run with the rust/.venv interpreter:
  & "F:\\flashlfq-rust\\rust\\.venv\\Scripts\\python.exe" rust/flashlfq-py/pep_smoke_test.py
"""

import os
import sys
import tempfile

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import flashlfq_py as fl
from pep_model import compute_pep

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
BASE_COLS = 18


def roc_auc(labels, scores):
    """AUC of `scores` predicting `labels` (1 = positive). Rank-based (Mann-Whitney)."""
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=np.float64)
    ranks[order] = np.arange(1, len(scores) + 1)
    # average ranks for ties
    s = np.asarray(scores)[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j > i:
            avg = (i + 1 + j + 1) / 2.0
            ranks[order[i:j + 1]] = avg
        i = j + 1
    pos = np.asarray(labels) == 1
    n_pos = pos.sum()
    n_neg = (~pos).sum()
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    return (ranks[pos].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def main():
    # ---- In-memory: feature table -> Python PEP model -> Rust write-back -> 19-col table -------
    batch = fl.quant(PSMTSV, RAW, output_path=None, match_between_runs=True, pep_model=compute_pep)
    assert isinstance(batch, pa.RecordBatch), f"expected RecordBatch, got {type(batch)}"
    table = pa.Table.from_batches([batch])

    assert table.num_rows == EXPECTED_ROWS, f"expected {EXPECTED_ROWS} rows, got {table.num_rows}"
    assert table.num_columns == BASE_COLS + 1, f"expected {BASE_COLS + 1} cols, got {table.num_columns}"
    assert table.column_names[-1] == "mbr_pep", f"last column should be mbr_pep, got {table.column_names[-1]}"
    print(f"[1] Round trip OK: 19-column table, {table.num_rows} rows")

    peps = np.asarray(table.column("mbr_pep").to_pylist(), dtype=np.float64)
    assert np.all(np.isfinite(peps)), "non-finite PEP returned"
    assert np.all((peps >= 0.0) & (peps <= 1.0)), "PEP out of [0, 1]"
    print(f"[2] All {len(peps)} PEPs finite and in [0, 1]")

    # ---- The model separates targets from decoys --------------------------------------------
    random_rt = np.asarray(table.column("random_rt").to_pylist())
    target_pep = peps[~random_rt]
    decoy_pep = peps[random_rt]
    assert len(target_pep) and len(decoy_pep), "need both targets and decoys"
    assert decoy_pep.mean() > target_pep.mean(), (
        f"decoys should score worse: mean target PEP {target_pep.mean():.4f} "
        f"vs decoy {decoy_pep.mean():.4f}"
    )
    # label 1 = target (low PEP good), so use (1 - PEP) as the target score.
    auc = roc_auc((~random_rt).astype(int), 1.0 - peps)
    assert auc > 0.6, f"ROC-AUC too low ({auc:.3f}); model is not learning"
    print(f"[3] Separation OK: mean PEP target={target_pep.mean():.4f} "
          f"decoy={decoy_pep.mean():.4f}, ROC-AUC={auc:.3f}")

    # ---- Determinism: re-running the model on the 18-col table reproduces the same PEPs -------
    base_batch = fl.quant(PSMTSV, RAW, output_path=None, match_between_runs=True)
    base_table = pa.Table.from_batches([base_batch])
    assert base_table.num_columns == BASE_COLS, "no-model path must stay 18 columns"
    recomputed = np.asarray(compute_pep(base_table), dtype=np.float64)
    assert np.allclose(recomputed, peps, atol=1e-9), (
        "Rust-returned mbr_pep disagrees with a direct compute_pep on the same table"
    )
    print("[4] Deterministic round trip: Rust mbr_pep == compute_pep on the 18-col table")

    # ---- Parquet path agrees with in-memory --------------------------------------------------
    out_path = os.path.join(tempfile.gettempdir(), "flashlfq_pep_smoke.parquet")
    returned = fl.quant(PSMTSV, RAW, output_path=out_path, match_between_runs=True,
                        pep_model=compute_pep)
    assert returned == out_path and os.path.exists(out_path), "Parquet not written"
    parquet_table = pq.read_table(out_path)
    assert parquet_table.equals(table), "Parquet and in-memory tables disagree"
    os.remove(out_path)
    print("[5] Parquet == in-memory (with mbr_pep column)")

    # ---- No-model MBR path is unchanged (P3.2e) ----------------------------------------------
    assert "mbr_pep" not in base_table.column_names, "no-model path must not add mbr_pep"
    print(f"[6] No-model MBR path unaffected ({base_table.num_rows} rows, {base_table.num_columns} cols)")

    print("ALL PEP SMOKE CHECKS PASSED")


if __name__ == "__main__":
    main()
