"""P0.5 wheel smoke test — exercises the P0.2–P0.4 binding spikes end-to-end.

Run inside the env where `maturin develop` installed `flashlfq` (the `flashlfq_py`
extension module). Proves the compiled wheel imports and that the three FFI paths
(version string, NumPy move-return, Arrow C Data Interface) work from Python.

Usage: python smoke_test.py   (exits non-zero on any failure)
"""

import sys

import numpy as np
import pyarrow as pa

import flashlfq_py


def check(label, cond):
    status = "ok" if cond else "FAIL"
    print(f"  [{status}] {label}")
    if not cond:
        raise AssertionError(label)


def main():
    print(f"python   {sys.version.split()[0]}")
    print(f"numpy    {np.__version__}")
    print(f"pyarrow  {pa.__version__}")

    # P0.2 surface: core version string reaches Python through the bindings.
    print("core_version():")
    ver = flashlfq_py.core_version()
    check(f"returns a str: {ver!r}", isinstance(ver, str) and len(ver) > 0)

    # P0.3: iota_array(n) -> numpy.ndarray of f64, length n, moved (no copy).
    print("iota_array(5):")
    arr = flashlfq_py.iota_array(5)
    check("is numpy.ndarray", isinstance(arr, np.ndarray))
    check("dtype is float64", arr.dtype == np.float64)
    check("length == 5", arr.shape == (5,))
    check("values are [0,1,2,3,4]", np.array_equal(arr, np.arange(5, dtype=np.float64)))

    # P0.4: demo_table() -> pyarrow.RecordBatch, 2 cols (mz, intensity), 3 rows.
    print("demo_table():")
    batch = flashlfq_py.demo_table()
    check("is pyarrow.RecordBatch", isinstance(batch, pa.RecordBatch))
    check("2 columns", batch.num_columns == 2)
    check("3 rows", batch.num_rows == 3)
    check("column names == ['mz', 'intensity']",
          batch.schema.names == ["mz", "intensity"])
    check("both columns are float64",
          all(pa.types.is_float64(f.type) for f in batch.schema))

    print("\nALL SMOKE CHECKS PASSED")


if __name__ == "__main__":
    main()
