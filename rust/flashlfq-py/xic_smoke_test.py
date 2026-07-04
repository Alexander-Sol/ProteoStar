"""P1.19 acceptance: build a PeakIndex from an mzML and plot one peptide's XIC.

Exercises the IndexingEngine/XIC Python surface added in P1.19:
  * `PeakIndex.from_mzml(path)`            — build/load the binned m/z index
  * `index.get_indexed_peak(mz, scan, ppm)` — single point query
  * `index.get_xic(mz, rt, ppm, ...)`       — trace a peptide's XIC → RT/intensity NumPy
                                              arrays + integration bounds

It reads the first usable identification from the live `AllPSMs.psmtsv`, builds the index
for that id's K562 mzML, computes the peptide's monoisotopic m/z, traces its XIC, and writes
a plot (`xic_smoke_plot.png`) if matplotlib is available. Prints `ALL XIC SMOKE CHECKS
PASSED` on success.

Run from `rust/flashlfq-py/` with `rust/.venv` active after `maturin develop`:
    & "F:\\flashlfq-rust\\rust\\.venv\\Scripts\\python.exe" rust/flashlfq-py/xic_smoke_test.py
"""

import os
import sys

import numpy as np

import flashlfq_py

PROTON_MASS = 1.007276466879  # mzLib Chemistry/Constants.cs

REPO = os.environ.get("MZLIB_DIR", r"F:\mzLib")  # external mzLib checkout (see README)
TESTDATA = os.path.join(REPO, "mzLib", "Test", "FlashLFQ", "TestData")
PSMTSV = os.path.join(TESTDATA, "AllPSMs.psmtsv")


def read_first_id_for_available_file():
    """Return (bare_file_name, mzml_path, mono_mass, charge, rt) for the first psmtsv row
    whose spectra file exists in TestData as an mzML."""
    with open(PSMTSV, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        fn_i = col["File Name"]
        mono_i = col["Peptide Monoisotopic Mass"]
        charge_i = col["Precursor Charge"]
        rt_i = col["Scan Retention Time"]
        full_i = col["Full Sequence"]
        for line in fh:
            cells = line.rstrip("\n").split("\t")
            bare = os.path.splitext(os.path.basename(cells[fn_i].strip().strip('"')))[0]
            mzml = os.path.join(TESTDATA, bare + ".mzML")
            if not os.path.exists(mzml):
                continue
            mono_token = cells[mono_i].strip().strip('"').split("|")[0]
            try:
                mono = float(mono_token)
                charge = int(float(cells[charge_i].strip().strip('"')))
                rt = float(cells[rt_i].strip().strip('"'))
            except ValueError:
                continue
            return bare, mzml, mono, charge, rt, cells[full_i].strip().strip('"')
    raise RuntimeError("no psmtsv row referenced an available mzML in TestData")


def main():
    bare, mzml, mono, charge, rt, full_seq = read_first_id_for_available_file()
    print(f"Identification: {full_seq}  file={bare}  charge={charge}  rt={rt:.4f}  mono={mono:.5f}")

    # Build the index once (the engine is not cheap to rebuild — P1.19 holds it on the object).
    index = flashlfq_py.PeakIndex.from_mzml(mzml)
    n_scans = index.num_ms1_scans
    print(f"Indexed {n_scans} MS1 scans from {os.path.basename(mzml)}")
    assert n_scans > 0, "expected MS1 scans"

    # Monoisotopic m/z of the precursor.
    mz = (mono + charge * PROTON_MASS) / charge
    print(f"Monoisotopic m/z (z={charge}): {mz:.5f}")

    # Trace the XIC at the engine's peak-finding tolerance (20 ppm).
    xic = index.get_xic(mz, rt, ppm=20.0)
    rts = xic.retention_times
    ints = xic.intensities
    mzs = xic.mzs
    print(repr(xic))
    print(f"XIC length: {len(xic)}")

    # --- Checks on the returned surface ---
    assert isinstance(rts, np.ndarray) and rts.dtype == np.float64, "retention_times must be f64 ndarray"
    assert isinstance(ints, np.ndarray) and ints.dtype == np.float64, "intensities must be f64 ndarray"
    assert isinstance(mzs, np.ndarray) and mzs.dtype == np.float64, "mzs must be f64 ndarray"
    assert len(rts) == len(ints) == len(mzs) == len(xic), "parallel arrays must share length"
    assert len(xic) > 0, "expected a non-empty XIC for the first identification"

    # RT-ascending invariant (core sorts the trace by RT).
    assert np.all(np.diff(rts) >= 0), "retention times must be ascending"

    # Integration bounds are consistent with the arrays.
    assert xic.start_rt == rts.min(), "start_rt must be the earliest RT"
    assert xic.end_rt == rts.max(), "end_rt must be the latest RT"
    assert xic.apex_intensity == ints.max(), "apex_intensity must be the max intensity"
    apex_pos = int(np.argmax(ints))
    assert xic.apex_rt == rts[apex_pos], "apex_rt must be the RT of the max-intensity peak"
    assert xic.start_rt <= xic.apex_rt <= xic.end_rt, "apex RT must lie within the bounds"
    print(
        f"Bounds: start_rt={xic.start_rt:.4f}  apex_rt={xic.apex_rt:.4f}  "
        f"end_rt={xic.end_rt:.4f}  apex_intensity={xic.apex_intensity:.1f}  "
        f"apex_scan_index={xic.apex_scan_index}"
    )

    # Point query: the apex peak must be findable in its own scan.
    found = index.get_indexed_peak(float(mzs[apex_pos]), int(xic.apex_scan_index), ppm=5.0)
    assert found is not None, "apex peak must be queryable via get_indexed_peak"
    fmz, fint, frt, fscan = found
    assert fscan == xic.apex_scan_index, "get_indexed_peak must return the apex scan"
    assert abs(frt - xic.apex_rt) < 1e-6, "get_indexed_peak RT must match the apex RT"
    print(f"get_indexed_peak(apex): mz={fmz:.5f}  intensity={fint:.1f}  rt={frt:.4f}  scan={fscan}")

    # Absent m/z → empty XIC with NaN bounds.
    empty = index.get_xic(50.0, rt, ppm=5.0)
    assert len(empty) == 0, "an absent m/z must trace to an empty XIC"
    assert np.isnan(empty.start_rt) and np.isnan(empty.apex_rt), "empty XIC bounds are NaN"
    assert empty.apex_scan_index == -1, "empty XIC apex scan index is -1"

    # Plot the XIC if matplotlib is present (the literal P1.19 acceptance: "plots one
    # peptide's XIC"). Optional so the script passes in a headless venv without matplotlib.
    out_png = os.path.join(os.path.dirname(__file__), "xic_smoke_plot.png")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.figure(figsize=(7, 4))
        plt.plot(rts, ints, marker="o", ms=3, lw=1)
        plt.axvline(xic.apex_rt, color="r", ls="--", lw=0.8, label="apex")
        plt.axvspan(xic.start_rt, xic.end_rt, color="0.85", alpha=0.4, label="integration bounds")
        plt.xlabel("Retention time (min)")
        plt.ylabel("Intensity")
        plt.title(f"XIC: {full_seq} (z={charge}) in {bare}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_png, dpi=120)
        print(f"Wrote plot to {out_png}")
    except ImportError:
        print("matplotlib not installed — skipping the PNG (arrays above are plot-ready)")

    print("ALL XIC SMOKE CHECKS PASSED")


if __name__ == "__main__":
    sys.exit(main())
