#!/usr/bin/env python3
"""L4 golden dumper — per-MS1-scan (m/z, intensity) peak lists, decoded straight from the
mzML file's binary data arrays.

P1.12 is a **report-only reader-calibration** step, not a pass/fail gate. The goal: quantify
how far the Rust `mzdata` reader's per-scan peak lists drift from the file's actual decoded
contents, so later (L5/L6) parity tolerances are set with eyes open.

This script is the independent reference reader. With no .NET SDK in this environment we cannot
run mzLib's `MsDataFileReader.GetDataFile` + `MassSpectrum.XArray/YArray` directly, but mzLib's
mzML reader does exactly what this script does: for each spectrum it base64-decodes each
`<binaryDataArray>`, zlib-decompresses it (when `zlib compression` / MS:1000574 is present),
and reinterprets the bytes as 32- or 64-bit little-endian floats (MS:1000521 / MS:1000523).
No re-centroiding, no sorting, no zero-trimming — a faithful decode. So the values this script
emits ARE the file's ground-truth peaks, i.e. what any faithful reader (mzLib included) yields.

The Rust side (`read_ms1_scans` over `mzdata`) is diffed against this golden by the report-only
test `tests/l4_reader_calibration.rs`, which records per-scan count deltas and the worst
|Δm/z| / relative |Δintensity| into PLAN.md.

Emits `golden/L4_reader_peaks.tsv`:
  S<TAB>zero_based_scan_index<TAB>one_based_scan_number<TAB>rt_minutes<TAB>peak_count
  P<TAB>zero_based_scan_index<TAB>peak_index<TAB>mz<TAB>intensity      (repr() round-trip floats)

Only MS1 (survey) spectra are emitted, in file order, mirroring `read_ms1_scans`.

Run from anywhere:  python rust/flashlfq-core/parity/dump_l4_reader_peaks.py
"""

import base64
import os
import struct
import zlib
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
CORE = os.path.dirname(HERE)                       # rust/flashlfq-core
REPO = os.environ.get("MZLIB_DIR", r"F:\mzLib")    # external mzLib checkout (see README)
MZML = os.path.join(REPO, "mzLib", "Test", "FlashLFQ", "TestData", "sliced-mzml.mzML")
GOLDEN_TSV = os.path.join(HERE, "golden", "L4_reader_peaks.tsv")

MZML_NS = "http://psi.hupo.org/ms/mzml"

# cvParam accessions we care about.
ACC_MS_LEVEL = "MS:1000511"
ACC_32_BIT = "MS:1000521"
ACC_64_BIT = "MS:1000523"
ACC_ZLIB = "MS:1000574"
ACC_NO_COMP = "MS:1000576"
ACC_MZ_ARRAY = "MS:1000514"
ACC_INTENSITY_ARRAY = "MS:1000515"
ACC_SCAN_START_TIME = "MS:1000016"
# scan start time unit accessions (minute vs second)
ACC_UNIT_MINUTE = "MS:1000038"
ACC_UNIT_SECOND = "MS:1000031"


def q(tag):
    return f"{{{MZML_NS}}}{tag}"


def cvparams(elem):
    """Direct-child cvParams of `elem` as a list of (accession, value, unitAccession)."""
    out = []
    for cv in elem.findall(q("cvParam")):
        out.append((cv.get("accession"), cv.get("value"), cv.get("unitAccession")))
    return out


def decode_binary_array(bda):
    """Decode one <binaryDataArray> to a list of Python floats (matching the file precision)."""
    accs = {a for (a, _v, _u) in cvparams(bda)}
    if ACC_64_BIT in accs:
        fmt_char, width = "d", 8
    elif ACC_32_BIT in accs:
        fmt_char, width = "f", 4
    else:
        raise ValueError("binaryDataArray missing 32/64-bit float cvParam")

    binary_el = bda.find(q("binary"))
    text = (binary_el.text or "").strip()
    raw = base64.b64decode(text)
    if ACC_ZLIB in accs:
        raw = zlib.decompress(raw)
    # else: ACC_NO_COMP or unspecified -> raw bytes as-is

    n = len(raw) // width
    return list(struct.unpack(f"<{n}{fmt_char}", raw[: n * width]))


def array_kind(bda):
    accs = {a for (a, _v, _u) in cvparams(bda)}
    if ACC_MZ_ARRAY in accs:
        return "mz"
    if ACC_INTENSITY_ARRAY in accs:
        return "intensity"
    return None


def scan_rt_minutes(spectrum):
    """Scan start time in minutes (None if absent)."""
    scan_list = spectrum.find(q("scanList"))
    if scan_list is None:
        return None
    for scan in scan_list.findall(q("scan")):
        for (acc, val, unit) in cvparams(scan):
            if acc == ACC_SCAN_START_TIME:
                t = float(val)
                if unit == ACC_UNIT_SECOND:
                    t /= 60.0
                return t
    return None


def main():
    tree = ET.parse(MZML)
    root = tree.getroot()

    # Spectra live under run/spectrumList/spectrum.
    spectra = root.findall(f".//{q('spectrumList')}/{q('spectrum')}")

    rows = []
    zero_based_scan_index = 0
    n_ms1 = 0
    total_peaks = 0
    for spectrum in spectra:
        accs = {a for (a, _v, _u) in cvparams(spectrum)}
        # ms level cvParam value is on the spectrum's own cvParams.
        ms_level = None
        for (acc, val, _u) in cvparams(spectrum):
            if acc == ACC_MS_LEVEL:
                ms_level = int(val)
        if ms_level != 1:
            continue

        one_based = zero_based_scan_index + 1
        rt = scan_rt_minutes(spectrum)

        mzs, intensities = None, None
        bdal = spectrum.find(q("binaryDataArrayList"))
        if bdal is not None:
            for bda in bdal.findall(q("binaryDataArray")):
                kind = array_kind(bda)
                if kind == "mz":
                    mzs = decode_binary_array(bda)
                elif kind == "intensity":
                    intensities = decode_binary_array(bda)
        mzs = mzs or []
        intensities = intensities or []
        if len(mzs) != len(intensities):
            raise ValueError(
                f"scan {zero_based_scan_index}: m/z len {len(mzs)} != intensity len {len(intensities)}"
            )

        rows.append(
            ("S", zero_based_scan_index, one_based, repr(rt if rt is not None else float("nan")), len(mzs))
        )
        for j, (mz, inten) in enumerate(zip(mzs, intensities)):
            rows.append(("P", zero_based_scan_index, j, repr(mz), repr(inten)))

        n_ms1 += 1
        total_peaks += len(mzs)
        zero_based_scan_index += 1

    os.makedirs(os.path.dirname(GOLDEN_TSV), exist_ok=True)
    with open(GOLDEN_TSV, "w", newline="\n", encoding="utf-8") as f:
        f.write("# L4 reader peaks (raw mzML decode) for sliced-mzml.mzML\n")
        f.write("# S\tscan_index\tone_based_scan_number\trt_minutes\tpeak_count\n")
        f.write("# P\tscan_index\tpeak_index\tmz\tintensity\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"wrote {GOLDEN_TSV}")
    print(f"MS1 scans: {n_ms1}; total MS1 peaks: {total_peaks}")


if __name__ == "__main__":
    main()
