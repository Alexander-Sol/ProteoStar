#!/usr/bin/env python
"""Ground-truth builders for the top-down recall benchmark.

Emits the score_any.py schema (one row per expected peak):
    mono_mass   monoisotopic neutral mass (Da)
    mz          monoisotopic m/z at `charge`
    charge      identifying precursor charge
    rt          retention time (min)  -- MS2/precursor RT, so score with a LENIENT window
    intensity   (blank)
    detection   source tag / confidence tier

JURKAT (top-down, 02-18-20_jurkat_td_rep1_fract6.raw):
  Union of the Classic and IsoDec AllProteoforms tables (two deconvolution algorithms
  over the same raw). Targets only (Decoy != 'Y'), QValue <= QVAL. Dedupe the union to
  unique (Full Sequence, Precursor Charge); mono_mass = Monoisotopic Mass; rt = Scan
  Retention Time (the MS2 scan RT, not a chromatographic apex -> lenient RT window).

GOLDEN (CD_FDR MSV000082367, golden.raw):
  golden_v4.csv carries Precursor Mass (= monoisotopic neutral mass, verified:
  (Precursor Monoisotopic MZ - proton)*Charge == Precursor Mass), Charge, Retention
  Time, and a Quality tier in {L, M, H}. Emit one combined table plus per-tier tables so
  recall can be reported by confidence.
"""
import csv, os

PROTON = 1.007276466812
QVAL = 0.01

JURKAT_CLASSIC = r"D:\JurkatTopdown\MultipleDeconTest\Classic_AllProteoforms.psmtsv"
JURKAT_ISODEC = r"D:\JurkatTopdown\MultipleDeconTest\IsoDec_AllProteoforms.psmtsv"
GOLDEN_V4 = r"D:\CD_FDR_MSV000082367\other\Manual Validation\golden_v4.csv"

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def mono_mz(mass, z):
    return (mass + z * PROTON) / z


def _read_proteoforms(path):
    """(FullSeq, charge) -> (mono_mass, charge, rt) for target rows at QValue <= QVAL."""
    rows = {}
    n_raw = n_target = 0
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            n_raw += 1
            if (row.get("Decoy") or "").strip().upper() == "Y":
                continue
            try:
                q = float(row["QValue"])
                if q > QVAL:
                    continue
                mass = float(row["Monoisotopic Mass"])
                z = int(float(row["Precursor Charge"]))
                rt = float(row["Scan Retention Time"])
            except (ValueError, KeyError):
                continue
            n_target += 1
            key = (row["Full Sequence"], z)
            if key not in rows:  # keep first-seen RT
                rows[key] = (mass, z, rt)
    return rows, n_raw, n_target


def build_jurkat():
    classic, cr, ct = _read_proteoforms(JURKAT_CLASSIC)
    isodec, ir, it = _read_proteoforms(JURKAT_ISODEC)
    union = dict(classic)
    added = 0
    for k, v in isodec.items():
        if k not in union:
            union[k] = v
            added += 1
    out = os.path.join(OUT_DIR, "gt_jurkat_td.tsv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mono_mass", "mz", "charge", "rt", "intensity", "detection"])
        for (mass, z, rt) in sorted(union.values()):
            w.writerow([f"{mass:.5f}", f"{mono_mz(mass, z):.5f}", z, f"{rt:.4f}", "", "PSM"])
    print(f"JURKAT Classic: {cr} rows -> {len(classic)} unique target proteoforms")
    print(f"JURKAT IsoDec:  {ir} rows -> {len(isodec)} unique target proteoforms")
    print(f"JURKAT union:   {len(union)} unique (FullSeq,charge)  (+{added} unique to IsoDec) -> {out}")
    return len(union)


def build_jurkat_intersection():
    """High-confidence GT: Classic proteoforms confirmed by an IsoDec proteoform agreeing on mass
    (<=10 ppm) AND RT (<=1.0 min). The union carries phantom targets from each engine's mis-assignments
    (which deflate recall); the intersection is the reliable subset both deconvolution engines agree on."""
    classic, _, _ = _read_proteoforms(JURKAT_CLASSIC)
    isodec, _, _ = _read_proteoforms(JURKAT_ISODEC)
    iso = sorted(isodec.values())  # (mass, z, rt), sorted by mass
    import bisect
    iso_masses = [m for (m, _, _) in iso]
    kept = []
    for (mass, z, rt) in classic.values():
        tol = mass * 10 / 1e6
        lo = bisect.bisect_left(iso_masses, mass - tol)
        hi = bisect.bisect_right(iso_masses, mass + tol)
        if any(abs(iso[k][2] - rt) <= 1.0 for k in range(lo, hi)):
            kept.append((mass, z, rt))
    out = os.path.join(OUT_DIR, "gt_jurkat_intersect.tsv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mono_mass", "mz", "charge", "rt", "intensity", "detection"])
        for (mass, z, rt) in sorted(kept):
            w.writerow([f"{mass:.5f}", f"{mono_mz(mass, z):.5f}", z, f"{rt:.4f}", "", "BOTH"])
    print(f"JURKAT intersection: {len(classic)} Classic vs IsoDec (mass 10ppm + RT 1min) -> "
          f"{len(kept)} confirmed proteoforms -> {out}")
    return len(kept)


def build_golden():
    tiers = {"L": [], "M": [], "H": []}
    allrows = []
    n = 0
    with open(GOLDEN_V4, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        for row in r:
            try:
                mass = float(row["Precursor Mass"])
                z = int(float(row["Charge"]))
                rt = float(row["Retention Time"])
            except (ValueError, KeyError):
                continue
            q = (row.get("Quality") or "").strip().upper()
            rec = (f"{mass:.5f}", f"{mono_mz(mass, z):.5f}", z, f"{rt:.4f}", "", q or "?")
            n += 1
            allrows.append(rec)
            if q in tiers:
                tiers[q].append(rec)

    def _write(name, data):
        out = os.path.join(OUT_DIR, name)
        with open(out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["mono_mass", "mz", "charge", "rt", "intensity", "detection"])
            for rec in data:
                w.writerow(rec)
        return out

    _write("gt_golden_all.tsv", allrows)
    for t in ("L", "M", "H"):
        _write(f"gt_golden_{t}.tsv", tiers[t])
    print(f"GOLDEN: {n} rows -> all={len(allrows)}  L={len(tiers['L'])}  "
          f"M={len(tiers['M'])}  H={len(tiers['H'])} -> {OUT_DIR}")
    return n


if __name__ == "__main__":
    build_jurkat()
    build_jurkat_intersection()
    build_golden()
