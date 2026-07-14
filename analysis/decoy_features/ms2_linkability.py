#!/usr/bin/env python3
"""How many *unmatched* target features above the 95%-decoy-reject threshold could plausibly be
linked to a real MS2 acquisition?

For each unmatched target feature that clears the threshold, we ask whether the raw file contains an
MS2 scan that:
  (1) fired while the feature was eluting or immediately afterwards
        ms2_rt in [RT Start, RT End + RT_TAIL_MIN]
  (2) whose isolation window contains at least one isotope m/z of the feature at its detected charge
        any isotope-k m/z (rel. abundance >= ISO_FLOOR) in [iso_lo, iso_hi]
  (3) BONUS: whose isolation-window center (iso_target) sits on the feature's most abundant isotope
        |iso_target - most_abundant_isotope_mz| <= CENTER_PPM

Controls, scored with the identical logic:
  - ID-matched targets above thr  -> should be ~fully linkable (they were identified, so an MS2 exists)
  - weird-averagine decoys above thr -> coincidental-link background (fake features, real MS2 list)

Threshold per file = 95th percentile of the weird-averagine decoy score2d (rejects 95% of decoys),
matching precision_threshold.py.

Feature isotope envelopes use an averagine model (Senko) so we know which isotopes carry intensity and
which one is the base peak. m/z ladder: mono_mz + k * 1.0033548 / z.

Usage:  python ms2_linkability.py
"""
import os, bisect, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = [
    ("10", r"F:\flashlfq-rust\analysis\dinosaur_bench\gt_10min.tsv", 0.5, "10-min (CA/Lumos)"),
    ("65", r"F:\flashlfq-rust\analysis\longer_gradients\medium_ground_truth.tsv", 1.0, "65-min (glyco)"),
    ("2h", r"F:\flashlfq-rust\analysis\longer_gradients\long_ground_truth_msms.tsv", 0.5, "2-hr (IonStar)"),
]
MASS_PPM = 20.0          # target<->GT mass match
DECOY_REJECT = 0.95      # threshold = 95th pct of weird decoy scores
PROTON = 1.0072764668
C13 = 1.0033548          # 13C - 12C
RT_TAIL_MIN = 0.20       # "immediately afterwards" grace past RT End
ISO_FLOOR = 0.05         # isotopes >=5% of base peak are considered isolable
CENTER_PPM = 20.0        # bonus: isolation center on most-abundant isotope, tight

# ---- averagine isotope envelope -------------------------------------------------
_ISO = {
    "C": [(0, 0.9893), (1, 0.0107)],
    "H": [(0, 0.999885), (1, 0.000115)],
    "N": [(0, 0.99636), (1, 0.00364)],
    "O": [(0, 0.99757), (1, 0.00038), (2, 0.00205)],
    "S": [(0, 0.9499), (1, 0.0075), (2, 0.0425), (4, 0.0001)],
}
_AVG = {"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417}
_UNIT = 111.0543
_MAXISO = 12
_env_cache = {}

def _elem_vec(e):
    v = np.zeros(_MAXISO + 1)
    for d, a in _ISO[e]:
        if d <= _MAXISO:
            v[d] += a
    return v

_EVEC = {e: _elem_vec(e) for e in _ISO}

def envelope(mass):
    """Averagine isotope distribution (by nucleon shift) for a neutral mass. Cached per 25-Da bin."""
    key = int(mass / 25.0)
    hit = _env_cache.get(key)
    if hit is not None:
        return hit
    u = (key * 25.0 + 12.5) / _UNIT
    dist = np.zeros(_MAXISO + 1); dist[0] = 1.0
    for e in _ISO:
        n = int(round(u * _AVG[e]))
        if n <= 0:
            continue
        acc = np.zeros(_MAXISO + 1); acc[0] = 1.0
        b = _EVEC[e].copy(); p = n
        while p > 0:
            if p & 1:
                acc = np.convolve(acc, b)[: _MAXISO + 1]
            p >>= 1
            if p > 0:
                b = np.convolve(b, b)[: _MAXISO + 1]
        dist = np.convolve(dist, acc)[: _MAXISO + 1]
    dist /= dist.sum()
    _env_cache[key] = dist
    return dist

# ---- IO -------------------------------------------------------------------------
def load(fk, model):
    return np.atleast_1d(np.genfromtxt(os.path.join(HERE, f"s{fk}_{model}.tsv"), delimiter="\t", names=True))

def load_detected(fk, model):
    return np.atleast_1d(np.genfromtxt(os.path.join(HERE, f"t{fk}_{model}.detected.tsv"),
                                       delimiter="\t", names=True))

def load_gt(path):
    gt = np.atleast_1d(np.genfromtxt(path, delimiter="\t", names=True))
    o = np.argsort(gt["mono_mass"]); return gt["mono_mass"][o], gt["rt"][o]

def load_ms2(fk):
    m = np.atleast_1d(np.genfromtxt(os.path.join(HERE, f"ms2_{fk}.tsv"), delimiter="\t", names=True))
    o = np.argsort(m["rt_min"])
    return (m["rt_min"][o], m["iso_lo"][o], m["iso_hi"][o], m["iso_target"][o], m["charge"][o].astype(int))

def match_mask(mono, rt, gm, gr, rt_tol):
    out = np.zeros(len(mono), dtype=bool)
    for i in range(len(mono)):
        m = mono[i]; tol = m * MASS_PPM * 1e-6
        lo = bisect.bisect_left(gm, m - tol); hi = bisect.bisect_right(gm, m + tol)
        if lo < hi and np.any(np.abs(gr[lo:hi] - rt[i]) <= rt_tol):
            out[i] = True
    return out

# ---- linkability scan -----------------------------------------------------------
def linkability(mono, z, rt_start, rt_end, ms2):
    """Vectorised per-feature scan. Three nested criteria, weakest to strongest:
        link     - some MS2 window (during elution) overlaps any isolable isotope m/z
        targeted - some MS2 window is *centred* (within CENTER_PPM) on any isotope of the feature
        center   - some MS2 window is centred on the feature's MOST-ABUNDANT isotope (the bonus test)
    Returns three boolean arrays."""
    m_rt, m_lo, m_hi, m_tgt, m_ch = ms2
    n = len(mono)
    link = np.zeros(n, dtype=bool)
    targeted = np.zeros(n, dtype=bool)
    center = np.zeros(n, dtype=bool)
    for i in range(n):
        zz = int(z[i])
        if zz <= 0:
            continue
        mono_mz = (mono[i] + zz * PROTON) / zz
        dist = envelope(mono[i])
        base_k = int(np.argmax(dist))
        ks = np.where(dist >= ISO_FLOOR)[0]
        iso_mz = mono_mz + ks * (C13 / zz)
        base_mz = mono_mz + base_k * (C13 / zz)
        lo = bisect.bisect_left(m_rt, rt_start[i])
        hi = bisect.bisect_right(m_rt, rt_end[i] + RT_TAIL_MIN)
        if lo >= hi:
            continue
        wlo = m_lo[lo:hi]; whi = m_hi[lo:hi]; wtgt = m_tgt[lo:hi]
        hit = np.zeros(hi - lo, dtype=bool)
        for mz in iso_mz:
            hit |= (wlo <= mz) & (mz <= whi)
        if not hit.any():
            continue
        link[i] = True
        tgt_hit = wtgt[hit]
        # centred on any isotope?
        for mz in iso_mz:
            if np.any(np.abs(tgt_hit - mz) <= mz * CENTER_PPM * 1e-6):
                targeted[i] = True
                break
        # centred specifically on the most-abundant isotope?
        if np.any(np.abs(tgt_hit - base_mz) <= base_mz * CENTER_PPM * 1e-6):
            center[i] = True
    return link, targeted, center

def summarize(name, link, targeted, center):
    n = len(link)
    if n == 0:
        return dict(set=name, n=0, linkable=0, linkable_pct=0.0, targeted=0, targeted_pct=0.0,
                    center=0, center_pct=0.0)
    return dict(set=name, n=int(n),
                linkable=int(link.sum()), linkable_pct=round(100 * link.mean(), 1),
                targeted=int(targeted.sum()), targeted_pct=round(100 * targeted.mean(), 1),
                center=int(center.sum()), center_pct=round(100 * center.mean(), 1))

def main():
    all_rows = []
    for fk, gt_path, rt_tol, title in FILES:
        tgt = load(fk, "target"); tdet = load_detected(fk, "target")
        weird = load(fk, "weird"); wdet = load_detected(fk, "weird")
        gm, gr = load_gt(gt_path)
        ms2 = load_ms2(fk)

        matched = match_mask(tgt["mono"], tgt["rt"], gm, gr, rt_tol)
        thr = float(np.percentile(weird["score2d"], DECOY_REJECT * 100))

        um = (~matched) & (tgt["score2d"] >= thr)     # unmatched, above threshold  (the question)
        mm = (matched) & (tgt["score2d"] >= thr)       # matched, above threshold    (positive control)
        wm = (weird["score2d"] >= thr)                 # weird decoys above threshold (null)

        print(f"\n=== {title}  thr={thr:.4f}  MS2 scans={len(ms2[0])} ===")
        print(f"  unmatched>=thr {int(um.sum())}  matched>=thr {int(mm.sum())}  weird>=thr {int(wm.sum())}")

        sets = [
            ("unmatched_target", um, tgt, tdet),
            ("matched_target",   mm, tgt, tdet),
            ("weird_decoy",      wm, weird, wdet),
        ]
        for sname, mask, sf, det in sets:
            idx = np.where(mask)[0]
            link, targeted, center = linkability(sf["mono"][idx], sf["z"][idx],
                                                 det["RT_Start"][idx], det["RT_End"][idx], ms2)
            row = summarize(sname, link, targeted, center)
            row["file"] = fk; row["title"] = title; row["thr"] = thr
            all_rows.append(row)
            print(f"    {sname:<18} n={row['n']:>7}  overlap={row['linkable']:>7} ({row['linkable_pct']:>5.1f}%)"
                  f"  targeted={row['targeted']:>7} ({row['targeted_pct']:>5.1f}%)"
                  f"  center={row['center']:>7} ({row['center_pct']:>5.1f}%)")

    with open(os.path.join(HERE, "ms2_linkability.json"), "w") as fh:
        json.dump(all_rows, fh, indent=2)
    print(f"\nwrote {os.path.join(HERE, 'ms2_linkability.json')}")

if __name__ == "__main__":
    main()
