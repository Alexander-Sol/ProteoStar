#!/usr/bin/env python
"""Diagnose whether the intensity-banded analysis fails on large files because there is no signal,
or because the label-free method loses it.

Per intensity band, reports:
  n_match          : PSM-matched targets in the band
  AUC(decon|PSM)   : raw Decon Score, matched-vs-unmatched TARGETS  -> is there raw per-band signal?
  AUC(iso|PSM)     : raw IsoCorr Top3, matched-vs-unmatched targets
  AUC(decon|T-vs-D): Decon, TARGET vs NOISE-DECOY  -> is the decoy a useful negative in this band?
                     (~0.5 => decoy shape == target shape here => decoy is a useless negative)

If AUC(decon|PSM) > 0.55 but AUC(decon|T-vs-D) ~ 0.5, the signal EXISTS but the decoy can't teach it —
so the label-free within-band training failing is a decoy/method problem, not absence of signal.

Usage: band_signal_check.py <target.tsv> <decoy.tsv> <gt.tsv> <rt_delta> [nbands]
"""
import csv, sys, bisect
import numpy as np

MASS_PPM = 20.0


def load(path, want_iso):
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        try:
            d = {"li": np.log1p(float(r["Summed Intensity"])), "decon": float(r["Decon Score"])}
            if want_iso:
                d["iso"] = float(r.get("IsoCorr Top3", -2.0))
                d["mass"] = float(r["Monoisotopic Mass"]); d["rt"] = float(r["RT Apex"])
            out.append(d)
        except (KeyError, ValueError):
            continue
    return out


def load_refs(path):
    refs = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        try:
            refs.append((float(r["mono_mass"]), float(r["rt"])))
        except (ValueError, KeyError):
            continue
    return refs


def auc(pos, neg):
    pos = np.asarray(pos); neg = np.asarray(neg)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    _, inv, counts = np.unique(allv, return_inverse=True, return_counts=True)
    r = np.argsort(np.argsort(allv)) + 1.0
    sums = np.zeros(len(counts)); np.add.at(sums, inv, r); r = (sums / counts)[inv]
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def main():
    T = load(sys.argv[1], True); D = load(sys.argv[2], False)
    refs = load_refs(sys.argv[3]); rt = float(sys.argv[4]); nb = int(sys.argv[5]) if len(sys.argv) > 5 else 10
    order = np.argsort([t["mass"] for t in T]); masses = [T[i]["mass"] for i in order]
    y = np.zeros(len(T), bool)
    for (rm, rr) in refs:
        tol = rm * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rm - tol); hi = bisect.bisect_right(masses, rm + tol)
        for j in range(lo, hi):
            if abs(T[order[j]]["rt"] - rr) <= rt:
                y[order[j]] = True
                break
    li_t = np.array([t["li"] for t in T]); li_d = np.array([d["li"] for d in D])
    dec_t = np.array([t["decon"] for t in T]); dec_d = np.array([d["decon"] for d in D])
    iso_t = np.array([t["iso"] for t in T])
    edges = np.quantile(li_t, np.linspace(0, 1, nb + 1)[1:-1])
    bt = np.digitize(li_t, edges); bd = np.digitize(li_d, edges)
    print(f"targets {len(T)}  decoys {len(D)}  matched {int(y.sum())}  bands {nb}")
    print(f"  {'band':>4} {'n_tgt':>7} {'n_match':>7} {'AUC(decon|PSM)':>15} {'AUC(iso|PSM)':>13} {'AUC(decon|T-vs-D)':>18}")
    for b in range(nb):
        mt = bt == b; md = bd == b
        a_dec = auc(dec_t[mt & y], dec_t[mt & ~y])
        a_iso = auc(iso_t[mt & y], iso_t[mt & ~y])
        a_td = auc(dec_t[mt], dec_d[md])
        print(f"  Q{b+1:<3} {mt.sum():7d} {int(y[mt].sum()):7d} {a_dec:15.3f} {a_iso:13.3f} {a_td:18.3f}")


if __name__ == "__main__":
    main()
