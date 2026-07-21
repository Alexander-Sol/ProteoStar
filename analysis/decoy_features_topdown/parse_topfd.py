#!/usr/bin/env python3
"""Parse a TopFD `*_feature.xml` into a refined-feature TSV our scorer can read.

TopFD reports one <frac_feature> per proteoform (mono_mass, apex_time in SECONDS, a charge range) with a
<single_charge_feature_list> of the charges it was seen at. We emit one row per (proteoform, charge) —
analogous to our per-charge refined features — with the proteoform's mono mass and apex RT (converted to
minutes). decoy_score_export then lays the averagine envelope at each (mono, charge) and computes the
same 2-D / 3-D envelope-fit cosine we use for our own features.

Usage: parse_topfd.py <feature.xml> <out.refined.tsv>
"""
import sys, xml.etree.ElementTree as ET

def main():
    xml_path, out_path = sys.argv[1], sys.argv[2]
    n_feat = n_row = 0
    monos, rts, charges = [], [], []
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("Refined Monoisotopic Mass\tCharge\tApex RT\n")
        # iterparse to keep memory bounded on the large file
        for event, elem in ET.iterparse(xml_path, events=("end",)):
            if elem.tag != "frac_feature":
                continue
            n_feat += 1
            mono = float(elem.findtext("mono_mass"))
            apex_min = float(elem.findtext("apex_time")) / 60.0   # TopFD apex_time is in seconds
            scl = elem.find("single_charge_feature_list")
            seen = set()
            if scl is not None:
                for sc in scl.findall("single_charge_feature"):
                    z = int(float(sc.findtext("charge")))
                    if z in seen:
                        continue
                    seen.add(z)
                    out.write(f"{mono:.5f}\t{z}\t{apex_min:.4f}\n")
                    n_row += 1
                    monos.append(mono); rts.append(apex_min); charges.append(z)
            elem.clear()
    print(f"frac_features: {n_feat}   per-charge rows: {n_row}")
    if monos:
        import statistics
        print(f"mono mass  min/med/max: {min(monos):.0f} / {statistics.median(monos):.0f} / {max(monos):.0f} Da")
        print(f"apex RT    min/med/max: {min(rts):.2f} / {statistics.median(rts):.2f} / {max(rts):.2f} min")
        print(f"charge     min/max:     {min(charges)} / {max(charges)}")

if __name__ == "__main__":
    main()
