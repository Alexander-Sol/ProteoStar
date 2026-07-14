import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(HERE, "ms2_linkability.json")))
def get(fk, s): return next(r for r in rows if r["file"]==fk and r["set"]==s)
print(f"{'file':<6}{'n_unm':>9}{'center%':>9}{'null%':>8}{'exp_coinc':>10}{'net':>9}{'enrich':>8}")
for fk in ["10","65","2h"]:
    u=get(fk,"unmatched_target"); d=get(fk,"weird_decoy")
    nullrate=d["center"]/d["n"]
    exp=nullrate*u["n"]; net=u["center"]-exp; enr=(u["center"]/u["n"])/nullrate
    print(f"{fk:<6}{u['n']:>9}{u['center_pct']:>9.1f}{100*nullrate:>8.2f}{exp:>10.0f}{net:>9.0f}{enr:>8.1f}")
