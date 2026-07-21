import os, bisect, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
C13=1.0033548; KMAX=3; MASS_PPM=20.0
def load(ds,m):
    p=os.path.join(HERE,f"s_{ds}_{m}.tsv")
    return np.atleast_1d(np.genfromtxt(p,delimiter="\t",names=True)) if os.path.exists(p) else None
def load_gt(p):
    gt=np.atleast_1d(np.genfromtxt(p,delimiter="\t",names=True)); o=np.argsort(gt["mono_mass"]); return gt["mono_mass"][o],gt["rt"][o]
def match(feat,gm,gr):
    out=np.zeros(len(feat),bool)
    for i in range(len(feat)):
        m=feat["mono"][i]; rt=feat["rt"][i]
        for k in range(-KMAX,KMAX+1):
            c=m+k*C13; t=c*MASS_PPM*1e-6
            lo=bisect.bisect_left(gm,c-t); hi=bisect.bisect_right(gm,c+t)
            if lo<hi and np.any(np.abs(gr[lo:hi]-rt)<=1.0): out[i]=True; break
    return out
DS=[("golden",os.path.join(HERE,"..","topdown_bench","gt_golden_all.tsv")),
    ("jurkat",os.path.join(HERE,"..","topdown_bench","gt_jurkat_intersect.tsv"))]
for ds,gtp in DS:
    gm,gr=load_gt(gtp)
    print(f"\n=== {ds} ===")
    print(f"{'run':<10}{'n':>9}{'ID-match%':>10}{'med s2d':>9}{'mean s2d':>10}")
    for m in ["target","shifted","weird","shuffled"]:
        f=load(ds,m)
        if f is None: continue
        mm=match(f,gm,gr)
        print(f"{m:<10}{len(f):>9}{100*mm.mean():>9.1f}%{np.median(f['score2d']):>9.3f}{f['score2d'].mean():>10.3f}")
