r"""Builds `isotope_envelope_explorer.ipynb` — a self-contained notebook for overlaying
isotopic envelopes: an **averagine** at a given mass, a **real protein** at that same mass,
and the **protein + its deamidated form** combination, all on one plot in different colors.

The isotope machinery is pure Python/NumPy (element isotope-pattern convolution), the same
family of model used by the Rust `deconvolution::averagine_*` code and the existing
`analysis/build_envelope_notebook.py` in flashlfq-rust, extended here with:
  - an exact elemental-formula isotope calculator (protein sequence -> C/H/N/O/S -> envelope),
  - deamidation as a formula delta (N,Q amide -> acid: -N -H +O, +0.98402 Da),
  - a binomial partial-deamidation mixture over n sites, and
  - Gaussian peak-shape rendering so the three envelopes overlay as smooth colored curves.

Run:  F:\ProteoStar\.venv\Scripts\python.exe analysis\build_isotope_envelope_notebook.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "isotope_envelope_explorer.ipynb")

_ID = [0]


def _next_id():
    _ID[0] += 1
    return f"iso-{_ID[0]:02d}"


def md(text):
    return {"cell_type": "markdown", "id": _next_id(), "metadata": {},
            "source": text.splitlines(keepends=True)}


def code(text):
    return {"cell_type": "code", "id": _next_id(), "metadata": {}, "execution_count": None,
            "outputs": [], "source": text.strip("\n").splitlines(keepends=True)}


cells = []

cells.append(md(r'''# Isotopic envelope explorer

Overlay and compare isotopic envelopes on a single plot, in different colors:

- **Averagine** at a given mass (green) — the mass-only statistical model.
- **Protein** at that same mass (blue) — the *exact* envelope from a real amino-acid sequence's
  elemental formula.
- **Protein + deamidate** (red) — the combination of the native protein and its deamidated
  form(s). Deamidation (Asn/Gln amide -> acid) adds **+0.98402 Da** per site, which is *not* a
  whole C-13 spacing (1.00336 Da), so the composite envelope is shifted/broadened and its
  isotope-peak ratios are distorted relative to the pure protein.

Everything is pure Python + NumPy — no data files, no repo build needed. Change the sequence,
charge, resolution, and deamidation parameters in the interactive cell or call `plot_envelopes(...)`
directly.
'''))

# ---------------------------------------------------------------------------
cells.append(md(r'''## 1 — Isotope engine (elements, amino acids, formula convolution)

Element isotope abundances match the Rust averagine model. A protein's exact envelope is the
convolution of each element's per-atom isotope pattern raised to its atom count (done by
binary exponentiation so even large proteins are fast).'''))

cells.append(code(r'''
import numpy as np
import matplotlib.pyplot as plt
from math import comb

PROTON = 1.007276466
C13    = 1.0033548381          # C-13 minus C-12; the nominal isotope spacing

# Element -> (mass of the lightest isotope, isotope abundance vector indexed by neutron offset
# from that lightest isotope, zeros filling the gaps). CHNOS drive proteins/averagine; the rest
# are here so you can bolt arbitrary elements onto a formula (see the "weird averagines" section).
EL = {
    "H":  (1.0078250319,  np.array([0.999885, 0.000115])),
    "B":  (10.0129370,    np.array([0.199, 0.801])),
    "C":  (12.0,          np.array([0.9893, 0.0107])),
    "N":  (14.0030740052, np.array([0.99636, 0.00364])),
    "O":  (15.9949146221, np.array([0.99757, 0.00038, 0.00205])),
    "F":  (18.99840322,   np.array([1.0])),
    "Na": (22.98976928,   np.array([1.0])),
    "Mg": (23.985041697,  np.array([0.7899, 0.1000, 0.1101])),
    "Si": (27.97692653,   np.array([0.92223, 0.04685, 0.03092])),
    "P":  (30.97376199,   np.array([1.0])),
    "S":  (31.97207069,   np.array([0.9499, 0.0075, 0.0425, 0.0, 0.0001])),
    "Cl": (34.96885268,   np.array([0.7576, 0.0, 0.2424])),
    "K":  (38.96370649,   np.array([0.932581, 0.000117, 0.067302])),
    "Ca": (39.96259086,   np.array([0.96941, 0.0, 0.00647, 0.00135, 0.02086, 0.0, 0.00004, 0.0, 0.00187])),
    "Mn": (54.93804391,   np.array([1.0])),
    "Fe": (53.9396147,    np.array([0.05845, 0.0, 0.91754, 0.02119, 0.00282])),
    "Co": (58.93319429,   np.array([1.0])),
    "Ni": (57.93534241,   np.array([0.68077, 0.0, 0.26223, 0.011399, 0.036346, 0.0, 0.009255])),
    "Cu": (62.92959772,   np.array([0.6915, 0.0, 0.3085])),
    "Zn": (63.92914201,   np.array([0.4917, 0.0, 0.2773, 0.0404, 0.1845, 0.0, 0.0061])),
    "Se": (73.9224759,    np.array([0.0089, 0.0, 0.0937, 0.0763, 0.2377, 0.0, 0.4961, 0.0, 0.0873])),
    "Br": (78.9183376,    np.array([0.5069, 0.0, 0.4931])),
    "I":  (126.904473,    np.array([1.0])),
}
ELEMENTS = ("C", "H", "N", "O", "S")   # protein/averagine backbone elements only

# Amino-acid *residue* composition (after peptide-bond water loss), order = (C, H, N, O, S).
AA = {
    "G": (2, 3, 1, 1, 0),  "A": (3, 5, 1, 1, 0),  "S": (3, 5, 1, 2, 0),  "P": (5, 7, 1, 1, 0),
    "V": (5, 9, 1, 1, 0),  "T": (4, 7, 1, 2, 0),  "C": (3, 5, 1, 1, 1),  "L": (6, 11, 1, 1, 0),
    "I": (6, 11, 1, 1, 0), "N": (4, 6, 2, 2, 0),  "D": (4, 5, 1, 3, 0),  "Q": (5, 8, 2, 2, 0),
    "K": (6, 12, 2, 1, 0), "E": (5, 7, 1, 3, 0),  "M": (5, 9, 1, 1, 1),  "H": (6, 7, 3, 1, 0),
    "F": (9, 9, 1, 1, 0),  "R": (6, 12, 4, 1, 0), "Y": (9, 9, 1, 2, 0),  "W": (11, 10, 2, 1, 0),
}
WATER = (0, 2, 0, 1, 0)        # one H2O for the intact chain
# Deamidation: amide -> carboxyl (N->D, Q->E). Per site: -N -H +O  ==  +0.98402 Da.
DEAM_DELTA = {"N": -1, "H": -1, "O": +1}


def formula_from_sequence(seq):
    """Sequence string -> {element: count}. Ignores whitespace; upper-cases; skips unknown chars."""
    counts = {e: 0 for e in ELEMENTS}
    seq = "".join(seq.split()).upper()
    unknown = set()
    for aa in seq:
        comp = AA.get(aa)
        if comp is None:
            unknown.add(aa)
            continue
        for e, n in zip(ELEMENTS, comp):
            counts[e] += n
    for e, n in zip(ELEMENTS, WATER):
        counts[e] += n
    if unknown:
        print("  (ignored non-standard residues: %s)" % ", ".join(sorted(unknown)))
    return counts


def apply_delta(formula, delta, times=1):
    """Return a copy of `formula` with `delta` (element->count) applied `times` times."""
    out = dict(formula)
    for e, n in delta.items():
        out[e] = out.get(e, 0) + n * times
    return out


def formula_mono_mass(formula):
    return sum(EL[e][0] * n for e, n in formula.items() if n)


def _poly_pow(base, n, maxlen):
    """(n-fold self-convolution of `base`), truncated to `maxlen`, via binary exponentiation."""
    result = np.array([1.0])
    b = np.asarray(base, dtype=float)
    while n > 0:
        if n & 1:
            result = np.convolve(result, b)[:maxlen]
        n >>= 1
        if n:
            b = np.convolve(b, b)[:maxlen]
    return result


def formula_isotopes(formula, maxlen=80):
    """Exact aggregated isotope distribution of an elemental formula (any elements in EL).

    Returns a probability array indexed by neutron count above the monoisotope (sums ~1)."""
    dist = np.array([1.0])
    for e, n in formula.items():
        if n <= 0:
            continue
        if e not in EL:
            raise KeyError(f"no isotope data for element {e!r}; add it to EL to use it")
        dist = np.convolve(dist, _poly_pow(EL[e][1], n, maxlen))[:maxlen]
    return dist / dist.sum()
'''))

# ---------------------------------------------------------------------------
cells.append(md(r'''## 2 — Averagine model, and envelope/profile builders

`averagine_weights(mass)` reproduces the Rust `deconvolution::averagine_*` shape (round the
averagine residue counts to the target mass, then convolve). `sticks_*` turn a distribution into
`(m/z, intensity)` peaks at charge *z*; `gaussian_profile` renders sticks as a smooth curve at a
given resolving power so envelopes overlay cleanly.'''))

cells.append(code(r'''
# Averagine residue (Senko): average elemental content per unit residue mass.
_AVG = {"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417}
_AVG_RES_MASS = sum(_AVG[e] * EL[e][0] for e in _AVG)


def averagine_formula(mass):
    mult = max(mass / _AVG_RES_MASS, 0.5)
    return {e: round(_AVG[e] * mult) for e in _AVG}


def averagine_weights(mass, maxlen=80):
    """Averagine aggregated isotope distribution at `mass` (probabilities, sum ~1)."""
    return formula_isotopes(averagine_formula(mass), maxlen)


def sticks_from_dist(mono_mass, dist, z, min_rel=1e-4, weight=1.0):
    """Distribution -> [(m/z, intensity)] peaks, anchored at `mono_mass`, spaced by C13/z."""
    base_mz = mono_mass / z + PROTON
    peak = dist.max() or 1.0
    return [(base_mz + k * C13 / z, weight * dist[k])
            for k in range(len(dist)) if dist[k] / peak >= min_rel]


def gaussian_profile(sticks, resolution, oversample=8, pad_mz=0.15):
    """Sum per-peak Gaussians (FWHM = m/z / resolution) into a smooth (x, y) curve."""
    if not sticks:
        return np.array([]), np.array([])
    mzs = np.array([p[0] for p in sticks])
    its = np.array([p[1] for p in sticks])
    sig = mzs / resolution / 2.35482                       # per-peak sigma
    lo, hi = mzs.min() - pad_mz, mzs.max() + pad_mz
    step = max(sig.min() / oversample, (hi - lo) / 4000.0)
    x = np.arange(lo, hi, step)
    y = np.zeros_like(x)
    for mz, it, s in zip(mzs, its, sig):
        y += it * np.exp(-0.5 * ((x - mz) / s) ** 2)
    return x, y


def combination_sticks(formula, z, n_deam, deam_frac, min_rel=1e-4):
    """Native + deamidated mixture as one stick list.

    Models `n_deam` independent sites each deamidated with probability `deam_frac`, i.e. a binomial
    mixture over j = 0..n_deam deamidations (j=0 is native). Each component is weighted by its
    binomial probability and placed at its own monoisotope (+ j * 0.98402 Da)."""
    sticks = []
    for j in range(n_deam + 1):
        w = comb(n_deam, j) * (deam_frac ** j) * ((1.0 - deam_frac) ** (n_deam - j))
        if w <= 0:
            continue
        f = apply_delta(formula, DEAM_DELTA, times=j)
        sticks += sticks_from_dist(formula_mono_mass(f), formula_isotopes(f), z,
                                   min_rel=min_rel, weight=w)
    return sticks
'''))

# ---------------------------------------------------------------------------
cells.append(md(r'''## 3 — The overlay plot

`plot_envelopes(seq, z, n_deam, deam_frac, resolution, ...)` overlays the three envelopes.
By default the averagine is taken at the protein's *own* monoisotopic mass (so "same mass" holds);
pass `avg_mass=` to decouple it. Each curve is normalized to unit height so the **shape**
differences (protein vs averagine) and the **deamidation distortion** (combination vs protein)
are what you read off the plot. Set `show_sticks=True` to also see the underlying peaks.'''))

cells.append(code(r'''
C_AVG, C_PROT, C_COMB = "#2ca02c", "#1f77b4", "#d62728"   # green / blue / red


def _norm_profile(sticks, resolution):
    x, y = gaussian_profile(sticks, resolution)
    if y.size and y.max() > 0:
        y = y / y.max()
    return x, y


def plot_envelopes(seq, z=10, n_deam=1, deam_frac=0.5, resolution=60000,
                   avg_mass=None, show_sticks=False, show_averagine=True,
                   show_protein=True, show_combination=True, mz_pad=None, ax=None):
    formula = formula_from_sequence(seq)
    mono = formula_mono_mass(formula)
    M = float(avg_mass) if avg_mass else mono

    layers = []   # (enabled, label, sticks, color)
    if show_averagine:
        layers.append((True, f"averagine  @ {M:,.2f} Da",
                       sticks_from_dist(M, averagine_weights(M), z), C_AVG))
    if show_protein:
        layers.append((True, f"protein  @ {mono:,.2f} Da",
                       sticks_from_dist(mono, formula_isotopes(formula), z), C_PROT))
    if show_combination:
        layers.append((True, f"protein + deamidate  ({n_deam} site x {deam_frac:.0%})",
                       combination_sticks(formula, z, n_deam, deam_frac), C_COMB))
    layers = [l for l in layers if l[3]]

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 5.2))

    for _, label, sticks, color in layers:
        x, y = _norm_profile(sticks, resolution)
        if x.size:
            ax.plot(x, y, "-", lw=1.8, color=color, label=label, alpha=0.9, zorder=3)
        if show_sticks and sticks:
            peak = max(i for _, i in sticks) or 1.0
            for mz, it in sticks:
                ax.vlines(mz, 0, it / peak, color=color, lw=0.8, alpha=0.35, zorder=2)

    ax.set_xlabel("m/z"); ax.set_ylabel("relative intensity")
    ax.set_ylim(0, 1.08)
    ax.set_title(f"{len(''.join(seq.split()))}-residue protein   mono {mono:,.2f} Da   "
                 f"z = {z}+   (deamidate +{0.98402:.3f} Da/site)", fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.15)
    return ax
'''))

# ---------------------------------------------------------------------------
cells.append(md(r'''## 4 — Example: ubiquitin

Human ubiquitin (76 residues, ~8.56 kDa). At high charge the deamidate's +0.984 Da offset falls
just short of a full isotope spacing, so the red composite leans to the low-mass side of the blue
protein envelope and its peak ratios shift — the classic isotope-ratio signature used to quantify
deamidation.'''))

cells.append(code(r'''
UBIQUITIN = ("MQIFVKTLTG KTITLEVEPS DTIENVKAKI QDKEGIPPDQ QRLIFAGKQL "
             "EDGRTLSDYN IQKESTLHLV LRLRGG")

plot_envelopes(UBIQUITIN, z=10, n_deam=1, deam_frac=0.5, resolution=60000, show_sticks=True)
plt.tight_layout()
'''))

cells.append(md("Zoom on the first few isotope peaks to see the per-peak distortion clearly:"))

cells.append(code(r'''
ax = plot_envelopes(UBIQUITIN, z=10, n_deam=2, deam_frac=0.5, resolution=120000, show_sticks=True)
mono = formula_mono_mass(formula_from_sequence(UBIQUITIN))
base = mono / 10 + PROTON
ax.set_xlim(base - 0.2, base + 6 * C13 / 10 + 0.2)
plt.tight_layout()
'''))

# ---------------------------------------------------------------------------
cells.append(md(r'''## 5 — Interactive explorer

Edit the sequence, then drag the sliders. `charge` sets the plotted charge state; `deam sites` is
the number of deamidation-competent sites; `deam frac` is the per-site deamidation probability;
`resolution` is the rendering resolving power (m/z ÷ FWHM). Toggle the three curves independently.

(Needs `ipywidgets` from the repo `.venv`. If unavailable, call `plot_envelopes(...)` directly.)'''))

cells.append(code(r'''
try:
    from ipywidgets import interact, Textarea, IntSlider, FloatSlider, Checkbox, fixed

    def _explore(seq, charge, deam_sites, deam_frac, resolution, sticks,
                 averagine, protein, combination):
        fig, ax = plt.subplots(figsize=(12, 5.4))
        plot_envelopes(seq, z=charge, n_deam=deam_sites, deam_frac=deam_frac,
                       resolution=resolution, show_sticks=sticks, show_averagine=averagine,
                       show_protein=protein, show_combination=combination, ax=ax)
        fig.tight_layout(); plt.show()

    interact(
        _explore,
        seq=Textarea(value=UBIQUITIN, description="sequence",
                     layout={"width": "90%", "height": "90px"}),
        charge=IntSlider(value=10, min=1, max=40, step=1, description="charge"),
        deam_sites=IntSlider(value=1, min=0, max=6, step=1, description="deam sites"),
        deam_frac=FloatSlider(value=0.5, min=0.0, max=1.0, step=0.05, description="deam frac"),
        resolution=IntSlider(value=60000, min=5000, max=300000, step=5000, description="resolution"),
        sticks=Checkbox(value=False, description="show sticks"),
        averagine=Checkbox(value=True, description="averagine"),
        protein=Checkbox(value=True, description="protein"),
        combination=Checkbox(value=True, description="protein+deamidate"),
    )
except Exception as e:
    print("ipywidgets not available (%s) -- call plot_envelopes(...) directly." % e)
    plot_envelopes(UBIQUITIN)
'''))

# ---------------------------------------------------------------------------
cells.append(md(r'''## 6 — Free exploration

A few ready patterns — copy and tweak:

```python
# Averagine at an arbitrary mass, decoupled from the protein:
plot_envelopes(UBIQUITIN, z=8, avg_mass=8600.0)

# Heavier partial deamidation across several sites (broadens the composite the most):
plot_envelopes(UBIQUITIN, z=12, n_deam=4, deam_frac=0.4, resolution=80000)

# Just protein vs its fully-deamidated form (no averagine):
plot_envelopes(UBIQUITIN, z=10, n_deam=1, deam_frac=1.0,
               show_averagine=False, show_sticks=True)

# Your own sequence:
plot_envelopes("MKVLAT...", z=15)
```

The model assumptions: peaks are anchored at the exact monoisotopic mass of each elemental formula
and spaced by the C-13 increment (1.00336 Da / z) — the same convention as the detector's averagine
model — with relative intensities from exact element-pattern convolution. Deamidation is applied as
the exact formula delta (−N −H +O per site, +0.98402 Da), so both the mass offset and the (tiny)
change in isotope ratios are captured.'''))

cells.append(code(r'''
plot_envelopes(UBIQUITIN, z=12, n_deam=4, deam_frac=0.4, resolution=80000, show_sticks=True)
plt.tight_layout()
'''))

# ---------------------------------------------------------------------------
cells.append(md(r'''## 7 — Arbitrary elements & mods: weird averagines

Bolt any element in the extended isotope table onto a formula and see what it does to the envelope.
The table now carries C H N O S plus **F Na Mg Si P Cl K Ca Mn Fe Co Ni Cu Zn Se Br I B** — enough
to make thoroughly non-peptide-like patterns. Chlorine/bromine add M+2 humps; iron/zinc/selenium
smear the envelope with several nearly-equal isotopes; halogen or metal clusters can flip which peak
is tallest so the monoisotope is no longer the base peak.

Three building blocks feed one overlay plotter:

- `averagine_of(mass)` / `protein_of(seq)` — formula dicts you already know.
- `molecule_of("FeCl3")` — parse a formula string (`Fe Cl3`, `FeCl3`, `C6H12O6` all work).
- `add_mod(formula, "Fe Cl3")` — add a mod (a formula string or dict; negative counts subtract).

Then `plot_formulas([(label, formula), ...], z=..., resolution=...)` overlays them, each normalized
to unit height, labelled with its monoisotopic mass.'''))

cells.append(code(r'''
import re

def molecule_of(s):
    """'Fe2 Cl3', 'FeCl3', 'C6H12O6' -> {element: count}. A bare symbol counts as 1; negatives ok."""
    if isinstance(s, dict):
        return {k: v for k, v in s.items()}
    counts = {}
    for sym, num in re.findall(r"([A-Z][a-z]?)(-?\d*)", s):
        if not sym:
            continue
        counts[sym] = counts.get(sym, 0) + (int(num) if num not in ("", "-") else 1)
    bad = [e for e in counts if e not in EL]
    if bad:
        raise KeyError(f"no isotope data for {bad}; add them to EL to use them")
    return counts

def averagine_of(mass):
    return averagine_formula(mass)

def protein_of(seq):
    return formula_from_sequence(seq)

def add_mod(base_formula, mod):
    """base_formula + mod (formula string or dict). Returns a new dict; counts may go negative."""
    out = dict(base_formula)
    for e, n in molecule_of(mod).items():
        out[e] = out.get(e, 0) + n
    return out


def plot_formulas(specs, z=1, resolution=60000, show_sticks=True, min_rel=1e-4, ax=None):
    """Overlay isotope envelopes for a list of (label, formula_dict[, color]).

    Each curve is a Gaussian profile normalized to unit height and labelled with its mono mass."""
    if ax is None:
        _, ax = plt.subplots(figsize=(11, 5.2))
    cyc = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i, spec in enumerate(specs):
        label, formula = spec[0], spec[1]
        color = spec[2] if len(spec) > 2 else cyc[i % len(cyc)]
        mono = formula_mono_mass(formula)
        sticks = sticks_from_dist(mono, formula_isotopes(formula), z, min_rel=min_rel)
        x, y = gaussian_profile(sticks, resolution)
        if y.size and y.max() > 0:
            y = y / y.max()
        ax.plot(x, y, "-", lw=1.8, color=color, alpha=0.9,
                label=f"{label}   {mono:,.2f} Da")
        if show_sticks and sticks:
            peak = max(v for _, v in sticks) or 1.0
            for mz, it in sticks:
                ax.vlines(mz, 0, it / peak, color=color, lw=0.8, alpha=0.3, zorder=2)
    ax.set_xlabel("m/z"); ax.set_ylabel("relative intensity"); ax.set_ylim(0, 1.08)
    ax.set_title(f"envelope overlay   z = {z}+   (resolution {resolution:,})", fontsize=10)
    ax.legend(fontsize=9, loc="upper right"); ax.grid(True, alpha=0.15)
    return ax

print("elements available:", ", ".join(EL))
'''))

cells.append(md(r'''### A protein with a metal/halogen adduct

Ubiquitin as-is versus the same protein carrying an FeCl3 adduct and, separately, six chlorines —
watch the M+2 / M+4 shoulders fill in and the base peak drift up the envelope.'''))

cells.append(code(r'''
base = protein_of(UBIQUITIN)
plot_formulas([
    ("ubiquitin",        base),
    ("+ FeCl3",          add_mod(base, "Fe Cl3")),
    ("+ Cl6 (halogen)",  add_mod(base, "Cl6")),
], z=10, resolution=60000, show_sticks=False)
plt.tight_layout()
'''))

cells.append(md(r'''### Weird averagines

Start from an averagine at a chosen mass, then stuff it with heavy multi-isotope elements. At z = 1
and modest resolution the distortions are obvious: Cl/Br add spaced M+2 humps, Fe and Se broaden and
shift the whole cluster, and enough of them moves the base peak well off the monoisotope.'''))

cells.append(code(r'''
m = 2000.0
plot_formulas([
    ("averagine",        averagine_of(m)),
    ("+ Cl8",            add_mod(averagine_of(m), "Cl8")),
    ("+ Fe3",            add_mod(averagine_of(m), "Fe3")),
    ("+ Se2 Br2",        add_mod(averagine_of(m), "Se2 Br2")),
], z=1, resolution=15000, show_sticks=True)
plt.tight_layout()
'''))

cells.append(md(r'''### Free-form: mix and match

`plot_formulas` takes any formulas you build. A few starting points:

```python
# Pure inorganic molecules, no protein at all:
plot_formulas([("FeCl3", molecule_of("FeCl3")),
               ("ZnBr2", molecule_of("ZnBr2")),
               ("SeCl4", molecule_of("SeCl4"))], z=1, resolution=8000)

# One weird species at several charge states (plot_formulas uses one z, so call it per z):
w = add_mod(averagine_of(5000), "Fe4 Cl6")
fig, axes = plt.subplots(1, 3, figsize=(16, 4))
for ax, z in zip(axes, (1, 3, 5)):
    plot_formulas([(f"Fe4Cl6-averagine", w)], z=z, resolution=20000, ax=ax)

# Subtract atoms too (negative counts) — e.g. water loss:
plot_formulas([("peptide", protein_of("PEPTIDE")),
               ("- H2O",   add_mod(protein_of("PEPTIDE"), "H-2 O-1"))], z=1, resolution=30000)
```'''))

cells.append(code(r'''
plot_formulas([
    ("FeCl3", molecule_of("FeCl3")),
    ("ZnBr2", molecule_of("ZnBr2")),
    ("SeCl4", molecule_of("SeCl4")),
], z=1, resolution=8000, show_sticks=True)
plt.tight_layout()
'''))

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "flashlfq (.venv)", "language": "python", "name": "flashlfq"},
        "language_info": {"name": "python", "version": "3"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

os.makedirs(HERE, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)
print(f"wrote {OUT}  ({len(cells)} cells)")
