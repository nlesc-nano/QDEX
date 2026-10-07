"""Bulk PBE gaps at the experimental and at the PBE equilibrium lattice constant (make_lattice_scan.py runs).

    python analyse_lattice_scan.py RUNS_DIR BASIS_MOLOPT_UZH GTH_SOC_POTENTIALS [OUT_JSON]

For each material: the total energies of the five scaled lattices are fitted with a third-order
Birch-Murnaghan equation of state (minimum = PBE lattice constant a_PBE at the level of theory of the
dots); the gap is the fundamental spin-free gap along the band path, or, for zinc blende with an
inverted band order (HgX, InAs, InSb), E(Gamma1) - E(Gamma15) from the real-space matrices
(qdex.bulk_soc), as in MATERIAL_DB. The gap at a_PBE comes from a quadratic fit of gap(a) over the
scan; the gap at a_exp from the run at the experimental lattice.
"""
import glob
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from make_lattice_scan import NAMES, BAND_ORDER  # noqa: E402
from qdex.bulk_bands import parse_cp2k_bs  # noqa: E402
from qdex.constants import valence_electrons  # noqa: E402
HA = 27.211386245988


def total_energy(out):
    m = re.findall(r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[(?:a\.u\.|hartree)\]:?\s+([-\d.E+]+)", open(out).read())
    return float(m[-1]) if m else None


def cell(inp):
    s = open(inp).read()
    return np.array([np.array(re.search(r"\n\s+%s ([-\d. ]+)\n" % k, s).group(1).split(), float) for k in "ABC"])


def conventional(cellm, name):
    """Conventional lattice constants (a, c) of the cell: fcc primitive |A| = a/sqrt2, cubic and wurtzite |A| = a."""
    la = np.linalg.norm(cellm, axis=1)
    if name.endswith("_wz"):
        return float(la[0]), float(la[2])
    return (float(la[0] * np.sqrt(2.0)) if not name.endswith("_cubic") else float(la[0])), None


def gap_of_run(d, name, basis, gth):
    """Spin-free gap along the path, or the Gamma s-p band order when inverted (zinc blende)."""
    syms = re.findall(r"\n\s+([A-Z][a-z]?) [-\d.]+ [-\d.]+ [-\d.]+", open(f"{d}/cp2k_job.in").read().split("&COORD")[1])
    n_occ = sum(valence_electrons[e] for e in syms) // 2
    bs = parse_cp2k_bs(glob.glob(f"{d}/*.bs")[0])
    E = bs["bands"]
    gap = float(E[:, n_occ].min() - E[:, n_occ - 1].max())
    if name not in BAND_ORDER:
        return gap, None
    from qdex.bulk_soc import soc_band_structure
    res = soc_band_structure(d, basis, gth, zincblende=True)
    bo = res["info"].get("gamma_band_order") or {}
    return (bo["sf"] if bo.get("inverted") else res["info"]["sf"]["gap"]), bo


def birch_murnaghan(V, E):
    """Minimum (V0, E0, B0 in Ha/A^3) of a third-order Birch-Murnaghan fit of E(V)."""
    from scipy.optimize import curve_fit
    V = np.asarray(V)
    def bm(V, E0, V0, B0, Bp):
        x = (V0 / V) ** (2.0 / 3.0)
        return E0 + 9.0 * V0 * B0 / 16.0 * ((x - 1) ** 3 * Bp + (x - 1) ** 2 * (6 - 4 * x))
    c = np.polyfit(V, E, 2)
    V0 = -c[1] / (2 * c[0])
    p, _ = curve_fit(bm, V, E, p0=[min(E), V0, 2 * c[0] * V0, 4.5], maxfev=20000)
    return p[1], p[0], p[2]


if __name__ == "__main__":
    runs, basis, gth = sys.argv[1], sys.argv[2], sys.argv[3]
    sys.path.insert(0, HERE)
    out_json = sys.argv[4] if len(sys.argv) > 4 else os.path.join(HERE, "lattice_scan.json")
    from qdex.hardness import MATERIAL_DB
    from make_lattice_scan import EXTRA
    table = {}
    for key, name in {**NAMES, **EXTRA}.items():
        cif_cell = cell(f"{HERE}/inputs/{name}/cp2k_job.in")
        V_cif = abs(np.linalg.det(cif_cell))
        a_cif, c_cif = conventional(cif_cell, name)
        pts = []
        for d in sorted(glob.glob(f"{runs}/{name}/s*")):
            sc = float(os.path.basename(d)[1:])
            pts.append((sc, total_energy(f"{d}/cp2k_job.out"), gap_of_run(d, name, basis, gth)[0]))
        sc, E, g = map(np.array, zip(*pts))
        V0, E0, B0 = birch_murnaghan(V_cif * sc ** 3, E)
        s0 = float((V0 / V_cif) ** (1.0 / 3.0))
        g_pbe = float(np.polyval(np.polyfit(sc, g, 2), s0))
        t = dict(name=name, a_cif=round(a_cif, 4), gap_cif_lattice=round(float(np.polyval(np.polyfit(sc, g, 2), 1.0)), 4),
                 scale_pbe=round(s0, 5), a_pbe=round(s0 * a_cif, 4), gap_pbe_lattice=round(g_pbe, 4),
                 bulk_modulus_gpa=round(float(B0) * 4359.744, 1),               # Ha/A^3 -> GPa
                 scan=dict(scale=sc.tolist(), energy_ha=E.tolist(), gap=g.round(4).tolist()))
        if c_cif:
            t.update(c_cif=round(c_cif, 4), c_pbe=round(s0 * c_cif, 4))
        d_exp = f"{runs}/{name}/exp"
        if os.path.isdir(d_exp):
            g_exp, bo_exp = gap_of_run(d_exp, name, basis, gth)
            s_exp = (abs(np.linalg.det(cell(f"{d_exp}/cp2k_job.in"))) / V_cif) ** (1.0 / 3.0)
            t.update(a_exp=round(s_exp * a_cif, 4), gap_exp_lattice=round(g_exp, 4),
                     gap_exp_lattice_from_scan=round(float(np.polyval(np.polyfit(sc, g, 2), s_exp)), 4),
                     inverted=bool(bo_exp and bo_exp.get("inverted")))
        table[key] = t
        print(f"{name:14s} a_CIF {a_cif:.3f} | a_PBE {t['a_pbe']:.3f} ({100 * (s0 - 1):+.2f}% vs CIF) gap {g_pbe:7.3f} | "
              + (f"a_exp {t['a_exp']:.3f} gap {t['gap_exp_lattice']:7.3f} (scan check {t['gap_exp_lattice_from_scan'] - t['gap_exp_lattice']:+.3f})" if 'a_exp' in t else "no exp entry")
              + f" | B0 {t['bulk_modulus_gpa']:5.1f} GPa{'  inverted' if t.get('inverted') else ''}")
    json.dump(table, open(out_json, "w"), indent=1)
    print("wrote", out_json)
