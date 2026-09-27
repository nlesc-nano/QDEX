"""Size of a nanocrystal from its coordinates, defined as in small-angle X-ray scattering (SAXS).

SAXS measures the electron-density contrast with the solvent. Organic ligands have about the
solvent's density and are invisible; the inorganic part (core and inorganic surface atoms such as
halides) scatters. The reported diameter is the one a SAXS experiment would give:

* Debye intensity of the inorganic atoms, I(q) = sum_ij Z_i Z_j sin(q r_ij)/(q r_ij), fitted with the
  form factor of a homogeneous sphere over the first lobe (d_saxs);
* Guinier radius of gyration of the electron density, D = 2 sqrt(5/3) R_g (d_guinier).

The formula-unit volume diameter and the convex-hull radius of the core are reported for comparison.
"""
import numpy as np
from scipy.optimize import least_squares

# Atoms treated as organic (invisible to SAXS) unless system.inorganic_elements says otherwise.
ORGANIC_ELEMENTS = {"H", "C", "N", "O", "P", "B", "SI", "F"}

_Z = {
    "H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "NE": 10, "NA": 11,
    "MG": 12, "AL": 13, "SI": 14, "P": 15, "S": 16, "CL": 17, "AR": 18, "K": 19, "CA": 20, "SC": 21,
    "TI": 22, "V": 23, "CR": 24, "MN": 25, "FE": 26, "CO": 27, "NI": 28, "CU": 29, "ZN": 30, "GA": 31,
    "GE": 32, "AS": 33, "SE": 34, "BR": 35, "KR": 36, "RB": 37, "SR": 38, "Y": 39, "ZR": 40, "NB": 41,
    "MO": 42, "TC": 43, "RU": 44, "RH": 45, "PD": 46, "AG": 47, "CD": 48, "IN": 49, "SN": 50, "SB": 51,
    "TE": 52, "I": 53, "XE": 54, "CS": 55, "BA": 56, "LA": 57, "HF": 72, "TA": 73, "W": 74, "RE": 75,
    "OS": 76, "IR": 77, "PT": 78, "AU": 79, "HG": 80, "TL": 81, "PB": 82, "BI": 83,
}


def _sphere_form_factor(q, radius):
    x = np.maximum(q * radius, 1e-8)
    return (3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3) ** 2


def debye_intensity(coords, weights, q, bin_width=0.02):
    """Debye intensity from a histogram of weighted pair distances (one O(N^2) pass, O(n_bins) per q)."""
    from scipy.spatial.distance import cdist
    coords = np.asarray(coords, dtype=float)
    w = np.asarray(weights, dtype=float)
    c = coords.mean(axis=0)
    r_max = 2.0 * float(np.linalg.norm(coords - c, axis=1).max()) + 2.0 * bin_width
    n_bins = int(r_max / bin_width) + 2
    hist = np.zeros(n_bins)
    block = max(1, int(2e7 // max(1, len(coords))))
    for i0 in range(0, len(coords), block):
        d = cdist(coords[i0:i0 + block], coords)
        idx = (d / bin_width).astype(np.int64).ravel()
        hist += np.bincount(idx, weights=(w[i0:i0 + block, None] * w[None, :]).ravel(), minlength=n_bins)[:n_bins]
    centers = (np.arange(n_bins) + 0.5) * bin_width
    return (hist[None, :] * np.sinc(np.outer(q, centers) / np.pi)).sum(axis=1)


def cluster_size(coords, symbols, material=None, inorganic_elements=None):
    """Size metrics of a cluster; diameters in nm."""
    from qdex.hardness import MATERIAL_DB, MATERIAL_ELEMENTS, get_cluster_size_metrics

    coords = np.asarray(coords, dtype=float)
    syms = [str(s).capitalize() for s in symbols]
    up = [s.upper() for s in syms]
    if inorganic_elements:
        inorg_set = {str(e).upper() for e in inorganic_elements}
        inorg = np.array([s in inorg_set for s in up])
    else:
        inorg = np.array([s not in ORGANIC_ELEMENTS for s in up])
    counts = {}
    for s in syms:
        counts[s] = counts.get(s, 0) + 1
    out = {"n_atoms": len(syms), "element_counts": counts,
           "inorganic_elements": sorted({s for s, k in zip(syms, inorg) if k}),
           "n_inorganic": int(inorg.sum())}

    if inorg.sum() >= 2:
        r = coords[inorg]
        z = np.array([_Z.get(s, 20) for s, k in zip(up, inorg) if k], dtype=float)
        c = (z[:, None] * r).sum(0) / z.sum()
        rg = float(np.sqrt((z * ((r - c) ** 2).sum(1)).sum() / z.sum()))
        d_guinier = 2.0 * np.sqrt(5.0 / 3.0) * rg
        # Debye intensity and sphere fit over the first lobe (q R up to ~ 1.3 x first minimum)
        r0 = 0.5 * d_guinier
        q = np.linspace(0.01, 1.3 * 4.493 / max(r0, 1.0), 240)
        I = debye_intensity(r, z, q)
        res = least_squares(lambda p: np.log(p[0] * _sphere_form_factor(q, p[1]) + 1e-30) - np.log(np.abs(I) + 1e-30),
                            [I[0], r0], bounds=([0.0, 0.5], [np.inf, 10.0 * r0 + 10.0]))
        out.update(radius_of_gyration_ang=rg, d_guinier_nm=d_guinier / 10.0, d_saxs_nm=2.0 * res.x[1] / 10.0)

    m = str(material).upper() if material else None
    if m and m in MATERIAL_ELEMENTS and m in MATERIAL_DB and len(MATERIAL_ELEMENTS[m]) == 2:
        els = [e.upper() for e in MATERIAL_ELEMENTS[m]]
        n_core = sum(1 for s in up if s in els)
        v_fu = float(MATERIAL_DB[m][2]) ** 3 / 4.0
        out["d_formula_unit_nm"] = 0.1 * (6.0 * 0.5 * n_core * v_fu / np.pi) ** (1.0 / 3.0)
        out["core_stoichiometry"] = "".join(f"{e.capitalize()}{counts.get(e.capitalize(), 0)}" for e in els)
    hull = get_cluster_size_metrics(coords, symbols, material)
    out["d_hull_nm"] = hull["diameter_hull"] / 10.0
    out["hull_radius_ang"] = hull["R_eff_hull"]
    out["principal_extents_nm"] = [x / 10.0 for x in hull["principal_extents_ang"]]
    out["anisotropy_ratio"] = hull["anisotropy_ratio"]
    return out


def format_cluster_size(size, radius_used=None):
    """Printable block."""
    lines = ["", "--- Cluster size ---"]
    counts = ", ".join(f"{k} {v}" for k, v in size["element_counts"].items())
    lines.append(f"  Atoms              : {size['n_atoms']} ({counts})")
    if "core_stoichiometry" in size:
        lines.append(f"  Core               : {size['core_stoichiometry']}")
    lines.append(f"  SAXS-visible atoms : {size['n_inorganic']} ({', '.join(size['inorganic_elements'])})")
    if "d_saxs_nm" in size:
        lines.append(f"  Diameter (SAXS)    : {size['d_saxs_nm']:.2f} nm   (Debye I(q), sphere fit; reference size)")
        lines.append(f"  Diameter (Guinier) : {size['d_guinier_nm']:.2f} nm   (R_g = {size['radius_of_gyration_ang']:.2f} A)")
    if "d_formula_unit_nm" in size:
        lines.append(f"  Diameter (volume)  : {size['d_formula_unit_nm']:.2f} nm   (core formula units x bulk volume)")
    lines.append(f"  Diameter (hull)    : {size['d_hull_nm']:.2f} nm   (core convex hull + 1.25 A)")
    ext = size["principal_extents_nm"]
    lines.append(f"  Principal extents  : {ext[0]:.2f} x {ext[1]:.2f} x {ext[2]:.2f} nm (anisotropy {size['anisotropy_ratio']:.2f})")
    if radius_used:
        lines.append(f"  Used by the model  : {radius_used}")
    else:
        lines.append("  Used by the model  : none (the size enters only through the orbitals)")
    return "\n".join(lines)
