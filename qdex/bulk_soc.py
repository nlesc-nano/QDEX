"""Spin-orbit band structure of a bulk crystal from a CP2K k-point calculation.

CP2K writes the Kohn-Sham and overlap matrices of the periodic images
(``&PRINT &KS_CSR_WRITE`` / ``&S_CSR_WRITE`` with ``REAL_SPACE .TRUE.``), the cell
and the atoms (``&PRINT &TREXIO``) and the spin-free bands (``&PRINT &BAND_STRUCTURE``).
From these, at every k-point of the band path:

    H(k) = sum_R e^{ik.R} <chi_mu(0)|H|chi_nu(R)>,   S(k) likewise,
    V_SO(k) = sum_kappa sigma_kappa (x) H_kappa(k),
    H_kappa(k) = 1/2 A(k) (k_ij (x) L_kappa) A(k)^dagger,   A_mu,p(k) = sum_L e^{ik.L} <chi_mu(0)|p(L)>,

with the same GTH projectors, L matrices and SOC constants as the QDEX spinor
Hamiltonian of the quantum dots (qdex.soc_utils). The two-component problem
[[H + Hz, Hx - iHy], [Hx + iHy, H - Hz]] c = E [[S, 0], [0, S]] c is solved in the
full AO basis, and the spinor bands are written in the CP2K .bs format.

Usage::

    python -m qdex.bulk_soc RUN_DIR --basis BASIS_MOLOPT_UZH --gth GTH_SOC_POTENTIALS \\
        [--basis-name DZVP-MOLOPT-PBE-GTH] [--functional PBE] [--name CdSe_zb] [--cif CdSe_zb.cif] [-d DEST]

writes <name>.bs (spin-free, rebuilt from H(R), S(R)), <name>_soc.bs and <name>.json (band edges, gaps,
Gamma band characters, s-p band order, SOC splittings and the semicore manifold used to place the bulk
bands on a dot's energy axis).
"""
import argparse
import glob
import logging
import os
import re

import numpy as np
from scipy.linalg import eigh

logger = logging.getLogger(__name__)

BOHR_PER_ANG = 1.8897259886


def read_image_cells(out_file, kind="KS"):
    """{image number: (n1, n2, n3)} from the 'CSR write| N periodic images' table of a CP2K output."""
    lines = open(out_file).read().splitlines()
    tag = "KS CSR write|" if kind == "KS" else "S CSR write|"
    i = next(n for n, l in enumerate(lines) if l.strip().startswith(tag) and "periodic images" in l)
    n_img = int(lines[i].split("|")[1].split()[0])
    cells = {}
    for l in lines[i + 2:i + 2 + n_img]:
        a, x, y, z = (int(v) for v in l.split())
        cells[a] = (x, y, z)
    return cells


def read_real_space_matrices(run_dir, kind, n_ao, out_file):
    """{(n1, n2, n3): (n_ao, n_ao) matrix} from CP2K real-space CSR text files (Hartree / unitless)."""
    cells = read_image_cells(out_file, kind)
    mats = {}
    for f in glob.glob(os.path.join(run_dir, f"*-{kind}_SPIN_1_R_*-1_0.csr")):
        img = int(re.search(r"_R_(\d+)-1_0\.csr$", os.path.basename(f)).group(1))
        dat = np.loadtxt(f, ndmin=2)
        M = np.zeros((n_ao, n_ao))
        M[dat[:, 0].astype(int) - 1, dat[:, 1].astype(int) - 1] = dat[:, 2]
        mats[cells[img]] = M
    if not mats:
        raise FileNotFoundError(f"No {kind} real-space CSR files in {run_dir}")
    return mats


def read_trexio_structure(path):
    """Cell vectors (rows, Angstrom), atom symbols and Cartesian coordinates (Angstrom)."""
    import h5py
    with h5py.File(path, "r") as f:
        cell = np.array([f["cell/cell_a"][:], f["cell/cell_b"][:], f["cell/cell_c"][:]]) / BOHR_PER_ANG
        coords = np.asarray(f["nucleus/nucleus_coord"][:]) / BOHR_PER_ANG
        labels = [x.decode() if isinstance(x, bytes) else str(x) for x in f["nucleus/nucleus_label"][:]]
    syms = [re.match(r"[A-Z][a-z]?", l).group(0) for l in labels]
    return cell, syms, coords


def _bloch(mats, kfrac):
    n = next(iter(mats.values())).shape[0]
    out = np.zeros((n, n), dtype=complex)
    for cell, M in mats.items():
        out += np.exp(2j * np.pi * np.dot(kfrac, cell)) * M
    return 0.5 * (out + out.conj().T)


def _projector_matrices(shells, syms, coords, cell, soc_tbl, radius_ang):
    """A_L = <chi_mu(0)|p(L)> for the lattice vectors L within radius_ang, and K_kappa = k_ij (x) L_kappa."""
    import libint_cpp
    from qdex.soc_utils import _build_soc_projectors, get_angular_momentum_matrices

    proj0, groups = _build_soc_projectors(syms, coords, soc_tbl)
    n_max = int(np.ceil(radius_ang / np.min(np.linalg.norm(cell, axis=1)))) + 2
    rng = range(-n_max, n_max + 1)
    A_L = {}
    for c in ((i, j, k) for i in rng for j in rng for k in rng):
        L = np.dot(c, cell)
        if np.linalg.norm(L) > radius_ang:
            continue
        projs = [{**p, "center": np.asarray(p["center"]) + L * BOHR_PER_ANG} for p in proj0]
        M = libint_cpp.compute_hgh_overlaps(shells, projs, 4)
        if np.abs(M).max() > 1e-12:
            A_L[c] = M

    pos, col = {}, 0
    for p in proj0:                                # column order of compute_hgh_overlaps
        pos[(p["atom_idx"], p["l"], p["i"])] = col
        col += 2 * p["l"] + 1
    K = [np.zeros((col, col), dtype=complex) for _ in range(3)]
    for (atom, l), grp in groups.items():
        nprj = grp["nprj"]
        kmat = np.zeros((nprj, nprj))
        idx = 0
        for i in range(nprj):
            for j in range(i, nprj):
                kmat[i, j] = kmat[j, i] = grp["k_coeffs"][idx]
                idx += 1
        Ls = get_angular_momentum_matrices(l)
        n = 2 * l + 1
        for i in range(nprj):
            for j in range(nprj):
                ci, cj = pos[(atom, l, i + 1)], pos[(atom, l, j + 1)]
                for kap in range(3):
                    K[kap][ci:ci + n, cj:cj + n] = kmat[i, j] * Ls[kap]
    return A_L, K


class BulkModel:
    """H(k), S(k) and the spinor Hamiltonian of a CP2K bulk run at any k (fractional, reciprocal basis)."""

    def __init__(self, run_dir, basis_file, gth_file, basis_name="DZVP-MOLOPT-PBE-GTH", functional="PBE",
                 trexio_file=None, out_file=None, projector_radius_ang=16.0):
        from qdex.constants import valence_electrons
        from qdex.io_utils import build_shell_dicts, parse_basis, parse_gth_soc_potentials
        import libint_cpp

        trexio_file = trexio_file or glob.glob(os.path.join(run_dir, "*.h5"))[0]
        out_file = out_file or os.path.join(run_dir, "cp2k_job.out")
        self.cell, self.syms, coords = read_trexio_structure(trexio_file)
        basis = parse_basis(basis_file, basis_name, required_elements=set(self.syms))
        n_ao = sum(2 * s["l"] + 1 for s in build_shell_dicts(self.syms, coords, basis))
        self.H_R = read_real_space_matrices(run_dir, "KS", n_ao, out_file)
        self.S_R = read_real_space_matrices(run_dir, "S", n_ao, out_file)
        # CP2K assigns the periodic images with the atoms in scaled coordinates [-0.5, 0.5): place them
        # there (an atom on the +-0.5 boundary goes to whichever image reproduces CP2K's S(0))
        frac = coords @ np.linalg.inv(self.cell)
        frac -= np.floor(frac + 0.5)
        self.coords = self._match_images(frac, basis, build_shell_dicts, libint_cpp)
        self.shells = [{**s, "pure": True} for s in build_shell_dicts(self.syms, self.coords, basis)]
        self.n_ao = n_ao
        self.ao_label = []
        for s in self.shells:
            self.ao_label += [f"{s['sym']}-{'spdfg'[s['l']]}"] * (2 * s["l"] + 1)
        dS = np.abs(libint_cpp.overlap(self.shells, 1) - self.S_R[(0, 0, 0)]).max()
        logger.info(f"  [Bulk SOC] {len(self.H_R)} real-space images; S(0) libint vs CP2K: {dS:.1e}")
        if dS > 1e-3:
            raise ValueError(f"AO basis does not match the CP2K run (max |dS| = {dS:.2e}); check basis file/name.")
        self.soc_tbl = parse_gth_soc_potentials(
            gth_file, {s: valence_electrons.get(s) for s in set(self.syms)}, functional=functional)
        self.A_L, self.K = _projector_matrices(self.shells, self.syms, self.coords, self.cell, self.soc_tbl,
                                               projector_radius_ang)

    def _match_images(self, frac, basis, build_shell_dicts, libint_cpp):
        """Cartesian coordinates (A) with boundary atoms (scaled coordinate ~ +-0.5) in CP2K's image."""
        import itertools
        edge = [(a, c) for a in range(len(frac)) for c in range(3) if abs(abs(frac[a, c]) - 0.5) < 1e-4]
        best, best_d = frac, np.inf
        for signs in itertools.product((0.5, -0.5), repeat=len(edge)):
            f = frac.copy()
            for (a, c), v in zip(edge, signs):
                f[a, c] = v
            xyz = f @ self.cell
            sh = [{**s, "pure": True} for s in build_shell_dicts(self.syms, xyz, basis)]
            d = np.abs(libint_cpp.overlap(sh, 1) - self.S_R[(0, 0, 0)]).max()
            if d < best_d:
                best, best_d = f, d
            if d < 1e-6:
                break
        return best @ self.cell

    def sf(self, kf, vectors=False):
        Hk, Sk = _bloch(self.H_R, kf), _bloch(self.S_R, kf)
        return eigh(Hk, Sk) if vectors else eigh(Hk, Sk, eigvals_only=True)

    def soc(self, kf, vectors=False):
        Hk, Sk = _bloch(self.H_R, kf), _bloch(self.S_R, kf)
        Ak = sum(np.exp(2j * np.pi * np.dot(kf, c)) * M for c, M in self.A_L.items())
        Hx, Hy, Hz = (0.5 * Ak @ Kk @ Ak.conj().T for Kk in self.K)
        Z = np.zeros_like(Sk)
        Htot = np.block([[Hk + Hz, Hx - 1j * Hy], [Hx + 1j * Hy, Hk - Hz]])
        Stot = np.block([[Sk, Z], [Z, Sk]])
        Htot = 0.5 * (Htot + Htot.conj().T)
        return eigh(Htot, Stot) if vectors else eigh(Htot, Stot, eigvals_only=True)

    def characters(self, kf, soc=False):
        """Energies (eV) and Mulliken fractions {element-l: weight} of every band at kf."""
        Sk = _bloch(self.S_R, kf)
        E, C = (self.soc if soc else self.sf)(kf, vectors=True)
        n = self.n_ao
        if soc:
            q = np.real(C[:n].conj() * (Sk @ C[:n])) + np.real(C[n:].conj() * (Sk @ C[n:]))
        else:
            q = np.real(C.conj() * (Sk @ C))
        labels = sorted(set(self.ao_label))
        idx = {l: [i for i, a in enumerate(self.ao_label) if a == l] for l in labels}
        frac = np.array([q[idx[l]].sum(axis=0) for l in labels]).T          # (n_bands, n_labels)
        frac /= frac.sum(axis=1, keepdims=True)
        return E * HA_EV, labels, frac, C


HA_EV = 27.211386245988

# Semicore manifolds used to place the bulk bands on a dot's energy axis, in order of preference:
# a narrow, chemically inert manifold below the valence band (cation d, Cs 5p, else anion s).
SEMICORE_PRIORITY = [("Zn", 2), ("Cd", 2), ("Hg", 2), ("Ga", 2), ("In", 2), ("Cs", 1),
                     ("S", 0), ("Se", 0), ("Te", 0), ("P", 0), ("As", 0), ("Sb", 0),
                     ("Cl", 0), ("Br", 0), ("I", 0)]


def choose_semicore(E_gamma, labels, frac, counts, vbm, min_depth_ev=2.0):
    """Semicore manifolds at Gamma (bands as 0-based [lo, hi)), in order of preference."""
    found = []
    for el, l in SEMICORE_PRIORITY:
        lab = f"{el}-{'spdf'[l]}"
        if lab not in labels or el not in counts:
            continue
        n = counts[el] * (2 * l + 1)
        w = frac[:, labels.index(lab)]
        occ = np.where(E_gamma < vbm - min_depth_ev)[0]
        if occ.size < n:
            continue
        best = np.sort(occ[np.argsort(-w[occ])[:n]])
        if best[-1] - best[0] + 1 == n and w[best].min() > 0.5:
            found.append(dict(element=el, l=l, label=lab, bands_sf=[int(best[0]), int(best[-1]) + 1],
                              min_fraction=float(w[best].min())))
    return found


def _triplet_split(model, kf, triplet, E_soc, C_soc):
    """Spin-orbit splitting of a spin-free triplet: the six spinors with the largest weight on it."""
    Sk = _bloch(model.S_R, kf)
    n = model.n_ao
    _, C = model.sf(kf, vectors=True)
    T = C[:, triplet]
    w = np.sum(np.abs(T.conj().T @ Sk @ C_soc[:n]) ** 2 + np.abs(T.conj().T @ Sk @ C_soc[n:]) ** 2, axis=0)
    six = np.sort(np.argsort(-w)[:6])
    e = np.sort(E_soc[six])
    return float(e[2:].mean() - e[:2].mean()), float(e[:2].mean()), float(e[2:].mean())


def analyse(model, bs, kfrac_segments, sf_segs, soc_segs, n_occ, zincblende=False):
    """Gaps, band edges, Gamma band characters, s-p band order and SOC splittings."""
    allk = np.concatenate(kfrac_segments)
    sf_all = np.concatenate(sf_segs)
    soc_all = np.concatenate(soc_segs)

    def edges(E, nv):
        iv, ic = int(np.argmax(E[:, nv - 1])), int(np.argmin(E[:, nv]))
        return dict(vbm=float(E[iv, nv - 1]), cbm=float(E[ic, nv]), gap=float(E[ic, nv] - E[iv, nv - 1]),
                    vbm_k=allk[iv].round(5).tolist(), cbm_k=allk[ic].round(5).tolist())

    out = dict(sf=edges(sf_all, n_occ), soc=edges(soc_all, 2 * n_occ))
    counts = {e: model.syms.count(e) for e in set(model.syms)}
    out["atoms"] = counts

    g = np.zeros(3)
    E, labels, frac, _ = model.characters(g)
    out["gamma_characters_sf"] = [dict(band=i, energy=round(float(E[i]), 4),
                                       character={labels[j]: round(float(frac[i, j]), 3)
                                                  for j in np.argsort(-frac[i])[:3] if frac[i, j] > 0.05})
                                  for i in range(min(len(E), n_occ + 8))]
    cands = choose_semicore(E, labels, frac, counts, out["sf"]["vbm"])
    for semi in cands:
        lo, hi = semi["bands_sf"]
        semi["bands_soc"] = [2 * lo, 2 * hi]
        semi["level_sf"] = float(sf_all[:, lo:hi].mean())
        semi["level_soc"] = float(soc_all[:, 2 * lo:2 * hi].mean())
        semi["depth_below_vbm_sf"] = float(out["sf"]["vbm"] - semi["level_sf"])
    out["semicore"] = cands[0] if cands else None
    out["semicore_candidates"] = cands

    # Zinc blende at Gamma: the p-like triplet Gamma15 near the Fermi level (partly filled in the
    # inverted semimetals HgX, InAs, InSb), its spin-orbit splitting Gamma8/Gamma7, and the s-p band
    # order E(Gamma1/Gamma6) - E(Gamma15/Gamma8) (negative = inverted)
    if zincblende:
        Es, lab_s, fr_s, C_soc = model.characters(g, soc=True)
        s_w = fr_s[:, [i for i, l in enumerate(lab_s) if l.endswith("-s")]].sum(axis=1)
        p_sf = frac[:, [i for i, l in enumerate(labels) if l.endswith("-p")]].sum(axis=1)
        s_sf = frac[:, [i for i, l in enumerate(labels) if l.endswith("-s")]].sum(axis=1)
        ef = 0.5 * (E[n_occ - 1] + E[n_occ])
        cands = [list(range(i, i + 3)) for i in range(max(0, n_occ - 5), min(len(E) - 2, n_occ + 2))
                 if np.ptp(E[i:i + 3]) < 0.01 and p_sf[i:i + 3].min() > 0.5 and abs(E[i:i + 3].mean() - ef) < 2.5]
        if cands:
            triplet = min(cands, key=lambda t: abs(E[t].mean() - ef))
            e15 = float(E[triplet].mean())
            d_so, e7, e8 = _triplet_split(model, g, triplet, Es, C_soc)
            out["gamma_vb_triplet"] = dict(sf=e15, delta_so=d_so, soc_lower=e7, soc_upper=e8,
                                           filled=bool(triplet[-1] < n_occ))
            near = [i for i in range(len(E)) if i not in triplet and abs(E[i] - e15) < 3.5]
            i_s = max(near, key=lambda i: s_sf[i])
            # Gamma6: the most s-like spinor pair within 3.5 eV, outside the six triplet spinors
            Sk = _bloch(model.S_R, g)
            T = model.sf(g, vectors=True)[1][:, triplet]
            n = model.n_ao
            w_t = np.sum(np.abs(T.conj().T @ Sk @ C_soc[:n]) ** 2 + np.abs(T.conj().T @ Sk @ C_soc[n:]) ** 2, axis=0)
            six = set(np.argsort(-w_t)[:6])
            win = [j for j in range(len(Es)) if j not in six and abs(Es[j] - e8) < 3.5]
            j_s = max(win, key=lambda j: s_w[j])
            out["gamma_band_order"] = dict(sf=float(E[i_s] - e15), soc=float(Es[j_s] - e8),
                                           inverted=bool(E[i_s] < e15),
                                           s_character_sf=float(s_sf[i_s]), s_character_soc=float(s_w[j_s]))
    # conduction-band triplet at the CBM (e.g. Pb 6p at R in the cubic perovskites)
    kc = np.array(out["sf"]["cbm_k"])
    Ec, Cc = model.sf(kc, vectors=True)
    Ec = Ec * HA_EV
    trip_c = list(range(n_occ, n_occ + 3))
    if np.ptp(Ec[trip_c]) < 0.01:
        Esc, Csc = model.soc(kc, vectors=True)
        d_so, lower, upper = _triplet_split(model, kc, trip_c, Esc * HA_EV, Csc)
        out["cbm_triplet"] = dict(k=kc.tolist(), sf=float(Ec[trip_c].mean()), delta_so=d_so,
                                  soc_lower=lower, soc_upper=upper)
    return out


def soc_band_structure(run_dir, basis_file, gth_file, basis_name="DZVP-MOLOPT-PBE-GTH", functional="PBE",
                       bs_file=None, trexio_file=None, out_file=None, projector_radius_ang=16.0, zincblende=False):
    """Spin-free and spin-orbit bands along the CP2K band path, with an analysis of edges and splittings."""
    from qdex.bulk_bands import parse_cp2k_bs

    bs_file = bs_file or glob.glob(os.path.join(run_dir, "*.bs"))[0]
    model = BulkModel(run_dir, basis_file, gth_file, basis_name, functional, trexio_file, out_file,
                      projector_radius_ang)
    bs = parse_cp2k_bs(bs_file)
    sf_segs, soc_segs = [], []
    for seg in bs["segments"]:
        sf_segs.append(np.array([model.sf(kf) for kf in seg["kfrac"]]) * HA_EV)
        soc_segs.append(np.array([model.soc(kf) for kf in seg["kfrac"]]) * HA_EV)

    dev = max(np.abs(s[:, :bs["n_bands"]] - seg["bands"]).max() for s, seg in zip(sf_segs, bs["segments"]))
    # occupied bands from the valence charges (CP2K occupations are fractional in semimetals/smeared runs)
    from qdex.constants import valence_electrons
    n_occ = sum(valence_electrons[e] for e in model.syms) // 2
    kfrac = [s["kfrac"] for s in bs["segments"]]
    info = analyse(model, bs, kfrac, sf_segs, soc_segs, n_occ, zincblende=zincblende)
    info["sf_vs_cp2k_max_dev_ev"] = float(dev)
    logger.info(f"  [Bulk SOC] spin-free bands vs CP2K: max deviation {1000 * dev:.3f} meV; "
                f"gap SF {info['sf']['gap']:.4f} eV, SOC {info['soc']['gap']:.4f} eV")
    return dict(kfrac=kfrac, sf=sf_segs, soc=soc_segs, n_occ=n_occ, n_bands=bs["n_bands"],
                sf_deviation_ev=dev, info=info, model=model)


def write_bs(path, kfrac_segments, band_segments, n_occupied, occupation=1.0, n_bands=None):
    """Bands in the CP2K .bs text format (one set per segment), readable by parse_cp2k_bs; gzipped
    when the path ends in .gz."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as f:
        for s, (kf, bands) in enumerate(zip(kfrac_segments, band_segments), start=1):
            nb = bands.shape[1] if n_bands is None else min(n_bands, bands.shape[1])
            f.write(f"# Set {s}: 2 special points, {len(kf)} k-points, {nb} bands\n")
            f.write(f"#  Special point 1  {kf[0][0]:14.8f} {kf[0][1]:14.8f} {kf[0][2]:14.8f}  not specified\n")
            f.write(f"#  Special point 2  {kf[-1][0]:14.8f} {kf[-1][1]:14.8f} {kf[-1][2]:14.8f}  not specified\n")
            w = 1.0 / len(kf)
            for p, (k, e) in enumerate(zip(kf, bands), start=1):
                f.write(f"#  Point {p}      Spin 1:  {k[0]:14.8f} {k[1]:14.8f} {k[2]:14.8f}  {w:14.8f}\n")
                f.write("#   Band    Energy [eV]     Occupation\n")
                for b in range(nb):
                    f.write(f"  {b + 1:6d} {e[b]:16.8f} {occupation if b < n_occupied else 0.0:14.8f}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m qdex.bulk_soc", description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir", help="CP2K run with *-KS/S_SPIN_1_R_*.csr, the TREXIO .h5, the .bs and cp2k_job.out")
    ap.add_argument("--basis", required=True, help="basis set file of the CP2K run (e.g. BASIS_MOLOPT_UZH)")
    ap.add_argument("--basis-name", default="DZVP-MOLOPT-PBE-GTH")
    ap.add_argument("--gth", required=True, help="GTH_SOC_POTENTIALS file")
    ap.add_argument("--functional", default="PBE", help="functional of the GTH SOC constants (GTH-<functional>-q<n>)")
    ap.add_argument("--out-file", default=None, help="CP2K output with the CSR image table (default RUN_DIR/cp2k_job.out)")
    ap.add_argument("--name", default=None, help="output name (default: the .bs file name), e.g. CdSe_zb")
    ap.add_argument("-d", "--dest", default=None, help="folder for <name>.bs, <name>_soc.bs and <name>.json (default RUN_DIR)")
    ap.add_argument("--cif", default=None, help="CIF of the crystal (labels the band edges, records the space group)")
    ap.add_argument("--n-bands", type=int, default=None,
                    help="spin-free bands to keep (two spinor bands each; default: all of the .bs file)")
    ap.add_argument("--gzip", action="store_true", help="write <name>.bs.gz and <name>_soc.bs.gz")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    import json
    import shutil
    zincblende = False
    if args.cif:
        from pymatgen.core import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        zincblende = SpacegroupAnalyzer(Structure.from_file(args.cif)).get_space_group_symbol() == "F-43m"
    res = soc_band_structure(args.run_dir, args.basis, args.gth, args.basis_name, args.functional,
                             out_file=args.out_file, zincblende=zincblende)
    bs_in = glob.glob(os.path.join(args.run_dir, "*.bs"))[0]
    name = args.name or os.path.basename(bs_in)[:-3]
    dest = args.dest or args.run_dir
    os.makedirs(dest, exist_ok=True)
    n_sf = min(args.n_bands or res["n_bands"], res["n_bands"])
    ext = ".bs.gz" if args.gzip else ".bs"
    write_bs(os.path.join(dest, f"{name}{ext}"), res["kfrac"], res["sf"], n_occupied=res["n_occ"],
             occupation=2.0, n_bands=n_sf)
    write_bs(os.path.join(dest, f"{name}_soc{ext}"), res["kfrac"], res["soc"], n_occupied=2 * res["n_occ"],
             occupation=1.0, n_bands=2 * n_sf)

    info = res["info"]
    model = res["model"]
    info.update(name=name, functional=args.functional, basis=args.basis_name,
                potentials={e: f"GTH-{args.functional}-q{(model.soc_tbl.get(e) or {}).get('name', '').split('-q')[-1]}"
                            for e in sorted(set(model.syms))},
                cell_ang=np.asarray(model.cell).round(6).tolist(), n_bands_sf=n_sf)
    if args.cif:
        from pymatgen.core import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        from pymatgen.symmetry.bandstructure import HighSymmKpath
        st = Structure.from_file(args.cif)
        sga = SpacegroupAnalyzer(st)
        info["spacegroup"] = sga.get_space_group_symbol()
        info["formula"] = st.composition.reduced_formula
        info["cif"] = os.path.basename(args.cif)
        special = HighSymmKpath(sga.get_primitive_standard_structure()).kpath["kpoints"]
        for key in ("sf", "soc"):
            for edge in ("vbm", "cbm"):
                kk = np.array(info[key][f"{edge}_k"])
                lab = [l.replace("\\Gamma", "Gamma") for l, v in special.items() if np.allclose(v, kk, atol=1e-4)]
                info[key][f"{edge}_label"] = lab[0] if lab else "path"
    with open(os.path.join(dest, f"{name}.json"), "w") as f:
        json.dump(info, f, indent=1)
    logger.info(f"  [Bulk SOC] wrote {name}{ext}, {name}_soc{ext}, {name}.json in {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
