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
        [--basis-name DZVP-MOLOPT-PBE-GTH] [--functional PBE] [-o CdSe_bulk_soc.bs]
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
    i = next(n for n, l in enumerate(lines) if tag in l and "periodic images" in l)
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


def soc_band_structure(run_dir, basis_file, gth_file, basis_name="DZVP-MOLOPT-PBE-GTH", functional="PBE",
                       bs_file=None, trexio_file=None, out_file=None, projector_radius_ang=16.0):
    """Spin-free and spin-orbit bands along the CP2K band path. Returns a dict (eV)."""
    from qdex.bulk_bands import parse_cp2k_bs
    from qdex.constants import valence_electrons
    from qdex.io_utils import build_shell_dicts, parse_basis, parse_gth_soc_potentials

    bs_file = bs_file or glob.glob(os.path.join(run_dir, "*.bs"))[0]
    trexio_file = trexio_file or glob.glob(os.path.join(run_dir, "*.h5"))[0]
    out_file = out_file or os.path.join(run_dir, "cp2k_job.out")

    cell, syms, coords = read_trexio_structure(trexio_file)
    basis = parse_basis(basis_file, basis_name, required_elements=set(syms))
    shells = [{**s, "pure": True} for s in build_shell_dicts(syms, coords, basis)]
    n_ao = sum(2 * s["l"] + 1 for s in shells)
    H_R = read_real_space_matrices(run_dir, "KS", n_ao, out_file)
    S_R = read_real_space_matrices(run_dir, "S", n_ao, out_file)

    import libint_cpp
    dS = np.abs(libint_cpp.overlap(shells, 1) - S_R[(0, 0, 0)]).max()
    logger.info(f"  [Bulk SOC] {len(H_R)} real-space images; S(0) libint vs CP2K: {dS:.1e}")
    if dS > 1e-5:
        raise ValueError(f"AO basis does not match the CP2K run (max |dS| = {dS:.2e}); check basis file/name.")

    soc_tbl = parse_gth_soc_potentials(gth_file, {s: valence_electrons.get(s) for s in set(syms)}, functional=functional)
    A_L, K = _projector_matrices(shells, syms, coords, cell, soc_tbl, projector_radius_ang)

    HA = 27.211386245988
    bs = parse_cp2k_bs(bs_file)
    sf_segs, soc_segs = [], []
    for seg in bs["segments"]:
        sf, so = [], []
        for kf in seg["kfrac"]:
            Hk, Sk = _bloch(H_R, kf), _bloch(S_R, kf)
            sf.append(eigh(Hk, Sk, eigvals_only=True))
            Ak = sum(np.exp(2j * np.pi * np.dot(kf, c)) * M for c, M in A_L.items())
            Hx, Hy, Hz = (0.5 * Ak @ Kk @ Ak.conj().T for Kk in K)
            Z = np.zeros_like(Sk)
            Htot = np.block([[Hk + Hz, Hx - 1j * Hy], [Hx + 1j * Hy, Hk - Hz]])
            so.append(eigh(0.5 * (Htot + Htot.conj().T), np.block([[Sk, Z], [Z, Sk]]), eigvals_only=True))
        sf_segs.append(np.array(sf) * HA)
        soc_segs.append(np.array(so) * HA)

    dev = max(np.abs(s[:, :bs["n_bands"]] - seg["bands"]).max() for s, seg in zip(sf_segs, bs["segments"]))
    n_occ = bs["vbm_band"] + 1
    soc_G = soc_segs[0][0]
    d_so = soc_G[2 * n_occ - 1] - soc_G[2 * n_occ - 6]
    logger.info(f"  [Bulk SOC] spin-free bands vs CP2K: max deviation {1000 * dev:.3f} meV; "
                f"Delta_so at the first k-point {d_so:.4f} eV")
    return dict(kfrac=[s["kfrac"] for s in bs["segments"]], sf=sf_segs, soc=soc_segs, n_occ=n_occ,
                n_bands=bs["n_bands"], sf_deviation_ev=dev, delta_so_first_k=d_so)


def write_bs(path, kfrac_segments, band_segments, n_occupied, occupation=1.0, n_bands=None):
    """Bands in the CP2K .bs text format (one set per segment), readable by parse_cp2k_bs."""
    with open(path, "w") as f:
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
    ap.add_argument("-o", "--output", default=None, help="SOC .bs file to write (default RUN_DIR/<name>_soc.bs)")
    ap.add_argument("--n-bands", type=int, default=None,
                    help="spin-free bands to keep (two spinor bands each; default: all of the .bs file)")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    res = soc_band_structure(args.run_dir, args.basis, args.gth, args.basis_name, args.functional,
                             out_file=args.out_file)
    bs_in = glob.glob(os.path.join(args.run_dir, "*.bs"))[0]
    out = args.output or bs_in.replace(".bs", "_soc.bs")
    # two spinor bands per spin-free band
    n_sf = args.n_bands or res["n_bands"]
    write_bs(out, res["kfrac"], res["soc"], n_occupied=2 * res["n_occ"], occupation=1.0, n_bands=2 * n_sf)
    logger.info(f"  [Bulk SOC] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
