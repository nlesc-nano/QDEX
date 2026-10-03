"""
xtb molden file -> per-atom basis file + MOs.mbse in the QDEX spherical AO convention.

xtb writes the basis per atom (the g-xTB q-vSZP contractions depend on the atomic charge)
and the MO coefficients over Cartesian d/f functions. The conversion:

1. rescales every contraction to unit norm, because libint renormalizes the shells it builds;
2. maps the Cartesian components onto libint's Cartesian order and normalization;
3. projects onto spherical functions with libint overlaps, C_sph = S_ss^-1 S_sc C_cart,
   which is exact when the orbitals contain no s/p contamination of the Cartesian shells.

Orthonormality, C^T S C = 1 in the libint metric, is the acceptance test of the conversion.
"""
import os
import json

import numpy as np

from qdex.constants import BOHR_PER_ANG
from qdex.io_utils import write_mos_mbse, write_xyz
import libint_cpp
import logging

logger = logging.getLogger(__name__)

L_OF = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}
L_NAME = "spdfg"

# Molden component order of Cartesian shells, as (lx, ly, lz)
MOLDEN_CART = {
    0: [(0, 0, 0)],
    1: [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    2: [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    3: [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0), (2, 0, 1),
        (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)],
}
# Molden order of spherical shells, as m; libint (standard) runs m = -l..l
MOLDEN_PURE_M = {2: [0, 1, -1, 2, -2], 3: [0, 1, -1, 2, -2, 3, -3]}


def libint_cart(l):
    """libint Cartesian component order: xx, xy, xz, yy, yz, zz for d."""
    return [(lx, l - lx - lz, lz) for lx in range(l, -1, -1) for lz in range(0, l - lx + 1)]


def _dfact(n):
    return 1.0 if n <= 0 else float(np.prod(np.arange(n, 0, -2)))


def cart_rel_norm(lx, ly, lz):
    """<x^lx y^ly z^lz|same> relative to the axis component x^l of the same shell."""
    l = lx + ly + lz
    return _dfact(2 * lx - 1) * _dfact(2 * ly - 1) * _dfact(2 * lz - 1) / _dfact(2 * l - 1)


def contraction_norm(l, exps, coefs):
    """Norm of sum_i c_i g_i over normalized primitives g_i of angular momentum l."""
    a = np.asarray(exps)
    c = np.asarray(coefs)
    s = (2.0 * np.sqrt(np.outer(a, a)) / np.add.outer(a, a)) ** (l + 1.5)
    return float(np.sqrt(c @ s @ c))


def primitive_norm(l, exps):
    """Normalization constant of the axis Cartesian primitive x^l exp(-a r^2)."""
    a = np.asarray(exps, dtype=float)
    return np.sqrt((2.0 * a / np.pi) ** 1.5 * (4.0 * a) ** l / _dfact(2 * l - 1))


def to_normalized_primitives(atom_shells, convention):
    """
    Contraction coefficients over normalized primitives.

    ``convention`` is how the file gives them: 'normalized' (molden standard) or 'raw'
    (over unnormalized primitives x^l exp(-a r^2), as xtb writes them).
    """
    if convention == "normalized":
        return atom_shells
    return [[(l, e, c / primitive_norm(l, e)) for l, e, c in shells] for shells in atom_shells]


def parse_molden(path):
    """
    Parse a molden file.

    Returns dict with ``syms``, ``coords_bohr`` (natoms x 3), ``atom_shells`` (per atom a list of
    (l, exps, coefs)), ``pure`` ({l: bool}), ``eps`` (Eh), ``occ`` and ``C`` (n_bf x n_mo, molden order).
    """
    with open(path) as f:
        lines = f.read().splitlines()
    text_flags = {l.strip().upper() for l in lines if l.strip().startswith("[")}
    pure = {0: False, 1: False,
            2: bool(text_flags & {"[5D]", "[5D7F]", "[5D10F]"}),
            3: bool(text_flags & {"[7F]", "[5D7F]"}),
            4: "[9G]" in text_flags}

    def section(name):
        for i, l in enumerate(lines):
            if l.strip().upper().startswith(name):
                return i
        raise ValueError(f"{path}: no {name} section")

    i_at = section("[ATOMS]")
    unit = BOHR_PER_ANG if "ANG" in lines[i_at].upper() else 1.0
    syms, coords = [], []
    i = i_at + 1
    while i < len(lines) and not lines[i].strip().startswith("["):
        t = lines[i].split()
        if len(t) >= 6:
            syms.append(t[0].capitalize())
            coords.append([float(v) * unit for v in t[3:6]])
        i += 1

    atom_shells = [[] for _ in syms]
    i = section("[GTO]") + 1
    while i < len(lines) and not lines[i].strip().startswith("["):
        t = lines[i].split()
        if len(t) == 2 and t[0].isdigit():
            ia = int(t[0]) - 1
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith("["):
                h = lines[i].split()
                lname, nprim = h[0].lower(), int(h[1])
                if lname not in L_OF:
                    raise ValueError(f"{path}: unsupported shell type '{h[0]}'")
                scale = float(h[2].replace("D", "E")) if len(h) > 2 else 1.0
                prim = np.array([[float(v.replace("D", "E")) for v in lines[i + 1 + k].split()[:2]]
                                 for k in range(nprim)])
                atom_shells[ia].append((L_OF[lname], prim[:, 0] * scale ** 2, prim[:, 1]))
                i += nprim + 1
        else:
            i += 1

    n_bf = sum((2 * l + 1) if pure[l] else (l + 1) * (l + 2) // 2 for sh in atom_shells for l, _, _ in sh)
    eps, occ, cols, spins = [], [], [], []
    i = section("[MO]") + 1
    current = None
    while i < len(lines) and not lines[i].strip().startswith("["):
        s = lines[i].strip()
        if "=" in s:
            key, val = [x.strip() for x in s.split("=", 1)]
            if key.lower() in ("sym", "ene") and (current is None or current["coef_seen"]):
                current = {"coef_seen": False, "c": np.zeros(n_bf)}
                cols.append(current)
            if key.lower() == "ene":
                eps.append(float(val.replace("D", "E")))
            elif key.lower() == "occup":
                occ.append(float(val.replace("D", "E")))
            elif key.lower() == "spin":
                spins.append(val.lower())
        elif s:
            t = s.split()
            current["c"][int(t[0]) - 1] = float(t[1].replace("D", "E"))
            current["coef_seen"] = True
        i += 1
    C = np.column_stack([c["c"] for c in cols])
    eps, occ = np.array(eps), np.array(occ)
    if not (len(eps) == len(occ) == C.shape[1]):
        raise ValueError(f"{path}: inconsistent MO block ({len(eps)} energies, {len(occ)} occupations, "
                         f"{C.shape[1]} vectors)")
    spins = np.array(spins) if len(spins) == len(eps) else np.array(["alpha"] * len(eps))
    if np.any(spins == "beta"):
        # Spin-unrestricted output of a closed shell: identical alpha and beta sets are merged.
        a, b = spins == "alpha", spins == "beta"
        same = (a.sum() == b.sum() and np.allclose(eps[a], eps[b], atol=1e-6)
                and np.allclose(occ[a], occ[b], atol=1e-6))
        if not same:
            raise ValueError(f"{path}: open-shell molden files are not supported")
        eps, occ, C = eps[a], occ[a] + occ[b], C[:, a]
    return dict(syms=syms, coords_bohr=np.array(coords), atom_shells=atom_shells, pure=pure,
                eps=eps, occ=occ, C=C)


def shells_from_atoms(syms, coords_bohr, atom_shells, pure=True):
    """QDEX/libint shell dicts from per-atom (l, exps, coefs) lists."""
    # libint_cpp reads the raw buffers: the arrays must be contiguous, not strided views
    shells = []
    for ia, (sym, xyz) in enumerate(zip(syms, coords_bohr)):
        for l, exps, coefs in atom_shells[ia]:
            shells.append(dict(sym=sym, atom_idx=ia, l=int(l), exps=np.ascontiguousarray(exps, dtype=float),
                               coefs=np.ascontiguousarray(coefs, dtype=float),
                               center=np.ascontiguousarray(xyz, dtype=float), pure=bool(pure)))
    return shells


def _cart_to_libint(mol, convention):
    """Molden Cartesian coefficients -> libint Cartesian basis (unit contractions, libint order and norm)."""
    C = mol["C"]
    rows, scale = [], []
    i0 = 0
    for shells in mol["atom_shells"]:
        for l, exps, coefs in shells:
            if mol["pure"][l]:
                raise ValueError("mixed spherical/Cartesian molden shells are not supported here")
            order = MOLDEN_CART[l]
            n_c = contraction_norm(l, exps, coefs)
            pos = {comp: k for k, comp in enumerate(order)}
            for comp in libint_cart(l):
                f = n_c / np.sqrt(cart_rel_norm(*comp)) if convention == "component" else n_c
                rows.append(i0 + pos[comp])
                scale.append(f)
            i0 += len(order)
    return C[rows, :] * np.asarray(scale)[:, None]


def _pure_to_libint(mol):
    """Molden spherical coefficients -> libint spherical order (m = -l..l), unit contractions."""
    C = mol["C"]
    rows, scale = [], []
    i0 = 0
    for shells in mol["atom_shells"]:
        for l, exps, coefs in shells:
            n_c = contraction_norm(l, exps, coefs)
            if l == 1:
                order = [1, -1, 0]                   # molden x, y, z -> m = +1, -1, 0
            elif l == 0:
                order = [0]
            else:
                order = MOLDEN_PURE_M[l]
            pos = {m: k for k, m in enumerate(order)}
            for m in range(-l, l + 1):
                rows.append(i0 + pos[m])
                scale.append(n_c)
            i0 += len(order)
    return C[rows, :] * np.asarray(scale)[:, None]


def molden_to_spherical(mol, nthreads=1):
    """
    MO coefficients of a parsed molden file in the QDEX spherical AO basis.

    The primitive convention ('normalized' or 'raw') and, for Cartesian files, the component
    normalization ('component' or 'axis') are detected as the combination that makes the
    orbitals orthonormal: xtb's GFN2 writer uses raw/component, the g-xTB (tblite) writer raw/axis.

    Returns (atom_shells, shells_sph, C_sph, info). ``atom_shells`` are over normalized
    primitives; ``info`` holds the detected conventions and orthonormality errors.
    """
    syms, X = mol["syms"], mol["coords_bohr"]
    cartesian = not all(mol["pure"][l] for sh in mol["atom_shells"] for l, _, _ in sh if l >= 2)
    best = None
    for prim in ("raw", "normalized"):
        atom_shells = to_normalized_primitives(mol["atom_shells"], prim)
        m = dict(mol, atom_shells=atom_shells)
        if cartesian:
            shells = shells_from_atoms(syms, X, atom_shells, pure=False)
            candidates = [(conv, _cart_to_libint(m, conv)) for conv in ("component", "axis")]
        else:
            shells = shells_from_atoms(syms, X, atom_shells, pure=True)
            candidates = [(None, _pure_to_libint(m))]
        S = libint_cpp.cross_overlap_geometries(shells, shells, nthreads)
        for conv, C in candidates:
            err = float(np.max(np.abs(C.T @ S @ C - np.eye(C.shape[1]))))
            if best is None or err < best["err"]:
                best = dict(prim=prim, conv=conv, err=err, C=C, atom_shells=atom_shells, shells=shells)

    atom_shells = best["atom_shells"]
    shells_sph = shells_from_atoms(syms, X, atom_shells, pure=True)
    S_ss = libint_cpp.cross_overlap_geometries(shells_sph, shells_sph, nthreads)
    if cartesian:
        S_sc = libint_cpp.cross_overlap_geometries(shells_sph, best["shells"], nthreads)
        C_sph = np.linalg.solve(S_ss, S_sc @ best["C"])
    else:
        C_sph = best["C"]
    M = C_sph.T @ S_ss @ C_sph
    info = dict(primitive_convention=best["prim"], cartesian_convention=best["conv"],
                orthonormality_error_input=best["err"],
                orthonormality_error=float(np.max(np.abs(M - np.eye(M.shape[0])))),
                min_norm=float(np.min(np.diag(M))))
    return atom_shells, shells_sph, C_sph, info


def write_atom_basis(path, syms, atom_shells, title="g-xTB basis (per atom, from molden)"):
    """
    CP2K-style basis file with one block per atom, named ``GXTB-A<index>`` (1-based).

    Coefficients refer to normalized primitives; the contraction need not be normalized
    (libint renormalizes it). Every shell is its own set.
    """
    with open(path, "w") as f:
        f.write(f"# {title}\n")
        for ia, (sym, shells) in enumerate(zip(syms, atom_shells)):
            f.write(f"{sym} GXTB-A{ia + 1}\n  {len(shells)}\n")
            for l, exps, coefs in shells:
                f.write(f"  1 {l} {l} {len(exps)} 1\n")
                for a, c in zip(exps, coefs):
                    f.write(f"    {a:24.16e} {c:24.16e}\n")


def read_atom_basis(path):
    """Inverse of ``write_atom_basis``: (syms, atom_shells)."""
    with open(path) as f:
        lines = [l.split() for l in f if l.strip() and not l.lstrip().startswith("#")]
    syms, atom_shells = [], []
    i = 0
    while i < len(lines):
        sym, name = lines[i][0], lines[i][1]
        if not name.startswith("GXTB-A") or int(name[6:]) != len(syms) + 1:
            raise ValueError(f"{path}: expected block GXTB-A{len(syms) + 1}, found '{name}'")
        nset = int(lines[i + 1][0])
        i += 2
        shells = []
        for _ in range(nset):
            l, nexp = int(lines[i][1]), int(lines[i][3])
            prim = np.array([[float(v) for v in lines[i + 1 + k][:2]] for k in range(nexp)])
            shells.append((l, prim[:, 0], prim[:, 1]))
            i += nexp + 1
        syms.append(sym)
        atom_shells.append(shells)
    return syms, atom_shells


def build_frame_shells(basis_path, syms, coords_ang):
    """Shell dicts of one frame from its per-atom basis file (the build_shell_dicts equivalent)."""
    b_syms, atom_shells = read_atom_basis(basis_path)
    if list(b_syms) != list(syms):
        raise ValueError(f"{basis_path}: atoms do not match the frame geometry")
    return shells_from_atoms(syms, np.asarray(coords_ang) * BOHR_PER_ANG, atom_shells, pure=True)


def convert_molden(molden_path, out_dir, basis_name="BASIS_GXTB", mo_name="MOs.mbse", xyz_name="frame.xyz",
                   nthreads=1, tol=1e-4, comment=""):
    """
    Write ``frame.xyz``, the per-atom basis file and ``MOs.mbse`` (spherical, energies in Eh) to ``out_dir``.

    Raises if the converted orbitals are not orthonormal to ``tol``. Returns the diagnostics,
    which are also saved as ``conversion.json``.
    """
    mol = parse_molden(molden_path)
    atom_shells, _shells, C_sph, info = molden_to_spherical(mol, nthreads=nthreads)
    if info["orthonormality_error"] > tol:
        raise ValueError(f"{molden_path}: converted MOs are not orthonormal "
                         f"(max |C^T S C - 1| = {info['orthonormality_error']:.2e} > {tol:g})")
    os.makedirs(out_dir, exist_ok=True)
    write_xyz(os.path.join(out_dir, xyz_name), mol["syms"], mol["coords_bohr"] / BOHR_PER_ANG, comment)
    write_atom_basis(os.path.join(out_dir, basis_name), mol["syms"], atom_shells)
    write_mos_mbse(os.path.join(out_dir, mo_name), C_sph, mol["eps"], mol["occ"])
    n_occ = int(np.sum(mol["occ"] > 0.5))
    info.update(n_ao=int(C_sph.shape[0]), n_mo=int(C_sph.shape[1]), n_occ=n_occ,
                homo_ev=float(mol["eps"][n_occ - 1] * 27.211386245988),
                lumo_ev=float(mol["eps"][n_occ] * 27.211386245988))
    info["gap_ev"] = info["lumo_ev"] - info["homo_ev"]
    with open(os.path.join(out_dir, "conversion.json"), "w") as f:
        json.dump(info, f, indent=1)
    return info
