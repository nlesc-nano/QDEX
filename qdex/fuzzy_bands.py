import json
import numpy as np
import time
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from scipy.spatial import cKDTree
import logging

logger = logging.getLogger(__name__)


def _kabsch(ref, target):
    """Proper rotation R minimizing sum |R ref_i - target_i|^2 (rows are vectors)."""
    H = ref.T @ target
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)) or 1.0
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def _bond_match(units, R, cif_dirs):
    """Mean cosine between each QD bond direction and its closest rotated CIF bond direction."""
    return float(np.mean(np.max(units @ (cif_dirs @ R.T).T, axis=1)))


def _refine_rotation(units, R, cif_dirs, n_iter=8):
    for _ in range(n_iter):
        assign = np.argmax(units @ (cif_dirs @ R.T).T, axis=1)
        R = _kabsch(cif_dirs[assign], units)
    return R


def cif_bond_star(prim_struct):
    """Nearest-neighbour bond directions of the crystal (both directions of every bond),
    the bond length, the coordination of each element and the bonded element pairs."""
    dmin = min(n.nn_distance for site in prim_struct for n in prim_struct.get_neighbors(site, 6.0))
    dirs, coordination, pairs = [], {}, set()
    for site in prim_struct:
        nbrs = prim_struct.get_neighbors(site, 1.15 * dmin)
        coordination[site.specie.symbol] = len(nbrs)
        for n in nbrs:
            v = n.coords - site.coords
            dirs.append(v / np.linalg.norm(v))
            pairs.add((site.specie.symbol, n.specie.symbol))
    dirs = np.unique(np.round(np.array(dirs), 6), axis=0)
    return dirs, float(dmin), coordination, pairs


def _sublattice_element(prim_struct, coordination):
    """Most coordinated element and its nearest same-element distance in the crystal.

    Its sublattice carries the lattice constant: the B site of a perovskite (Pb-Pb = a_pc),
    the cation of zinc blende, wurtzite or rock salt. Octahedral tilts shorten the B-B distance
    but not the B-X bond, so the bond length overestimates the lattice of a tilted dot.
    """
    el = max(sorted(coordination), key=lambda e: coordination[e])
    d = min(n.nn_distance for site in prim_struct if site.specie.symbol == el
            for n in prim_struct.get_neighbors(site, 12.0) if n.specie.symbol == el)
    return el, float(d)


def fit_lattice_orientation(prim_struct, coords_ang, syms=None, interior_fraction=0.7):
    """Rotation R (crystal frame -> dot frame) and isotropic scale of the dot's lattice.

    The interior nearest-neighbour bond directions of the dot are matched to the crystal's
    bond star (Kabsch fit with re-assignment, started from every pair of crystal bonds whose
    angle matches a pair of dot bonds). R is fixed up to the symmetry of the bond star, which
    leaves the fuzzy weights ``|<phi|k>|^2`` unchanged (crystal point group plus k -> -k).
    """
    coords = np.asarray(coords_ang, dtype=float)
    cif_dirs, cif_bond, coordination, pairs = cif_bond_star(prim_struct)
    elements = set(coordination)
    mask = np.ones(len(coords), bool) if syms is None else np.array([s in elements for s in syms])
    X = coords[mask]
    xs = None if syms is None else [s for s, m in zip(syms, mask) if m]
    info = dict(rotation=np.eye(3), scale=1.0, match=float("nan"), n_bonds=0, cif_bond_ang=cif_bond,
                qd_bond_ang=float("nan"), coordination=coordination, elements=sorted(elements))
    if len(X) < 2:
        logger.warning("  [Fuzzy] No atoms of the crystal's elements in the dot: k-path left in the CIF frame.")
        return info

    r = np.linalg.norm(X - X.mean(axis=0), axis=1)
    interior = np.where(r <= interior_fraction * r.max())[0]
    if len(interior) < 4:
        interior = np.arange(len(X))
    tree = cKDTree(X)
    vecs, ref_atom = [], None
    for i in interior[np.argsort(r[interior])]:
        nbrs = [j for j in tree.query_ball_point(X[i], 1.3 * cif_bond)
                if j != i and (xs is None or (xs[i], xs[j]) in pairs)]
        if ref_atom is None and len(nbrs) >= 2:
            ref_atom = len(vecs)
        vecs.extend(X[j] - X[i] for j in nbrs)
    if len(vecs) < 2 or ref_atom is None:
        logger.warning("  [Fuzzy] No interior bonds found: k-path left in the CIF frame.")
        return info
    vecs = np.array(vecs)
    dists = np.linalg.norm(vecs, axis=1)
    units = vecs / dists[:, None]
    qd_bond = float(np.median(dists))

    # Seeds: map every pair of crystal bonds onto two bonds of the most central atom
    b1, b2 = units[ref_atom], units[ref_atom + 1]
    ang_b = np.arccos(np.clip(b1 @ b2, -1, 1))
    seeds = [np.eye(3)]
    for i, c1 in enumerate(cif_dirs):
        for j, c2 in enumerate(cif_dirs):
            if i != j and abs(np.arccos(np.clip(c1 @ c2, -1, 1)) - ang_b) < np.radians(15):
                seeds.append(_kabsch(np.array([c1, c2]), np.array([b1, b2])))
    best_R, best = np.eye(3), -np.inf
    for R0 in seeds:
        R = _refine_rotation(units, R0, cif_dirs)
        m = _bond_match(units, R, cif_dirs)
        if m > best + 1e-9:
            best_R, best = R, m

    # Lattice scale from the sublattice of the most coordinated element (B-B for perovskites)
    scale, how = qd_bond / cif_bond, f"median bond {qd_bond:.4f} A vs CIF {cif_bond:.4f} A"
    if xs is not None:
        el, d_cif = _sublattice_element(prim_struct, coordination)
        sub = np.array([k for k, e in enumerate(xs) if e == el])
        if len(sub) >= 2:
            inner = sub[r[sub] <= interior_fraction * r.max()]
            # every sublattice neighbour (six B-B in a perovskite), not just the nearest: in a relaxed
            # dot the B-B distances spread by several per cent and the nearest one is biased short
            k = min(13, len(sub))
            dd, _ = cKDTree(X[sub]).query(X[inner] if len(inner) else X[sub], k=k)
            dd = dd[:, 1:].ravel()
            dd = dd[(dd > 0.7 * d_cif) & (dd < 1.25 * d_cif)]
            if len(dd) >= 1:
                d_qd = float(np.median(dd))
                scale, how = d_qd / d_cif, f"median {el}-{el} {d_qd:.4f} A vs CIF {d_cif:.4f} A"
                info.update(sublattice_element=el, qd_sublattice_ang=d_qd, cif_sublattice_ang=d_cif)
    info.update(rotation=best_R, scale=scale, match=best, n_bonds=len(units), qd_bond_ang=qd_bond, bond_pairs=pairs)
    angle = np.degrees(np.arccos(np.clip((np.trace(best_R) - 1) / 2, -1, 1)))
    logger.info(f"  [Fuzzy] Lattice orientation: {len(units)} interior bonds fitted to the crystal bond star "
                f"(mean cos {best:.4f}, rotation {angle:.1f} deg from the CIF frame).")
    logger.info(f"  [Fuzzy] Lattice scale: {how} -> k-points scaled by 1/{scale:.4f} "
                f"(median bond {qd_bond:.4f} A vs CIF {cif_bond:.4f} A).")
    if abs(qd_bond / cif_bond - 1) > 0.08:
        logger.warning(f"  [Fuzzy] The dot's bonds are {100 * (qd_bond / cif_bond - 1):+.1f}% off the CIF: "
                       f"check that the CIF is the dot's material.")
    if best < 0.95:
        logger.warning(f"  [Fuzzy] Weak lattice match (mean cos {best:.3f}): the dot has no well-ordered core "
                       f"of the CIF structure, so the k-path directions are uncertain.")
    return info


def generate_automated_kpath(cif_path, coords_ang, line_density=50, return_reciprocal=False, syms=None,
                             return_info=False):
    """High-symmetry k-path of the CIF crystal, expressed in the dot's Cartesian frame (1/A).

    The path and the reciprocal lattice are rotated onto the dot's lattice orientation and scaled
    to its bond length (fit_lattice_orientation). With return_info the fit and the fractional
    k-points (in the primitive reciprocal basis) are returned as well.
    """
    logger.info(f"  [Fuzzy] Loading CIF: {cif_path}")
    struct = Structure.from_file(cif_path)
    sga = SpacegroupAnalyzer(struct)
    prim_struct = sga.get_primitive_standard_structure()
    kpath = HighSymmKpath(prim_struct)
    kpts_frac, labels = kpath.get_kpoints(line_density=line_density, coords_are_cartesian=False)
    kpts_frac = np.asarray(kpts_frac, dtype=float)

    info = fit_lattice_orientation(prim_struct, coords_ang, syms=syms)
    R, scale = info["rotation"], info["scale"]
    reciprocal_matrix = prim_struct.lattice.reciprocal_lattice.matrix @ R.T / scale
    kpts_cart = kpts_frac @ reciprocal_matrix
    info["kpts_frac"] = kpts_frac

    logger.info(f"  [Fuzzy] Generated {len(kpts_cart)} k-points for spacegroup {sga.get_space_group_symbol()}.")
    out = (kpts_cart, labels)
    if return_reciprocal:
        out += (reciprocal_matrix,)
    if return_info:
        out += (info,)
    return out


def make_reciprocal_replicas(reciprocal_matrix, g_shell):
    coeffs = np.arange(-g_shell, g_shell + 1, dtype=int)
    hkl = np.array(np.meshgrid(coeffs, coeffs, coeffs, indexing="ij")).reshape(3, -1).T
    return hkl @ reciprocal_matrix


def fuzzy_energy_mask(eps_dft, ewin, sigma_ev, qp_energies=None):
    margin = 4.0 * sigma_ev
    mask = (eps_dft >= ewin[0] - margin) & (eps_dft <= ewin[1] + margin)
    if qp_energies is not None:
        mask |= (qp_energies >= ewin[0] - margin) & (qp_energies <= ewin[1] + margin)
    return mask


def fuzzy_energy_indices(eps_dft, ewin, sigma_ev, homo_index=None, qp_energies=None):
    mask = fuzzy_energy_mask(eps_dft, ewin, sigma_ev, qp_energies=qp_energies)
    if homo_index is not None:
        if 0 <= homo_index < len(mask):
            mask[homo_index] = True
        if 0 <= homo_index + 1 < len(mask):
            mask[homo_index + 1] = True
    return np.where(mask)[0].astype(int)


def folded_plane_wave_weights(shells, kpts_cart, G_vecs, nthreads, coeff_blocks, spinor_parts=None):
    """Folded plane-wave weights ``sum_G |<k+G|psi>|^2`` along kpts_cart (1/A).

    coeff_blocks: real MO coefficients (n_ao, n_b), one array per set of orbitals.
    spinor_parts: optional [(block, rows, U_part)], U_part (len(rows), n_spinor): the spin
    components ``sum_r U_part[r, n] phi_{block, rows[r]}`` of two-component spinors. The
    components are orthogonal in spin, so their weights add: ``|U_a^T F|^2 + |U_b^T F|^2``.

    The G replicas are transformed in batches (n_AO x n_batch n_k plane waves, real and imaginary
    parts side by side) and projected onto every block with one real GEMM per batch: one replica
    at a time gives GEMMs only n_k (~200) columns wide, which ran at a fraction of the BLAS peak
    (Cs324Pb216Br756: 60 s). All k+G points at once would need n_AO x n_k x n_G complex numbers
    (~25 GB for 2,000 atoms at g_shell 2); a batch is kept below ~1.5e8 numbers per array.
    Returns ([W_block (n_b, n_k)], W_spinor (n_spinor, n_k) or None).
    """
    import libint_cpp

    n_k = len(kpts_cart)
    G_vecs = np.atleast_2d(G_vecs)
    CT = [np.ascontiguousarray(np.asarray(C, dtype=float).T) for C in coeff_blocks]
    W = [np.zeros((C.shape[0], n_k)) for C in CT]
    parts = [(b, np.asarray(rows, dtype=int), np.ascontiguousarray(np.asarray(U).T))
             for b, rows, U in (spinor_parts or [])]
    W_spin = np.zeros((parts[0][2].shape[0], n_k)) if parts else None
    need_amp = {b for b, _, _ in parts}
    n_rows = max([CT[0].shape[1]] + [C.shape[0] for C in CT] + [p[2].shape[0] for p in parts])
    n_batch = int(np.clip(1.5e8 // (2 * n_rows * n_k), 1, len(G_vecs)))
    for g0 in range(0, len(G_vecs), n_batch):
        Gb = G_vecs[g0:g0 + n_batch]
        nb = len(Gb)
        kq = (kpts_cart[None, :, :] + Gb[:, None, :]).reshape(-1, 3)
        F_ao = libint_cpp.ao_ft_complex(shells, kq / 1.8897259886, nthreads)     # (n_AO, nb n_k)
        F_ri = np.empty((F_ao.shape[0], 2 * nb * n_k))
        F_ri[:, :nb * n_k] = F_ao.real
        F_ri[:, nb * n_k:] = F_ao.imag
        del F_ao
        amp = {}
        for b, C in enumerate(CT):
            RI = C @ F_ri
            re, im = RI[:, :nb * n_k], RI[:, nb * n_k:]
            W[b] += (re * re + im * im).reshape(-1, nb, n_k).sum(axis=1)
            if b in need_amp:
                amp[b] = re + 1j * im
            del RI
        if parts:
            for b, rows, UT in parts:
                A = UT @ amp[b][rows]           # <k+G|spin component of each spinor>
                W_spin += (A.real ** 2 + A.imag ** 2).reshape(-1, nb, n_k).sum(axis=1)
    return W, W_spin


def compute_fuzzy_intensity(C_dense, shells, kpts_cart, nthreads, fold_to_bz=False, g_shell=0, reciprocal_matrix=None, mo_indices=None):
    if g_shell < 0:
        raise ValueError("g_shell must be >= 0")

    C_project = C_dense if mo_indices is None else C_dense[:, mo_indices]
    if not fold_to_bz:
        G_vecs = np.zeros((1, 3))
    else:
        if reciprocal_matrix is None:
            raise ValueError("reciprocal_matrix is required when fold_to_bz=True")
        # BZ-folded spectral projection of finite QD molecular orbitals: sums reciprocal
        # replicas of the same finite-MO Fourier amplitude (not a true Bloch band structure).
        G_vecs = make_reciprocal_replicas(reciprocal_matrix, g_shell)
    W, _ = folded_plane_wave_weights(shells, kpts_cart, G_vecs, nthreads, [C_project])
    return W[0]


def classify_edges(energies, occupied, labels):
    """Band edges from the angular-coverage classes (qdex.angular: S-like, P-like, D-like, Facet, Localized), the
    same rule as the population table. Counting inward from the gap, the first S/P/D-like state is the delocalized
    HOMO (LUMO); the states passed on the way are flagged 1 (Localized) or 4 (Facet).
    Returns dict(flag, homo, lumo, deloc_homo, deloc_lumo, localized), flag 0 other, 1 localized, 2 delocalized
    HOMO, 3 delocalized LUMO, 4 facet; indices into ``energies`` (None if not found); localized marks every state
    that is not S/P/D-like."""
    from qdex.angular import BAND_CLASSES
    E = np.asarray(energies, dtype=float)
    occ = np.asarray(occupied, dtype=bool)
    lab = np.asarray([str(x) for x in labels])
    band = np.isin(lab, BAND_CLASSES)
    res = {"homo": None, "lumo": None, "deloc_homo": None, "deloc_lumo": None}
    flag = np.zeros(len(E), dtype=np.int8)
    for side, key, order in ((occ, "homo", np.argsort(-E)), (~occ, "lumo", np.argsort(E))):
        idx = [int(i) for i in order if side[i]]
        if not idx:
            continue
        res[key] = idx[0]
        for i in idx:
            if band[i]:
                res["deloc_" + key] = i
                flag[i] = 2 if key == "homo" else 3
                break
            flag[i] = 1 if lab[i] == "Localized" else 4
    res["flag"] = flag
    res["localized"] = ~band
    return res


def edge_summary(prefix, energies, res):
    """Log line [Traps:SF|SOC] of the delocalized band edges and the states between them and the nominal edges."""
    E = np.asarray(energies, dtype=float)
    if res["deloc_homo"] is None or res["deloc_lumo"] is None:
        logger.info(f"  [Traps:{prefix}] No S/P/D-like state on one side of the gap in the fuzzy window.")
        return
    f = res["flag"]
    def count(hi, lo):
        m = (E <= hi + 1e-9) & (E >= lo - 1e-9)
        return int(np.sum(m & (f == 4))), int(np.sum(m & (f == 1)))
    fh, lh = count(E[res["homo"]], E[res["deloc_homo"]])
    fe, le = count(E[res["deloc_lumo"]], E[res["lumo"]])
    logger.info(f"  [Traps:{prefix}] Delocalized HOMO {E[res['deloc_homo']]:+.3f} eV ({fh} Facet + {lh} Localized states "
                f"above, nominal HOMO {E[res['homo']]:+.3f}); delocalized LUMO {E[res['deloc_lumo']]:+.3f} eV ({fe} + {le} "
                f"below, nominal LUMO {E[res['lumo']]:+.3f}); gap {E[res['lumo']] - E[res['homo']]:.3f} -> "
                f"{E[res['deloc_lumo']] - E[res['deloc_homo']]:.3f} eV")


def state_classes(args, pops, shells, coords_ang, syms):
    """Angular-coverage class, IPR and Omega of the states whose Mulliken AO populations are the columns of
    ``pops`` (qdex.angular.states_omega), and the AngularCoverage with the class boundaries."""
    from qdex.angular import states_omega
    om, cls, ipr, cov = states_omega(pops, shells, coords_ang, syms,
                                     core_elements=getattr(args, "centroid_core_elements", None))
    return cls, ipr, om, cov


def spinor_ao_pops(C_dense, S_dense, active_indices, U, beta_C=None, beta_indices=None):
    """Mulliken AO populations of the spinors U (columns, rows: alpha MOs then beta MOs of the active space):
    P_mu = sum_sigma Re[c_sigma,mu* (S c_sigma)_mu], c_alpha = C[:, active] U_alpha, c_beta likewise."""
    n_a = len(active_indices)
    Cb = C_dense if beta_C is None else beta_C
    ib = active_indices if beta_indices is None else beta_indices
    Ca, Cbb = C_dense[:, active_indices], Cb[:, ib]
    # S c = (S C_active) U: one real S product over the active MOs instead of one complex product per spinor
    SCa = S_dense @ Ca
    SCb = SCa if (beta_C is None and beta_indices is None) else S_dense @ Cbb
    Ua, Ub = U[:n_a, :], U[n_a:, :]
    ca, cb = Ca @ Ua, Cbb @ Ub
    return np.real(np.conj(ca) * (SCa @ Ua)) + np.real(np.conj(cb) * (SCb @ Ub))


def build_qp_energies(eps_dft, homo_index, scissor_ev=None, sigma_occ=None, sigma_virt=None):
    if sigma_occ is None and sigma_virt is None and scissor_ev is None:
        return None

    eps_qp = np.array(eps_dft, dtype=float, copy=True)
    occ_slice = slice(0, homo_index + 1)
    virt_slice = slice(homo_index + 1, len(eps_qp))

    if sigma_occ is not None:
        eps_qp[occ_slice] += np.asarray(sigma_occ, dtype=float)
    if sigma_virt is not None:
        eps_qp[virt_slice] += np.asarray(sigma_virt, dtype=float)
    elif scissor_ev is not None:
        eps_qp[virt_slice] += float(scissor_ev)

    return eps_qp


def build_qp_energies_vacuum(eps_abs, homo_index, qp_homo=None, qp_lumo=None, occ_shift=None, virt_shift=None):
    eps_qp = np.array(eps_abs, dtype=float, copy=True)
    if occ_shift is None:
        occ_shift = float(qp_homo) - float(eps_abs[homo_index])
    if virt_shift is None:
        virt_shift = float(qp_lumo) - float(eps_abs[homo_index + 1])
    eps_qp[:homo_index + 1] += occ_shift
    eps_qp[homo_index + 1:] += virt_shift
    return eps_qp


def build_soc_qp_energies_vacuum(soc_E_abs, soc_U, eps_abs_spin, eps_qp_abs_spin):
    delta_qp = np.asarray(eps_qp_abs_spin, dtype=float) - np.asarray(eps_abs_spin, dtype=float)
    return np.asarray(soc_E_abs, dtype=float) + (np.abs(soc_U) ** 2).T @ delta_qp


def auto_energy_window(*energy_arrays, sigma_ev=0.03):
    vals = []
    for energies in energy_arrays:
        if energies is None:
            continue
        arr = np.asarray(energies, dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size:
            vals.append(arr)
    if not vals:
        return np.array([-5.0, 5.0], dtype=float)

    all_e = np.concatenate(vals)
    e_min = float(np.min(all_e))
    e_max = float(np.max(all_e))
    if np.isclose(e_min, e_max):
        pad = max(10.0 * float(sigma_ev), 0.1)
    else:
        pad = max(4.0 * float(sigma_ev), 0.02 * (e_max - e_min), 0.05)
    return np.array([e_min - pad, e_max + pad], dtype=float)


def build_smeared_fuzzy(intensity, eps_plot, ewin, sigma_ev):
    window_mask = (eps_plot >= ewin[0] - 4 * sigma_ev) & (eps_plot <= ewin[1] + 4 * sigma_ev)
    E_w = eps_plot[window_mask]
    I_w = intensity[window_mask, :]

    dE = max(0.5 * sigma_ev, 0.01)
    edges = np.arange(ewin[0], ewin[1] + dE, dE)
    centres = 0.5 * (edges[:-1] + edges[1:])

    Z = np.zeros((centres.size, I_w.shape[1]), dtype=float)
    for En, Ik in zip(E_w, I_w):
        w = np.exp(-0.5 * ((centres - En) / sigma_ev) ** 2)
        Z += np.outer(w, Ik)
    return centres, Z


def state_weights(intensity):
    """Weight of each state along the path normalised to a mean of 1 (1 = spread evenly over the path).

    The raw ``|<phi|k>|^2`` of a state summed over the folded replicas grows with its localisation in
    real space (a cation d state carries ~10x the weight of a band-edge state), so on the raw map
    the semicore and surface states dominate the colour scale; per state the k-profile is what
    carries the band information."""
    intensity = np.asarray(intensity, float)
    mean = intensity.mean(axis=1, keepdims=True)
    return np.divide(intensity, mean, out=np.zeros_like(intensity), where=mean > 0)


def _path_segments(n_k, kpts_frac=None):
    """Index ranges of the continuous pieces of the k-path (split where the path jumps, e.g. K|U)."""
    if kpts_frac is None or len(kpts_frac) != n_k or n_k < 3:
        return [(0, n_k)]
    step = np.linalg.norm(np.diff(np.asarray(kpts_frac, float), axis=0), axis=1)
    jumps = np.flatnonzero(step > 3.0 * np.median(step[step > 1e-9])) + 1
    edges = [0, *jumps.tolist(), n_k]
    return [(a, b) for a, b in zip(edges[:-1], edges[1:]) if b > a]


def fuzzy_state_peaks(P, energies, ewin, kpts_frac=None, min_weight=1.5, rel_prominence=0.25, smooth=1.0):
    """Dominant k of every state: the maxima of its normalised k-profile (state_weights).

    A peak is kept when it is above `min_weight` (times the mean of the state) and stands out by
    `rel_prominence` of the state's maximum. Returns (k index, energy, weight, state index), one
    entry per peak. These are the sharp counterpart of the energy-smeared map: one marker per state
    at its own energy, at the wavevectors where its envelope is concentrated."""
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
    P = np.asarray(P, float)
    energies = np.asarray(energies, float)
    segs = _path_segments(P.shape[1], kpts_frac)
    ks, es, ws, ns = [], [], [], []
    for n in np.flatnonzero((energies >= ewin[0]) & (energies <= ewin[1])):
        p = np.concatenate([gaussian_filter1d(P[n, a:b], smooth, mode="nearest") if b - a > 2 else P[n, a:b]
                            for a, b in segs])
        if p.max() < min_weight:
            continue
        pk = np.concatenate([a + find_peaks(np.r_[0.0, p[a:b], 0.0], height=min_weight,
                                            prominence=rel_prominence * p.max())[0] - 1 for a, b in segs])
        ks.append(pk); ws.append(p[pk])
        es.append(np.full(pk.size, energies[n])); ns.append(np.full(pk.size, n))
    if not ks:
        return np.zeros(0, int), np.zeros(0), np.zeros(0), np.zeros(0, int)
    return np.concatenate(ks), np.concatenate(es), np.concatenate(ws), np.concatenate(ns)


def spinor_soc_energy(soc_E, soc_U, eps_spin):
    """<psi_n|V_SOC|psi_n> of every spinor of the window, in the units of soc_E.

    In the basis (phi_i alpha, phi_i beta) the spinor Hamiltonian is diag(eps) + V_SOC, so
    ``<V>_n = E_n - sum_i |U_in|^2 eps_i``. V_SOC has a zero trace (L is purely imaginary in a real
    basis), so a constant offset between the frames of soc_E and eps_spin is removed by setting the
    mean over the window to zero. For p-like states <V> = +Delta_so/3 for j = 3/2 and -2 Delta_so/3
    for j = 1/2: the sign tells the heavy/light-hole states from the split-off ones."""
    soc_U = np.asarray(soc_U)
    v = np.asarray(soc_E, float) - (np.abs(soc_U) ** 2 * np.asarray(eps_spin, float)[:, None]).sum(axis=0)
    return v - v.mean()


def deloc_edge_mos(cls, energies, mo_indices, n_homos, n_lumos):
    """MO indices of the first n_homos band states from the delocalized HOMO down and the first n_lumos from the
    delocalized LUMO up (S/P/D-like states of classify_edges), as {mo: label}, labels
    dHOMO, dHOMO-1, ..., dLUMO, dLUMO+1, ..."""
    E = np.asarray(energies, dtype=float)
    mo = np.asarray(mo_indices)
    band = ~np.asarray(cls["localized"], dtype=bool)
    out = {}
    for key, n, sign in (("homo", n_homos, -1), ("lumo", n_lumos, +1)):
        i0 = cls.get("deloc_" + key)
        if i0 is None or n <= 0:
            continue
        occ_side = mo <= mo[i0] if key == "homo" else mo >= mo[i0]
        cand = [int(i) for i in np.argsort(sign * E) if band[i] and occ_side[i] and sign * (E[i] - E[i0]) >= -1e-9]
        name = "dHOMO" if key == "homo" else "dLUMO"
        for k, i in enumerate(cand[:n]):
            out[int(mo[i])] = name if k == 0 else f"{name}{'-' if key == 'homo' else '+'}{k}"
    return out


def write_deloc_cubes(args, cls, energies, mo_indices, homo_index, C_dense, shells, syms, coords_ang):
    """Cubes of the delocalized band edges (output.cube_nhomos_deloc / cube_nlumos_deloc): spatial_MO_<nominal
    label>_<dHOMO|dLUMO...>.cube. MOs already written as nominal cubes (cube / mo_cubes, cube_nhomos and
    cube_nlumos around the gap) are not written again."""
    n_h, n_l = int(getattr(args, "cube_nhomos_deloc", 0) or 0), int(getattr(args, "cube_nlumos_deloc", 0) or 0)
    if n_h <= 0 and n_l <= 0:
        return []
    wanted = deloc_edge_mos(cls, energies, mo_indices, n_h, n_l)
    nominal = set()
    if getattr(args, "cube", False) or getattr(args, "mo_cubes", False):
        nh0, nl0 = getattr(args, "cube_nhomos", 2), getattr(args, "cube_nlumos", 2)
        nominal = {homo_index - i for i in range(nh0)} | {homo_index + 1 + i for i in range(nl0)}
    skipped = {m: lbl for m, lbl in wanted.items() if m in nominal}
    todo = {m: lbl for m, lbl in wanted.items() if m not in nominal}
    if skipped:
        logger.info("  [Cube] Delocalized edges already among the nominal cubes: "
                    + ", ".join(f"{lbl} = MO {m}" for m, lbl in sorted(skipped.items())))
    if not todo:
        return []
    from qdex.exciton_cube import generate_cubes
    from types import SimpleNamespace
    spacing = getattr(args, "mo_cube_spacing", 0.8) if (getattr(args, "mo_cubes", False) or not getattr(args, "cube", False)) \
        else getattr(args, "cube_spacing", 0.5)
    logger.info(f"\n--- Delocalized band-edge cubes ({len(todo)} spin-free MOs, {spacing} A grid): "
                + ", ".join(f"{lbl} = MO {m}" for m, lbl in sorted(todo.items())) + " ---")
    mo_list = sorted(todo)
    generate_cubes(solver=SimpleNamespace(C=C_dense, homo_index=homo_index), bse_states_dict={}, mo_list=mo_list,
                   spinor_list=[], soc_U=None, shells=shells, symbols=syms, coords=coords_ang, spacing_ang=spacing,
                   nthreads=args.nthreads, use_cpp=not getattr(args, "disable_cpp_cube", False), mo_suffix=todo)
    return mo_list


def smear_and_export_fuzzy(intensity, eps_plot, labels, ewin, sigma_ev, prefix="sf", export=True, kpts_frac=None,
                           soc_energy=None, trap=None):
    """Write fuzzy_data_<prefix>.npz for the dashboard: the energy-smeared map of the raw weights
    (intensity), of the per-state normalised weights (intensity_norm), the k-peaks of every state
    (peak_*) and, for spinors, their spin-orbit energy <V_SOC> (soc_energy: per state, NaN where
    no SOC was applied) as a weighted map (soc_energy_map) and on the peaks. ``trap`` (rows as
    intensity: flag from classify_edges, ipr, omega; plus thresholds and L of qdex.angular) adds state_energy,
    state_flag, state_ipr, state_omega, omega_thresholds (Facet, D, P, S lower bounds) and omega_L, from which the
    dashboard draws the angular-coverage panel and the delocalized band edges."""
    if not export:
        return
    t0 = time.time()
    centres, Z = build_smeared_fuzzy(intensity, eps_plot, ewin, sigma_ev)
    P = state_weights(intensity)
    _, Zn = build_smeared_fuzzy(P, eps_plot, ewin, sigma_ev)
    pk_k, pk_e, pk_w, pk_n = fuzzy_state_peaks(P, eps_plot, ewin, kpts_frac=kpts_frac)
    soc_extra = {}
    if soc_energy is not None:
        v = np.asarray(soc_energy, float)
        has = np.isfinite(v)
        _, Zv = build_smeared_fuzzy(P[has] * v[has, None], eps_plot[has], ewin, sigma_ev)
        _, Zh = build_smeared_fuzzy(P[has], eps_plot[has], ewin, sigma_ev)
        vmap = np.divide(Zv, Zh, out=np.full_like(Zv, np.nan), where=Zh > 1e-3 * max(Zh.max(), 1e-30))
        soc_extra = dict(soc_energy=v.astype(np.float32), soc_energy_map=vmap.astype(np.float32),
                         peak_soc_energy=v[pk_n].astype(np.float32))

    # ====================================================================
    # SCIENTIFIC FIX: CLEAN K-PATH LABELS & MERGE PATH BREAKS
    # ====================================================================
    valid_idx = []
    valid_labels = []
    for i, lbl in enumerate(labels):
        if lbl:
            # 1. Clean up LaTeX and Pymatgen formatting (convert to Unicode)
            clean_lbl = lbl.replace("\\Gamma", "Γ").replace("GAMMA", "Γ").replace("$", "")
            
            # 2. Handle Path Breaks (Adjacent duplicate or different labels)
            if valid_idx and (i - valid_idx[-1] <= 1):
                # Combine labels with a vertical bar if they are different
                if clean_lbl not in valid_labels[-1].split(" | "):
                    valid_labels[-1] = f"{valid_labels[-1]} | {clean_lbl}"
                valid_idx[-1] = i  # Snap to the exact current index
            else:
                valid_idx.append(i)
                valid_labels.append(clean_lbl)

    out_name = f"fuzzy_data_{prefix}.npz"
    extra = {} if kpts_frac is None else {"kpath_frac": np.asarray(kpts_frac, dtype=np.float64)}
    if trap is not None and len(eps_plot):
        extra.update(state_energy=np.asarray(eps_plot, dtype=np.float32),
                     state_flag=np.asarray(trap["flag"], dtype=np.int8),
                     state_ipr=np.asarray(trap["ipr"], dtype=np.float32),
                     state_omega=np.asarray(trap["omega"], dtype=np.float32))
        if trap.get("thresholds"):
            t = trap["thresholds"]
            extra.update(omega_thresholds=np.array([t["Facet"], t["D"], t["P"], t["S"]], dtype=np.float32),
                         omega_L=np.int32(trap.get("L", 0)))
    np.savez_compressed(
        out_name,
        **extra,
        **soc_extra,
        centres=centres.astype(np.float32),
        intensity=Z.astype(np.float32),
        intensity_norm=Zn.astype(np.float32),
        peak_k=pk_k.astype(np.int32), peak_energy=pk_e.astype(np.float32),
        peak_weight=pk_w.astype(np.float32), peak_state=pk_n.astype(np.int32),
        tick_positions=np.array(valid_idx, dtype=np.float32),
        tick_labels=np.array(valid_labels, dtype=object),
        ewin=np.array(ewin, dtype=np.float32),
        extent=np.array([0.0, float(Z.shape[1] - 1), float(ewin[0]), float(ewin[1])])
    )
    logger.debug(f"  [Fuzzy] Exported {out_name} in {time.time()-t0:.2f} s")


def smear_and_export_spin_fuzzy(intensity_alpha, eps_alpha, intensity_beta, eps_beta, labels, ewin, sigma_ev, prefix="uks", kpts_frac=None):
    t0 = time.time()
    centres, Z_a = build_smeared_fuzzy(intensity_alpha, eps_alpha, ewin, sigma_ev)
    _, Z_b = build_smeared_fuzzy(intensity_beta, eps_beta, ewin, sigma_ev)
    Z_total = Z_a + Z_b
    P_ab = state_weights(np.vstack([intensity_alpha, intensity_beta]))
    E_ab = np.concatenate([eps_alpha, eps_beta])
    _, Zn = build_smeared_fuzzy(P_ab, E_ab, ewin, sigma_ev)
    pk_k, pk_e, pk_w, pk_n = fuzzy_state_peaks(P_ab, E_ab, ewin, kpts_frac=kpts_frac)
    spinpol = np.divide(Z_a - Z_b, Z_total, out=np.zeros_like(Z_total), where=Z_total > 1e-14)

    valid_idx = []
    valid_labels = []
    for i, lbl in enumerate(labels):
        if lbl:
            clean_lbl = lbl.replace("\\Gamma", "Γ").replace("GAMMA", "Γ").replace("$", "")
            if valid_idx and (i - valid_idx[-1] <= 1):
                if clean_lbl not in valid_labels[-1].split(" | "):
                    valid_labels[-1] = f"{valid_labels[-1]} | {clean_lbl}"
                valid_idx[-1] = i
            else:
                valid_idx.append(i)
                valid_labels.append(clean_lbl)

    common = dict(
        centres=centres.astype(np.float32),
        tick_positions=np.array(valid_idx, dtype=np.float32),
        tick_labels=np.array(valid_labels, dtype=object),
        ewin=np.array(ewin, dtype=np.float32),
        extent=np.array([0.0, float(Z_total.shape[1] - 1), float(ewin[0]), float(ewin[1])])
    )
    if kpts_frac is not None:
        common["kpath_frac"] = np.asarray(kpts_frac, dtype=np.float64)
    np.savez_compressed(f"fuzzy_data_{prefix}.npz", intensity=Z_total.astype(np.float32), spinpol=spinpol.astype(np.float32),
                        intensity_norm=Zn.astype(np.float32), peak_k=pk_k.astype(np.int32),
                        peak_energy=pk_e.astype(np.float32), peak_weight=pk_w.astype(np.float32),
                        peak_state=pk_n.astype(np.int32), **common)
    np.savez_compressed(f"fuzzy_data_{prefix}_alpha.npz", intensity=Z_a.astype(np.float32), **common)
    np.savez_compressed(f"fuzzy_data_{prefix}_beta.npz", intensity=Z_b.astype(np.float32), **common)
    logger.debug(f"  [Fuzzy-UKS] Exported fuzzy_data_{prefix}.npz with spin polarization overlay in {time.time()-t0:.2f} s")

def _matmul_real_matrix(A_real, B, device="numpy"):
    from qdex.device_utils import is_gpu
    if not np.iscomplexobj(B):
        if is_gpu(device):
            import torch
            dev = torch.device(device)
            A_t = torch.from_numpy(A_real).to(device=dev, dtype=torch.float64)
            B_t = torch.from_numpy(np.ascontiguousarray(B)).to(device=dev, dtype=torch.float64)
            return torch.matmul(A_t, B_t).cpu().numpy()
        return A_real @ B
    if is_gpu(device):
        import torch
        dev = torch.device(device)
        A_t = torch.from_numpy(A_real).to(device=dev, dtype=torch.float64)
        B_r = torch.from_numpy(np.ascontiguousarray(np.real(B))).to(device=dev, dtype=torch.float64)
        B_i = torch.from_numpy(np.ascontiguousarray(np.imag(B))).to(device=dev, dtype=torch.float64)
        res_r = torch.matmul(A_t, B_r).cpu().numpy()
        res_i = torch.matmul(A_t, B_i).cpu().numpy()
        return res_r + 1j * res_i
    return (A_real @ np.ascontiguousarray(np.real(B))) + 1j * (A_real @ np.ascontiguousarray(np.imag(B)))

def run_fuzzy_bands_and_pdos(args, C_dense, S_dense, eps_shifted, occ, homo_index, e_homo, e_lumo, e_fermi_raw, syms, coords_ang, shells, pops_sf, soc_active_indices=None, soc_E_act=None, soc_U_act=None, spinor_homo_idx=None, qp_energies=None, eps_abs=None, qp_energies_abs=None, soc_E_abs_act=None, C_beta_dense=None, eps_beta_shifted=None, eps_beta_abs=None, homo_index_beta=None, qp_energies_beta=None, qp_energies_beta_abs=None, soc_active_indices_beta=None):
    import time
    import numpy as np
    from qdex.device_utils import is_gpu
    from qdex.pdos_coop import compute_pdos_and_coop, export_pdos_coop_data
    
    logger.info("\n===================================================")
    logger.info(" [ FUZZY BANDS & PDOS ]")
    logger.info("===================================================")
    
    fold_to_bz = bool(getattr(args, 'fold_to_bz', True))
    g_shell = int(getattr(args, 'g_shell', 2))
    kpts_cart, labels, reciprocal_matrix, kpath_info = generate_automated_kpath(
        args.cif, np.array(coords_ang), line_density=50, return_reciprocal=True, syms=syms, return_info=True)
    kpts_frac = kpath_info["kpts_frac"]
    if not fold_to_bz:
        reciprocal_matrix = None
    
    # --- 1. SPIN-FREE CALCULATION ---
    n_occ = homo_index + 1
    n_virt = len(eps_shifted) - n_occ
    logger.info(f"\n  [Fuzzy] --- Spin-Free MO Statistics ---")
    logger.info(f"  [Fuzzy] Total MOs: {len(eps_shifted)} ({n_occ} Occupied, {n_virt} Virtual)")
    logger.info(f"  [Fuzzy] MO HOMO (Idx {homo_index}): {e_homo:8.4f} eV")
    logger.info(f"  [Fuzzy] MO LUMO (Idx {homo_index + 1}): {e_lumo:8.4f} eV")
    logger.info(f"  [Fuzzy] Fermi Level (raw shifted to 0.0): {e_fermi_raw:8.4f} eV")
    logger.info(f"  [Fuzzy] -------------------------------")

    sigma_use = getattr(args, 'fuzzy_sigma', 0.03)
    pdos_sigma_use = getattr(args, 'pdos_sigma', 0.10)
    # Files for the HTML dashboards (which read them back) or for csv output
    export_files = bool(getattr(args, 'html', True) or getattr(args, 'write_csv', False))
    store = getattr(args, 'qdex_store', None)
    dashboard_energy_mode = str(getattr(args, 'dashboard_energy_mode', 'dft')).lower()
    if dashboard_energy_mode not in ("dft", "qp", "both"):
        raise ValueError("dashboard_energy_mode must be one of: dft, qp, both")
    qp_energy_reference = str(getattr(args, "qp_energy_reference", "vacuum")).lower()
    if qp_energy_reference not in ("vacuum", "fermi"):
        raise ValueError("qp_energy_reference must be one of: vacuum, fermi")
    qp_plot_energies = qp_energies_abs if qp_energy_reference == "vacuum" else qp_energies
    qp_plot_energies_beta = qp_energies_beta_abs if qp_energy_reference == "vacuum" else qp_energies_beta
    dft_ewin = list(getattr(args, "ewin", [-5.0, 5.0]))
    qp_ewin = dft_ewin if qp_energy_reference == "fermi" else (auto_energy_window(qp_plot_energies, sigma_ev=sigma_use) if qp_plot_energies is not None else None)
    qp_mask_energies = qp_plot_energies if qp_energy_reference == "fermi" else None
    fuzzy_indices = fuzzy_energy_indices(
        eps_shifted, dft_ewin, sigma_use, homo_index=homo_index, qp_energies=qp_mask_energies
    )
    eps_fuzzy = eps_shifted[fuzzy_indices]
    qp_fuzzy = qp_plot_energies[fuzzy_indices] if qp_plot_energies is not None else None

    logger.info(f"  [Fuzzy] Projecting {len(fuzzy_indices)} / {len(eps_shifted)} MOs in ewin [{dft_ewin[0]:.3f}, {dft_ewin[1]:.3f}] eV relative to mid-gap.")
    if qp_ewin is not None:
        logger.info(f"  [Fuzzy] QP plot window: [{qp_ewin[0]:.3f}, {qp_ewin[1]:.3f}] eV")

    is_uks = C_beta_dense is not None and eps_beta_shifted is not None and homo_index_beta is not None
    if is_uks:
        qp_mask_energies_beta = qp_plot_energies_beta if qp_energy_reference == "fermi" else None
        fuzzy_indices_b = fuzzy_energy_indices(
            eps_beta_shifted, dft_ewin, sigma_use, homo_index=homo_index_beta, qp_energies=qp_mask_energies_beta
        )
        eps_fuzzy_b = eps_beta_shifted[fuzzy_indices_b]

    # --- SOC states of the map: spin-free MOs below/above the SOC window and the active spinors ---
    do_soc = bool(args.soc_flag and soc_active_indices is not None)
    soc_uks = do_soc and is_uks and soc_active_indices_beta is not None
    if do_soc:
        n_act_mo = len(soc_active_indices)
        n_act_occ = np.sum(soc_active_indices <= homo_index)
        n_act_virt = n_act_mo - n_act_occ

        logger.info(f"\n  [Fuzzy-SOC] Applying Precomputed Unified SOC Projection...")
        logger.info(f"  [Fuzzy-SOC] --- Dual-Window SOC Statistics ---")
        logger.info(f"  [Fuzzy-SOC] Active Space Spatial MOs: {n_act_mo} ({n_act_occ} Occ, {n_act_virt} Virt)")
        logger.info(f"  [Fuzzy-SOC] Full spinor basis available: {len(eps_shifted) * 2} states")

        alpha_plot_indices = fuzzy_energy_indices(eps_shifted, dft_ewin, sigma_use, homo_index=homo_index)
        core_idx = alpha_plot_indices[alpha_plot_indices < soc_active_indices[0]]
        virt_idx = alpha_plot_indices[alpha_plot_indices > soc_active_indices[-1]]
        if soc_uks:
            beta_plot_indices = fuzzy_energy_indices(eps_beta_shifted, dft_ewin, sigma_use, homo_index=homo_index_beta)
            core_idx_b = beta_plot_indices[beta_plot_indices < soc_active_indices_beta[0]]
            virt_idx_b = beta_plot_indices[beta_plot_indices > soc_active_indices_beta[-1]]

        if soc_uks:
            E_core = np.concatenate([eps_shifted[core_idx], eps_beta_shifted[core_idx_b]])
            E_virt = np.concatenate([eps_shifted[virt_idx], eps_beta_shifted[virt_idx_b]])
            if eps_abs is not None and eps_beta_abs is not None and qp_energies_abs is not None and qp_energies_beta_abs is not None and soc_E_abs_act is not None:
                eps_abs_spin_act = np.concatenate([eps_abs[soc_active_indices], eps_beta_abs[soc_active_indices_beta]])
                eps_qp_abs_spin_act = np.concatenate([qp_energies_abs[soc_active_indices], qp_energies_beta_abs[soc_active_indices_beta]])
                soc_E_qp_act = build_soc_qp_energies_vacuum(soc_E_abs_act, soc_U_act, eps_abs_spin_act, eps_qp_abs_spin_act)
                E_core_qp = np.concatenate([qp_energies_abs[core_idx], qp_energies_beta_abs[core_idx_b]])
                E_virt_qp = np.concatenate([qp_energies_abs[virt_idx], qp_energies_beta_abs[virt_idx_b]])
            else:
                soc_E_qp_act = E_core_qp = E_virt_qp = None
        else:
            E_core = np.concatenate([eps_shifted[core_idx], eps_shifted[core_idx]])
            E_virt = np.concatenate([eps_shifted[virt_idx], eps_shifted[virt_idx]])
            if eps_abs is not None and qp_energies_abs is not None and soc_E_abs_act is not None:
                eps_abs_spin_act = np.concatenate([eps_abs[soc_active_indices], eps_abs[soc_active_indices]])
                eps_qp_abs_spin_act = np.concatenate([qp_energies_abs[soc_active_indices], qp_energies_abs[soc_active_indices]])
                soc_E_qp_act = build_soc_qp_energies_vacuum(soc_E_abs_act, soc_U_act, eps_abs_spin_act, eps_qp_abs_spin_act)
                E_core_qp = np.concatenate([qp_energies_abs[core_idx], qp_energies_abs[core_idx]])
                E_virt_qp = np.concatenate([qp_energies_abs[virt_idx], qp_energies_abs[virt_idx]])
            else:
                soc_E_qp_act = E_core_qp = E_virt_qp = None

        # soc_E_act is centred on the spinor mid-gap; put the spin-free core and virtual
        # states of this map (and the bulk anchor) on that same axis
        soc_axis_offset = 0.0
        if soc_E_abs_act is not None:
            soc_axis_offset = float(np.mean(np.asarray(soc_E_act) - (np.asarray(soc_E_abs_act) - e_fermi_raw)))
            E_core = E_core + soc_axis_offset
            E_virt = E_virt + soc_axis_offset
        logger.info(f"  [Fuzzy-SOC] Spinor mid-gap is {-soc_axis_offset:+.4f} eV from the spin-free mid-gap.")

        eps_soc_unsorted = np.concatenate([E_core, soc_E_act, E_virt])
        plot_keep = fuzzy_energy_mask(eps_soc_unsorted, dft_ewin, sigma_use)
        below_zero = np.where(eps_soc_unsorted <= 0.0)[0]
        above_zero = np.where(eps_soc_unsorted > 0.0)[0]
        if below_zero.size:
            plot_keep[below_zero[np.argmax(eps_soc_unsorted[below_zero])]] = True
        if above_zero.size:
            plot_keep[above_zero[np.argmin(eps_soc_unsorted[above_zero])]] = True
        n_core_rows, n_act = len(E_core), len(soc_E_act)

    # --- Folded plane-wave weights of every plotted state: one pass over the G replicas ---
    if fold_to_bz:
        G_vecs = make_reciprocal_replicas(reciprocal_matrix, g_shell)
        logger.info(f"  [Fuzzy] BZ folding enabled: g_shell={g_shell} ({len(G_vecs)} reciprocal replicas).")
    else:
        G_vecs = np.zeros((1, 3))
    alpha_set = fuzzy_indices
    if do_soc:
        alpha_set = np.unique(np.concatenate([fuzzy_indices, core_idx, soc_active_indices, virt_idx])).astype(int)
    pos_a = {int(i): p for p, i in enumerate(alpha_set)}
    blocks = [C_dense[:, alpha_set]]
    if is_uks:
        beta_set = fuzzy_indices_b
        if soc_uks:
            beta_set = np.unique(np.concatenate([fuzzy_indices_b, core_idx_b, soc_active_indices_beta, virt_idx_b])).astype(int)
        pos_b = {int(i): p for p, i in enumerate(beta_set)}
        blocks.append(C_beta_dense[:, beta_set])
    spinor_parts = None
    if do_soc:
        # active spinors kept in the plot: psi_n = sum_i U[i, n] phi_i alpha + U[n_a + i, n] phi'_i beta
        U_kept = soc_U_act[:, plot_keep[n_core_rows:n_core_rows + n_act]]
        n_a = len(soc_active_indices)
        rows_a = [pos_a[int(i)] for i in soc_active_indices]
        if soc_uks:
            rows_b = [pos_b[int(i)] for i in soc_active_indices_beta]
            spinor_parts = [(0, rows_a, U_kept[:n_a]), (1, rows_b, U_kept[n_a:])]
        else:
            spinor_parts = [(0, rows_a, U_kept[:n_a]), (0, rows_a, U_kept[n_a:])]

    logger.info("  [Fuzzy] Computing Analytic AO-FT via C++ ...")
    t_ft = time.time()
    W_blocks, W_spinor = folded_plane_wave_weights(shells, kpts_cart, G_vecs, args.nthreads, blocks, spinor_parts)
    logger.debug(f"  [Fuzzy] Folded plane-wave weights of {sum(b.shape[1] for b in blocks)} MOs"
                 f"{f' and {W_spinor.shape[0]} spinors' if W_spinor is not None else ''} in {time.time() - t_ft:.2f} s")

    def _rows(W, pos, indices):
        return W[[pos[int(i)] for i in indices]] if len(indices) else np.zeros((0, W.shape[1]))

    intensity_sf = _rows(W_blocks[0], pos_a, fuzzy_indices)

    if fold_to_bz and g_shell == 0:
        intensity_ref = compute_fuzzy_intensity(C_dense, shells, kpts_cart, args.nthreads, fold_to_bz=False, mo_indices=fuzzy_indices)
        _, Z_ref = build_smeared_fuzzy(intensity_ref, eps_fuzzy, dft_ewin, sigma_use)
        _, Z_fold = build_smeared_fuzzy(intensity_sf, eps_fuzzy, dft_ewin, sigma_use)
        logger.info(
            "  [Fuzzy] g_shell=0 folding diagnostic: "
            f"max |dW|={np.max(np.abs(intensity_sf - intensity_ref)):.3e}, "
            f"max |dA|={np.max(np.abs(Z_fold - Z_ref)):.3e}"
        )
    
    # Semicore level of the interior atoms: places the bulk bands on this energy axis
    semicore = None
    try:
        from qdex.bulk_bands import qd_semicore_level
        if np.isfinite(kpath_info["qd_bond_ang"]):
            semicore = qd_semicore_level(
                args.material, C_dense, S_dense, shells, syms, coords_ang, eps_shifted, homo_index,
                coordination=kpath_info["coordination"], bond_ang=kpath_info["qd_bond_ang"],
                bs_path=getattr(args, "bulk_bs", None), cif=args.cif,
                anchor=getattr(args, "bulk_anchor", "auto"), bond_pairs=kpath_info.get("bond_pairs"))
    except Exception as exc:
        logger.warning(f"  [Bulk Bands] Semicore level of the dot not computed: {exc}")
    semicore_rel = {"sf": semicore["level_ev"]} if semicore is not None else {}
    semicore_label = semicore.get("label") if semicore is not None else None
    if semicore is None:
        logger.info("  [Bulk Bands] No semicore level measurable for this dot: bulk bands aligned at mid-gap.")
    if store is not None and semicore is not None:
        store.put("electronic", "sf/bulk_anchor/semicore_level_ev", semicore["level_ev"])
        store.attr("electronic", "sf/bulk_anchor", element=semicore["element"], label=semicore_label,
                   n_atoms=semicore["n_atoms"],
                   spread_ev=semicore["spread_ev"], all_atoms_level_ev=semicore["all_atoms_level_ev"],
                   note="semicore level of the interior bulk-like atoms on the fuzzy energy axis")

    # surface / core labels of the population analysis (one rule for table, log, dashboard and cubes)
    occ_f = np.asarray(fuzzy_indices) <= homo_index
    lab_sf, ipr_sf, om_sf, cov = state_classes(args, np.asarray(pops_sf)[:, fuzzy_indices], shells, coords_ang, syms)
    cls_sf = classify_edges(eps_fuzzy, occ_f, lab_sf)
    edge_summary("SF", eps_fuzzy, cls_sf)
    trap_sf = dict(flag=cls_sf["flag"], ipr=ipr_sf, omega=om_sf, thresholds=cov.thresholds, L=cov.L)
    write_deloc_cubes(args, cls_sf, eps_fuzzy, fuzzy_indices, homo_index, C_dense, shells, syms, coords_ang)

    smear_and_export_fuzzy(intensity_sf, eps_fuzzy, labels, dft_ewin, sigma_use, prefix="sf", export=export_files, kpts_frac=kpts_frac,
                           trap=trap_sf)
    if store is not None:
        from qdex.store import put_fuzzy
        put_fuzzy(store, "sf", kpts_cart, labels, eps_fuzzy, intensity_sf, sigma_use, dft_ewin, indices=fuzzy_indices,
                  kpts_frac=kpts_frac, cif=args.cif)
        store.put("electronic", "sf/fuzzy/trap_flag", np.asarray(trap_sf["flag"], np.int8))
        store.put("electronic", "sf/fuzzy/ipr", np.asarray(trap_sf["ipr"], float))
        store.put("electronic", "sf/fuzzy/omega", np.asarray(trap_sf["omega"], float))
        store.attr("electronic", "sf/fuzzy", omega_L=int(cov.L), omega_thresholds=json.dumps(cov.thresholds))

    pdos_analysis_sf = None
    if getattr(args, 'pdos_atoms', None) and getattr(args, 'coop_pairs', None):
        logger.info("  [PDOS/COOP] Computing Spin-Free population analysis...")
        pdos_analysis_sf = compute_pdos_and_coop(
            C_dense, S_dense, eps_shifted, shells, args.pdos_atoms, args.coop_pairs, dft_ewin,
            sigma=pdos_sigma_use, is_soc=False, prefix="sf", pops=pops_sf,
            population_bars=getattr(args, "population_bars", None),
            device=getattr(args, "device", "numpy"), export=export_files
        )
        if store is not None:
            pairs = [p for p in args.coop_pairs if p in pdos_analysis_sf["coop_results"]]
            if pairs:
                store.put("electronic", "sf/mo/coop_pairs", pairs)
                store.put("electronic", "sf/mo/coop", np.stack([pdos_analysis_sf["coop_results"][p] for p in pairs], axis=1))

    if dashboard_energy_mode in ("qp", "both"):
        if qp_plot_energies is None:
            msg = "QP dashboard requested, but QP-corrected orbital energies were not found."
            if dashboard_energy_mode == "qp":
                raise ValueError(msg)
            logger.warning(f"  [Warning] {msg} Skipping QP dashboard.")
        else:
            smear_and_export_fuzzy(intensity_sf, qp_fuzzy, labels, qp_ewin, sigma_use, prefix="sf_qp", kpts_frac=kpts_frac,
                                   trap=trap_sf)
            if pdos_analysis_sf is not None:
                export_pdos_coop_data(
                    pdos_analysis_sf, qp_plot_energies, args.pdos_atoms, args.coop_pairs, qp_ewin,
                    sigma=pdos_sigma_use, is_soc=False, prefix="sf_qp",
                    population_bars=getattr(args, "population_bars", None)
                )

    if is_uks:
        logger.info(f"\n  [Fuzzy-UKS] Computing alpha/beta fuzzy channels with total intensity and spin polarization...")
        dft_uks_ewin = dft_ewin
        qp_uks_ewin = auto_energy_window(qp_plot_energies, qp_plot_energies_beta, sigma_ev=sigma_use) if qp_plot_energies is not None and qp_plot_energies_beta is not None else None
        intensity_a = intensity_sf
        intensity_b = _rows(W_blocks[1], pos_b, fuzzy_indices_b)
        smear_and_export_spin_fuzzy(intensity_a, eps_fuzzy, intensity_b, eps_fuzzy_b, labels, dft_uks_ewin, sigma_use, prefix="uks", kpts_frac=kpts_frac)

        if dashboard_energy_mode in ("qp", "both") and qp_plot_energies is not None and qp_plot_energies_beta is not None:
            smear_and_export_spin_fuzzy(
                intensity_a, qp_plot_energies[fuzzy_indices],
                intensity_b, qp_plot_energies_beta[fuzzy_indices_b],
                labels, qp_uks_ewin, sigma_use, prefix="uks_qp", kpts_frac=kpts_frac
            )

    # --- 2. SOC CALCULATION ---
    if do_soc:
        if semicore is not None:
            semicore_rel["soc"] = semicore["level_ev"] + soc_axis_offset
            if store is not None:
                store.put("electronic", "soc/bulk_anchor/semicore_level_ev", semicore_rel["soc"])

        # folded weights of the plotted states, in the order of eps_soc_unsorted[plot_keep]
        if soc_uks:
            I_core_rows = np.vstack([_rows(W_blocks[0], pos_a, core_idx), _rows(W_blocks[1], pos_b, core_idx_b)])
            I_virt_rows = np.vstack([_rows(W_blocks[0], pos_a, virt_idx), _rows(W_blocks[1], pos_b, virt_idx_b)])
        else:
            I_core_one = _rows(W_blocks[0], pos_a, core_idx)
            I_virt_one = _rows(W_blocks[0], pos_a, virt_idx)
            I_core_rows = np.vstack([I_core_one, I_core_one])
            I_virt_rows = np.vstack([I_virt_one, I_virt_one])
        I_plot_unsorted = np.vstack([
            I_core_rows[plot_keep[:n_core_rows]],
            W_spinor,
            I_virt_rows[plot_keep[n_core_rows + n_act:]],
        ])
        eps_soc_plot_unsorted = eps_soc_unsorted[plot_keep]
        # <V_SOC> of the active spinors (j = 3/2 vs 1/2 character); NaN for the spin-free core/virtual rows
        if soc_uks:
            eps_spin_act = np.concatenate([eps_shifted[soc_active_indices], eps_beta_shifted[soc_active_indices_beta]])
        else:
            eps_spin_act = np.concatenate([eps_shifted[soc_active_indices], eps_shifted[soc_active_indices]])
        v_act = spinor_soc_energy(soc_E_act, soc_U_act, eps_spin_act)
        v_soc_plot_unsorted = np.concatenate([
            np.full(int(plot_keep[:n_core_rows].sum()), np.nan),
            v_act[plot_keep[n_core_rows:n_core_rows + n_act]],
            np.full(int(plot_keep[n_core_rows + n_act:].sum()), np.nan),
        ])
        # surface / core labels of the plotted rows: spin-free MOs (core, virtual) and the SOC spinors
        P_sf_all = np.asarray(pops_sf)
        U_kept_act = soc_U_act[:, plot_keep[n_core_rows:n_core_rows + n_act]]
        if soc_uks:
            P_beta_all = np.real(C_beta_dense * (S_dense @ C_beta_dense))
            P_core = np.hstack([P_sf_all[:, core_idx], P_beta_all[:, core_idx_b]])
            P_virt = np.hstack([P_sf_all[:, virt_idx], P_beta_all[:, virt_idx_b]])
            P_act = spinor_ao_pops(C_dense, S_dense, soc_active_indices, U_kept_act, C_beta_dense, soc_active_indices_beta)
        else:
            P_core = np.hstack([P_sf_all[:, core_idx]] * 2)
            P_virt = np.hstack([P_sf_all[:, virt_idx]] * 2)
            P_act = spinor_ao_pops(C_dense, S_dense, soc_active_indices, U_kept_act)
        P_plot_unsorted = np.hstack([P_core[:, plot_keep[:n_core_rows]], P_act,
                                     P_virt[:, plot_keep[n_core_rows + n_act:]]])
        lab_u, ipr_u, om_u, cov = state_classes(args, P_plot_unsorted, shells, coords_ang, syms)
        cls_u = classify_edges(eps_soc_unsorted[plot_keep], eps_soc_unsorted[plot_keep] <= 0.0, lab_u)
        edge_summary("SOC", eps_soc_unsorted[plot_keep], cls_u)
        trap_soc_u = dict(flag=cls_u["flag"], ipr=ipr_u, omega=om_u)
        soc_meta = dict(thresholds=cov.thresholds, L=cov.L)
        sort_idx = np.argsort(eps_soc_plot_unsorted)
        eps_soc = eps_soc_plot_unsorted[sort_idx]
        intensity_soc = I_plot_unsorted[sort_idx, :]
        soc_energy_plot = v_soc_plot_unsorted[sort_idx]
        soc_ewin = dft_ewin
        occupied_plot = np.where(eps_soc <= 0.0)[0]
        global_spinor_homo_idx = int(occupied_plot[-1]) if occupied_plot.size else 0
 
        logger.info(f"  [Fuzzy-SOC] Plotting {len(eps_soc)} spinors in ewin [{soc_ewin[0]:.3f}, {soc_ewin[1]:.3f}] eV relative to mid-gap.")
        logger.info(f"  [Fuzzy-SOC] Spinor HOMO (Idx {global_spinor_homo_idx}): {eps_soc[global_spinor_homo_idx]:8.4f} eV")
        logger.info(f"  [Fuzzy-SOC] Spinor LUMO (Idx {global_spinor_homo_idx + 1}): {eps_soc[global_spinor_homo_idx + 1]:8.4f} eV")
        logger.info(f"  [Fuzzy-SOC] ----------------------------------") 
        
        smear_and_export_fuzzy(intensity_soc, eps_soc, labels, soc_ewin, sigma_use, prefix="soc", export=export_files,
                               kpts_frac=kpts_frac, soc_energy=soc_energy_plot,
                               trap={**{k: v[sort_idx] for k, v in trap_soc_u.items()}, **soc_meta})
        if store is not None:
            from qdex.store import put_fuzzy
            put_fuzzy(store, "soc", kpts_cart, labels, eps_soc, intensity_soc, sigma_use, soc_ewin, kpts_frac=kpts_frac,
                      cif=args.cif, soc_energy=soc_energy_plot)
            store.put("electronic", "soc/fuzzy/trap_flag", np.asarray(trap_soc_u["flag"][sort_idx], np.int8))
            store.put("electronic", "soc/fuzzy/ipr", np.asarray(trap_soc_u["ipr"][sort_idx], float))
            store.put("electronic", "soc/fuzzy/omega", np.asarray(trap_soc_u["omega"][sort_idx], float))
            store.attr("electronic", "soc/fuzzy", omega_L=int(soc_meta["L"]), omega_thresholds=json.dumps(soc_meta["thresholds"]))

        eps_soc_qp = None
        sort_idx_qp = None
        soc_qp_ewin = None
        if dashboard_energy_mode in ("qp", "both") and soc_E_qp_act is not None:
            eps_soc_qp_abs_unsorted = np.concatenate([E_core_qp, soc_E_qp_act, E_virt_qp])[plot_keep]
            if qp_energy_reference == "fermi":
                eps_soc_qp_abs_sorted = np.sort(eps_soc_qp_abs_unsorted)
                qp_midgap = 0.5 * (eps_soc_qp_abs_sorted[global_spinor_homo_idx] + eps_soc_qp_abs_sorted[global_spinor_homo_idx + 1])
                eps_soc_qp_unsorted = eps_soc_qp_abs_unsorted - qp_midgap
                soc_qp_ewin = dft_ewin
            else:
                eps_soc_qp_unsorted = eps_soc_qp_abs_unsorted
                soc_qp_ewin = auto_energy_window(eps_soc_qp_unsorted, sigma_ev=sigma_use)
            sort_idx_qp = np.argsort(eps_soc_qp_unsorted)
            eps_soc_qp = eps_soc_qp_unsorted[sort_idx_qp]
            intensity_soc_qp = I_plot_unsorted[sort_idx_qp, :]
            smear_and_export_fuzzy(intensity_soc_qp, eps_soc_qp, labels, soc_qp_ewin, sigma_use, prefix="soc_qp", kpts_frac=kpts_frac,
                                   soc_energy=v_soc_plot_unsorted[sort_idx_qp],
                                   trap={**{k: v[sort_idx_qp] for k, v in trap_soc_u.items()}, **soc_meta})
        
        if getattr(args, 'pdos_atoms', None) and getattr(args, 'coop_pairs', None):
            logger.info("  [PDOS/COOP] Computing SOC Spinor population analysis...")
            t_pop = time.time()
            n_ao = S_dense.shape[0]
            device = getattr(args, "device", "numpy")

            if is_uks and soc_active_indices_beta is not None:
                n_ca = len(core_idx)
                n_cb = len(core_idx_b)
                n_act = soc_U_act.shape[1]
                n_va = len(virt_idx)
                n_vb = len(virt_idx_b)
                n_alpha_act = len(soc_active_indices)
                n_beta_act = len(soc_active_indices_beta)

                b_core_a = np.zeros((2 * n_ao, 0), dtype=complex)
                b_core_b = np.zeros((2 * n_ao, 0), dtype=complex)
                b_act = np.zeros((2 * n_ao, 0), dtype=complex)
                b_virt_a = np.zeros((2 * n_ao, 0), dtype=complex)
                b_virt_b = np.zeros((2 * n_ao, 0), dtype=complex)

                m_ca = plot_keep[:n_ca]
                if np.any(m_ca):
                    b_core_a = np.zeros((2 * n_ao, np.sum(m_ca)), dtype=complex)
                    b_core_a[:n_ao, :] = C_dense[:, core_idx[m_ca]]

                m_cb = plot_keep[n_ca : n_ca + n_cb]
                if np.any(m_cb):
                    b_core_b = np.zeros((2 * n_ao, np.sum(m_cb)), dtype=complex)
                    b_core_b[n_ao:, :] = C_beta_dense[:, core_idx_b[m_cb]]

                m_act = plot_keep[n_ca + n_cb : n_ca + n_cb + n_act]
                if np.any(m_act):
                    U_act_kept = soc_U_act[:, m_act]
                    C_act_a = C_dense[:, soc_active_indices]
                    C_act_b = C_beta_dense[:, soc_active_indices_beta]
                    C_act_top = _matmul_real_matrix(C_act_a, U_act_kept[:n_alpha_act, :], device=device)
                    C_act_bot = _matmul_real_matrix(C_act_b, U_act_kept[n_alpha_act:n_alpha_act + n_beta_act, :], device=device)
                    b_act = np.vstack([C_act_top, C_act_bot])

                m_va = plot_keep[n_ca + n_cb + n_act : n_ca + n_cb + n_act + n_va]
                if np.any(m_va):
                    b_virt_a = np.zeros((2 * n_ao, np.sum(m_va)), dtype=complex)
                    b_virt_a[:n_ao, :] = C_dense[:, virt_idx[m_va]]

                m_vb = plot_keep[n_ca + n_cb + n_act + n_va :]
                if np.any(m_vb):
                    b_virt_b = np.zeros((2 * n_ao, np.sum(m_vb)), dtype=complex)
                    b_virt_b[n_ao:, :] = C_beta_dense[:, virt_idx_b[m_vb]]

                blocks = [b_core_a, b_core_b, b_act, b_virt_a, b_virt_b]
                C_spinor_ao_plot = np.hstack([b for b in blocks if b.shape[1] > 0])
            else:
                n_c = len(core_idx)
                n_act = len(soc_E_act)
                n_v = len(virt_idx)
                n_act_half = len(soc_active_indices)

                b_core_a = np.zeros((2 * n_ao, 0), dtype=complex)
                b_core_b = np.zeros((2 * n_ao, 0), dtype=complex)
                b_act = np.zeros((2 * n_ao, 0), dtype=complex)
                b_virt_a = np.zeros((2 * n_ao, 0), dtype=complex)
                b_virt_b = np.zeros((2 * n_ao, 0), dtype=complex)

                m_ca = plot_keep[:n_c]
                if np.any(m_ca):
                    b_core_a = np.zeros((2 * n_ao, np.sum(m_ca)), dtype=complex)
                    b_core_a[:n_ao, :] = C_dense[:, core_idx[m_ca]]

                m_cb = plot_keep[n_c : 2 * n_c]
                if np.any(m_cb):
                    b_core_b = np.zeros((2 * n_ao, np.sum(m_cb)), dtype=complex)
                    b_core_b[n_ao:, :] = C_dense[:, core_idx[m_cb]]

                m_act = plot_keep[2 * n_c : 2 * n_c + n_act]
                if np.any(m_act):
                    U_act_kept = soc_U_act[:, m_act]
                    C_act = C_dense[:, soc_active_indices]
                    C_act_top = _matmul_real_matrix(C_act, U_act_kept[:n_act_half, :], device=device)
                    C_act_bot = _matmul_real_matrix(C_act, U_act_kept[n_act_half:, :], device=device)
                    b_act = np.vstack([C_act_top, C_act_bot])

                m_va = plot_keep[2 * n_c + n_act : 2 * n_c + n_act + n_v]
                if np.any(m_va):
                    b_virt_a = np.zeros((2 * n_ao, np.sum(m_va)), dtype=complex)
                    b_virt_a[:n_ao, :] = C_dense[:, virt_idx[m_va]]

                m_vb = plot_keep[2 * n_c + n_act + n_v :]
                if np.any(m_vb):
                    b_virt_b = np.zeros((2 * n_ao, np.sum(m_vb)), dtype=complex)
                    b_virt_b[n_ao:, :] = C_dense[:, virt_idx[m_vb]]

                blocks = [b_core_a, b_core_b, b_act, b_virt_a, b_virt_b]
                C_spinor_ao_plot = np.hstack([b for b in blocks if b.shape[1] > 0])

            SC_spinor_top = _matmul_real_matrix(S_dense, C_spinor_ao_plot[:n_ao, :], device=device)
            SC_spinor_bot = _matmul_real_matrix(S_dense, C_spinor_ao_plot[n_ao:, :], device=device)
            SC_spinor_ao_plot = np.vstack([SC_spinor_top, SC_spinor_bot])

            C_spinor_ao = C_spinor_ao_plot[:, sort_idx]
            SC_spinor_ao = SC_spinor_ao_plot[:, sort_idx]

            pops_soc_full = np.real(C_spinor_ao[:n_ao, :].conj() * SC_spinor_ao[:n_ao, :]) + \
                            np.real(C_spinor_ao[n_ao:, :].conj() * SC_spinor_ao[n_ao:, :])
            logger.debug(f"  [PDOS/COOP] Pre-computed populations in {time.time() - t_pop:.2f}s")
            
            pdos_analysis_soc = compute_pdos_and_coop(
                C_spinor_ao, S_dense, eps_soc, shells, args.pdos_atoms, args.coop_pairs, soc_ewin,
                sigma=pdos_sigma_use, is_soc=True, prefix="soc", pops=pops_soc_full,
                population_bars=getattr(args, "population_bars", None),
                device=device, export=export_files
            )
            if store is not None:
                from qdex.store import put_orbitals
                put_orbitals(store, "soc/spinor", eps_soc, (eps_soc <= 0.0).astype(float), pdos_analysis_soc, shells,
                             args.coop_pairs, soc_ewin, pdos_sigma_use, spin_factor=1.0)
                store.attr("electronic", "soc/spinor", energy_reference="QDEX fuzzy axis: mid-gap at 0 (DFT, eV)",
                           window="spinors plotted in the fuzzy window (ewin)")

            if dashboard_energy_mode in ("qp", "both") and eps_soc_qp is not None and sort_idx_qp is not None:
                C_spinor_ao_qp = C_spinor_ao_plot[:, sort_idx_qp]
                SC_spinor_ao_qp = SC_spinor_ao_plot[:, sort_idx_qp]
                pops_soc_qp = np.real(C_spinor_ao_qp[:n_ao, :].conj() * SC_spinor_ao_qp[:n_ao, :]) + \
                               np.real(C_spinor_ao_qp[n_ao:, :].conj() * SC_spinor_ao_qp[n_ao:, :])
                compute_pdos_and_coop(
                    C_spinor_ao_qp, S_dense, eps_soc_qp, shells, args.pdos_atoms, args.coop_pairs, soc_qp_ewin,
                    sigma=pdos_sigma_use, is_soc=True, prefix="soc_qp", pops=pops_soc_qp,
                    population_bars=getattr(args, "population_bars", None),
                    device=device
                )

    # --- 3. Generate Multi-Row Interactive Plotly HTML ---
    if getattr(args, 'html', True):
        from qdex.plot_fuzzy import generate_interactive_plot
        ef_dict = {"sf": 0.0}; homo_dict = {"sf": e_homo}; lumo_dict = {"sf": e_lumo}
        
        if args.soc_flag:
            ef_dict["soc"] = 0.0
            homo_dict["soc"] = eps_soc[global_spinor_homo_idx]
            lumo_dict["soc"] = eps_soc[global_spinor_homo_idx + 1]

        # bulk-band overlay of the DFT dashboards (fuzzy.bulk_overlay / bulk_alignment / bulk_unfolded)
        bulk_opts = dict(bulk_overlay=getattr(args, "bulk_overlay", True),
                         bulk_alignment=getattr(args, "bulk_alignment", "core_level"),
                         bulk_unfolded=getattr(args, "bulk_unfolded", "auto"))
        if dashboard_energy_mode in ("dft", "both"):
            generate_interactive_plot(
                prefix="sf",
                material=args.material,
                ef=ef_dict.get("sf", 0.0),
                e_homo=homo_dict.get("sf"),
                e_lumo=lumo_dict.get("sf"),
                normalize_coop=False,
                energy_label="DFT MO energy (eV)",
                output_html="fuzzy_dashboard_sf.html",
                bulk_semicore_rel=semicore_rel.get("sf"), bulk_cif=args.cif, bulk_semicore_label=semicore_label, **bulk_opts
            )

        if dashboard_energy_mode in ("qp", "both") and qp_plot_energies is not None:
            generate_interactive_plot(
                prefix="sf_qp",
                material=args.material,
                ef=0.0 if qp_energy_reference == "fermi" else None,
                e_homo=qp_plot_energies[homo_index],
                e_lumo=qp_plot_energies[homo_index + 1],
                normalize_coop=False,
                energy_label="QP energy vs vacuum (eV)" if qp_energy_reference == "vacuum" else "QP-corrected energy relative to Fermi (eV)",
                output_html="fuzzy_dashboard_sf_qp.html",
                bulk_overlay=False  # the bulk bands are DFT; no QP axis for them
            )

        if is_uks and dashboard_energy_mode in ("dft", "both"):
            generate_interactive_plot(
                prefix="uks",
                material=args.material,
                ef=0.0,
                e_homo=max(eps_shifted[homo_index], eps_beta_shifted[homo_index_beta]),
                e_lumo=min(eps_shifted[homo_index + 1], eps_beta_shifted[homo_index_beta + 1]),
                normalize_coop=False,
                energy_label="DFT MO energy (eV)",
                output_html="fuzzy_dashboard_uks.html",
                bulk_semicore_rel=semicore_rel.get("sf"), bulk_cif=args.cif, bulk_semicore_label=semicore_label, **bulk_opts
            )

        if is_uks and dashboard_energy_mode in ("qp", "both") and qp_plot_energies is not None and qp_plot_energies_beta is not None:
            generate_interactive_plot(
                prefix="uks_qp",
                material=args.material,
                ef=0.0 if qp_energy_reference == "fermi" else None,
                e_homo=max(qp_plot_energies[homo_index], qp_plot_energies_beta[homo_index_beta]),
                e_lumo=min(qp_plot_energies[homo_index + 1], qp_plot_energies_beta[homo_index_beta + 1]),
                normalize_coop=False,
                energy_label="QP energy vs vacuum (eV)" if qp_energy_reference == "vacuum" else "QP-corrected energy relative to Fermi (eV)",
                output_html="fuzzy_dashboard_uks_qp.html",
                bulk_overlay=False  # the bulk bands are DFT; no QP axis for them
            )
    
        # Generate SOC Dashboard (if requested)
        if args.soc_flag and dashboard_energy_mode in ("dft", "both"):
            generate_interactive_plot(
                prefix="soc", 
                material=args.material, 
                ef=ef_dict.get("soc", 0.0), 
                e_homo=homo_dict.get("soc"), 
                e_lumo=lumo_dict.get("soc"), 
                normalize_coop=False,
                energy_label="DFT MO energy (eV)",
                output_html="fuzzy_dashboard_soc.html",
                bulk_semicore_rel=semicore_rel.get("soc"), bulk_cif=args.cif, bulk_semicore_label=semicore_label, **bulk_opts
            )

        if args.soc_flag and dashboard_energy_mode in ("qp", "both") and eps_soc_qp is not None:
            generate_interactive_plot(
                prefix="soc_qp",
                material=args.material,
                ef=0.0 if qp_energy_reference == "fermi" else None,
                e_homo=eps_soc_qp[global_spinor_homo_idx],
                e_lumo=eps_soc_qp[global_spinor_homo_idx + 1],
                normalize_coop=False,
                energy_label="QP+SOC energy vs vacuum (eV)" if qp_energy_reference == "vacuum" else "QP+SOC energy relative to Fermi (eV)",
                output_html="fuzzy_dashboard_soc_qp.html",
                bulk_overlay=False  # the bulk bands are DFT; no QP axis for them
            )
