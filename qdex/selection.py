"""Perturbative selection of transitions (Grimme, J. Chem. Phys. 138, 244104 (2013)).

Follows ``ptselect`` of the std2 program:

* primary transitions P: diagonal element A_ii <= E_thr;
* a candidate u with A_uu > E_thr is kept when its second-order coupling to P,
  E_u = sum_{v in P} A_uv^2 / (A_uu - A_vv), exceeds t (1e-4 Eh in std2);
* the second-order contributions of the rejected candidates lower the diagonal of
  the primaries they couple to.

A_ia,jb = D_ia delta + k_x (ia|jb)_K - (ij|ab)_W, with k_x = 2 (singlet) or 0 (triplet).
Energies in eV.
"""
import numpy as np

HARTREE_EV = 27.211386245988


def perturbative_selection(D, q_ov, q_occ, W_virt, gamma_k, kx_factor, include_direct,
                           e_thr, t_pt_hartree=1e-4, base_mask=None):
    """Return (mask, diag_shift, info) over the (n_occ, n_virt) transition grid.

    D        : (n_occ, n_virt) transition energies (QP or DFT differences) in eV
    q_ov     : (n_occ, n_virt, m) transition charges
    q_occ    : (n_occ, n_occ, m) hole densities;  W_virt: (n_virt, n_virt, m) = q_virt @ W^T
    gamma_k  : (m, m) exchange interaction in eV
    """
    n_o, n_v = D.shape
    m = q_ov.shape[-1]
    base = np.ones((n_o, n_v), dtype=bool) if base_mask is None else np.asarray(base_mask, dtype=bool)
    qf = np.asarray(q_ov, dtype=float).reshape(n_o * n_v, m)
    pk = qf @ np.asarray(gamma_k, dtype=float)
    kx_diag = np.einsum("pm,pm->p", pk, qf).reshape(n_o, n_v) if kx_factor else np.zeros((n_o, n_v))
    if include_direct and W_virt is not None:
        q_ii = np.einsum("iim->im", q_occ)
        w_aa = np.einsum("aam->am", W_virt)
        kd_diag = q_ii @ w_aa.T
    else:
        kd_diag = np.zeros((n_o, n_v))
    A = D + kx_factor * kx_diag - kd_diag

    prim = base & (A <= e_thr)
    cand = base & ~prim
    ip, ap = np.nonzero(prim)
    ic, ac = np.nonzero(cand)
    n_p, n_c = len(ip), len(ic)
    shift = np.zeros((n_o, n_v))
    info = {"n_primary": int(n_p), "n_candidates": int(n_c), "n_added": 0, "pt2_mean_ev": 0.0, "pt2_max_ev": 0.0}
    if n_p == 0:
        raise ValueError(f"Perturbative selection: no transition below E_thr = {e_thr:.2f} eV; raise selection_energy.")
    if n_c == 0:
        return prim, shift, info

    fp = ip * n_v + ap
    fc = ic * n_v + ac
    coupling = np.zeros((n_c, n_p))
    if kx_factor:
        coupling += kx_factor * (pk[fc] @ qf[fp].T)
    if include_direct and W_virt is not None:
        for k in range(n_p):
            coupling[:, k] -= np.einsum("um,um->u", q_occ[ic, ip[k], :], W_virt[ac, ap[k], :])
    denom = A[ic, ac][:, None] - A[ip, ap][None, :] + 1e-10 * HARTREE_EV
    pt = coupling ** 2 / denom
    e_u = pt.sum(axis=1)
    keep = e_u > t_pt_hartree * HARTREE_EV

    mask = prim.copy()
    mask[ic[keep], ac[keep]] = True
    pt2 = pt[~keep].sum(axis=0)
    shift[ip, ap] = -pt2
    info.update(n_added=int(keep.sum()), pt2_mean_ev=float(pt2.mean()), pt2_max_ev=float(pt2.max(initial=0.0)))
    return mask, shift, info
