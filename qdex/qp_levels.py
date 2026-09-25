"""Orbital-resolved quasiparticle levels and the xs form of a shared W.

A Delta-W QP model defines an atom-pair screened interaction W_QD, its bulk
reference W_bulk and additive classical terms (solvent reaction field or
sphere image).  This module

* applies the same correction to every orbital p of the active window,

      occupied:  e_p -> e_p - f_b D_bulk       - Z_p sigma_p
      virtual:   e_p -> e_p + (1 - f_b) D_bulk + Z_p sigma_p
      sigma_p = 1/2 q_p^T Delta W q_p,

  with Z_p from the model's plasmon pole (``compute_dynamic_z``), and
* builds the xs (AO density-pair) representation of the same W, so that the
  QP correction and the BSE kernel use one W with either two-electron
  representation (mnok or xs).

In the atom (mnok) representation q_p are atomic populations (Mulliken or
Löwdin, following the BSE charge_type); in the xs representation they are
Löwdin AO populations and Delta W_{mu nu} =
Delta S_{AB} (mu mu | nu nu) with Delta S_{AB} = Delta W_{AB} / gamma_{AB}.
"""
import numpy as np

from qdex.hardness import compute_dynamic_z


def _owner(atom_ao_ranges, n_ao):
    owner = np.empty(n_ao, dtype=int)
    for A, (a0, a1) in enumerate(atom_ao_ranges):
        owner[a0:a1] = A
    return owner


def _screened_ao(ratio_atom, add_atom, gamma_ao, owner):
    """ratio[A(mu), B(nu)] * gamma_ao[mu, nu] + add[A(mu), B(nu)], built row block by row block.

    Peak memory is one n_ao x n_ao result plus a small block, which matters
    for 10^4 basis functions (one matrix is 8 n_ao^2 bytes).
    """
    n_ao = gamma_ao.shape[0]
    out = np.empty((n_ao, n_ao), dtype=float)
    block = 512
    for i0 in range(0, n_ao, block):
        i1 = min(n_ao, i0 + block)
        rows = owner[i0:i1]
        blk = ratio_atom[np.ix_(rows, owner)]
        blk *= gamma_ao[i0:i1]
        if add_atom is not None:
            blk += add_atom[np.ix_(rows, owner)]
        out[i0:i1] = blk
    return out


def xs_shared_w(parts, z_eh, gamma_ao, atom_ao_ranges):
    """AO kernel W_bulk + Z (W_QD - W_bulk) in the xs representation.

    Returns (W_kernel_ao, None, dW_ao_clipped) where dW_ao_clipped is the
    Delta W used for the QP levels (screening part clipped at zero, as in the
    atom-resolved models).
    """
    n_ao = gamma_ao.shape[0]
    owner = _owner(atom_ao_ranges, n_ao)
    gamma = np.asarray(parts["gamma"], dtype=float)
    w_bulk = np.asarray(parts["w_bulk"], dtype=float)
    w_qd = np.asarray(parts["w_qd"], dtype=float)
    add = parts.get("w_add")
    add = None if add is None else np.asarray(add, dtype=float)
    z = float(z_eh)
    ratio_kernel = (w_bulk + z * (w_qd - w_bulk)) / gamma
    W_kernel_ao = _screened_ao(ratio_kernel, None if add is None else z * add, gamma_ao, owner)
    ratio_qp = np.maximum(0.0, w_qd - w_bulk) / gamma
    dW_qp = _screened_ao(ratio_qp, add, gamma_ao, owner)
    return W_kernel_ao, None, dW_qp


def atom_delta_w(parts):
    """Delta W for the QP levels in the atom representation."""
    dW = np.maximum(0.0, parts["w_qd"] - parts["w_bulk"])
    if parts.get("w_add") is not None:
        dW = dW + parts["w_add"]
    return dW


def orbital_populations(C_cols, S, atom_ao_ranges, mode="mulliken", representation="atom"):
    """Populations q[feature, p] of the orbitals in the columns of C (AO basis).

    mode 'mulliken': q_mu = C_mu (S C)_mu  -- needs only S C, no diagonalization.
    mode 'lowdin'  : q_mu = (S^1/2 C)_mu^2 -- one cached eigh(S).
    representation 'atom' sums the AO populations per atom; 'ao' keeps them.
    Use the same partition as the BSE transition charges (charge_type).
    """
    C_cols = C_cols.toarray() if hasattr(C_cols, "toarray") else np.asarray(C_cols, dtype=float)
    if mode == "lowdin":
        from qdex.lowdin import lowdin_apply
        pop = np.abs(lowdin_apply(S, C_cols)) ** 2
    else:
        S_d = S.toarray() if hasattr(S, "toarray") else S
        pop = C_cols * (S_d @ C_cols)
    if representation == "atom":
        n_ao = pop.shape[0]
        owner = _owner(atom_ao_ranges, n_ao)
        q = np.zeros((len(atom_ao_ranges), pop.shape[1]))
        np.add.at(q, owner, pop)
        return q
    return pop


def orbital_sigma(q, dW):
    """sigma_p = 1/2 q_p^T dW q_p for each column p of the populations q."""
    return 0.5 * np.einsum("ap,ap->p", q, dW @ q)


def orbital_qp_energies(eps, q_occ, q_virt, occ_idx, virt_idx, dW,
                        bulk_shift, z_mode, z_fixed, eps_z, material,
                        homo_fraction_bulk=0.5, edge_shifts=None, representation="atom"):
    """Return (eps_qp, info) with orbital-resolved QP energies on the window.

    q_occ, q_virt are the populations of the window orbitals (see
    ``orbital_populations``), in the representation of dW (atom or AO).
    Orbitals outside the window get the shift of the nearest window edge.
    ``edge_shifts`` = (dH, dL) pins the HOMO and LUMO shifts (two-anchor model):
    every orbital then gets the edge shift plus its sigma relative to the edge
    orbital (Z = 1), so the anchor-calibrated frontier levels are kept.
    """
    eps = np.asarray(eps, dtype=float)
    sig_o = orbital_sigma(q_occ, dW)
    sig_v = orbital_sigma(q_virt, dW)
    if edge_shifts is not None:
        dH, dL = edge_shifts
        shift_o = -(dH + (sig_o - sig_o[-1]))
        shift_v = dL + (sig_v - sig_v[0])
        z_o = np.ones_like(sig_o)
        z_v = np.ones_like(sig_v)
    else:
        if z_mode == "derived":
            z_o = np.array([compute_dynamic_z(s, None, eps_z, material) for s in sig_o])
            z_v = np.array([compute_dynamic_z(s, None, eps_z, material) for s in sig_v])
        else:
            z_o = np.full_like(sig_o, float(z_fixed))
            z_v = np.full_like(sig_v, float(z_fixed))
        fb = float(homo_fraction_bulk)
        shift_o = -(fb * bulk_shift + z_o * sig_o)
        shift_v = (1.0 - fb) * bulk_shift + z_v * sig_v

    eps_qp = eps.copy()
    homo = int(occ_idx[-1])
    lumo = int(virt_idx[0])
    eps_qp[: occ_idx[0]] += shift_o[0]
    eps_qp[occ_idx] += shift_o
    eps_qp[homo + 1: lumo] += shift_o[-1]
    eps_qp[virt_idx] += shift_v
    eps_qp[virt_idx[-1] + 1:] += shift_v[-1]
    info = {
        "qp_levels": "orbital",
        "qp_levels_representation": representation,
        "qp_levels_window": [int(occ_idx[0]), int(virt_idx[-1])],
        "qp_homo_shift_ev": float(shift_o[-1]),
        "qp_lumo_shift_ev": float(shift_v[0]),
        "qp_shift_spread_occ_ev": float(np.ptp(shift_o)),
        "qp_shift_spread_virt_ev": float(np.ptp(shift_v)),
        "z_homo_orbital": float(z_o[-1]),
        "z_lumo_orbital": float(z_v[0]),
        "z_min_window": float(min(z_o.min(), z_v.min())),
        "z_max_window": float(max(z_o.max(), z_v.max())),
    }
    return eps_qp, info
