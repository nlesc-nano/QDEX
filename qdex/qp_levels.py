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

In the atom (mnok) representation q_p are Löwdin atomic populations; in the xs
representation they are Löwdin AO populations and Delta W_{mu nu} =
Delta S_{AB} (mu mu | nu nu) with Delta S_{AB} = Delta W_{AB} / gamma_{AB}.
"""
import numpy as np

from qdex.hardness import (
    ao_delta_w,
    ao_screening_ratio,
    compute_dynamic_z,
    expand_atom_to_ao,
)


def xs_shared_w(parts, z_eh, gamma_ao, atom_ao_ranges):
    """AO kernel W_bulk + Z (W_QD - W_bulk) in the xs representation.

    Returns (W_kernel_ao, W_bulk_ao, dW_ao_clipped) where dW_ao_clipped is the
    Delta W used for the QP levels (screening part clipped at zero, as in the
    atom-resolved models).
    """
    n_ao = gamma_ao.shape[0]
    W_bulk_ao = ao_screening_ratio(parts["w_bulk"], parts["gamma"], atom_ao_ranges, n_ao) * gamma_ao
    dW_kernel = ao_delta_w(parts["w_qd"] - parts["w_bulk"], parts["gamma"], parts.get("w_add"),
                           gamma_ao, atom_ao_ranges)
    dW_qp = ao_delta_w(np.maximum(0.0, parts["w_qd"] - parts["w_bulk"]), parts["gamma"], parts.get("w_add"),
                       gamma_ao, atom_ao_ranges)
    return W_bulk_ao + float(z_eh) * dW_kernel, W_bulk_ao, dW_qp


def atom_delta_w(parts):
    """Delta W for the QP levels in the atom representation."""
    dW = np.maximum(0.0, parts["w_qd"] - parts["w_bulk"])
    if parts.get("w_add") is not None:
        dW = dW + parts["w_add"]
    return dW


def orbital_sigma(C_low, dW, atom_ao_ranges=None, representation="atom"):
    """sigma_p = 1/2 q_p^T dW q_p for each column p of the Löwdin MOs C_low."""
    pop_ao = np.abs(C_low) ** 2
    if representation == "atom":
        q = np.array([pop_ao[a0:a1].sum(axis=0) for a0, a1 in atom_ao_ranges])
    else:
        q = pop_ao
    return 0.5 * np.einsum("ap,ab,bp->p", q, dW, q, optimize=True)


def orbital_qp_energies(eps, C_low_occ, C_low_virt, occ_idx, virt_idx, dW, atom_ao_ranges,
                        representation, bulk_shift, z_mode, z_fixed, eps_z, material,
                        homo_fraction_bulk=0.5, edge_shifts=None):
    """Return (eps_qp, info) with orbital-resolved QP energies on the window.

    Orbitals outside the window get the shift of the nearest window edge.
    ``edge_shifts`` = (dH, dL) pins the HOMO and LUMO shifts (two-anchor model):
    every orbital then gets the edge shift plus its sigma relative to the edge
    orbital (Z = 1), so the anchor-calibrated frontier levels are kept.
    """
    eps = np.asarray(eps, dtype=float)
    sig_o = orbital_sigma(C_low_occ, dW, atom_ao_ranges, representation)
    sig_v = orbital_sigma(C_low_virt, dW, atom_ao_ranges, representation)
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
