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


def orbital_populations(C_cols, S, atom_ao_ranges, mode="mulliken", representation="atom", SC=None):
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
        if SC is None:
            S_d = S.toarray() if hasattr(S, "toarray") else S
            SC = S_d @ C_cols
        pop = C_cols * SC
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
    Orbitals outside the window (if a window is used) get the shift of the nearest window edge.
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
            # vectorized compute_dynamic_z: Z = 1 / (1 + |sigma| / omega_tilde), clipped to [0.5, 1]
            from qdex.hardness import valence_plasmon_ev
            eps_val = max(1.01, float(eps_z) if eps_z is not None else 1.01)
            omega = valence_plasmon_ev(material) / np.sqrt(1.0 - 1.0 / eps_val)
            z_o = np.clip(1.0 / (1.0 + np.abs(sig_o) / omega), 0.5, 1.0)
            z_v = np.clip(1.0 / (1.0 + np.abs(sig_v) / omega), 0.5, 1.0)
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


def cohsex_diagonal(C, S, homo_index, atom_ao_ranges, dW_atom=None, dW_ao=None):
    """Static Delta-COHSEX diagonal <n|Sigma|n> for ALL orbitals (one-shot, no orbital update).

    Same form as the qsGW models, in the Löwdin basis (c = S^1/2 C):

        COH_n = 1/2 sum_mu c_{mu n}^2 dW_{mu mu}
        SEX_n = -1/2 sum_{mu nu} c_{mu n} c_{nu n} P_{mu nu} dW_{mu nu},   P = 2 c_occ c_occ^T

    dW is Delta W in the atom-block (mnok, ``dW_atom``) or AO (xs, ``dW_ao``)
    representation.  The classical charging term 1/2 q^T dW q is the limit in
    which the occupied states act as a complete set (SEX -> -dW(r, r)); the
    difference is the non-classical screened exchange.  Cost: one cached
    eigh(S) and three n^3 matrix products; memory about four n_ao^2 arrays.
    Returns (coh, sex), each of length n_mo.
    """
    from qdex.lowdin import lowdin_apply
    C = C.toarray() if hasattr(C, "toarray") else np.asarray(C, dtype=float)
    n_ao = C.shape[0]
    Cl = lowdin_apply(S, C)
    n_occ = homo_index + 1
    P = Cl[:, :n_occ] @ Cl[:, :n_occ].T
    P *= 2.0
    if dW_ao is None:
        owner = _owner(atom_ao_ranges, n_ao)
        dW_atom = np.asarray(dW_atom, dtype=float)
        block = 512
        for i0 in range(0, n_ao, block):
            i1 = min(n_ao, i0 + block)
            P[i0:i1] *= dW_atom[np.ix_(owner[i0:i1], owner)]
        diag_dW = dW_atom[owner, owner]
    else:
        P *= dW_ao
        diag_dW = np.diag(dW_ao).copy()
    MC = P @ Cl
    del P
    sex = -0.5 * np.einsum("mn,mn->n", Cl, MC, optimize=True)
    del MC
    coh = 0.5 * ((Cl * Cl).T @ diag_dW)
    return coh, sex


def plasmon_pole_z(sigma, eps_z, material):
    """Vectorized compute_dynamic_z: Z = 1 / (1 + |Sigma| / omega_tilde), clipped to [0.5, 1]."""
    from qdex.hardness import valence_plasmon_ev
    eps_val = max(1.01, float(eps_z) if eps_z is not None else 1.01)
    omega = valence_plasmon_ev(material) / np.sqrt(1.0 - 1.0 / eps_val)
    return np.clip(1.0 / (1.0 + np.abs(sigma) / omega), 0.5, 1.0)


def cohsex_qp_energies(eps, coh, sex, homo_index, bulk_shift, z_mode, z_fixed, eps_z, material,
                       homo_fraction_bulk=0.5):
    """QP energies of all orbitals from the one-shot Delta-COHSEX diagonal.

    occupied:  e_n - f_b D_bulk + Z_n (COH_n + SEX_n)
    virtual:   e_n + (1 - f_b) D_bulk + Z_n (COH_n + SEX_n)
    """
    eps = np.asarray(eps, dtype=float)
    sig = np.asarray(coh) + np.asarray(sex)
    z = plasmon_pole_z(sig, eps_z, material) if z_mode == "derived" else np.full_like(sig, float(z_fixed))
    fb = float(homo_fraction_bulk)
    bulk = np.where(np.arange(len(eps)) <= homo_index, -fb * bulk_shift, (1.0 - fb) * bulk_shift)
    eps_qp = eps + bulk + z * sig
    h, l = homo_index, homo_index + 1
    occ, vir = slice(0, h + 1), slice(l, len(eps))
    info = {
        "qp_levels": "orbital",
        "qp_selfenergy": "cohsex",
        "qp_homo_shift_ev": float(eps_qp[h] - eps[h]),
        "qp_lumo_shift_ev": float(eps_qp[l] - eps[l]),
        "cohsex_homo_coh_ev": float(coh[h]), "cohsex_homo_sex_ev": float(sex[h]),
        "cohsex_lumo_coh_ev": float(coh[l]), "cohsex_lumo_sex_ev": float(sex[l]),
        "qp_shift_spread_occ_ev": float(np.ptp((eps_qp - eps)[occ])),
        "qp_shift_spread_virt_ev": float(np.ptp((eps_qp - eps)[vir])),
        "z_homo_orbital": float(z[h]), "z_lumo_orbital": float(z[l]),
        "z_min_window": float(z.min()), "z_max_window": float(z.max()),
    }
    return eps_qp, info


# ---------------------------------------------------------------------------
# Anchor residuals of the Delta-W models (calibrated on the evGW anchor cluster)
# ---------------------------------------------------------------------------
import json as _json
import os as _os

ANCHOR_TABLE = _os.path.join(_os.path.dirname(__file__), "data", "dw_anchor_residuals.json")


def anchor_key(material, qp_model, selfenergy, representation, populations, z_label, solvent_term="sphere"):
    return "|".join(str(x).lower() for x in (material, qp_model, selfenergy, representation, populations, z_label,
                                             solvent_term))


def load_anchor_table(path=None):
    path = path or ANCHOR_TABLE
    if not _os.path.exists(path):
        return {}
    with open(path) as fh:
        return _json.load(fh)


def save_anchor_entry(key, entry, path=None):
    path = path or ANCHOR_TABLE
    table = load_anchor_table(path)
    table[key] = entry
    _os.makedirs(_os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        _json.dump(table, fh, indent=2, sort_keys=True)
    return path
