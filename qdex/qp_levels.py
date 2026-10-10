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
        n_atoms = len(atom_ao_ranges)
        owner = _owner(atom_ao_ranges, n_ao)
        import scipy.sparse as sp
        P_atom = sp.csr_matrix((np.ones(n_ao, dtype=pop.dtype), (owner, np.arange(n_ao))),
                               shape=(n_atoms, n_ao))
        return P_atom @ pop
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
        "bulk_shifts": np.where(np.arange(len(eps)) <= homo, -float(homo_fraction_bulk) * bulk_shift, (1.0 - float(homo_fraction_bulk)) * bulk_shift),
        "z_factors": np.concatenate([np.full(occ_idx[0], z_o[0]), z_o, np.full(lumo - homo - 1, z_o[-1]), z_v, np.full(len(eps) - virt_idx[-1] - 1, z_v[-1])]),
        "delta_sigmas": np.concatenate([np.full(occ_idx[0], sig_o[0]), sig_o, np.full(lumo - homo - 1, sig_o[-1]), sig_v, np.full(len(eps) - virt_idx[-1] - 1, sig_v[-1])]),
    }
    return eps_qp, info


def cohsex_diagonal(C, S, homo_index, atom_ao_ranges, dW_atom=None, dW_ao=None, eval_indices=None, occ_density=None):
    """Static Delta-COHSEX diagonal <n|Sigma|n> in the Löwdin basis (c = S^1/2 C):

        COH_n = 1/2 sum_mu c_{mu n}^2 dW_{mu mu}
        SEX_n = -1/2 sum_{mu nu} c_{mu n} c_{nu n} P_{mu nu} dW_{mu nu},   P = 2 c_occ c_occ^T

    dW is Delta W in the atom-block (mnok, ``dW_atom``) or AO (xs, ``dW_ao``)
    representation.  The classical charging term 1/2 q^T dW q is the limit in
    which the occupied states act as a complete set (SEX -> -dW(r, r)); the
    difference is the non-classical screened exchange.

    Parameters
    ----------
    eval_indices : array-like of int, optional
        If provided, evaluates the self-energy diagonal specifically for the
        requested orbital indices (e.g. active space around the Fermi level),
        returning (coh, sex) of length n_mo where eval_indices have exact values
        and outer indices are clamped to the window edges.
        If None, evaluates all n_mo orbitals in memory-efficient chunks.
    occ_density : array-like of int, optional
        Occupied orbitals that build the screened-exchange density matrix P (default: all occupied). The SEX term is
        -sum_i <ni|dW|ni> over these i; states far below the gap have transition densities phi_n phi_i of zero net
        charge, which a smooth dW sees little (quasiparticles.cohsex_occ_window).
    """
    from qdex.lowdin import lowdin_apply
    C = C.toarray() if hasattr(C, "toarray") else np.asarray(C, dtype=float)
    n_ao = C.shape[0]
    n_mo = C.shape[1]
    n_occ = homo_index + 1

    # Form the 1-RDM P = 2 c_occ c_occ^T in the Löwdin basis (only occupied orbitals are needed)
    if occ_density is not None and eval_indices is not None:
        # truncated density: the occupied orbitals of P plus the occupied orbitals to evaluate, Löwdin-transformed once
        occ_d = np.asarray(occ_density, dtype=int)
        eval_occ = np.asarray(eval_indices, dtype=int)
        cols = np.union1d(occ_d, eval_occ[eval_occ < n_occ])
        Cl_cols = lowdin_apply(S, C[:, cols])
        pos = {int(c): k for k, c in enumerate(cols)}
        Cl_d = Cl_cols[:, [pos[int(i)] for i in occ_d]]
        P = Cl_d @ Cl_d.T
        P *= 2.0
        del Cl_d
        Cl_occ = None
        Cl_occ_map = (Cl_cols, pos)
    else:
        Cl_occ = lowdin_apply(S, C[:, :n_occ])
        P = Cl_occ @ Cl_occ.T
        P *= 2.0
        Cl_occ_map = None

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

    if eval_indices is not None:
        eval_idx = np.asarray(eval_indices, dtype=int)
        occ_mask = eval_idx < n_occ
        virt_indices = eval_idx[~occ_mask]

        if len(virt_indices) > 0:
            Cl_virt = lowdin_apply(S, C[:, virt_indices])
        else:
            Cl_virt = np.empty((n_ao, 0), dtype=float)

        n_eval = len(eval_idx)
        Cl_eval = np.empty((n_ao, n_eval), dtype=float)
        if np.any(occ_mask) and Cl_occ_map is not None:
            Cl_cols, pos = Cl_occ_map
            Cl_eval[:, occ_mask] = Cl_cols[:, [pos[int(i)] for i in eval_idx[occ_mask]]]
        elif np.any(occ_mask):
            Cl_eval[:, occ_mask] = Cl_occ[:, eval_idx[occ_mask]]
        if len(virt_indices) > 0:
            Cl_eval[:, ~occ_mask] = Cl_virt

        del Cl_occ, Cl_occ_map
        del Cl_virt

        MC_eval = P @ Cl_eval
        del P
        sex_eval = -0.5 * np.sum(Cl_eval * MC_eval, axis=0)
        del MC_eval
        coh_eval = 0.5 * ((Cl_eval ** 2).T @ diag_dW)
        del Cl_eval

        coh = np.empty(n_mo, dtype=float)
        sex = np.empty(n_mo, dtype=float)
        coh[eval_idx] = coh_eval
        sex[eval_idx] = sex_eval

        # Edge clamping for indices outside the evaluated window
        eval_sorted = np.sort(eval_idx)
        occ_eval_sorted = eval_sorted[eval_sorted < n_occ]
        virt_eval_sorted = eval_sorted[eval_sorted >= n_occ]

        if len(occ_eval_sorted) > 0:
            min_occ = occ_eval_sorted[0]
            max_occ = occ_eval_sorted[-1]
            coh[:min_occ] = coh[min_occ]
            sex[:min_occ] = sex[min_occ]
            if max_occ < homo_index:
                coh[max_occ + 1: n_occ] = coh[max_occ]
                sex[max_occ + 1: n_occ] = sex[max_occ]

        if len(virt_eval_sorted) > 0:
            min_virt = virt_eval_sorted[0]
            max_virt = virt_eval_sorted[-1]
            if min_virt > n_occ:
                coh[n_occ: min_virt] = coh[min_virt]
                sex[n_occ: min_virt] = sex[min_virt]
            coh[max_virt + 1:] = coh[max_virt]
            sex[max_virt + 1:] = sex[max_virt]

        return coh, sex

    else:
        # Full evaluation: lowdin transform remaining columns
        if n_occ < n_mo:
            Cl_virt = lowdin_apply(S, C[:, n_occ:])
            Cl = np.hstack([Cl_occ, Cl_virt])
            del Cl_occ
            del Cl_virt
        else:
            Cl = Cl_occ

        coh = np.empty(n_mo, dtype=float)
        sex = np.empty(n_mo, dtype=float)
        chunk = 1024
        for j0 in range(0, n_mo, chunk):
            j1 = min(n_mo, j0 + chunk)
            Cl_blk = Cl[:, j0:j1]
            MC_blk = P @ Cl_blk
            sex[j0:j1] = -0.5 * np.sum(Cl_blk * MC_blk, axis=0)
            coh[j0:j1] = 0.5 * ((Cl_blk ** 2).T @ diag_dW)
        del P
        del Cl
        return coh, sex


def plasmon_pole_z(sigma, eps_z, material):
    """Vectorized compute_dynamic_z: Z = 1 / (1 + |Sigma| / omega_tilde), clipped to [0.5, 1]."""
    from qdex.hardness import valence_plasmon_ev
    eps_val = max(1.01, float(eps_z) if eps_z is not None else 1.01)
    omega = valence_plasmon_ev(material) / np.sqrt(1.0 - 1.0 / eps_val)
    return np.clip(1.0 / (1.0 + np.abs(sigma) / omega), 0.5, 1.0)


def cohsex_qp_energies(eps, coh, sex, homo_index, bulk_shift, z_mode, z_fixed, eps_z, material,
                       homo_fraction_bulk=0.5, occ_idx=None, virt_idx=None):
    """QP energies of all orbitals from the one-shot Delta-COHSEX diagonal.

    occupied:  e_n - f_b D_bulk + Z_n (COH_n + SEX_n)
    virtual:   e_n + (1 - f_b) D_bulk + Z_n (COH_n + SEX_n)
    """
    eps = np.asarray(eps, dtype=float)
    sig = np.asarray(coh) + np.asarray(sex)
    z = plasmon_pole_z(sig, eps_z, material) if z_mode == "derived" else np.full_like(sig, float(z_fixed))
    fb = float(homo_fraction_bulk)
    bulk = np.where(np.arange(len(eps)) <= homo_index, -fb * bulk_shift, (1.0 - fb) * bulk_shift)
    shifts = bulk + z * sig
    eps_qp = eps.copy()
    h, l = homo_index, homo_index + 1

    if occ_idx is not None and virt_idx is not None:
        occ_idx = np.asarray(occ_idx, dtype=int)
        virt_idx = np.asarray(virt_idx, dtype=int)
        shift_o = shifts[occ_idx]
        shift_v = shifts[virt_idx]
        eps_qp[: occ_idx[0]] += shift_o[0]
        eps_qp[occ_idx] += shift_o
        eps_qp[h + 1: virt_idx[0]] += shift_o[-1]
        eps_qp[virt_idx] += shift_v
        eps_qp[virt_idx[-1] + 1:] += shift_v[-1]
        spread_occ = float(np.ptp(shift_o))
        spread_virt = float(np.ptp(shift_v))
        z_min = float(min(z[occ_idx].min(), z[virt_idx].min()))
        z_max = float(max(z[occ_idx].max(), z[virt_idx].max()))
        w_range = [int(occ_idx[0]), int(virt_idx[-1])]
    else:
        eps_qp += shifts
        occ, vir = slice(0, h + 1), slice(l, len(eps))
        spread_occ = float(np.ptp(shifts[occ]))
        spread_virt = float(np.ptp(shifts[vir]))
        z_min = float(z.min())
        z_max = float(z.max())
        w_range = [0, len(eps) - 1]

    info = {
        "qp_levels": "orbital",
        "qp_selfenergy": "cohsex",
        "qp_levels_window": w_range,
        "qp_homo_shift_ev": float(eps_qp[h] - eps[h]),
        "qp_lumo_shift_ev": float(eps_qp[l] - eps[l]),
        "cohsex_homo_coh_ev": float(coh[h]), "cohsex_homo_sex_ev": float(sex[h]),
        "cohsex_lumo_coh_ev": float(coh[l]), "cohsex_lumo_sex_ev": float(sex[l]),
        "qp_shift_spread_occ_ev": spread_occ,
        "qp_shift_spread_virt_ev": spread_virt,
        "z_homo_orbital": float(z[h]), "z_lumo_orbital": float(z[l]),
        "z_min_window": z_min, "z_max_window": z_max,
        "bulk_shifts": bulk,
        "z_factors": z,
        "delta_sigmas": sig,
    }
    return eps_qp, info


def resta_edge_shifts(C, S, homo_index, atom_ao_ranges, coords, atom_symbols, material, eps_out, dft_gap,
                      bulk_shift, homo_fraction_bulk, alpha=1.0, solvent_term="sphere", selfenergy="cohsex",
                      z_mode="derived", z_fixed=0.8, populations="mulliken"):
    """HOMO and LUMO QP shifts of the sgw-resta model, for the absolute edges (IP/EA) of a gap-only model.

    The shifts are those of the orbital-resolved sgw-resta levels (one shot, Penn-scaled Resta W,
    sphere reaction field, atom representation), with the bulk correction of the calling model:

        HOMO:  -f_b D_bulk       + Z_H Sigma_H
        LUMO:  (1 - f_b) D_bulk  + Z_L Sigma_L

    Sigma_n = COH_n + SEX_n (``cohsex_diagonal``) or -/+ 1/2 q_n^T dW q_n (classical). The finite-size
    part is what the charged states (IP, EA) feel; in the optical gap it cancels against the
    electron-hole binding. Returns (d_homo, d_lumo, info).
    """
    from qdex.hardness import _ScreeningModel
    h, l = int(homo_index), int(homo_index) + 1
    model = _ScreeningModel("resta", np.asarray(coords, dtype=float), atom_symbols, material, eps_out,
                            solvent_term, alpha=alpha, penn_scaling=True)
    w_qd, eps_z = model.one_shot(float(dft_gap))
    dW = atom_delta_w(model.w_parts(w_qd, eps_z))
    if str(selfenergy).lower() == "cohsex":
        coh, sex = cohsex_diagonal(C, S, h, atom_ao_ranges, dW_atom=dW, eval_indices=np.array([h, l]))
        sig = np.array([coh[h] + sex[h], coh[l] + sex[l]])
    else:
        C_cols = C[:, [h, l]]
        q = orbital_populations(C_cols, S, atom_ao_ranges, populations, "atom")
        s = orbital_sigma(q, dW)
        sig = np.array([-s[0], s[1]])
    z = plasmon_pole_z(sig, eps_z, material) if z_mode == "derived" else np.full(2, float(z_fixed))
    fb = float(homo_fraction_bulk)
    d_h = -fb * float(bulk_shift) + float(z[0] * sig[0])
    d_l = (1.0 - fb) * float(bulk_shift) + float(z[1] * sig[1])
    info = {"ipea_edge_model": "sgw_resta_edges", "ipea_selfenergy": str(selfenergy).lower(),
            "ipea_sigma_homo_ev": float(sig[0]), "ipea_sigma_lumo_ev": float(sig[1]),
            "ipea_z_homo": float(z[0]), "ipea_z_lumo": float(z[1]), "ipea_eps_eff": float(eps_z),
            "ipea_homo_shift_ev": float(d_h), "ipea_lumo_shift_ev": float(d_l)}
    return d_h, d_l, info


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
