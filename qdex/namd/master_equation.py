import numpy as np
from scipy.linalg import expm
from qdex.namd.integrator import HBAR_EV_FS, KB_EV

def compute_rate_matrix(E_vec, d_mat, temp_k=300.0, tau_dec_fs=15.0):
    """
    Computes the transition rate matrix R for the Pauli Master Equation::

      k_{I -> J} = 2 * |d_IJ|^2 * [tau_dec / (1 + (dE * tau_dec / hbar)^2)] * B_{IJ}

    where B_{IJ} = min(1, exp(-max(0, E_J - E_I) / (kB * T))) enforces detailed balance.
    
    The rate matrix R satisfies::

      dP_I / dt = sum_J (k_{J -> I} * P_J - k_{I -> J} * P_I)
      dP / dt = R @ P

    where R[I, J] = k_{J -> I} for I != J, and R[I, I] = - sum_{J != I} k_{I -> J}.
    """
    n_states = len(E_vec)
    R = np.zeros((n_states, n_states), dtype=np.float64)
    beta = 1.0 / (KB_EV * max(temp_k, 1e-3))

    for I in range(n_states):
        for J in range(n_states):
            if I == J:
                continue
            d_val = np.abs(d_mat[I, J])
            if d_val < 1e-12:
                continue

            dE = E_vec[J] - E_vec[I]
            lorentzian = tau_dec_fs / (1.0 + (dE * tau_dec_fs / HBAR_EV_FS) ** 2)

            # Detailed balance: Boltzmann factor on upward transitions
            boltz = np.exp(-dE * beta) if dE > 0.0 else 1.0
            boltz = min(1.0, max(boltz, 0.0))

            k_IJ = 2.0 * (d_val ** 2) * lorentzian * boltz
            R[J, I] += k_IJ  # transition from I to J adds to P_J

    # Diagonal elements: loss from state I
    for I in range(n_states):
        R[I, I] = - np.sum(R[:, I]) + R[I, I]

    return R


def run_master_equation_step(P_current, R_matrix, dt_fs):
    """
    Propagates population vector P over dt_fs using matrix exponential::

      P(t + dt) = expm(R * dt) @ P(t)
    """
    prop = expm(R_matrix * dt_fs)
    P_next = prop @ P_current
    # Ensure positivity and normalization
    P_next = np.maximum(P_next, 0.0)
    norm = np.sum(P_next)
    if norm > 1e-12:
        P_next /= norm
    return P_next


def propagate_pme_tensor(
    P_mat,
    E_mat,
    d_occ,
    d_virt,
    dt_fs,
    temp_k=300.0,
    tau_dec_fs=None,
    n_substeps=20,
    eps_occ=None,
    eps_virt=None,
    n_atoms=775,
    k_loss=None,
    return_flux=False
):
    """
    Propagates exciton population matrix P_mat of shape (n_occ, n_virt) via the Pauli Master Equation.
    Fully vectorized using BLAS matrix operations: runs in ~0.2 s per step.
    Supports on-the-fly state-dependent EDC decoherence: ``tau_kj = hbar / |dE_kj| * (1 + C / E_kin)``.
    If return_flux=True, also returns instantaneous transition probability fluxes (flux_virt, flux_occ).
    """
    n_occ, n_virt = P_mat.shape
    dt = dt_fs / n_substeps

    beta = 1.0 / (KB_EV * max(temp_k, 1e-3))
    P = np.array(P_mat, dtype=np.float64, copy=True)

    d_virt_sq = 2.0 * (np.abs(d_virt) ** 2)
    d_occ_sq = 2.0 * (np.abs(d_occ) ** 2)

    # Nuclear kinetic energy for EDC (3/2 * N_atoms * kB * T in eV)
    E_kin_ev = 1.5 * max(n_atoms, 1) * KB_EV * max(temp_k, 1e-3)
    c_param_ev = 2.7211386  # 0.1 Hartree
    edc_factor = 1.0 + (c_param_ev / E_kin_ev)

    # 1. Electron channel rate matrix k_virt (from a to b)
    if eps_virt is not None:
        dE_virt = eps_virt[np.newaxis, :] - eps_virt[:, np.newaxis]
    else:
        E_row = np.mean(E_mat, axis=0)
        dE_virt = E_row[np.newaxis, :] - E_row[:, np.newaxis]

    if tau_dec_fs is None or str(tau_dec_fs).lower() in ("edc", "auto", "dynamic", "on_the_fly"):
        # On-the-fly state-dependent EDC decoherence
        abs_dE_v = np.abs(dE_virt)
        tau_v = np.minimum(HBAR_EV_FS / np.maximum(abs_dE_v, 1e-4) * edc_factor, 50.0)
        lor_virt = tau_v / (1.0 + (dE_virt * tau_v / HBAR_EV_FS) ** 2)
    else:
        tau_val = float(tau_dec_fs)
        lor_virt = tau_val / (1.0 + (dE_virt * tau_val / HBAR_EV_FS) ** 2)

    boltz_virt = np.exp(-np.maximum(dE_virt, 0.0) * beta)
    k_virt = d_virt_sq * lor_virt * boltz_virt
    np.fill_diagonal(k_virt, 0.0)
    loss_virt = np.sum(k_virt, axis=1)  # (n_virt,)

    # 2. Hole channel rate matrix k_occ (from i to j)
    # Exciton energy change: dE = E_j - E_i = eps_occ[i] - eps_occ[j]
    if eps_occ is not None:
        dE_occ = eps_occ[:, np.newaxis] - eps_occ[np.newaxis, :]
    else:
        E_col = np.mean(E_mat, axis=1)
        dE_occ = E_col[np.newaxis, :] - E_col[:, np.newaxis]

    if tau_dec_fs is None or str(tau_dec_fs).lower() in ("edc", "auto", "dynamic", "on_the_fly"):
        abs_dE_o = np.abs(dE_occ)
        tau_o = np.minimum(HBAR_EV_FS / np.maximum(abs_dE_o, 1e-4) * edc_factor, 50.0)
        lor_occ = tau_o / (1.0 + (dE_occ * tau_o / HBAR_EV_FS) ** 2)
    else:
        tau_val = float(tau_dec_fs)
        lor_occ = tau_val / (1.0 + (dE_occ * tau_val / HBAR_EV_FS) ** 2)

    boltz_occ = np.exp(-np.maximum(dE_occ, 0.0) * beta)
    k_occ = d_occ_sq * lor_occ * boltz_occ
    np.fill_diagonal(k_occ, 0.0)
    loss_occ = np.sum(k_occ, axis=1)  # (n_occ,)

    def _row_stochastic(k_from_to, dt_sub):
        """Transition matrix T[i, j] = probability of i -> j. Rows sum to 1."""
        n = k_from_to.shape[0]
        T = dt_sub * np.array(k_from_to, dtype=np.float64, copy=True)
        np.fill_diagonal(T, 0.0)
        loss = np.sum(T, axis=1)
        stay = 1.0 - loss
        overflow = stay < 0.0
        if np.any(overflow):
            scale = 1.0 / np.maximum(loss[overflow], 1e-30)
            T[overflow] *= scale[:, np.newaxis]
            stay[overflow] = 0.0
        T[np.arange(n), np.arange(n)] = np.maximum(stay, 0.0)
        return T

    # Positivity-preserving channel updates. Each row-stochastic step
    # conserves probability. Recombination is applied once afterwards.
    T_virt = _row_stochastic(k_virt, dt)
    T_occ = _row_stochastic(k_occ, dt)
    for step in range(n_substeps):
        P = P @ T_virt
        P = T_occ.T @ P

    flux_virt = None
    flux_occ = None
    if return_flux:
        p_virt_init = np.sum(P_mat, axis=0)
        p_occ_init = np.sum(P_mat, axis=1)

        def _extract_sparse_flux(p_vec, k_mat, max_transitions=50):
            pop_idx = np.where(p_vec > 1e-6)[0]
            if len(pop_idx) == 0:
                return []
            sub_flux = p_vec[pop_idx, np.newaxis] * (k_mat[pop_idx, :] * dt_fs)
            tot_flux = np.sum(sub_flux)
            if tot_flux < 1e-14:
                return []
            thresh = 1e-4 * tot_flux
            r, c = np.where(sub_flux >= thresh)
            orig_r = pop_idx[r]
            off_diag = orig_r != c
            orig_r = orig_r[off_diag]
            c = c[off_diag]
            v = sub_flux[r[off_diag], c]
            if len(v) > max_transitions:
                top_k = np.argpartition(v, -max_transitions)[-max_transitions:]
                orig_r, c, v = orig_r[top_k], c[top_k], v[top_k]
            return [(int(i), int(j), float(val)) for i, j, val in zip(orig_r, c, v)]

        flux_virt = _extract_sparse_flux(p_virt_init, k_virt, max_transitions=50)
        flux_occ = _extract_sparse_flux(p_occ_init, k_occ, max_transitions=50)

    if k_loss is not None:
        P *= np.exp(-k_loss * dt_fs)

    if return_flux:
        return P, flux_virt, flux_occ
    return P

