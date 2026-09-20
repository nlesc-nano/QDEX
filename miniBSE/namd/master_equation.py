import numpy as np
from scipy.linalg import expm
from miniBSE.namd.integrator import HBAR_EV_FS, KB_EV

def compute_rate_matrix(E_vec, d_mat, temp_k=300.0, tau_dec_fs=15.0):
    """
    Computes the transition rate matrix R for the Pauli Master Equation:
      k_{I -> J} = 2 * |d_IJ|^2 * [tau_dec / (1 + (dE * tau_dec / hbar)^2)] * B_{IJ}
    where B_{IJ} = min(1, exp(-max(0, E_J - E_I) / (kB * T))) enforces detailed balance.
    
    The rate matrix R satisfies:
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
    Propagates population vector P over dt_fs using matrix exponential:
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
    k_loss=None
):
    """
    Propagates exciton population matrix P_mat of shape (n_occ, n_virt) via the Pauli Master Equation.
    Fully vectorized using BLAS matrix operations: runs in ~0.2 s per step.
    Supports on-the-fly state-dependent EDC decoherence: tau_kj = hbar / |dE_kj| * (1 + C / E_kin).
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

    # Vectorized Substep Propagation via BLAS matrix multiplications
    for step in range(n_substeps):
        # Electron flux: gain from all virtuals + loss from current virtual
        flux_e = (P @ k_virt) - (P * loss_virt[np.newaxis, :])

        # Hole flux: gain from all holes + loss from current hole
        flux_h = (k_occ.T @ P) - (loss_occ[:, np.newaxis] * P)

        P += dt * (flux_e + flux_h)
        P = np.maximum(P, 0.0)
        norm = np.sum(P)
        if norm > 1e-12:
            P /= norm

    # Ground state recombination loss (radiative + non-radiative)
    if k_loss is not None:
        P *= np.exp(-k_loss * dt_fs)

    return P
