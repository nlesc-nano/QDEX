import numpy as np
from scipy.linalg import expm

# Constants
HBAR_EV_FS = 0.6582119569  # hbar in eV * fs
KB_EV = 8.617333262e-5      # Boltzmann constant in eV / K


def build_effective_hamiltonian(E_vec, d_mat):
    """
    Builds the effective electronic Hamiltonian in eV:
      H_eff = diag(E) - i * hbar * d_mat
    Since d_mat is anti-Hermitian (d_IJ = -d_JI*),
    -i * hbar * d_mat is Hermitian, making H_eff strictly Hermitian.
    """
    H = np.diag(E_vec).astype(np.complex128)
    if d_mat is not None:
        H -= 1j * HBAR_EV_FS * d_mat
    return H


def step_unitary_matrix_exp(c, H_eff, dt_elec):
    """
    Exact unitary step via diagonalization:
      c(t + dt) = V @ exp(-i * Lambda * dt / hbar) @ V^dagger @ c(t)
    Guarantees norm conservation to machine precision.
    """
    w, v = np.linalg.eigh(H_eff)
    phase = np.exp(-1j * w * (dt_elec / HBAR_EV_FS))
    return v @ (phase * (v.conj().T @ c))


def step_cayley(c, H_eff, dt_elec):
    """
    Crank-Nicolson / Cayley unitary step:
      c(t + dt) = (I + i H dt / 2 hbar)^(-1) @ (I - i H dt / 2 hbar) @ c(t)
    """
    n = len(c)
    I = np.eye(n, dtype=np.complex128)
    A = 0.5 * (dt_elec / HBAR_EV_FS) * H_eff
    lhs = I + 1j * A
    rhs = (I - 1j * A) @ c
    return np.linalg.solve(lhs, rhs)


def step_rk4(c, H_eff, dt_elec):
    """
    Standard 4th-order Runge-Kutta step for dc/dt = -i/hbar * H_eff * c.
    """
    fac = -1j / HBAR_EV_FS
    k1 = fac * (H_eff @ c)
    k2 = fac * (H_eff @ (c + 0.5 * dt_elec * k1))
    k3 = fac * (H_eff @ (c + 0.5 * dt_elec * k2))
    k4 = fac * (H_eff @ (c + dt_elec * k3))
    c_next = c + (dt_elec / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    # Renormalize to correct RK4 truncation drift
    norm = np.linalg.norm(c_next)
    if norm > 1e-12:
        c_next /= norm
    return c_next


def apply_edc_decoherence(c, active_surface, E_vec, dt_elec, c_param=0.1):
    """
    Energy-based Decoherence Correction (EDC / Granucci-Persico):
    Damp off-diagonal coefficients J != K:
      c_J -> c_J * exp(-dt_elec / tau_KJ)
      tau_KJ = hbar / |E_K - E_J| * (1 + C / E_kin)
    Then renormalize active state c_K.
    """
    if active_surface is None or active_surface < 0:
        return c

    K = active_surface
    n = len(c)
    dE = np.abs(E_vec - E_vec[K])

    for J in range(n):
        if J == K:
            continue
        if dE[J] > 1e-4:
            tau_KJ = HBAR_EV_FS / dE[J]
            damp = np.exp(-dt_elec / tau_KJ)
            c[J] *= damp

    # Renormalize active state K so total norm is 1
    sum_other_sq = np.sum(np.abs(c)**2) - np.abs(c[K])**2
    if sum_other_sq < 1.0:
        new_cK_mag = np.sqrt(1.0 - sum_other_sq)
        phase_k = np.angle(c[K]) if np.abs(c[K]) > 1e-12 else 0.0
        c[K] = new_cK_mag * np.exp(1j * phase_k)
    else:
        norm = np.linalg.norm(c)
        if norm > 1e-12:
            c /= norm

    return c


def propagate_electronic_substeps(
    c_init,
    E_k,
    E_kplus1,
    d_mat,
    dt_nuc_fs,
    n_substeps=500,
    active_surface=None,
    decoherence="edc",
    tau_dec_const_fs=15.0,
    integrator="unitary_matrix_exp"
):
    """
    Propagates electronic state coefficients over nuclear interval [t_k, t_k+1]
    using sub-stepping dt_elec = dt_nuc / n_substeps.

    Returns:
      c_final: coefficients at t_{k+1}
      rho_avg: average density matrix over the interval
      flux_KJ: cumulative flux for Tully hopping from active surface K
    """
    n_states = len(c_init)
    dt_elec = dt_nuc_fs / n_substeps
    c = np.array(c_init, dtype=np.complex128, copy=True)

    rho_accum = np.zeros((n_states, n_states), dtype=np.complex128)
    flux_KJ = np.zeros(n_states, dtype=np.float64)

    step_func = {
        "unitary_matrix_exp": step_unitary_matrix_exp,
        "cayley": step_cayley,
        "rk4": step_rk4,
    }.get(integrator, step_unitary_matrix_exp)

    for step in range(n_substeps):
        # Linear interpolation of energy
        tau_mid = (step + 0.5) * dt_elec
        s = tau_mid / dt_nuc_fs
        E_mid = E_k + s * (E_kplus1 - E_k)

        H_eff = build_effective_hamiltonian(E_mid, d_mat)

        # Electronic step
        c = step_func(c, H_eff, dt_elec)

        # Apply decoherence
        if decoherence == "edc" and active_surface is not None:
            c = apply_edc_decoherence(c, active_surface, E_mid, dt_elec)
        elif decoherence == "constant" and active_surface is not None:
            K = active_surface
            for J in range(n_states):
                if J != K:
                    c[J] *= np.exp(-dt_elec / tau_dec_const_fs)
            sum_other = np.sum(np.abs(c)**2) - np.abs(c[K])**2
            if sum_other < 1.0:
                c[K] = np.sqrt(1.0 - sum_other) * np.exp(1j * np.angle(c[K]))

        # Density matrix rho_IJ = c_I * conj(c_J)
        rho = np.outer(c, c.conj())
        rho_accum += rho

        # If running surface hopping, accumulate Tully flux from active surface K
        if active_surface is not None and d_mat is not None:
            K = active_surface
            rho_KK = max(rho[K, K].real, 1e-12)
            # 2 * dt_elec / rho_KK * Im(rho_KJ* * d_KJ)
            for J in range(n_states):
                if J != K:
                    val = 2.0 * dt_elec * (rho[K, J].conj() * d_mat[K, J]).imag / rho_KK
                    if val > 0.0:
                        flux_KJ[J] += val

    rho_avg = rho_accum / n_substeps
    return c, rho_avg, flux_KJ


def propagate_channel_rk4(c, E_k, E_kplus1, d_mat, dt_nuc_fs, n_substeps=50):
    """
    Propagates 1D state vector c (electron or hole channel) over [t_k, t_k+1] via RK4:
      dc/dt = -i/hbar * (E_mid * c) - d_mat @ c
    """
    dt_elec = dt_nuc_fs / n_substeps
    fac = -1j / HBAR_EV_FS
    c = np.array(c, dtype=np.complex128, copy=True)

    for step in range(n_substeps):
        s = (step + 0.5) / n_substeps
        E_mid = E_k + s * (E_kplus1 - E_k)

        def H_apply(c_in):
            return fac * (E_mid * c_in) - (d_mat @ c_in)

        k1 = H_apply(c)
        k2 = H_apply(c + 0.5 * dt_elec * k1)
        k3 = H_apply(c + 0.5 * dt_elec * k2)
        k4 = H_apply(c + dt_elec * k3)
        c = c + (dt_elec / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    norm = np.linalg.norm(c)
    if norm > 1e-12:
        c /= norm
    return c


def propagate_channel_batch_rk4(C, E_k_batch, E_kplus1_batch, d_mat, dt_nuc_fs, n_substeps=50):
    """
    Batched RK4 propagator for an ensemble of trajectories:
      dC/dt = -i/hbar * (E_mid * C) - d_mat @ C
    where:
      C: state amplitudes of shape (n_states, n_trajectories)
      E_k_batch, E_kplus1_batch: energies of shape (n_states, n_trajectories)
      d_mat: non-adiabatic coupling matrix of shape (n_states, n_states)
    """
    dt_elec = dt_nuc_fs / n_substeps
    fac = -1j / HBAR_EV_FS
    C = np.array(C, dtype=np.complex128, copy=True)
    dE = E_kplus1_batch - E_k_batch

    for step in range(n_substeps):
        s = (step + 0.5) / n_substeps
        E_mid = E_k_batch + s * dE

        k1 = fac * (E_mid * C) - (d_mat @ C)
        C2 = C + 0.5 * dt_elec * k1
        k2 = fac * (E_mid * C2) - (d_mat @ C2)
        C3 = C + 0.5 * dt_elec * k2
        k3 = fac * (E_mid * C3) - (d_mat @ C3)
        C4 = C + dt_elec * k3
        k4 = fac * (E_mid * C4) - (d_mat @ C4)
        C += (dt_elec / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    norms = np.linalg.norm(C, axis=0, keepdims=True)
    norms = np.maximum(norms, 1e-12)
    C /= norms
    return C


def propagate_channel_batch_strang(C, E_k_batch, E_kplus1_batch, d_mat, dt_nuc_fs, n_substeps=2, device="auto"):
    """
    Second-order Strang splitting (Trotter) unitary propagator for an ensemble of trajectories:
      U(tau) = exp(-i E tau / 2hbar) * exp(-d tau) * exp(-i E tau / 2hbar)
    where:
      C: state amplitudes of shape (n_states, n_trajectories)
      E_k_batch, E_kplus1_batch: energies of shape (n_states, n_trajectories)
      d_mat: skew-symmetric non-adiabatic coupling matrix of shape (n_states, n_states)
    """
    from scipy.linalg import expm
    dt_sub = dt_nuc_fs / n_substeps
    n_states, n_tr = C.shape
    dE = E_kplus1_batch - E_k_batch

    # expm(-d * dt_sub) is a real orthogonal matrix (O^T O = I)
    O_mat = expm(-d_mat * dt_sub)

    use_gpu = False
    if device in ("gpu", "mps", "cuda") or (device == "auto"):
        try:
            import torch
            dev = None
            if torch.backends.mps.is_available():
                dev = torch.device("mps")
            elif torch.cuda.is_available():
                dev = torch.device("cuda")

            if dev is not None:
                use_gpu = True
                O_t = torch.as_tensor(O_mat, dtype=torch.complex64, device=dev)
                C_t = torch.as_tensor(C, dtype=torch.complex64, device=dev)
                E_k_t = torch.as_tensor(E_k_batch, dtype=torch.float32, device=dev)
                dE_t = torch.as_tensor(dE, dtype=torch.float32, device=dev)

                for step in range(n_substeps):
                    s = (step + 0.5) / n_substeps
                    E_mid_t = E_k_t + s * dE_t
                    D_half = torch.exp(-1j * E_mid_t * (dt_sub / (2.0 * HBAR_EV_FS)))
                    C_t = D_half * (O_t @ (D_half * C_t))

                if dev.type == "mps":
                    torch.mps.synchronize()
                elif dev.type == "cuda":
                    torch.cuda.synchronize()
                return C_t.cpu().numpy().astype(np.complex128)
        except Exception:
            use_gpu = False

    # CPU path with NumPy
    C = np.array(C, dtype=np.complex128, copy=True)
    for step in range(n_substeps):
        s = (step + 0.5) / n_substeps
        E_mid = E_k_batch + s * dE
        D_half = np.exp(-1j * E_mid * (dt_sub / (2.0 * HBAR_EV_FS)))
        C = D_half * (O_mat @ (D_half * C))

    return C
