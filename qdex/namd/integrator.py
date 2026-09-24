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


def apply_edc_decoherence(c, active_surface, E_vec, dt_elec, c_param=0.1, n_atoms=1, temp_k=300.0, decay_type="exponential"):
    """
    Energy-based Decoherence Correction (EDC / Granucci-Persico):
    Damp off-diagonal coefficients J != K:
      c_J -> c_J * exp(-dt_elec / tau_KJ)       [exponential / EDC]
      c_J -> c_J * exp(-0.5 * (dt_elec / tau_KJ)^2)  [gaussian / GDC]
      tau_KJ = hbar / |E_K - E_J| * (1 + C / E_kin)
    C is c_param in Hartree. E_kin is the classical nuclear kinetic energy.
    Then renormalize active state c_K.
    """
    if active_surface is None or active_surface < 0:
        return c

    K = active_surface
    n = len(c)
    dE = np.abs(E_vec - E_vec[K])
    # 0.1 Ha -> eV. E_kin = (3/2) N_atoms kT, same convention as the master equation.
    HA_TO_EV = 27.211386245988
    E_kin = 1.5 * max(int(n_atoms), 1) * KB_EV * max(float(temp_k), 1e-3)
    edc_factor = 1.0 + (float(c_param) * HA_TO_EV) / max(E_kin, 1e-8)

    is_gaussian = str(decay_type).lower() in ("gaussian", "gdc")

    for J in range(n):
        if J == K:
            continue
        if dE[J] > 1e-4:
            tau_KJ = (HBAR_EV_FS / dE[J]) * edc_factor
            if is_gaussian:
                damp = np.exp(-0.5 * (dt_elec / max(tau_KJ, 1e-6)) ** 2)
            else:
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


def apply_gdc_decoherence(c, active_surface, E_vec, dt_elec, c_param=0.1, n_atoms=1, temp_k=300.0):
    """
    Gaussian Decoherence Correction (GDC / Granucci-Persico-Zoccante 2010):
    Damps off-diagonal coefficients via Gaussian factor exp(-0.5 * (dt / tau)^2).
    """
    return apply_edc_decoherence(
        c, active_surface, E_vec, dt_elec,
        c_param=c_param, n_atoms=n_atoms, temp_k=temp_k,
        decay_type="gaussian"
    )


def apply_edc_decoherence_batch(
    C,
    active_surfaces,
    E_batch,
    dt_fs,
    tau_mat=None,
    c_param=0.1,
    n_atoms=775,
    temp_k=300.0,
    min_tau_fs=1.0,
    max_tau_fs=500.0,
    decay_type="exponential"
):
    """
    Vectorized Energy-based Decoherence Correction (EDC / Granucci-Persico)
    or Gaussian Decoherence Correction (GDC / Granucci-Persico-Zoccante 2010)
    for an ensemble of trajectories C of shape (n_states, n_trajectories).

    For each trajectory tr with active state K = active_surfaces[tr]:
      c_J -> c_J * exp(-dt_fs / tau_KJ)             [decay_type="exponential"]
      c_J -> c_J * exp(-0.5 * (dt_fs / tau_KJ)^2)   [decay_type="gaussian"]
      for all J != K, then renormalize c_K so that sum_J |c_J|^2 = 1.

    If tau_mat of shape (n_states, n_states) is given, tau_KJ is taken from tau_mat[K, J].
    Otherwise, tau_KJ is computed via instantaneous EDC:
      tau_KJ = (hbar / |E_K - E_J|) * (1 + C / E_kin).
    """
    n_states, n_traj = C.shape
    active_surfaces = np.asarray(active_surfaces, dtype=int)
    tr_idx = np.arange(n_traj)
    is_gaussian = str(decay_type).lower() in ("gaussian", "gdc")

    if tau_mat is not None:
        tau_KJ = np.maximum(tau_mat[active_surfaces, :].T, min_tau_fs)  # shape: (n_states, n_traj)
        if is_gaussian:
            damp = np.exp(-0.5 * (dt_fs / tau_KJ) ** 2)
        else:
            damp = np.exp(-dt_fs / tau_KJ)
    else:
        E_act = E_batch[active_surfaces, tr_idx]  # shape: (n_traj,)
        dE = np.abs(E_batch - E_act[np.newaxis, :])  # shape: (n_states, n_traj)
        HA_TO_EV = 27.211386245988
        E_kin = 1.5 * max(int(n_atoms), 1) * KB_EV * max(float(temp_k), 1e-3)
        edc_factor = 1.0 + (float(c_param) * HA_TO_EV) / max(E_kin, 1e-8)
        tau_KJ = (HBAR_EV_FS / np.maximum(dE, 1e-4)) * edc_factor
        tau_KJ = np.clip(tau_KJ, min_tau_fs, max_tau_fs)
        if is_gaussian:
            damp = np.exp(-0.5 * (dt_fs / tau_KJ) ** 2)
        else:
            damp = np.exp(-dt_fs / tau_KJ)

    # Apply damping to non-active states
    C_damped = C * damp
    # Preserve active amplitude
    C_damped[active_surfaces, tr_idx] = C[active_surfaces, tr_idx]

    # Renormalize active state so sum_J |c_J|^2 = 1
    other_pop = np.sum(np.abs(C_damped)**2, axis=0) - np.abs(C[active_surfaces, tr_idx])**2
    other_pop = np.maximum(other_pop, 0.0)
    valid = other_pop < 1.0
    new_mag = np.sqrt(np.maximum(1.0 - other_pop, 0.0))
    orig_phase = np.angle(C[active_surfaces, tr_idx])
    C_damped[active_surfaces[valid], np.where(valid)[0]] = new_mag[valid] * np.exp(1j * orig_phase[valid])

    invalid = ~valid
    if np.any(invalid):
        norms = np.linalg.norm(C_damped[:, invalid], axis=0)
        C_damped[:, invalid] /= np.maximum(norms, 1e-12)[np.newaxis, :]

    return C_damped


def apply_gdc_decoherence_batch(
    C,
    active_surfaces,
    E_batch,
    dt_fs,
    tau_mat=None,
    c_param=0.1,
    n_atoms=775,
    temp_k=300.0,
    min_tau_fs=1.0,
    max_tau_fs=500.0
):
    """
    Vectorized Gaussian Decoherence Correction (GDC / Granucci-Persico-Zoccante 2010; Prezhdo & Rossky 1997)
    for an ensemble of trajectories C of shape (n_states, n_trajectories).
    Damping is Gaussian: c_J -> c_J * exp(-0.5 * (dt_fs / tau_KJ)^2).
    """
    return apply_edc_decoherence_batch(
        C,
        active_surfaces,
        E_batch,
        dt_fs,
        tau_mat=tau_mat,
        c_param=c_param,
        n_atoms=n_atoms,
        temp_k=temp_k,
        min_tau_fs=min_tau_fs,
        max_tau_fs=max_tau_fs,
        decay_type="gaussian"
    )



def step_dish_batch(
    C,
    active_surfaces,
    E_batch,
    dt_fs,
    tau_mat=None,
    beta=None,
    detailed_balance=True,
    min_tau_fs=1.0,
):
    """
    Decoherence-Induced Surface Hopping (DISH, Jaeger, Fischer, Prezhdo, JCP 2012)
    step for an ensemble of trajectories C of shape (n_states, n_trajectories).

    For each trajectory tr in active state K = active_surfaces[tr]:
      1. For each non-active state J != K:
         A dephasing event occurs with probability P_dec = 1 - exp(-dt_fs / tau_KJ).
      2. If a dephasing event occurs for J:
         A hop K -> J is attempted with probability:
           P_hop = |c_J|^2 * min(1, exp(-beta * max(E_J - E_K, 0)))
         - If accepted: active surface switches to J, and wavepacket collapses to J:
           c_J = 1.0, c_{L != J} = 0.
         - If rejected: state J is quenched (c_J -> 0).
      3. If no hop occurred, renormalize the remaining non-zero amplitudes of trajectory tr.

    Returns:
      C_new : updated complex ndarray of shape (n_states, n_trajectories)
      active_surfaces_new : updated int ndarray of shape (n_trajectories,)
      hops_occurred : bool ndarray of shape (n_trajectories,) indicating which hopped
    """
    n_states, n_traj = C.shape
    active_surfaces = np.asarray(active_surfaces, dtype=int).copy()
    C_out = np.array(C, dtype=np.complex128, copy=True)
    hops_occurred = np.zeros(n_traj, dtype=bool)

    if tau_mat is not None:
        tau_KJ = tau_mat[active_surfaces, :].T  # shape: (n_states, n_traj)
    else:
        tau_KJ = np.full((n_states, n_traj), 20.0)

    P_dec = 1.0 - np.exp(-dt_fs / np.maximum(tau_KJ, min_tau_fs))
    P_dec[active_surfaces, np.arange(n_traj)] = 0.0

    R1 = np.random.rand(n_states, n_traj)
    dec_events = R1 < P_dec

    E_act = E_batch[active_surfaces, np.arange(n_traj)]
    dE = E_batch - E_act[np.newaxis, :]
    if detailed_balance and beta is not None and beta > 0:
        boltz = np.exp(-np.maximum(dE, 0.0) * beta)
    else:
        boltz = np.ones((n_states, n_traj), dtype=np.float64)

    pop = np.abs(C_out)**2
    P_hop = pop * boltz

    R2 = np.random.rand(n_states, n_traj)
    hop_cands = dec_events & (R2 < P_hop)

    for tr in range(n_traj):
        K = active_surfaces[tr]
        cands = np.where(hop_cands[:, tr])[0]
        if len(cands) > 0:
            if len(cands) == 1:
                new_K = cands[0]
            else:
                p_cands = P_hop[cands, tr]
                sum_p = np.sum(p_cands)
                if sum_p > 0:
                    new_K = np.random.choice(cands, p=p_cands / sum_p)
                else:
                    new_K = np.random.choice(cands)
            active_surfaces[tr] = new_K
            C_out[:, tr] = 0.0
            C_out[new_K, tr] = 1.0
            hops_occurred[tr] = True
        else:
            decs = np.where(dec_events[:, tr])[0]
            if len(decs) > 0:
                C_out[decs, tr] = 0.0
                norm = np.linalg.norm(C_out[:, tr])
                if norm > 1e-12:
                    C_out[:, tr] /= norm
                else:
                    C_out[K, tr] = 1.0

    return C_out, active_surfaces, hops_occurred


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
            # g_{K->J} = max(0, 2 dt Re(c_K* c_J d_KJ) / |c_K|^2), d_KJ = <K|d/dt|J>
            for J in range(n_states):
                if J != K:
                    val = 2.0 * dt_elec * np.real(np.conj(c[K]) * c[J] * d_mat[K, J]) / rho_KK
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


def accumulate_fssh_flux(C, active_idx, d_mat, dt_sub, flux):
    """
    Add one substep of fewest-switches probability onto flux.

    g_{a->j} = max(0, 2 dt Re(c_a* c_j d_aj) / |c_a|^2)
    with d_aj = <a | d/dt | j>. C has shape (n_states, n_traj).
    """
    n_tr = C.shape[1]
    tr = np.arange(n_tr)
    active_idx = np.asarray(active_idx, dtype=int)
    c_a = C[active_idx, tr]
    pop = np.maximum(np.abs(c_a) ** 2, 1e-12)
    d_rows = d_mat[active_idx, :]
    term = np.conj(c_a)[np.newaxis, :] * C * d_rows.T
    incr = 2.0 * dt_sub * np.real(term) / pop[np.newaxis, :]
    np.maximum(incr, 0.0, out=incr)
    incr[active_idx, tr] = 0.0
    flux += incr


def propagate_channel_batch_rk4(C, E_k_batch, E_kplus1_batch, d_mat, dt_nuc_fs, n_substeps=50, active_idx=None):
    """
    Batched RK4 propagator for an ensemble of trajectories:
      dC/dt = -i/hbar * (E_mid * C) - d_mat @ C
    where:
      C: state amplitudes of shape (n_states, n_trajectories)
      E_k_batch, E_kplus1_batch: energies of shape (n_states, n_trajectories)
      d_mat: non-adiabatic coupling matrix of shape (n_states, n_states)

    If active_idx is given, also return the substep-integrated hop flux
    of shape (n_states, n_trajectories).
    """
    dt_elec = dt_nuc_fs / n_substeps
    fac = -1j / HBAR_EV_FS
    C = np.array(C, dtype=np.complex128, copy=True)
    dE = E_kplus1_batch - E_k_batch
    flux = None
    if active_idx is not None:
        flux = np.zeros(C.shape, dtype=np.float64)

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
        if flux is not None:
            accumulate_fssh_flux(C, active_idx, d_mat, dt_elec, flux)

    norms = np.linalg.norm(C, axis=0, keepdims=True)
    norms = np.maximum(norms, 1e-12)
    C /= norms
    if flux is None:
        return C
    return C, flux


def propagate_channel_batch_strang(C, E_k_batch, E_kplus1_batch, d_mat, dt_nuc_fs, n_substeps=2, device="auto", active_idx=None):
    """
    Second-order Strang splitting (Trotter) unitary propagator for an ensemble of trajectories:
      U(tau) = exp(-i E tau / 2hbar) * exp(-d tau) * exp(-i E tau / 2hbar)
    where:
      C: state amplitudes of shape (n_states, n_trajectories)
      E_k_batch, E_kplus1_batch: energies of shape (n_states, n_trajectories)
      d_mat: skew-symmetric non-adiabatic coupling matrix of shape (n_states, n_states)

    If active_idx is given, the fewest-switches flux is integrated in float64
    and (C, flux) is returned. GPU propagation is skipped so the flux reads
    the same amplitudes that are propagated.
    """
    from scipy.linalg import expm
    dt_sub = dt_nuc_fs / n_substeps
    n_states, n_tr = C.shape
    dE = E_kplus1_batch - E_k_batch

    # expm(-d * dt_sub) is a real orthogonal matrix (O^T O = I)
    O_mat = expm(-d_mat * dt_sub)

    use_gpu = False
    if active_idx is not None:
        device = "cpu"
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
    flux = np.zeros(C.shape, dtype=np.float64) if active_idx is not None else None
    for step in range(n_substeps):
        s = (step + 0.5) / n_substeps
        E_mid = E_k_batch + s * dE
        D_half = np.exp(-1j * E_mid * (dt_sub / (2.0 * HBAR_EV_FS)))
        C = D_half * (O_mat @ (D_half * C))
        if flux is not None:
            accumulate_fssh_flux(C, active_idx, d_mat, dt_sub, flux)

    if flux is None:
        return C
    return C, flux
