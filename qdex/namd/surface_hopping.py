import os
import time
import numpy as np

from qdex.namd.integrator import (
    propagate_electronic_substeps,
    propagate_channel_rk4,
    propagate_channel_batch_rk4,
    propagate_channel_batch_strang,
    KB_EV
)
from qdex.namd.initial_condition import sample_initial_states
from qdex.namd.master_equation import (
    compute_rate_matrix,
    run_master_equation_step,
    propagate_pme_tensor
)
from qdex.namd.analysis import analyze_and_plot_namd_results


def build_exciton_coupling_matrix(i_pairs, a_pairs, S_occ, S_virt, dt_nuc_fs):
    """
    Constructs the exciton non-adiabatic coupling matrix d_IJ between pair states:
      I = (i_I -> a_I), J = (j_J -> b_J)
    
    Rules in diagonal BSE:
      - Electron transition (i_I == j_J, a_I != b_J):
          S_IJ = S_virt[a_I, b_J]
      - Hole transition (a_I == b_J, i_I != j_J):
          S_IJ = S_occ[i_I, j_J]
      - Same state (i_I == j_J, a_I == b_J):
          S_IJ = S_virt[a_I, a_I] * S_occ[i_I, i_I] ~ 1
      - Two-particle change:
          S_IJ = 0
    
    d_IJ = (S_IJ - S_JI*) / (2 * dt_nuc_fs)
    """
    n_pairs = len(i_pairs)
    d_mat = np.zeros((n_pairs, n_pairs), dtype=np.complex128)

    # Electron transitions: same hole (i_I == j_J)
    # Group pairs by hole index i
    occ_groups = {}
    for idx, i_val in enumerate(i_pairs):
        occ_groups.setdefault(i_val, []).append(idx)

    for i_val, pair_indices in occ_groups.items():
        if len(pair_indices) > 1:
            pair_arr = np.array(pair_indices)
            a_idx = a_pairs[pair_arr]
            sub_S_virt = S_virt[np.ix_(a_idx, a_idx)]
            sub_d = (sub_S_virt - sub_S_virt.conj().T) / (2.0 * dt_nuc_fs)
            np.fill_diagonal(sub_d, 0.0)
            d_mat[np.ix_(pair_arr, pair_arr)] += sub_d

    # Hole transitions: same electron (a_I == b_J)
    # Group pairs by electron index a
    virt_groups = {}
    for idx, a_val in enumerate(a_pairs):
        virt_groups.setdefault(a_val, []).append(idx)

    for a_val, pair_indices in virt_groups.items():
        if len(pair_indices) > 1:
            pair_arr = np.array(pair_indices)
            i_idx = i_pairs[pair_arr]
            sub_S_occ = S_occ[np.ix_(i_idx, i_idx)]
            sub_d = (sub_S_occ - sub_S_occ.conj().T) / (2.0 * dt_nuc_fs)
            np.fill_diagonal(sub_d, 0.0)
            d_mat[np.ix_(pair_arr, pair_arr)] += sub_d

    return d_mat


def run_namd_dynamics(config):
    """
    Main driver for the NAMD Dynamics stage.
    Loads precomputed data from namd_precomputed/, samples hot exciton initial states
    (e.g. at 2 * Eg), and propagates non-adiabatic carrier cooling via CPA-DC-FSSH
    or the Pauli Master Equation.
    """
    t0_start = time.time()
    namd_cfg = config.get("namd", {})
    storage_cfg = namd_cfg.get("storage", {})
    dyn_cfg = namd_cfg.get("dynamics", {})
    init_cfg = namd_cfg.get("initial_excitation", {})
    integ_cfg = namd_cfg.get("integration", {})
    out_cfg = namd_cfg.get("output", {})

    precompute_dir = storage_cfg.get("precompute_dir", "namd_precomputed")
    meta_path = os.path.join(precompute_dir, "namd_metadata.npz")
    frame0_path = os.path.join(precompute_dir, "frame_00000.npz")

    if not os.path.exists(meta_path) or not os.path.exists(frame0_path):
        raise FileNotFoundError(
            f"Precomputed NAMD data not found in '{precompute_dir}'. "
            "Please run 'minibse --namd-precompute' first!"
        )

    meta = np.load(meta_path)
    n_frames = int(meta["n_frames"])
    dt_nuc_fs = float(meta["dt_nuc_fs"])

    frame0 = np.load(frame0_path)
    E0_pairs = frame0["E_pairs"]
    f0_pairs = frame0["f_pairs"]
    dft_gap = float(meta["mean_dft_gap"]) if "mean_dft_gap" in meta else float(frame0["dft_gap"])
    scissor = float(frame0["scissor"])
    if "mean_lowest_exc" in meta:
        qp_gap = float(meta["mean_lowest_exc"])
    elif "mean_qp_gap" in meta:
        qp_gap = float(meta["mean_qp_gap"])
    else:
        qp_gap = dft_gap + scissor

    # Initial condition setup
    mode = init_cfg.get("mode", "ratio_eg")
    if mode == "ratio_eg":
        ratio = float(init_cfg.get("ratio", 2.0))
        pump_energy_ev = ratio * qp_gap
    else:
        pump_energy_ev = float(init_cfg.get("energy_ev", 2.0 * qp_gap))

    pulse_fwhm_ev = float(init_cfg.get("pulse_fwhm_ev", 0.08))
    filter_dark = bool(init_cfg.get("filter_dark_states", True))

    method = dyn_cfg.get("method", "cpa_fssh").lower()
    temp_k = float(dyn_cfg.get("temperature_k", 300.0))
    detailed_balance = bool(dyn_cfg.get("detailed_balance", True))
    decoherence = dyn_cfg.get("decoherence", "edc")
    tau_dec_raw = dyn_cfg.get("tau_dec_fs", "edc")
    if str(decoherence).lower() in ("cumulant", "lowest_state", "lowest_exciton", "gap_cumulant"):
        tau_dec_raw = "cumulant"

    if tau_dec_raw is None or str(tau_dec_raw).lower() in ("edc", "none", "auto", "dynamic"):
        tau_dec_fs = "edc"
    elif str(tau_dec_raw).lower() in ("cumulant", "lowest_state", "lowest_exciton", "gap_cumulant"):
        from qdex.namd.analysis import compute_band_gap_dynamics_and_spectral_density
        res_cumulant = compute_band_gap_dynamics_and_spectral_density(precompute_dir, use_lowest_exciton=True)
        if res_cumulant is not None and not np.isnan(res_cumulant.get("tau_dec_fs", np.nan)):
            tau_dec_fs = float(res_cumulant["tau_dec_fs"])
            std_g = np.sqrt(res_cumulant["var_g"])
            print(f"  [NAMD] Ab initio decoherence time from cumulant expansion of lowest excited state: tau_dec = {tau_dec_fs:.2f} fs (gap fluctuation std: {std_g:.4f} eV)")
            decoherence = f"cumulant ({tau_dec_fs:.1f} fs)"
        else:
            print("  [NAMD:Warn] Cumulant expansion from lowest excited state failed; falling back to tau_dec = 14.0 fs")
            tau_dec_fs = 14.0
            decoherence = f"cumulant (fallback 14.0 fs)"
    else:
        try:
            tau_dec_fs = float(tau_dec_raw)
        except ValueError:
            tau_dec_fs = "edc"
    n_trajectories = int(dyn_cfg.get("n_trajectories", 1000))

    integrator_type = integ_cfg.get("integrator", "strang").lower()
    device_cfg = integ_cfg.get("device", config.get("system", {}).get("device", "auto"))
    if integrator_type in ("strang", "trotter", "unitary_strang", "split_operator", "unitary_matrix_exp"):
        n_substeps = int(integ_cfg.get("n_substeps", 2))
    else:
        n_substeps = int(integ_cfg.get("n_substeps", 50))

    i_pairs0 = frame0["i_pairs"]
    a_pairs0 = frame0["a_pairs"]
    n_occ = int(meta["n_occ"]) if "n_occ" in meta else (int(np.max(i_pairs0)) + 1)
    n_virt = int(meta["n_virt"]) if "n_virt" in meta else (int(np.max(a_pairs0)) + 1)

    # Active Dynamics Window: restrict active MOs to the carrier cooling cascade
    dyn_window = dyn_cfg.get("active_energy_window_ev", None)
    if dyn_window is None:
        # Default cooling window: [qp_gap - 0.2, pump_energy_ev + 0.3]
        dyn_e_min = max(0.0, qp_gap - 0.2)
        dyn_e_max = pump_energy_ev + 0.3
    elif isinstance(dyn_window, (list, tuple)) and len(dyn_window) == 2:
        dyn_e_min = float(dyn_window[0])
        dyn_e_max = float(dyn_window[1])
    else:
        # "all" or False: use full precomputed window
        dyn_e_min = float(np.min(E0_pairs))
        dyn_e_max = float(np.max(E0_pairs))

    dyn_mask = (E0_pairs >= dyn_e_min) & (E0_pairs <= dyn_e_max)
    if np.sum(dyn_mask) == 0:
        dyn_mask = np.ones(len(E0_pairs), dtype=bool)

    dyn_occ_active = np.sort(np.unique(i_pairs0[dyn_mask]))
    dyn_virt_active = np.sort(np.unique(a_pairs0[dyn_mask]))
    n_occ_dyn = len(dyn_occ_active)
    n_virt_dyn = len(dyn_virt_active)

    occ_to_sub = np.full(n_occ, -1, dtype=int)
    occ_to_sub[dyn_occ_active] = np.arange(n_occ_dyn)
    virt_to_sub = np.full(n_virt, -1, dtype=int)
    virt_to_sub[dyn_virt_active] = np.arange(n_virt_dyn)

    pair_lookup = np.full((n_occ, n_virt), -1, dtype=int)
    pair_lookup[i_pairs0, a_pairs0] = np.arange(len(i_pairs0))

    # Recombination & Photoluminescence Rates (Radiative & Non-Radiative)
    from qdex.hardness import (
        get_refractive_index, compute_radiative_rates, compute_energy_gap_law_rate, MATERIAL_DB
    )
    from qdex.namd.integrator import KB_EV

    material = config.get("system", {}).get("material", "DEFAULT")
    n_refr = float(dyn_cfg.get("refractive_index", get_refractive_index(material)))
    k_rad_s, k_rad_fs = compute_radiative_rates(E0_pairs, f0_pairs, refractive_index=n_refr)

    # Non-radiative rate setup (trap-assisted or Energy Gap Law)
    if "tau_nr_ns" in dyn_cfg:
        tau_nr_ns = float(dyn_cfg["tau_nr_ns"])
        k_nr_s = 1e9 / max(tau_nr_ns, 1e-12)
        k_nr_fs = k_nr_s * 1e-15
    elif "k_nr_s" in dyn_cfg:
        k_nr_s = float(dyn_cfg["k_nr_s"])
        k_nr_fs = k_nr_s * 1e-15
        tau_nr_ns = 1e9 / max(k_nr_s, 1e-12)
    else:
        mat_entry = MATERIAL_DB.get(str(material).upper(), MATERIAL_DB["DEFAULT"])
        e_lo = mat_entry[6] if len(mat_entry) > 6 else 0.018
        k_nr_s, k_nr_fs = compute_energy_gap_law_rate(qp_gap, E_LO_ev=e_lo)
        tau_nr_ns = 1e9 / max(k_nr_s, 1e-30) if k_nr_s > 0 else np.inf

    k_loss_mat = np.zeros((n_occ, n_virt), dtype=np.float64)
    k_loss_mat[i_pairs0, a_pairs0] = k_rad_fs + k_nr_fs

    # Thermalized radiative rate and PLQY at band edge
    beta = 1.0 / (KB_EV * max(temp_k, 1e-3))
    delta_E = E0_pairs - np.min(E0_pairs)
    boltz_weights = np.exp(-delta_E * beta)
    sum_w = np.sum(boltz_weights)
    k_rad_therm_s = float(np.sum(k_rad_s * boltz_weights) / sum_w) if sum_w > 0 else 0.0
    tau_rad_therm_ns = (1e9 / k_rad_therm_s) if k_rad_therm_s > 0 else np.inf
    tau_rad_1_ns = (1e9 / k_rad_s[0]) if (len(k_rad_s) > 0 and k_rad_s[0] > 0) else np.inf
    plqy = (k_rad_therm_s / max(k_rad_therm_s + k_nr_s, 1e-30)) * 100.0

    recombination_info = {
        "material": str(material).upper(),
        "refractive_index": n_refr,
        "k_rad_s": k_rad_s,
        "k_rad_fs": k_rad_fs,
        "k_nr_s": k_nr_s,
        "k_nr_fs": k_nr_fs,
        "tau_rad_1_ns": tau_rad_1_ns,
        "tau_rad_therm_ns": tau_rad_therm_ns,
        "tau_nr_ns": tau_nr_ns,
        "k_rad_therm_s": k_rad_therm_s,
        "plqy_percent": plqy
    }

    print("=" * 65)
    print(f" QDEX - NAMD Carrier Cooling Simulation ({method.upper()})")
    print("=" * 65)
    print(f"  Precomputed Data     : {precompute_dir} ({n_frames} frames, dt = {dt_nuc_fs} fs)")
    print(f"  Band Gap (Eg)        : DFT = {dft_gap:.3f} eV, QP = {qp_gap:.3f} eV")
    print(f"  Photoexcitation Pump : {pump_energy_ev:.3f} eV ({pump_energy_ev/qp_gap:.2f} * Eg)")
    print(f"  Laser Pulse FWHM     : {pulse_fwhm_ev:.3f} eV (filter dark states: {filter_dark})")
    print(f"  Temperature          : {temp_k} K (detailed balance: {detailed_balance})")
    print(f"  Precomputed Pairs    : {len(E0_pairs)} (n_occ={n_occ}, n_virt={n_virt})")
    print(f"  Dynamics Window      : [{dyn_e_min:.2f}, {dyn_e_max:.2f}] eV ({np.sum(dyn_mask)} active pairs)")
    print(f"  Active Dynamics Space: {n_occ_dyn} occ (HOMO-{n_occ - 1 - dyn_occ_active[0]} .. HOMO-{n_occ - 1 - dyn_occ_active[-1]}), {n_virt_dyn} virt (LUMO+{dyn_virt_active[0]} .. LUMO+{dyn_virt_active[-1]})")
    if method == "cpa_fssh":
        print(f"  Trajectories         : {n_trajectories}")
        print(f"  Decoherence          : {decoherence} (tau_dec = {tau_dec_fs} fs)")
        print(f"  Electronic Sub-steps : {n_substeps} (dt_elec = {dt_nuc_fs/n_substeps:.5f} fs)")
        print(f"  Integrator           : {integrator_type} (device: {device_cfg})")
    print("=" * 65 + "\n")

    # Sample initial state distribution
    sampled_states, P_init = sample_initial_states(
        energies=E0_pairs,
        f_osc=f0_pairs,
        pump_energy_ev=pump_energy_ev,
        pulse_fwhm_ev=pulse_fwhm_ev,
        n_trajectories=n_trajectories,
        filter_dark_states=filter_dark
    )

    times = np.arange(n_frames) * dt_nuc_fs
    mean_energy = np.zeros(n_frames)
    mean_excess_e = np.zeros(n_frames)
    mean_excess_h = np.zeros(n_frames)
    populations = np.zeros((n_frames, len(E0_pairs)))
    trajectory_energies = np.zeros((n_frames, n_trajectories)) if method == "cpa_fssh" else None
    all_energies = np.zeros((n_frames, len(E0_pairs)))
    all_energies[0, :] = E0_pairs

    eps_occ_0 = frame0["eps_occ"] if "eps_occ" in frame0 else None
    eps_virt_0 = frame0["eps_virt"] if "eps_virt" in frame0 else None

    # -------------------------------------------------------------
    # SCHEME A: CPA-DC-FSSH (Fewest Switches Surface Hopping)
    # -------------------------------------------------------------
    if method == "cpa_fssh":
        active_pairs = [(int(i_pairs0[s]), int(a_pairs0[s])) for s in sampled_states]
        active_surfaces = np.array(sampled_states, dtype=int)

        mean_energy[0] = np.mean(E0_pairs[active_surfaces])
        trajectory_energies[0, :] = E0_pairs[active_surfaces]
        if eps_occ_0 is not None and eps_virt_0 is not None:
            mean_excess_e[0] = np.mean(eps_virt_0[a_pairs0[active_surfaces]] - eps_virt_0[0])
            mean_excess_h[0] = np.mean(eps_occ_0[-1] - eps_occ_0[i_pairs0[active_surfaces]])
        else:
            mean_excess_e[0] = 0.5 * (mean_energy[0] - qp_gap)
            mean_excess_h[0] = 0.5 * (mean_energy[0] - qp_gap)

        for s in active_surfaces:
            populations[0, s] += 1.0 / n_trajectories

        print(f"  [NAMD:FSSH] Starting batched propagation across {n_frames} frames ({n_trajectories} trajectories)...")
        beta = 1.0 / (KB_EV * max(temp_k, 1e-3))

        for k in range(n_frames - 1):
            t_step_start = time.time()
            step_file = os.path.join(precompute_dir, f"step_{k:05d}_to_{k+1:05d}.npz")
            if not os.path.exists(step_file):
                break

            step_data = np.load(step_file)
            E_k = step_data["E_prev"]
            E_kplus1 = step_data["E_curr"]
            all_energies[k + 1, :] = E_kplus1
            i_pairs = step_data["i_pairs"] if "i_pairs" in step_data else step_data["i_pairs_curr"]
            a_pairs = step_data["a_pairs"] if "a_pairs" in step_data else step_data["a_pairs_curr"]
            S_occ = step_data["S_occ"]
            S_virt = step_data["S_virt"]
            eps_occ_curr = step_data["eps_occ_curr"] if "eps_occ_curr" in step_data else eps_occ_0
            eps_virt_curr = step_data["eps_virt_curr"] if "eps_virt_curr" in step_data else eps_virt_0
            eps_occ_prev = step_data["eps_occ_prev"] if "eps_occ_prev" in step_data else eps_occ_0
            eps_virt_prev = step_data["eps_virt_prev"] if "eps_virt_prev" in step_data else eps_virt_0

            d_occ = (S_occ - S_occ.conj().T) / (2.0 * dt_nuc_fs)
            d_virt = (S_virt - S_virt.conj().T) / (2.0 * dt_nuc_fs)

            # Build exact full-space pair energy matrices from single-particle energies
            # This ensures no zeros or missing entries for states in the active dynamics window
            if eps_occ_prev is not None and eps_virt_prev is not None:
                E_mat_k = eps_virt_prev[np.newaxis, :] - eps_occ_prev[:, np.newaxis]
            else:
                E_mat_k = np.zeros((n_occ, n_virt), dtype=np.float64)
                E_mat_k[i_pairs, a_pairs] = E_k

            if eps_occ_curr is not None and eps_virt_curr is not None:
                E_mat_kplus1 = eps_virt_curr[np.newaxis, :] - eps_occ_curr[:, np.newaxis]
            else:
                E_mat_kplus1 = np.zeros((n_occ, n_virt), dtype=np.float64)
                E_mat_kplus1[i_pairs, a_pairs] = E_kplus1

            d_occ_dyn = d_occ[np.ix_(dyn_occ_active, dyn_occ_active)]
            d_virt_dyn = d_virt[np.ix_(dyn_virt_active, dyn_virt_active)]

            curr_i = np.array([p[0] for p in active_pairs], dtype=int)
            curr_a = np.array([p[1] for p in active_pairs], dtype=int)
            curr_i_sub = occ_to_sub[curr_i]
            curr_a_sub = virt_to_sub[curr_a]

            print(f"  [Step {k+1}/{n_frames-1}] Propagating {n_trajectories} trajectories (t = {k*dt_nuc_fs:.1f} -> {(k+1)*dt_nuc_fs:.1f} fs)...", flush=True)

            # 1. Batched Electron channel propagation
            t_e0 = time.time()
            C_e = np.zeros((n_virt_dyn, n_trajectories), dtype=np.complex128)
            C_e[curr_a_sub, np.arange(n_trajectories)] = 1.0
            E_e_k_batch = E_mat_k[curr_i, :][:, dyn_virt_active].T
            E_e_kplus1_batch = E_mat_kplus1[curr_i, :][:, dyn_virt_active].T

            if integrator_type in ("strang", "trotter", "unitary_strang", "split_operator", "unitary_matrix_exp"):
                C_e = propagate_channel_batch_strang(C_e, E_e_k_batch, E_e_kplus1_batch, d_virt_dyn, dt_nuc_fs, n_substeps, device=device_cfg)
            else:
                C_e = propagate_channel_batch_rk4(C_e, E_e_k_batch, E_e_kplus1_batch, d_virt_dyn, dt_nuc_fs, n_substeps)
            t_e_rk4 = time.time() - t_e0

            # 2. Batched Hole channel propagation
            t_h0 = time.time()
            C_h = np.zeros((n_occ_dyn, n_trajectories), dtype=np.complex128)
            C_h[curr_i_sub, np.arange(n_trajectories)] = 1.0
            E_h_k_batch = E_mat_k[:, curr_a][dyn_occ_active, :]
            E_h_kplus1_batch = E_mat_kplus1[:, curr_a][dyn_occ_active, :]

            if integrator_type in ("strang", "trotter", "unitary_strang", "split_operator", "unitary_matrix_exp"):
                C_h = propagate_channel_batch_strang(C_h, E_h_k_batch, E_h_kplus1_batch, d_occ_dyn, dt_nuc_fs, n_substeps, device=device_cfg)
            else:
                C_h = propagate_channel_batch_rk4(C_h, E_h_k_batch, E_h_kplus1_batch, d_occ_dyn, dt_nuc_fs, n_substeps)
            t_h_rk4 = time.time() - t_h0

            # 3. Surface Hopping evaluation
            t_hop0 = time.time()
            n_hops_e = 0
            n_hops_h = 0

            for tr in range(n_trajectories):
                i_c = curr_i[tr]
                a_c = curr_a[tr]
                a_sub = curr_a_sub[tr]
                i_sub = curr_i_sub[tr]
                E_curr = E_mat_kplus1[i_c, a_c]

                # Electron hop probabilities
                rho_curr_e = max(np.abs(C_e[a_sub, tr]) ** 2, 1e-12)
                flux_e = 2.0 * dt_nuc_fs * np.real(np.conj(C_e[a_sub, tr]) * C_e[:, tr] * d_virt_dyn[a_sub, :]) / rho_curr_e
                probs_b = np.maximum(flux_e, 0.0)
                probs_b[a_sub] = 0.0
                valid_b = pair_lookup[i_c, dyn_virt_active] >= 0
                probs_b[~valid_b] = 0.0

                if eps_virt_curr is not None:
                    dE_b = eps_virt_curr[dyn_virt_active] - eps_virt_curr[a_c]
                else:
                    dE_b = E_mat_kplus1[i_c, dyn_virt_active] - E_curr
                if detailed_balance:
                    probs_b *= np.exp(-np.maximum(dE_b, 0.0) * beta)

                # Hole hop probabilities
                rho_curr_h = max(np.abs(C_h[i_sub, tr]) ** 2, 1e-12)
                flux_h = 2.0 * dt_nuc_fs * np.real(np.conj(C_h[i_sub, tr]) * C_h[:, tr] * d_occ_dyn[i_sub, :]) / rho_curr_h
                probs_j = np.maximum(flux_h, 0.0)
                probs_j[i_sub] = 0.0
                valid_j = pair_lookup[dyn_occ_active, a_c] >= 0
                probs_j[~valid_j] = 0.0

                if eps_occ_curr is not None:
                    dE_j = eps_occ_curr[i_c] - eps_occ_curr[dyn_occ_active]
                else:
                    dE_j = E_mat_kplus1[dyn_occ_active, a_c] - E_curr
                if detailed_balance:
                    probs_j *= np.exp(-np.maximum(dE_j, 0.0) * beta)

                # Combine hops
                total_b = np.sum(probs_b)
                total_j = np.sum(probs_j)
                total_hop = total_b + total_j

                if total_hop > 0.0:
                    zeta = np.random.rand()
                    if zeta < total_b:
                        p_norm = probs_b / total_b
                        new_a_sub = np.random.choice(n_virt_dyn, p=p_norm)
                        active_pairs[tr] = (i_c, dyn_virt_active[new_a_sub])
                        n_hops_e += 1
                    elif zeta < total_hop:
                        p_norm = probs_j / total_j
                        new_i_sub = np.random.choice(n_occ_dyn, p=p_norm)
                        active_pairs[tr] = (dyn_occ_active[new_i_sub], a_c)
                        n_hops_h += 1

                new_idx = pair_lookup[active_pairs[tr][0], active_pairs[tr][1]]
                active_surfaces[tr] = new_idx if new_idx >= 0 else active_surfaces[tr]

            t_hop = time.time() - t_hop0

            # Record observables at step k + 1
            mean_energy[k + 1] = np.mean(E_kplus1[active_surfaces])
            trajectory_energies[k + 1, :] = E_kplus1[active_surfaces]
            for s in active_surfaces:
                populations[k + 1, s] += 1.0 / n_trajectories

            curr_a_new = np.array([p[1] for p in active_pairs])
            curr_i_new = np.array([p[0] for p in active_pairs])
            if eps_occ_curr is not None and eps_virt_curr is not None:
                mean_excess_e[k + 1] = np.mean(eps_virt_curr[curr_a_new] - eps_virt_curr[0])
                mean_excess_h[k + 1] = np.mean(eps_occ_curr[-1] - eps_occ_curr[curr_i_new])
            else:
                mean_excess_e[k + 1] = 0.5 * (mean_energy[k + 1] - qp_gap)
                mean_excess_h[k + 1] = 0.5 * (mean_energy[k + 1] - qp_gap)

            t_step = time.time() - t_step_start
            min_e_act = np.min(E_kplus1[active_surfaces])
            max_e_act = np.max(E_kplus1[active_surfaces])
            print(f"    -> e- RK4: {t_e_rk4:.2f} s | h+ RK4: {t_h_rk4:.2f} s | hops: {n_hops_e} e-, {n_hops_h} h+ ({t_hop:.2f} s)")
            print(f"    -> Step {k+1}/{n_frames-1} completed in {t_step:.2f} s | <E_exc> = {mean_energy[k+1]:.4f} eV | Active range: [{min_e_act:.3f}, {max_e_act:.3f}] eV\n", flush=True)

    # -------------------------------------------------------------
    # SCHEME B: Pauli Master Equation (Deterministic Kinetics)
    # -------------------------------------------------------------
    elif method == "master_equation":
        i_pairs0 = frame0["i_pairs"]
        a_pairs0 = frame0["a_pairs"]
        n_occ = int(np.max(i_pairs0)) + 1
        n_virt = int(np.max(a_pairs0)) + 1

        P_mat = np.zeros((n_occ, n_virt), dtype=np.float64)
        P_mat[i_pairs0, a_pairs0] = P_init
        mean_energy[0] = np.sum(P_init * E0_pairs)
        populations[0] = P_init

        p_virt0 = np.sum(P_mat, axis=0)
        p_occ0 = np.sum(P_mat, axis=1)
        if eps_occ_0 is not None and eps_virt_0 is not None:
            mean_excess_e[0] = np.sum(p_virt0 * (eps_virt_0 - eps_virt_0[0]))
            mean_excess_h[0] = np.sum(p_occ0 * (eps_occ_0[-1] - eps_occ_0))
        else:
            mean_excess_e[0] = 0.5 * (mean_energy[0] - qp_gap)
            mean_excess_h[0] = 0.5 * (mean_energy[0] - qp_gap)

        print(f"  [NAMD:PME] Starting deterministic Master Equation propagation (n_occ={n_occ}, n_virt={n_virt})...")

        for k in range(n_frames - 1):
            t_step_start = time.time()
            step_file = os.path.join(precompute_dir, f"step_{k:05d}_to_{k+1:05d}.npz")
            if not os.path.exists(step_file):
                break

            step_data = np.load(step_file)
            E_k = step_data["E_curr"]
            all_energies[k + 1, :] = E_k
            i_pairs = step_data["i_pairs"] if "i_pairs" in step_data else (step_data["i_pairs_curr"] if "i_pairs_curr" in step_data else i_pairs0)
            a_pairs = step_data["a_pairs"] if "a_pairs" in step_data else (step_data["a_pairs_curr"] if "a_pairs_curr" in step_data else a_pairs0)
            S_occ = step_data["S_occ"]
            S_virt = step_data["S_virt"]
            eps_occ_curr = step_data["eps_occ_curr"] if "eps_occ_curr" in step_data else eps_occ_0
            eps_virt_curr = step_data["eps_virt_curr"] if "eps_virt_curr" in step_data else eps_virt_0

            d_occ = (S_occ - S_occ.conj().T) / (2.0 * dt_nuc_fs)
            d_virt = (S_virt - S_virt.conj().T) / (2.0 * dt_nuc_fs)

            E_mat_k = np.zeros((n_occ, n_virt), dtype=np.float64)
            E_mat_k[i_pairs, a_pairs] = E_k

            P_mat = propagate_pme_tensor(
                P_mat=P_mat,
                E_mat=E_mat_k,
                d_occ=d_occ,
                d_virt=d_virt,
                dt_fs=dt_nuc_fs,
                temp_k=temp_k,
                tau_dec_fs=tau_dec_fs,
                eps_occ=eps_occ_curr,
                eps_virt=eps_virt_curr,
                k_loss=k_loss_mat
            )

            P_vec = P_mat[i_pairs, a_pairs]
            mean_energy[k + 1] = np.sum(P_vec * E_k)
            populations[k + 1] = P_vec

            p_virt = np.sum(P_mat, axis=0)
            p_occ = np.sum(P_mat, axis=1)
            if eps_occ_curr is not None and eps_virt_curr is not None:
                mean_excess_e[k + 1] = np.sum(p_virt * (eps_virt_curr - eps_virt_curr[0]))
                mean_excess_h[k + 1] = np.sum(p_occ * (eps_occ_curr[-1] - eps_occ_curr))
            else:
                mean_excess_e[k + 1] = 0.5 * (mean_energy[k + 1] - qp_gap)
                mean_excess_h[k + 1] = 0.5 * (mean_energy[k + 1] - qp_gap)

            t_step = time.time() - t_step_start
            print(f"    Step {k+1}/{n_frames-1} (t = {(k+1)*dt_nuc_fs:.1f} fs) in {t_step:.3f} s | <E_exc> = {mean_energy[k+1]:.4f} eV")

    total_sim_time = time.time() - t0_start
    print(f"\n[NAMD Dynamics] Completed in {total_sim_time:.2f} s")

    # Analysis, CSV writing, and plotting
    analyze_and_plot_namd_results(
        times_fs=times,
        mean_energies_ev=mean_energy,
        populations=populations,
        out_cfg=out_cfg,
        qp_gap_ev=qp_gap,
        pump_energy_ev=pump_energy_ev,
        mean_excess_e=mean_excess_e,
        mean_excess_h=mean_excess_h,
        trajectory_energies=trajectory_energies,
        all_energies=all_energies,
        precompute_dir=precompute_dir,
        recombination_info=recombination_info
    )
