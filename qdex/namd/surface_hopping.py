import os
import time
import numpy as np

from qdex.constants import HA_TO_EV
from qdex.namd.integrator import (
    propagate_electronic_substeps,
    propagate_channel_rk4,
    propagate_channel_batch_rk4,
    propagate_channel_batch_strang,
    apply_edc_decoherence_batch,
    step_dish_batch,
    KB_EV,
    HBAR_EV_FS
)
from qdex.namd.initial_condition import sample_initial_states
from qdex.namd.master_equation import (
    compute_rate_matrix,
    run_master_equation_step,
    propagate_pme_tensor
)
from qdex.namd.analysis import analyze_and_plot_namd_results
from qdex.namd.ensemble import (
    auto_calibrate_ensemble_origins,
    load_origin_frame_data,
    sample_origin_initial_states,
    aggregate_multi_origin_results,
)



def pair_energy_matrix(eps_occ, eps_virt, i_pairs, a_pairs, E_pairs):
    """Quasiparticle pair energies with diagonal-BSE binding written on stored pairs."""
    E = eps_virt[np.newaxis, :] - eps_occ[:, np.newaxis]
    if E_pairs is not None and i_pairs is not None and len(E_pairs) == len(i_pairs):
        E = np.array(E, copy=True)
        E[np.asarray(i_pairs, dtype=int), np.asarray(a_pairs, dtype=int)] = np.asarray(E_pairs, dtype=np.float64)
    return E


def rescale_legacy_oscillator_strengths(f_pairs, soc=False, f_convention=None):
    """
    Older precomputes stored f = (2/3) * E_eV * mu^2.
    Closed-shell singlets need (4/3) * E_Ha * mu^2, a factor 2/HA_TO_EV.
    Spinors need (2/3) * E_Ha * mu^2, a factor 1/HA_TO_EV.
    """
    if f_convention in ("au_hartree", "hartree"):
        return np.asarray(f_pairs, dtype=np.float64)
    factor = (1.0 / HA_TO_EV) if soc else (2.0 / HA_TO_EV)
    return np.asarray(f_pairs, dtype=np.float64) * factor


def resolve_optin_auger_rate_fs(config, dyn_cfg):
    """
    Biexciton / ECSH rate. Never called from a plain cooling run.
    A user lifetime wins. Otherwise the frame-0 golden-rule k_XX is required.
    """
    if "k_auger_fs" in dyn_cfg:
        return float(dyn_cfg["k_auger_fs"]), "k_auger_fs from configuration"
    if "tau_auger_ps" in dyn_cfg:
        tau = max(float(dyn_cfg["tau_auger_ps"]), 1e-12)
        return (1.0 / tau) * 1e-3, f"tau_auger_ps = {tau} ps"
    if "tau_auger_ns" in dyn_cfg:
        tau_ps = max(float(dyn_cfg["tau_auger_ns"]), 1e-12) * 1.0e3
        return (1.0 / tau_ps) * 1e-3, f"tau_auger_ns = {dyn_cfg['tau_auger_ns']} ns"
    from qdex.auger import auger_rates_from_config
    res = auger_rates_from_config(config)
    if not np.isfinite(res.rate_biexciton_fs) or res.rate_biexciton_fs <= 0.0:
        raise RuntimeError(
            "ECSH or a biexciton clock was requested, but the frame-0 Auger "
            "golden-rule rate is zero. The NAMD orbital window does not contain "
            "the Auger continuum. Pass the full MO file or set "
            "namd.dynamics.tau_auger_ps. A cooling run does not need either."
        )
    return float(res.rate_biexciton_fs), (
        f"frame-0 golden rule k_XX = {res.rate_biexciton_ps:.4e} ps^-1 "
        f"(k_X- = {res.rate_eeh_ps:.4e}, k_X+ = {res.rate_hhe_ps:.4e} ps^-1)"
    )


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


def propagate_single_namd_origin(
    k0,
    n_win_steps,
    n_trajectories,
    precompute_dir,
    frame0,
    E0_pairs,
    f0_pairs,
    i_pairs0,
    a_pairs0,
    qp_gap_0,
    pump_energy_ev,
    pulse_fwhm_ev,
    filter_dark,
    dyn_mask,
    dyn_occ_active,
    dyn_virt_active,
    occ_to_sub,
    virt_to_sub,
    n_occ,
    n_virt,
    n_occ_dyn,
    n_virt_dyn,
    n_steps_base,
    dt_nuc_fs,
    method,
    decoherence,
    decoherence_decay_type,
    tau_dec_fs,
    integrator_type,
    n_substeps,
    device_cfg,
    temp_k,
    beta,
    detailed_balance,
    pair_lookup,
    tau_occ_dyn,
    tau_virt_dyn,
    is_surface_hopping,
    is_biexciton,
    ecsh_auger,
    k_auger_fs,
    ecsh_window_ev,
    k_loss_mat,
    run_tr_sd,
    hop_records=None,
    pme_flux_records=None,
    verbose=True,
    origin_idx=0,
    n_origins=1
):
    """
    Propagates NAMD carrier relaxation across a single trajectory window of length n_win_steps
    starting from frame k0.
    """
    if k0 == 0:
        E_init = E0_pairs
        f_init = f0_pairs
        eps_occ_init = frame0["eps_occ"] if "eps_occ" in frame0 else None
        eps_virt_init = frame0["eps_virt"] if "eps_virt" in frame0 else None
        qp_gap_init = qp_gap_0
    else:
        orig_data = load_origin_frame_data(precompute_dir, k0, frame0=frame0)
        E_init = orig_data["E_pairs"]
        f_init = orig_data["f_pairs"]
        eps_occ_init = orig_data["eps_occ"]
        eps_virt_init = orig_data["eps_virt"]
        qp_gap_init = orig_data["qp_gap"]

    # Preserve excess kinetic energy relative to instantaneous band edge
    pump_excess_ev = max(0.0, pump_energy_ev - qp_gap_0)
    pump_ev_k = qp_gap_init + pump_excess_ev

    sampled_states, P_init = sample_initial_states(
        energies=E_init,
        f_osc=f_init,
        pump_energy_ev=pump_ev_k,
        pulse_fwhm_ev=pulse_fwhm_ev,
        n_trajectories=n_trajectories,
        filter_dark_states=filter_dark,
        mask=dyn_mask
    )

    n_frames_win = n_win_steps + 1
    times = np.arange(n_frames_win) * dt_nuc_fs
    mean_energy = np.zeros(n_frames_win)
    mean_excess_e = np.zeros(n_frames_win)
    mean_excess_h = np.zeros(n_frames_win)
    populations = np.zeros((n_frames_win, len(E0_pairs)))
    trajectory_energies = np.zeros((n_frames_win, n_trajectories)) if is_surface_hopping else None
    all_energies = np.zeros((n_frames_win, len(E0_pairs)))
    all_energies[0, :] = E_init
    biexciton_pop = np.zeros(n_frames_win) if is_biexciton else None

    # -------------------------------------------------------------
    # SCHEME A: Surface Hopping (CPA-FSSH, CPA-FSSH-GDC, and DISH)
    # -------------------------------------------------------------
    if is_surface_hopping:
        active_pairs = [(int(i_pairs0[s]), int(a_pairs0[s])) for s in sampled_states]
        active_surfaces = np.array(sampled_states, dtype=int)

        curr_i_init = np.array([p[0] for p in active_pairs], dtype=int)
        curr_a_init = np.array([p[1] for p in active_pairs], dtype=int)
        curr_i_sub_init = occ_to_sub[curr_i_init]
        curr_a_sub_init = virt_to_sub[curr_a_init]

        C_e = np.zeros((n_virt_dyn, n_trajectories), dtype=np.complex128)
        ok_e_init = curr_a_sub_init >= 0
        if np.any(ok_e_init):
            C_e[curr_a_sub_init[ok_e_init], np.where(ok_e_init)[0]] = 1.0

        C_h = np.zeros((n_occ_dyn, n_trajectories), dtype=np.complex128)
        ok_h_init = curr_i_sub_init >= 0
        if np.any(ok_h_init):
            C_h[curr_i_sub_init[ok_h_init], np.where(ok_h_init)[0]] = 1.0

        biexciton_active = np.ones(n_trajectories, dtype=bool) if is_biexciton else None
        if is_biexciton:
            biexciton_pop[0] = 1.0
            mean_energy[0] = 2.0 * qp_gap_init
            if trajectory_energies is not None:
                trajectory_energies[0, :] = 2.0 * qp_gap_init
        else:
            mean_energy[0] = np.mean(E_init[active_surfaces])
            if trajectory_energies is not None:
                trajectory_energies[0, :] = E_init[active_surfaces]

        if eps_occ_init is not None and eps_virt_init is not None:
            mean_excess_e[0] = np.mean(eps_virt_init[a_pairs0[active_surfaces]] - eps_virt_init[0])
            mean_excess_h[0] = np.mean(eps_occ_init[-1] - eps_occ_init[i_pairs0[active_surfaces]])
        else:
            mean_excess_e[0] = 0.5 * (mean_energy[0] - qp_gap_init)
            mean_excess_h[0] = 0.5 * (mean_energy[0] - qp_gap_init)

        for s in active_surfaces:
            populations[0, s] += 1.0 / n_trajectories

        if verbose:
            scheme_label = "DISH" if method == "dish" else "FSSH-EDC"
            print(f"  [NAMD:{scheme_label}] Starting batched propagation across {n_win_steps} steps ({n_trajectories} trajectories)...")

        for step_idx in range(n_win_steps):
            k = (k0 + step_idx) % n_steps_base
            step_file = os.path.join(precompute_dir, f"step_{k:05d}_to_{k+1:05d}.npz")
            if not os.path.exists(step_file):
                break

            step_data = np.load(step_file)
            E_k = step_data["E_prev"] if "E_prev" in step_data else all_energies[step_idx, :]
            E_kplus1 = step_data["E_curr"]
            all_energies[step_idx + 1, :] = E_kplus1
            i_pairs = step_data["i_pairs"] if "i_pairs" in step_data else (step_data["i_pairs_curr"] if "i_pairs_curr" in step_data else i_pairs0)
            a_pairs = step_data["a_pairs"] if "a_pairs" in step_data else (step_data["a_pairs_curr"] if "a_pairs_curr" in step_data else a_pairs0)
            S_occ = step_data["S_occ"]
            S_virt = step_data["S_virt"]
            eps_occ_curr = step_data["eps_occ_curr"] if "eps_occ_curr" in step_data else eps_occ_init
            eps_virt_curr = step_data["eps_virt_curr"] if "eps_virt_curr" in step_data else eps_virt_init
            eps_occ_prev = step_data["eps_occ_prev"] if "eps_occ_prev" in step_data else eps_occ_init
            eps_virt_prev = step_data["eps_virt_prev"] if "eps_virt_prev" in step_data else eps_virt_init

            d_occ = (S_occ - S_occ.conj().T) / (2.0 * dt_nuc_fs)
            d_virt = (S_virt - S_virt.conj().T) / (2.0 * dt_nuc_fs)

            # Diagonal-BSE pair energies
            if eps_occ_prev is not None and eps_virt_prev is not None:
                E_mat_k = pair_energy_matrix(eps_occ_prev, eps_virt_prev, i_pairs, a_pairs, E_k)
            else:
                E_mat_k = np.zeros((n_occ, n_virt), dtype=np.float64)
                E_mat_k[i_pairs, a_pairs] = E_k

            if eps_occ_curr is not None and eps_virt_curr is not None:
                E_mat_kplus1 = pair_energy_matrix(eps_occ_curr, eps_virt_curr, i_pairs, a_pairs, E_kplus1)
            else:
                E_mat_kplus1 = np.zeros((n_occ, n_virt), dtype=np.float64)
                E_mat_kplus1[i_pairs, a_pairs] = E_kplus1

            d_occ_dyn = d_occ[np.ix_(dyn_occ_active, dyn_occ_active)]
            d_virt_dyn = d_virt[np.ix_(dyn_virt_active, dyn_virt_active)]

            curr_i = np.array([p[0] for p in active_pairs], dtype=int)
            curr_a = np.array([p[1] for p in active_pairs], dtype=int)
            curr_i_sub = occ_to_sub[curr_i]
            curr_a_sub = virt_to_sub[curr_a]

            if verbose and ((step_idx < 5) or ((step_idx + 1) % max(1, n_win_steps // 10) == 0) or (step_idx == n_win_steps - 1)):
                print(f"  [Step {step_idx+1}/{n_win_steps}] Propagating {n_trajectories} trajectories (t = {step_idx*dt_nuc_fs:.1f} -> {(step_idx+1)*dt_nuc_fs:.1f} fs)...", flush=True)
            elif not verbose and (step_idx + 1 == n_win_steps // 2 or step_idx == n_win_steps - 1):
                print(f"    Origin [{origin_idx+1}/{n_origins}] Progress: {step_idx+1}/{n_win_steps} steps (t = {(step_idx+1)*dt_nuc_fs:.1f} fs)", flush=True)

            # 1. Batched Electron channel propagation
            ok_e = curr_a_sub >= 0
            active_e = np.where(ok_e, curr_a_sub, 0)
            E_e_k_batch = E_mat_k[curr_i, :][:, dyn_virt_active].T
            E_e_kplus1_batch = E_mat_kplus1[curr_i, :][:, dyn_virt_active].T

            if method == "dish":
                if integrator_type in ("strang", "trotter", "unitary_strang", "split_operator"):
                    C_e = propagate_channel_batch_strang(
                        C_e, E_e_k_batch, E_e_kplus1_batch, d_virt_dyn, dt_nuc_fs,
                        n_substeps, device=device_cfg, active_idx=None,
                    )
                else:
                    C_e = propagate_channel_batch_rk4(
                        C_e, E_e_k_batch, E_e_kplus1_batch, d_virt_dyn, dt_nuc_fs,
                        n_substeps, active_idx=None,
                    )
                flux_e = None
            else:
                if integrator_type in ("strang", "trotter", "unitary_strang", "split_operator"):
                    C_e, flux_e = propagate_channel_batch_strang(
                        C_e, E_e_k_batch, E_e_kplus1_batch, d_virt_dyn, dt_nuc_fs,
                        n_substeps, device=device_cfg, active_idx=active_e,
                    )
                else:
                    C_e, flux_e = propagate_channel_batch_rk4(
                        C_e, E_e_k_batch, E_e_kplus1_batch, d_virt_dyn, dt_nuc_fs,
                        n_substeps, active_idx=active_e,
                    )

            # 2. Batched Hole channel propagation
            ok_h = curr_i_sub >= 0
            active_h = np.where(ok_h, curr_i_sub, 0)
            E_h_k_batch = E_mat_k[:, curr_a][dyn_occ_active, :]
            E_h_kplus1_batch = E_mat_kplus1[:, curr_a][dyn_occ_active, :]

            if method == "dish":
                if integrator_type in ("strang", "trotter", "unitary_strang", "split_operator"):
                    C_h = propagate_channel_batch_strang(
                        C_h, E_h_k_batch, E_h_kplus1_batch, d_occ_dyn, dt_nuc_fs,
                        n_substeps, device=device_cfg, active_idx=None,
                    )
                else:
                    C_h = propagate_channel_batch_rk4(
                        C_h, E_h_k_batch, E_h_kplus1_batch, d_occ_dyn, dt_nuc_fs,
                        n_substeps, active_idx=None,
                    )
                flux_h_mat = None
            else:
                if integrator_type in ("strang", "trotter", "unitary_strang", "split_operator"):
                    C_h, flux_h_mat = propagate_channel_batch_strang(
                        C_h, E_h_k_batch, E_h_kplus1_batch, d_occ_dyn, dt_nuc_fs,
                        n_substeps, device=device_cfg, active_idx=active_h,
                    )
                else:
                    C_h, flux_h_mat = propagate_channel_batch_rk4(
                        C_h, E_h_k_batch, E_h_kplus1_batch, d_occ_dyn, dt_nuc_fs,
                        n_substeps, active_idx=active_h,
                    )

            # 3. Surface Hopping (DISH vs CPA-FSSH)
            if method == "dish":
                C_e, new_curr_a_sub, hopped_e = step_dish_batch(
                    C_e, curr_a_sub, E_e_kplus1_batch, dt_nuc_fs,
                    tau_mat=tau_virt_dyn, beta=beta,
                    detailed_balance=detailed_balance,
                )
                C_h, new_curr_i_sub, hopped_h = step_dish_batch(
                    C_h, curr_i_sub, E_h_kplus1_batch, dt_nuc_fs,
                    tau_mat=tau_occ_dyn, beta=beta,
                    detailed_balance=detailed_balance,
                )

                if run_tr_sd and hop_records is not None:
                    t_curr_fs = (step_idx + 1) * dt_nuc_fs
                    if np.any(hopped_e):
                        for tr_h in np.where(hopped_e)[0]:
                            if ok_e[tr_h] and new_curr_a_sub[tr_h] >= 0:
                                hop_records.append({
                                    "time_fs": t_curr_fs,
                                    "channel": "electron",
                                    "from": int(curr_a[tr_h]),
                                    "to": int(dyn_virt_active[new_curr_a_sub[tr_h]]),
                                    "traj": int(origin_idx * n_trajectories + tr_h)
                                })
                    if np.any(hopped_h):
                        for tr_h in np.where(hopped_h)[0]:
                            if ok_h[tr_h] and new_curr_i_sub[tr_h] >= 0:
                                hop_records.append({
                                    "time_fs": t_curr_fs,
                                    "channel": "hole",
                                    "from": int(curr_i[tr_h]),
                                    "to": int(dyn_occ_active[new_curr_i_sub[tr_h]]),
                                    "traj": int(origin_idx * n_trajectories + tr_h)
                                })

                for tr in range(n_trajectories):
                    if is_biexciton and biexciton_active[tr]:
                        continue
                    new_i = dyn_occ_active[new_curr_i_sub[tr]] if ok_h[tr] and new_curr_i_sub[tr] >= 0 else curr_i[tr]
                    new_a = dyn_virt_active[new_curr_a_sub[tr]] if ok_e[tr] and new_curr_a_sub[tr] >= 0 else curr_a[tr]
                    active_pairs[tr] = (new_i, new_a)

                    if ecsh_auger and not hopped_e[tr] and not hopped_h[tr] and k_auger_fs is not None:
                        E_curr = E_mat_kplus1[new_i, new_a]
                        e_diff_all = np.abs(E_mat_kplus1[dyn_occ_active[:, None], dyn_virt_active[None, :]] - E_curr)
                        mask_2body = (dyn_occ_active[:, None] != new_i) & (dyn_virt_active[None, :] != new_a)
                        mask_resonant = (e_diff_all <= ecsh_window_ev) & mask_2body
                        cand_j, cand_b = np.where(mask_resonant)
                        n_cand = len(cand_j)
                        if n_cand > 0:
                            k_each = k_auger_fs / n_cand
                            p_each = 1.0 - np.exp(-k_each * dt_nuc_fs)
                            p_tot = min(1.0, p_each * n_cand)
                            if np.random.rand() < p_tot:
                                sel = np.random.randint(n_cand)
                                new_i = int(dyn_occ_active[cand_j[sel]])
                                new_a = int(dyn_virt_active[cand_b[sel]])
                                active_pairs[tr] = (new_i, new_a)
                                sub_a = virt_to_sub[new_a]
                                sub_i = occ_to_sub[new_i]
                                if sub_a >= 0:
                                    C_e[:, tr] = 0.0
                                    C_e[sub_a, tr] = 1.0
                                if sub_i >= 0:
                                    C_h[:, tr] = 0.0
                                    C_h[sub_i, tr] = 1.0

                    new_idx = pair_lookup[active_pairs[tr][0], active_pairs[tr][1]]
                    active_surfaces[tr] = new_idx if new_idx >= 0 else active_surfaces[tr]

            else:
                hopped_e_traj = np.zeros(n_trajectories, dtype=bool)
                hopped_h_traj = np.zeros(n_trajectories, dtype=bool)

                for tr in range(n_trajectories):
                    if is_biexciton and biexciton_active[tr]:
                        p_hop_auger = 1.0 - np.exp(-k_auger_fs * dt_nuc_fs)
                        if np.random.rand() < p_hop_auger:
                            target_e = 2.0 * qp_gap_init
                            pair_diffs = np.abs(E_mat_kplus1[i_pairs0, a_pairs0] - target_e)
                            valid_cands = np.where(pair_diffs <= ecsh_window_ev)[0]
                            if len(valid_cands) > 0:
                                chosen = int(np.random.choice(valid_cands))
                                active_pairs[tr] = (int(i_pairs0[chosen]), int(a_pairs0[chosen]))
                                active_surfaces[tr] = chosen
                                biexciton_active[tr] = False
                                new_a = virt_to_sub[active_pairs[tr][1]]
                                new_i = occ_to_sub[active_pairs[tr][0]]
                                if new_a >= 0:
                                    C_e[:, tr] = 0.0
                                    C_e[new_a, tr] = 1.0
                                if new_i >= 0:
                                    C_h[:, tr] = 0.0
                                    C_h[new_i, tr] = 1.0
                        continue

                    i_c = curr_i[tr]
                    a_c = curr_a[tr]
                    a_sub = curr_a_sub[tr]
                    i_sub = curr_i_sub[tr]
                    E_curr = E_mat_kplus1[i_c, a_c]

                    if a_sub >= 0:
                        probs_b = np.array(flux_e[:, tr], copy=True)
                        probs_b[a_sub] = 0.0
                    else:
                        probs_b = np.zeros(n_virt_dyn, dtype=np.float64)
                    valid_b = pair_lookup[i_c, dyn_virt_active] >= 0
                    probs_b[~valid_b] = 0.0
                    dE_b = E_mat_kplus1[i_c, dyn_virt_active] - E_curr
                    if detailed_balance:
                        probs_b *= np.exp(-np.maximum(dE_b, 0.0) * beta)

                    if i_sub >= 0:
                        probs_j = np.array(flux_h_mat[:, tr], copy=True)
                        probs_j[i_sub] = 0.0
                    else:
                        probs_j = np.zeros(n_occ_dyn, dtype=np.float64)
                    valid_j = pair_lookup[dyn_occ_active, a_c] >= 0
                    probs_j[~valid_j] = 0.0
                    dE_j = E_mat_kplus1[dyn_occ_active, a_c] - E_curr
                    if detailed_balance:
                        probs_j *= np.exp(-np.maximum(dE_j, 0.0) * beta)

                    total_b = np.sum(probs_b)
                    total_j = np.sum(probs_j)
                    total_hop = total_b + total_j

                    zeta = np.random.rand()
                    if zeta < total_b:
                        p_norm = probs_b / total_b
                        new_a_sub = np.random.choice(n_virt_dyn, p=p_norm)
                        new_a_val = int(dyn_virt_active[new_a_sub])
                        if run_tr_sd and hop_records is not None:
                            hop_records.append({
                                "time_fs": (step_idx + 1) * dt_nuc_fs,
                                "channel": "electron",
                                "from": int(a_c),
                                "to": new_a_val,
                                "traj": int(origin_idx * n_trajectories + tr)
                            })
                        active_pairs[tr] = (i_c, new_a_val)
                        hopped_e_traj[tr] = True
                        C_e[:, tr] = 0.0
                        C_e[new_a_sub, tr] = 1.0
                    elif zeta < total_hop:
                        p_norm = probs_j / total_j
                        new_i_sub = np.random.choice(n_occ_dyn, p=p_norm)
                        new_i_val = int(dyn_occ_active[new_i_sub])
                        if run_tr_sd and hop_records is not None:
                            hop_records.append({
                                "time_fs": (step_idx + 1) * dt_nuc_fs,
                                "channel": "hole",
                                "from": int(i_c),
                                "to": new_i_val,
                                "traj": int(origin_idx * n_trajectories + tr)
                            })
                        active_pairs[tr] = (new_i_val, a_c)
                        hopped_h_traj[tr] = True
                        C_h[:, tr] = 0.0
                        C_h[new_i_sub, tr] = 1.0
                    else:
                        if ecsh_auger and k_auger_fs is not None:
                            e_diff_all = np.abs(E_mat_kplus1[dyn_occ_active[:, None], dyn_virt_active[None, :]] - E_curr)
                            mask_2body = (dyn_occ_active[:, None] != i_c) & (dyn_virt_active[None, :] != a_c)
                            mask_resonant = (e_diff_all <= ecsh_window_ev) & mask_2body
                            cand_j, cand_b = np.where(mask_resonant)
                            n_cand = len(cand_j)
                            if n_cand > 0:
                                k_each = k_auger_fs / n_cand
                                p_each = 1.0 - np.exp(-k_each * dt_nuc_fs)
                                p_tot = min(1.0, p_each * n_cand)
                                if np.random.rand() < p_tot:
                                    sel = np.random.randint(n_cand)
                                    new_i = int(dyn_occ_active[cand_j[sel]])
                                    new_a = int(dyn_virt_active[cand_b[sel]])
                                    active_pairs[tr] = (new_i, new_a)
                                    sub_a = virt_to_sub[new_a]
                                    sub_i = occ_to_sub[new_i]
                                    if sub_a >= 0:
                                        C_e[:, tr] = 0.0
                                        C_e[sub_a, tr] = 1.0
                                        hopped_e_traj[tr] = True
                                    if sub_i >= 0:
                                        C_h[:, tr] = 0.0
                                        C_h[sub_i, tr] = 1.0
                                        hopped_h_traj[tr] = True

                    new_idx = pair_lookup[active_pairs[tr][0], active_pairs[tr][1]]
                    active_surfaces[tr] = new_idx if new_idx >= 0 else active_surfaces[tr]

                curr_a_sub_after = virt_to_sub[np.array([p[1] for p in active_pairs])]
                curr_i_sub_after = occ_to_sub[np.array([p[0] for p in active_pairs])]

                non_hopped_e = ~hopped_e_traj & ok_e
                if np.any(non_hopped_e):
                    idx_e = np.where(non_hopped_e)[0]
                    C_e[:, idx_e] = apply_edc_decoherence_batch(
                        C_e[:, idx_e],
                        curr_a_sub_after[idx_e],
                        E_e_kplus1_batch[:, idx_e],
                        dt_nuc_fs,
                        tau_mat=tau_virt_dyn,
                        temp_k=temp_k,
                        decay_type=decoherence_decay_type,
                    )

                non_hopped_h = ~hopped_h_traj & ok_h
                if np.any(non_hopped_h):
                    idx_h = np.where(non_hopped_h)[0]
                    C_h[:, idx_h] = apply_edc_decoherence_batch(
                        C_h[:, idx_h],
                        curr_i_sub_after[idx_h],
                        E_h_kplus1_batch[:, idx_h],
                        dt_nuc_fs,
                        tau_mat=tau_occ_dyn,
                        temp_k=temp_k,
                        decay_type=decoherence_decay_type,
                    )

            if is_biexciton:
                biexciton_pop[step_idx + 1] = np.mean(biexciton_active)
                e_tr = np.zeros(n_trajectories)
                for tr in range(n_trajectories):
                    if biexciton_active[tr]:
                        e_tr[tr] = 2.0 * qp_gap_init
                    else:
                        e_tr[tr] = E_kplus1[active_surfaces[tr]]
                mean_energy[step_idx + 1] = np.mean(e_tr)
                trajectory_energies[step_idx + 1, :] = e_tr
            else:
                mean_energy[step_idx + 1] = np.mean(E_kplus1[active_surfaces])
                trajectory_energies[step_idx + 1, :] = E_kplus1[active_surfaces]

            for s in active_surfaces:
                populations[step_idx + 1, s] += 1.0 / n_trajectories

            curr_a_new = np.array([p[1] for p in active_pairs])
            curr_i_new = np.array([p[0] for p in active_pairs])
            if eps_occ_curr is not None and eps_virt_curr is not None:
                mean_excess_e[step_idx + 1] = np.mean(eps_virt_curr[curr_a_new] - eps_virt_curr[0])
                mean_excess_h[step_idx + 1] = np.mean(eps_occ_curr[-1] - eps_occ_curr[curr_i_new])
            else:
                mean_excess_e[step_idx + 1] = 0.5 * (mean_energy[step_idx + 1] - qp_gap_init)
                mean_excess_h[step_idx + 1] = 0.5 * (mean_energy[step_idx + 1] - qp_gap_init)

            if verbose and ((step_idx < 5) or ((step_idx + 1) % max(1, n_win_steps // 10) == 0) or (step_idx == n_win_steps - 1)):
                if not is_biexciton:
                    min_e_act = np.min(E_kplus1[active_surfaces])
                    max_e_act = np.max(E_kplus1[active_surfaces])
                    print(f"    -> [Step {step_idx+1}/{n_win_steps}] <E_exc> = {mean_energy[step_idx+1]:.4f} eV | Active range: [{min_e_act:.3f}, {max_e_act:.3f}] eV", flush=True)

    # -------------------------------------------------------------
    # SCHEME B: Pauli Master Equation (Deterministic Kinetics)
    # -------------------------------------------------------------
    elif method in ("master_equation", "pme"):
        P_mat = np.zeros((n_occ, n_virt), dtype=np.float64)
        P_mat[i_pairs0, a_pairs0] = P_init
        mean_energy[0] = np.sum(P_init * E_init)
        populations[0] = P_init

        p_virt0 = np.sum(P_mat, axis=0)
        p_occ0 = np.sum(P_mat, axis=1)
        if eps_occ_init is not None and eps_virt_init is not None:
            mean_excess_e[0] = np.sum(p_virt0 * (eps_virt_init - eps_virt_init[0]))
            mean_excess_h[0] = np.sum(p_occ0 * (eps_occ_init[-1] - eps_occ_init))
        else:
            mean_excess_e[0] = 0.5 * (mean_energy[0] - qp_gap_init)
            mean_excess_h[0] = 0.5 * (mean_energy[0] - qp_gap_init)

        if verbose:
            print(f"  [NAMD:PME] Starting deterministic Master Equation propagation (n_occ={n_occ}, n_virt={n_virt})...")

        for step_idx in range(n_win_steps):
            t_step_start = time.time()
            k = (k0 + step_idx) % n_steps_base
            step_file = os.path.join(precompute_dir, f"step_{k:05d}_to_{k+1:05d}.npz")
            if not os.path.exists(step_file):
                break

            step_data = np.load(step_file)
            E_k = step_data["E_curr"]
            all_energies[step_idx + 1, :] = E_k
            i_pairs = step_data["i_pairs"] if "i_pairs" in step_data else (step_data["i_pairs_curr"] if "i_pairs_curr" in step_data else i_pairs0)
            a_pairs = step_data["a_pairs"] if "a_pairs" in step_data else (step_data["a_pairs_curr"] if "a_pairs_curr" in step_data else a_pairs0)
            S_occ = step_data["S_occ"]
            S_virt = step_data["S_virt"]
            eps_occ_curr = step_data["eps_occ_curr"] if "eps_occ_curr" in step_data else eps_occ_init
            eps_virt_curr = step_data["eps_virt_curr"] if "eps_virt_curr" in step_data else eps_virt_init

            d_occ = (S_occ - S_occ.conj().T) / (2.0 * dt_nuc_fs)
            d_virt = (S_virt - S_virt.conj().T) / (2.0 * dt_nuc_fs)

            if n_occ_dyn < n_occ or n_virt_dyn < n_virt:
                d_occ_pme = np.zeros_like(d_occ)
                d_occ_pme[np.ix_(dyn_occ_active, dyn_occ_active)] = d_occ[np.ix_(dyn_occ_active, dyn_occ_active)]
                d_virt_pme = np.zeros_like(d_virt)
                d_virt_pme[np.ix_(dyn_virt_active, dyn_virt_active)] = d_virt[np.ix_(dyn_virt_active, dyn_virt_active)]
            else:
                d_occ_pme = d_occ
                d_virt_pme = d_virt

            E_mat_k = np.zeros((n_occ, n_virt), dtype=np.float64)
            E_mat_k[i_pairs, a_pairs] = E_k

            if run_tr_sd and pme_flux_records is not None:
                P_mat, flux_virt, flux_occ = propagate_pme_tensor(
                    P_mat=P_mat,
                    E_mat=E_mat_k,
                    d_occ=d_occ_pme,
                    d_virt=d_virt_pme,
                    dt_fs=dt_nuc_fs,
                    temp_k=temp_k,
                    tau_dec_fs=tau_dec_fs,
                    eps_occ=eps_occ_curr,
                    eps_virt=eps_virt_curr,
                    k_loss=k_loss_mat,
                    return_flux=True,
                )
                pme_flux_records.append({
                    "time_fs": (step_idx + 1) * dt_nuc_fs,
                    "flux_virt": flux_virt,
                    "flux_occ": flux_occ,
                })
            else:
                P_mat = propagate_pme_tensor(
                    P_mat=P_mat,
                    E_mat=E_mat_k,
                    d_occ=d_occ_pme,
                    d_virt=d_virt_pme,
                    dt_fs=dt_nuc_fs,
                    temp_k=temp_k,
                    tau_dec_fs=tau_dec_fs,
                    eps_occ=eps_occ_curr,
                    eps_virt=eps_virt_curr,
                    k_loss=k_loss_mat,
                )

            P_vec = P_mat[i_pairs, a_pairs]
            mean_energy[step_idx + 1] = np.sum(P_vec * E_k)
            populations[step_idx + 1] = P_vec

            p_virt = np.sum(P_mat, axis=0)
            p_occ = np.sum(P_mat, axis=1)
            if eps_occ_curr is not None and eps_virt_curr is not None:
                mean_excess_e[step_idx + 1] = np.sum(p_virt * (eps_virt_curr - eps_virt_curr[0]))
                mean_excess_h[step_idx + 1] = np.sum(p_occ * (eps_occ_curr[-1] - eps_occ_curr))
            else:
                mean_excess_e[step_idx + 1] = 0.5 * (mean_energy[step_idx + 1] - qp_gap_init)
                mean_excess_h[step_idx + 1] = 0.5 * (mean_energy[step_idx + 1] - qp_gap_init)

            t_step = time.time() - t_step_start
            if verbose and ((step_idx < 5) or ((step_idx + 1) % max(1, n_win_steps // 10) == 0) or (step_idx == n_win_steps - 1)):
                print(f"    Step {step_idx+1}/{n_win_steps} (t = {(step_idx+1)*dt_nuc_fs:.1f} fs) in {t_step:.3f} s | <E_exc> = {mean_energy[step_idx+1]:.4f} eV")
            elif not verbose and (step_idx + 1 == n_win_steps // 2 or step_idx == n_win_steps - 1):
                print(f"    Origin [{origin_idx+1}/{n_origins}] Progress: {step_idx+1}/{n_win_steps} steps (t = {(step_idx+1)*dt_nuc_fs:.1f} fs) | <E_exc> = {mean_energy[step_idx+1]:.4f} eV", flush=True)

    return {
        "mean_energy": mean_energy,
        "mean_excess_e": mean_excess_e,
        "mean_excess_h": mean_excess_h,
        "populations": populations,
        "trajectory_energies": trajectory_energies,
        "all_energies": all_energies,
        "biexciton_pop": biexciton_pop,
        "times": times,
        "k0": k0
    }


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
        traj_dir = namd_cfg.get("trajectory", {}).get("dir", "")
        if traj_dir and os.path.isdir(traj_dir):
            cand = os.path.join(traj_dir, precompute_dir)
            if os.path.exists(os.path.join(cand, "namd_metadata.npz")) and os.path.exists(os.path.join(cand, "frame_00000.npz")):
                precompute_dir = cand
                meta_path = os.path.join(precompute_dir, "namd_metadata.npz")
                frame0_path = os.path.join(precompute_dir, "frame_00000.npz")

    if not os.path.exists(meta_path) or not os.path.exists(frame0_path):
        raise FileNotFoundError(
            f"Precomputed NAMD data not found in '{precompute_dir}'. "
            "Please run 'minibse --namd-precompute' first!"
        )

    meta = np.load(meta_path)
    n_frames = int(meta["n_frames"])
    dt_nuc_fs = float(meta["dt_nuc_fs"]) if "dt_nuc_fs" in meta else float(meta.get("dt_fs", 2.0))

    frame0 = np.load(frame0_path)
    E0_pairs = frame0["E_pairs"]
    f_convention = None
    if "f_convention" in frame0.files:
        f_convention = str(np.array(frame0["f_convention"]).reshape(-1)[0])
    elif "f_convention" in meta.files:
        f_convention = str(np.array(meta["f_convention"]).reshape(-1)[0])
    soc_meta = bool(meta["soc"]) if "soc" in meta.files else False
    f0_pairs = rescale_legacy_oscillator_strengths(
        frame0["f_pairs"], soc=soc_meta, f_convention=f_convention
    )
    dft_gap = float(meta["mean_dft_gap"]) if "mean_dft_gap" in meta else float(frame0.get("dft_gap", meta.get("dft_gap", 1.0)))
    scissor = float(frame0["scissor"]) if "scissor" in frame0 else 0.0

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
    if method in ("cpa_fssh_gdc", "fssh_gdc", "gdc") or str(decoherence).lower() in ("gdc", "gaussian"):
        decoherence_decay_type = "gaussian"
        if method in ("cpa_fssh_gdc", "fssh_gdc", "gdc"):
            method = "cpa_fssh_gdc"
            decoherence = "gdc"
    else:
        decoherence_decay_type = "exponential"
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
            print(
                f"  [NAMD] Cumulant dephasing time of the lowest exciton: "
                f"tau_dec = {tau_dec_fs:.2f} fs (gap fluctuation std: {std_g:.4f} eV). "
                "Diagnostic only: each nuclear step starts from the active orbital, "
                "so this time is not applied inside the step."
            )
            decoherence = f"cumulant diagnostic ({tau_dec_fs:.1f} fs)"
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
    if integrator_type == "unitary_matrix_exp":
        raise NotImplementedError(
            "integrator 'unitary_matrix_exp' is the single-wavefunction spectral step "
            "in step_unitary_matrix_exp. The ensemble path is Strang. "
            "Set namd.integration.integrator to 'strang'."
        )
    if integrator_type in ("strang", "trotter", "unitary_strang", "split_operator"):
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

    phys_cfg = config.get("physics", {})
    req_nhomos = phys_cfg.get("nhomos", None)
    req_nlumos = phys_cfg.get("nlumos", None)
    if req_nhomos is not None:
        scale_occ = 2 if soc_meta and n_occ == 2 * int(meta.get("nhomos", n_occ // 2)) else 1
        cutoff_i = n_occ - int(req_nhomos) * scale_occ
        if cutoff_i > 0:
            dyn_mask = dyn_mask & (i_pairs0 >= cutoff_i)
    if req_nlumos is not None:
        scale_virt = 2 if soc_meta and n_virt == 2 * int(meta.get("nlumos", n_virt // 2)) else 1
        cutoff_a = int(req_nlumos) * scale_virt
        if cutoff_a < n_virt:
            dyn_mask = dyn_mask & (a_pairs0 < cutoff_a)

    if np.sum(dyn_mask) == 0:
        dyn_mask = np.ones(len(E0_pairs), dtype=bool)
        if req_nhomos is not None:
            dyn_mask = dyn_mask & (i_pairs0 >= max(0, cutoff_i))
        if req_nlumos is not None:
            dyn_mask = dyn_mask & (a_pairs0 < min(n_virt, cutoff_a))
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

    # Load state-pair decoherence matrices if available
    dec_file = os.path.join(precompute_dir, "decoherence_times.npz")
    tau_occ_mat = None
    tau_virt_mat = None
    if os.path.exists(dec_file):
        try:
            d_dec = np.load(dec_file)
            tau_occ_mat = d_dec["tau_occ"]
            tau_virt_mat = d_dec["tau_virt"]
            print(f"  [NAMD] Loaded state-pair decoherence matrices from '{dec_file}':")
            print(f"         tau_occ : {tau_occ_mat.shape} (median = {np.median(tau_occ_mat):.2f} fs)")
            print(f"         tau_virt: {tau_virt_mat.shape} (median = {np.median(tau_virt_mat):.2f} fs)")
        except Exception as e:
            print(f"  [NAMD:Warn] Failed loading {dec_file}: {e}")

    if tau_occ_mat is not None and tau_occ_mat.shape[0] >= n_occ:
        tau_occ_dyn = tau_occ_mat[np.ix_(dyn_occ_active, dyn_occ_active)]
    else:
        tau_val = float(tau_dec_fs) if isinstance(tau_dec_fs, (int, float)) and tau_dec_fs > 0 else 20.0
        tau_occ_dyn = np.full((n_occ_dyn, n_occ_dyn), tau_val)
        np.fill_diagonal(tau_occ_dyn, 500.0)

    if tau_virt_mat is not None and tau_virt_mat.shape[0] >= n_virt:
        tau_virt_dyn = tau_virt_mat[np.ix_(dyn_virt_active, dyn_virt_active)]
    else:
        tau_val = float(tau_dec_fs) if isinstance(tau_dec_fs, (int, float)) and tau_dec_fs > 0 else 20.0
        tau_virt_dyn = np.full((n_virt_dyn, n_virt_dyn), tau_val)
        np.fill_diagonal(tau_virt_dyn, 500.0)

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
        S_hr = 1.0
        try:
            from qdex.namd.analysis import compute_band_gap_dynamics_and_spectral_density
            bg_for_gap = compute_band_gap_dynamics_and_spectral_density(
                precompute_dir, use_lowest_exciton=True, temp_k=temp_k
            )
            if bg_for_gap is not None and bg_for_gap.get("recomb_params"):
                e_lo = float(bg_for_gap["recomb_params"]["E_LO_ev"])
                S_hr = float(bg_for_gap["recomb_params"]["S_hr"])
        except Exception:
            pass
        A_nr = float(dyn_cfg.get("A_nr_s", dyn_cfg.get("A_nr", 1e13)))
        k_nr_s, k_nr_fs = compute_energy_gap_law_rate(
            qp_gap, E_LO_ev=e_lo, S_hr=S_hr, A_nr=A_nr
        )
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
    i_1s = int(np.argmin(E0_pairs)) if len(E0_pairs) else 0
    tau_rad_1_ns = (1e9 / k_rad_s[i_1s]) if (len(k_rad_s) > i_1s and k_rad_s[i_1s] > 0) else np.inf
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
    if method in ("cpa_fssh", "fssh", "dish", "cpa_fssh_gdc", "cpa_fssh_edc", "fssh_gdc", "fssh_edc"):
        if method == "dish":
            method_name = "DISH (Decoherence-Induced Surface Hopping)"
        elif method == "cpa_fssh_gdc" or decoherence_decay_type == "gaussian":
            method_name = "CPA-FSSH-GDC (Gaussian Decoherence Correction)"
        else:
            method_name = "CPA-FSSH-EDC (Energy-based Decoherence Correction)"
        print(f"  Trajectories         : {n_trajectories}")
        print(f"  Hopping Scheme       : {method_name}")
        print(f"  Decoherence Model    : {decoherence} (tau_dec = {tau_dec_fs} fs, decay = {decoherence_decay_type})")
        print(f"  Electronic Sub-steps : {n_substeps} (dt_elec = {dt_nuc_fs/n_substeps:.5f} fs)")
        print(f"  Integrator           : {integrator_type} (device: {device_cfg})")
    print("=" * 65 + "\n")

    trajectory_loops = int(dyn_cfg.get("trajectory_loops", 1))
    ecsh_auger = bool(dyn_cfg.get("ecsh_auger", False))
    initial_state_mode = str(dyn_cfg.get("initial_state", "ratio_eg")).lower()
    is_biexciton = (initial_state_mode == "biexciton")
    kT_ev = KB_EV * max(temp_k, 1e-3)
    if "ecsh_window_ev" in dyn_cfg:
        ecsh_window_ev = float(dyn_cfg["ecsh_window_ev"])
    elif ecsh_auger or is_biexciton:
        sigma_aug = float(config.get("auger", {}).get("sigma", kT_ev)) if config.get("auger") else kT_ev
        ecsh_window_ev = max(kT_ev, sigma_aug)
    else:
        ecsh_window_ev = kT_ev
    k_auger_fs = None
    if ecsh_auger or is_biexciton:
        k_auger_fs, auger_note = resolve_optin_auger_rate_fs(config, dyn_cfg)
        print("  Auger Dynamics Mode  : opt-in ECSH / biexciton (off for a plain cooling run)")
        print(f"  ECSH Resonance Window: {ecsh_window_ev*1e3:.2f} meV")
        print(f"  Auger rate           : {k_auger_fs:.4e} fs^-1 ({auger_note})")
        if is_biexciton:
            print("  Initial State        : Biexciton (XX) -> Auger annihilation")
    if trajectory_loops > 1:
        print(f"  Trajectory Looping   : {trajectory_loops} loops (simulating {(n_frames - 1) * trajectory_loops * dt_nuc_fs:.1f} fs / {(n_frames - 1) * trajectory_loops * dt_nuc_fs * 1e-3:.2f} ps)")

    n_steps_base = max(1, n_frames - 1)
    n_steps_total = n_steps_base * max(1, trajectory_loops)
    n_frames_total = n_steps_total + 1

    # Check initial conditions / ensemble mode
    init_cond_mode = str(dyn_cfg.get("initial_conditions", dyn_cfg.get("ensemble_mode", "single"))).lower()
    is_multi_origin = init_cond_mode in ("multiple", "multi", "ensemble", "many")

    # Time-resolved vibrational action spectrum / dynamical phonon spectrogram tracking
    tr_sd_cfg = namd_cfg.get("time_resolved_spectral_density", {})
    run_tr_sd = bool(tr_sd_cfg.get("run", False) or tr_sd_cfg.get("enabled", False))
    hop_records = [] if run_tr_sd else None
    pme_flux_records = [] if run_tr_sd else None

    is_surface_hopping = method in ("cpa_fssh", "fssh", "dish", "cpa_fssh_gdc", "cpa_fssh_edc", "fssh_gdc", "fssh_edc")

    if not is_multi_origin:
        n_origins = 1
        origin_frames = np.array([0], dtype=int)
        n_win_steps = n_steps_total
        n_traj_per_origin = n_trajectories
        calib = None
    else:
        calib = auto_calibrate_ensemble_origins(
            precompute_dir=precompute_dir,
            qp_gap_ev=qp_gap,
            pump_energy_ev=pump_energy_ev,
            dyn_cfg=dyn_cfg,
            config=config
        )
        n_origins = calib["n_origins"]
        origin_frames = calib["origin_frames"]
        n_win_steps = calib["n_win_steps"]
        if is_surface_hopping:
            n_traj_per_origin = max(1, int(np.ceil(n_trajectories / n_origins)))
        else:
            n_traj_per_origin = n_trajectories

        print("=" * 68)
        print("  [NAMD] Automated Multi-Origin Ensemble Setup")
        print("=" * 68)
        print(f"  Trajectory Length Available : {calib['t_md_total_fs']:.1f} fs ({n_frames} frames)")
        print(f"  Phonon Dephasing (tau_corr) : {calib['tau_corr_fs']:.1f} fs -> Safe spacing dt0 = {calib['dt0_fs']:.1f} fs")
        print(f"  Pilot Cooling (tau_cool)    : {calib['tau_cool_fs']:.1f} fs -> Simulation window = {calib['window_fs']:.1f} fs ({n_win_steps} steps)")
        if calib["is_fallback"]:
            print("  Notice                      : MD trajectory length is comparable to cooling window.")
            print("                                Executing single origin from t0 = 0 fs.")
        else:
            print(f"  Ensemble Origins Generated  : {n_origins} independent AIMD origins (t0 = {calib['origin_times_fs'][0]:.1f} .. {calib['origin_times_fs'][-1]:.1f} fs)")
            if is_surface_hopping:
                print(f"  Trajectory Allocation       : {n_traj_per_origin} traj/origin ({n_origins * n_traj_per_origin} total across ensemble)")
        print("=" * 68 + "\n")


    origin_results = []
    for m, k0 in enumerate(origin_frames):
        if is_multi_origin and not (calib and calib["is_fallback"]):
            print(f"\n  >>> Propagating Ensemble Origin [{m+1}/{n_origins}]: frame {k0:05d} (t0 = {k0 * dt_nuc_fs:.1f} fs) ...")
        res_m = propagate_single_namd_origin(
            k0=k0,
            n_win_steps=n_win_steps,
            n_trajectories=n_traj_per_origin,
            precompute_dir=precompute_dir,
            frame0=frame0,
            E0_pairs=E0_pairs,
            f0_pairs=f0_pairs,
            i_pairs0=i_pairs0,
            a_pairs0=a_pairs0,
            qp_gap_0=qp_gap,
            pump_energy_ev=pump_energy_ev,
            pulse_fwhm_ev=pulse_fwhm_ev,
            filter_dark=filter_dark,
            dyn_mask=dyn_mask,
            dyn_occ_active=dyn_occ_active,
            dyn_virt_active=dyn_virt_active,
            occ_to_sub=occ_to_sub,
            virt_to_sub=virt_to_sub,
            n_occ=n_occ,
            n_virt=n_virt,
            n_occ_dyn=n_occ_dyn,
            n_virt_dyn=n_virt_dyn,
            n_steps_base=n_steps_base,
            dt_nuc_fs=dt_nuc_fs,
            method=method,
            decoherence=decoherence,
            decoherence_decay_type=decoherence_decay_type,
            tau_dec_fs=tau_dec_fs,
            integrator_type=integrator_type,
            n_substeps=n_substeps,
            device_cfg=device_cfg,
            temp_k=temp_k,
            beta=beta,
            detailed_balance=detailed_balance,
            pair_lookup=pair_lookup,
            tau_occ_dyn=tau_occ_dyn,
            tau_virt_dyn=tau_virt_dyn,
            is_surface_hopping=is_surface_hopping,
            is_biexciton=is_biexciton,
            ecsh_auger=ecsh_auger,
            k_auger_fs=k_auger_fs,
            ecsh_window_ev=ecsh_window_ev,
            k_loss_mat=k_loss_mat,
            run_tr_sd=run_tr_sd,
            hop_records=hop_records,
            pme_flux_records=pme_flux_records,
            verbose=(not is_multi_origin or (calib and calib["is_fallback"]) or (m == 0)),
            origin_idx=m,
            n_origins=n_origins
        )
        origin_results.append(res_m)

    times = origin_results[0]["times"]
    agg = aggregate_multi_origin_results(origin_results, times_fs=times)

    mean_energy = agg["mean_energy"]
    std_energies_ev = agg["std_energy"] if (is_multi_origin and not (calib and calib["is_fallback"])) else None
    mean_excess_e = agg["mean_excess_e"]
    std_excess_e = agg["std_excess_e"] if (is_multi_origin and not (calib and calib["is_fallback"])) else None
    mean_excess_h = agg["mean_excess_h"]
    std_excess_h = agg["std_excess_h"] if (is_multi_origin and not (calib and calib["is_fallback"])) else None
    populations = agg["populations"]
    trajectory_energies = agg["trajectory_energies"]
    all_energies = agg["all_energies"]
    biexciton_pop = origin_results[0]["biexciton_pop"]

    total_sim_time = time.time() - t0_start
    print(f"\n[NAMD Dynamics] Completed in {total_sim_time:.2f} s")

    if is_biexciton and biexciton_pop is not None:
        biex_file = os.path.join(precompute_dir if os.path.isdir(precompute_dir) else ".", "biexciton_decay.csv")
        np.savetxt(
            biex_file,
            np.column_stack([times, biexciton_pop, mean_energy]),
            header="time_fs,P_biexciton,mean_energy_ev",
            delimiter=",",
            comments=""
        )
        print(f"  [NAMD:ECSH] Exported biexciton decay trace to: {biex_file}")

    # Flux averaging for PME in multi-origin mode
    if run_tr_sd and method in ("master_equation", "pme") and pme_flux_records and n_origins > 1:
        combined_flux = {}
        for r in pme_flux_records:
            t = r["time_fs"]
            if t not in combined_flux:
                combined_flux[t] = {
                    "flux_virt": np.zeros_like(r["flux_virt"]),
                    "flux_occ": np.zeros_like(r["flux_occ"]),
                }
            combined_flux[t]["flux_virt"] += r["flux_virt"] / n_origins
            combined_flux[t]["flux_occ"] += r["flux_occ"] / n_origins
        pme_flux_records = [
            {"time_fs": t, "flux_virt": v["flux_virt"], "flux_occ": v["flux_occ"]}
            for t, v in sorted(combined_flux.items())
        ]

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
        recombination_info=recombination_info,
        std_energies_ev=std_energies_ev,
        std_excess_e=std_excess_e,
        std_excess_h=std_excess_h,
        n_origins=n_origins,
    )

    # -------------------------------------------------------------

    # Transient Absorption (Pump-Probe) Spectroscopy Analysis
    # -------------------------------------------------------------
    ta_cfg = namd_cfg.get("transient_absorption", {})
    run_ta = bool(ta_cfg.get("run", False) or ta_cfg.get("enabled", False) or config.get("system", {}).get("namd_ta", False))
    if run_ta:
        from qdex.namd.transient_absorption import (
            compute_transient_absorption,
            plot_transient_absorption,
            export_transient_absorption_data,
        )
        print("\n" + "=" * 68)
        print("  [NAMD] Computing Ultrafast Pump-Probe Transient Absorption...")
        print("=" * 68)
        ta_sigma = float(ta_cfg.get("sigma", 0.03))
        ta_erange = ta_cfg.get("e_range", None)
        ta_n_e = int(ta_cfg.get("n_e_points", 300))
        ta_include_se = bool(ta_cfg.get("include_se", False))
        ta_degeneracy = 1.0 if soc_meta else float(ta_cfg.get("state_degeneracy", 2.0))

        ta_res = compute_transient_absorption(
            times_fs=times,
            populations=populations,
            E_pairs=E0_pairs,
            f_pairs=f0_pairs,
            i_pairs=i_pairs0,
            a_pairs=a_pairs0,
            sigma_ev=ta_sigma,
            e_range=ta_erange,
            n_e_points=ta_n_e,
            include_se=ta_include_se,
            state_degeneracy=ta_degeneracy,
            all_energies=all_energies if all_energies.shape[1] == populations.shape[1] else None,
        )

        fit = ta_res["fit_results"]
        print(f"  1S Band-Edge Energy          : {ta_res['e_1s_ev']:.3f} eV")
        if fit.get("success", False):
            print(f"  1S Bleach Rise Time (tau_C)  : {fit['tau_rise_fs']:.1f} fs ({fit['tau_rise_ps']:.3f} ps)")
            print(f"  Carrier Cooling Rate (k_C)   : {fit['k_cool_ps']:.2f} ps^-1")
        print("=" * 68 + "\n")

        # Plotting
        plot_ta = bool(ta_cfg.get("plot", True))
        if plot_ta:
            ta_plot_file = ta_cfg.get("plot_file", "transient_absorption_map.png")
            mat_name = config.get("system", {}).get("material", "CSPBBR3")
            plot_transient_absorption(ta_res, plot_file=ta_plot_file, material_name=mat_name)

        # Exporting
        csv_file = ta_cfg.get("csv_file", "ta_bleach_kinetics.csv")
        map_npz = ta_cfg.get("map_npz", "ta_2d_map.npz")
        export_transient_absorption_data(ta_res, kinetics_csv=csv_file, map_npz=map_npz)

    # -------------------------------------------------------------
    # Time-Resolved Vibrational Action Spectrum J(omega, t)
    # -------------------------------------------------------------
    if run_tr_sd:
        print("\n" + "=" * 68)
        print("  [NAMD] Computing Time-Resolved Vibrational Action Spectrum J(omega, t)...")
        print("=" * 68)
        from qdex.namd.analysis import (
            load_trajectory_orbital_energies,
            compute_time_resolved_spectral_density,
            plot_time_resolved_spectral_density,
            export_time_resolved_spectral_density,
            compute_band_gap_dynamics_and_spectral_density,
            compute_2d_vibronic_action_map,
            plot_2d_vibronic_action_map,
            export_2d_vibronic_action_map,
        )
        _, eps_o_traj, eps_v_traj = load_trajectory_orbital_energies(precompute_dir)
        bg_for_tr = compute_band_gap_dynamics_and_spectral_density(precompute_dir, use_lowest_exciton=True, temp_k=temp_k)

        w_max_cm = float(tr_sd_cfg.get("w_max_cm", 400.0))
        sigma_t = float(tr_sd_cfg.get("sigma_t_fs", 15.0))

        tr_sd_res = compute_time_resolved_spectral_density(
            times_fs=times,
            eps_occ_traj=eps_o_traj,
            eps_virt_traj=eps_v_traj,
            hop_records=hop_records if method not in ("master_equation", "pme") else None,
            pme_flux_records=pme_flux_records if method in ("master_equation", "pme") else None,
            n_trajectories=n_trajectories,
            sigma_t_fs=sigma_t,
            w_max_cm=w_max_cm,
            bg_data=bg_for_tr,
            method_name=method.upper(),
        )
        print(f"  Visited unique hopping/flux pairs: {tr_sd_res['n_unique_pairs']}")
        print(f"  Max action density amplitude    : {np.max(tr_sd_res['J_total_raw']):.4e}")

        plot_tr_sd = bool(tr_sd_cfg.get("plot", True))
        mat_name = config.get("system", {}).get("material", "CSPBBR3")
        if plot_tr_sd:
            plot_file = tr_sd_cfg.get("plot_file", "time_resolved_spectral_density.png")
            html_file = tr_sd_cfg.get("html_file", "time_resolved_spectral_density.html")
            plot_time_resolved_spectral_density(
                tr_sd_res,
                plot_file=plot_file,
                html_file=html_file,
                material_name=mat_name,
                method_name=method.upper(),
            )
        export_npz = tr_sd_cfg.get("npz_file", "time_resolved_spectral_density.npz")
        export_csv = tr_sd_cfg.get("csv_file", "time_resolved_spectral_density.csv")
        export_time_resolved_spectral_density(tr_sd_res, output_npz=export_npz, output_csv=export_csv)

        # 2D Non-Adiabatic Vibronic Action Map S(omega_acc, Omega_prom)
        compute_2d = bool(tr_sd_cfg.get("compute_2d_vibronic", True))
        if compute_2d:
            vib2d = compute_2d_vibronic_action_map(
                tr_sd_res,
                n_fft=int(tr_sd_cfg.get("n_fft_2d", 2048)),
                detrend_mode=str(tr_sd_cfg.get("detrend_2d", "linear")),
                wmax_acc=float(tr_sd_cfg.get("wmax_acc_2d", 200.0)),
                wmax_prom=float(tr_sd_cfg.get("wmax_prom_2d", 300.0)),
            )
            plot_2d = bool(tr_sd_cfg.get("plot_2d", True))
            if plot_2d:
                plot_file_2d = tr_sd_cfg.get("plot_file_2d", "2d_vibronic_action_map.png")
                html_file_2d = tr_sd_cfg.get("html_file_2d", "2d_vibronic_action_map.html")
                plot_2d_vibronic_action_map(
                    vib2d,
                    output_png=plot_file_2d,
                    output_html=html_file_2d,
                    material_name=mat_name,
                    method_name=method.upper()
                )
            npz_file_2d = tr_sd_cfg.get("npz_file_2d", "2d_vibronic_action_map.npz")
            csv_file_2d = tr_sd_cfg.get("csv_file_2d", "2d_vibronic_action_projections.csv")
            export_2d_vibronic_action_map(vib2d, output_npz=npz_file_2d, output_csv=csv_file_2d)

        print("=" * 68 + "\n")
