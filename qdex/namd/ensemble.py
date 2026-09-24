"""
Automated Multi-Origin Ensemble Sampling Module for QDEX NAMD.

Provides zero-configuration, physically grounded sampling of non-adiabatic
molecular dynamics trajectories from multiple statistically independent
initial nuclear configurations along an AIMD trajectory:
  1. Computes the phonon dephasing / correlation time (tau_corr) from band edge fluctuations.
  2. Runs a fast pilot calibration to determine the intrinsic cooling lifetime (tau_cool).
  3. Lays out an optimal simulation window T_window and a grid of origins separated by dt0 >= 2 * tau_corr.
  4. Dynamically references initial hot-carrier excitation energies to preserve excess energy Delta E_excess
     against thermal bandgap fluctuations.
  5. Aggregates multi-origin runs into grand ensemble means and thermal disorder confidence bands (sigma_thermal).
"""

import os
import glob
import time
import numpy as np

from qdex.namd.initial_condition import sample_initial_states
from qdex.namd.integrator import KB_EV
from qdex.namd.analysis import (
    load_trajectory_orbital_energies,
    fit_exponential_lifetime,
    natural_sort_key
)
from qdex.namd.master_equation import propagate_pme_tensor


def compute_energy_autocorrelation_time(times_fs, series):
    """
    Computes the normalized autocorrelation function (ACF) of an energy time series
    via zero-padded FFT and determines the 1/e dephasing/correlation time (tau_corr).

    Parameters:
      times_fs: 1D array of time stamps in fs
      series: 1D array of energy values (e.g. bandgap or orbital energy) in eV

    Returns:
      tau_corr_fs: float, correlation time in fs (fallback 200.0 fs if no decay)
    """
    series = np.asarray(series, dtype=np.float64)
    n = len(series)
    if n < 4:
        return 200.0

    dt_fs = float(times_fs[1] - times_fs[0]) if len(times_fs) > 1 else 1.0
    delta = series - np.mean(series)
    var = np.var(delta)
    if var <= 1e-14:
        return 200.0

    n_fft = 2 ** int(np.ceil(np.log2(2 * n - 1)))
    f_d = np.fft.rfft(delta, n=n_fft)
    acf = np.fft.irfft(f_d * np.conj(f_d), n=n_fft)[:n]
    normalization = np.arange(n, 0, -1, dtype=np.float64) * var
    acf /= np.maximum(normalization, 1e-14)

    # Find first crossing below 1/e
    decay_idx = np.where(acf <= (1.0 / np.e))[0]
    if len(decay_idx) > 0 and decay_idx[0] > 0:
        tau_corr = float(decay_idx[0] * dt_fs)
    else:
        # Fallback to half the trajectory length or 200.0 fs
        tau_corr = min(200.0, float(0.25 * times_fs[-1]))

    return max(tau_corr, 20.0)


def estimate_pilot_cooling_time(
    precompute_dir,
    qp_gap_ev,
    pump_energy_ev,
    dyn_cfg=None,
    temp_k=300.0,
    max_pilot_steps=120
):
    """
    Runs a fast pilot PME step (< 1-2 s) to measure the intrinsic carrier cooling
    lifetime (tau_cool), or retrieves user-specified tau_cool_fs.

    Returns:
      tau_cool_fs: float, estimated carrier cooling lifetime in fs
    """
    if dyn_cfg is None:
        dyn_cfg = {}

    if "tau_cool_fs" in dyn_cfg:
        return float(dyn_cfg["tau_cool_fs"])
    if "window_fs" in dyn_cfg:
        return float(dyn_cfg["window_fs"]) / 3.0

    frame0_path = os.path.join(precompute_dir, "frame_00000.npz")
    if not os.path.exists(frame0_path):
        return 1000.0

    step_files = sorted(glob.glob(os.path.join(precompute_dir, "step_*.npz")), key=natural_sort_key)
    if not step_files:
        return 1000.0

    try:
        f0 = np.load(frame0_path)
        E0 = f0["E_pairs"]
        i_p = f0["i_pairs"]
        a_p = f0["a_pairs"]
        n_occ = len(f0["eps_occ"])
        n_virt = len(f0["eps_virt"])

        # Active window restricted to cooling cascade
        mask = (E0 >= qp_gap_ev - 0.2) & (E0 <= pump_energy_ev + 0.3)
        if np.sum(mask) == 0:
            mask = np.ones(len(E0), dtype=bool)

        dyn_occ = np.sort(np.unique(i_p[mask]))
        dyn_virt = np.sort(np.unique(a_p[mask]))
        occ_to_sub = np.full(n_occ, -1, dtype=int)
        virt_to_sub = np.full(n_virt, -1, dtype=int)
        occ_to_sub[dyn_occ] = np.arange(len(dyn_occ))
        virt_to_sub[dyn_virt] = np.arange(len(dyn_virt))

        init_idx = np.argmin(np.abs(E0[mask] - pump_energy_ev))
        pair_sub = np.where(mask)[0][init_idx]
        i_init = occ_to_sub[i_p[pair_sub]]
        a_init = virt_to_sub[a_p[pair_sub]]

        P = np.zeros((len(dyn_occ), len(dyn_virt)), dtype=np.float64)
        if i_init >= 0 and a_init >= 0:
            P[i_init, a_init] = 1.0
        else:
            P[-1, -1] = 1.0

        n_test = min(len(step_files), max_pilot_steps)
        s0 = np.load(step_files[0])
        dt_fs = float(s0["time_curr_fs"] - s0["time_prev_fs"]) if "time_prev_fs" in s0 else 2.0

        mean_energies = []
        for k in range(n_test):
            s = np.load(step_files[k])
            d_occ = (s["S_occ"] - s["S_occ"].conj().T) / (2.0 * dt_fs)
            d_virt = (s["S_virt"] - s["S_virt"].conj().T) / (2.0 * dt_fs)
            d_occ_dyn = d_occ[np.ix_(dyn_occ, dyn_occ)]
            d_virt_dyn = d_virt[np.ix_(dyn_virt, dyn_virt)]

            eps_v = s["eps_virt_prev"][dyn_virt] if "eps_virt_prev" in s else s["eps_virt_curr"][dyn_virt]
            eps_o = s["eps_occ_prev"][dyn_occ] if "eps_occ_prev" in s else s["eps_occ_curr"][dyn_occ]
            E_mat = eps_v[np.newaxis, :] - eps_o[:, np.newaxis]

            mean_energies.append(np.sum(P * E_mat))
            P = propagate_pme_tensor(
                P, E_mat, d_occ_dyn, d_virt_dyn,
                dt_fs=dt_fs, temp_k=temp_k, n_substeps=2
            )

        mean_energies = np.array(mean_energies)
        excess = np.maximum(mean_energies - qp_gap_ev, 0.0)
        times_test = np.arange(len(mean_energies)) * dt_fs
        tau_pilot = fit_exponential_lifetime(times_test, excess)

        if np.isfinite(tau_pilot) and tau_pilot > 0.0:
            return float(tau_pilot)
        return 1000.0
    except Exception as exc:
        print(f"  [NAMD:Warn] Pilot cooling estimation encountered: {exc}; falling back to 1000.0 fs")
        return 1000.0


def auto_calibrate_ensemble_origins(
    precompute_dir,
    qp_gap_ev,
    pump_energy_ev,
    dyn_cfg=None,
    config=None
):
    """
    Automatically calibrates multi-origin ensemble sampling parameters:
      1. Determines tau_corr from band edge energy fluctuations.
      2. Determines tau_cool from fast pilot PME or config.
      3. Determines optimal simulation window T_window ~ 3 * tau_cool.
      4. Generates a grid of independent origins t_0,k separated by dt0 >= 2 * tau_corr.
      5. Gracefully falls back to a single origin if T_span is insufficient.

    Returns:
      calibration: dict containing:
        - tau_corr_fs: float
        - tau_cool_fs: float
        - dt0_fs: float
        - window_fs: float
        - n_win_steps: int
        - n_origins: int
        - origin_frames: 1D array of int
        - origin_times_fs: 1D array of float
        - t_md_total_fs: float
        - dt_nuc_fs: float
        - is_fallback: bool
    """
    if dyn_cfg is None:
        dyn_cfg = {}
    if config is None:
        config = {}

    times_fs, eps_occ_traj, eps_virt_traj = load_trajectory_orbital_energies(precompute_dir)
    n_frames = len(times_fs)
    dt_nuc_fs = float(times_fs[1] - times_fs[0]) if n_frames > 1 else 2.0
    t_md_total_fs = float(times_fs[-1])

    # 1. Phonon correlation / dephasing time
    gap_traj = eps_virt_traj[:, 0] - eps_occ_traj[:, -1]
    tau_corr_fs = compute_energy_autocorrelation_time(times_fs, gap_traj)

    # 2. Pilot cooling lifetime
    temp_k = float(dyn_cfg.get("temperature_k", 300.0))
    tau_cool_fs = estimate_pilot_cooling_time(
        precompute_dir=precompute_dir,
        qp_gap_ev=qp_gap_ev,
        pump_energy_ev=pump_energy_ev,
        dyn_cfg=dyn_cfg,
        temp_k=temp_k
    )

    # 3. Safe origin spacing
    user_dt0 = dyn_cfg.get("dt0_fs", None)
    if user_dt0 is not None:
        dt0_fs = max(float(user_dt0), dt_nuc_fs)
    else:
        dt0_fs = max(200.0, 2.0 * tau_corr_fs)

    # 4. Optimal simulation window
    user_win = dyn_cfg.get("window_fs", None)
    if user_win is not None:
        window_fs = min(t_md_total_fs, max(float(user_win), 2.0 * dt_nuc_fs))
    else:
        # Standard: 3.0 * tau_cool captures 95% of relaxation, capped at available MD
        rec_win = max(400.0, 3.0 * tau_cool_fs)
        window_fs = min(t_md_total_fs, rec_win)

    n_win_steps = max(1, int(np.round(window_fs / dt_nuc_fs)))
    window_fs = n_win_steps * dt_nuc_fs

    # 5. Feasible origin span
    t_span = t_md_total_fs - window_fs

    user_n_origins = dyn_cfg.get("n_origins", None)
    max_origins = int(dyn_cfg.get("max_origins", 25))

    if t_span < dt0_fs * 0.8:
        # MD trajectory too short to place multiple independent origins of length window_fs
        origin_frames = np.array([0], dtype=int)
        origin_times_fs = np.array([0.0], dtype=np.float64)
        n_origins = 1
        is_fallback = True
    else:
        if user_n_origins is not None:
            n_origins = max(1, min(int(user_n_origins), max_origins))
        else:
            calc_origins = int(np.floor(t_span / dt0_fs)) + 1
            n_origins = max(2, min(calc_origins, max_origins))

        step_frames = max(1, int(np.round(t_span / max(1, n_origins - 1) / dt_nuc_fs)))
        max_start_frame = max(0, n_frames - 1 - n_win_steps)

        origin_frames = []
        for m in range(n_origins):
            kf = min(m * step_frames, max_start_frame)
            if kf not in origin_frames:
                origin_frames.append(kf)
        origin_frames = np.array(sorted(origin_frames), dtype=int)
        origin_times_fs = origin_frames * dt_nuc_fs
        n_origins = len(origin_frames)
        is_fallback = (n_origins <= 1)

    return {
        "tau_corr_fs": tau_corr_fs,
        "tau_cool_fs": tau_cool_fs,
        "dt0_fs": dt0_fs,
        "window_fs": window_fs,
        "n_win_steps": n_win_steps,
        "n_origins": n_origins,
        "origin_frames": origin_frames,
        "origin_times_fs": origin_times_fs,
        "t_md_total_fs": t_md_total_fs,
        "dt_nuc_fs": dt_nuc_fs,
        "is_fallback": is_fallback,
    }


def load_origin_frame_data(precompute_dir, k0, frame0=None):
    """
    Loads initial pair energies, oscillator strengths, and orbital energies
    at frame index k0.

    Parameters:
      precompute_dir: path to precomputed directory
      k0: integer frame index
      frame0: optional pre-loaded frame_00000 dict

    Returns:
      data: dict with 'E_pairs', 'f_pairs', 'i_pairs', 'a_pairs', 'eps_occ', 'eps_virt', 'qp_gap'
    """
    if frame0 is None:
        frame0 = np.load(os.path.join(precompute_dir, "frame_00000.npz"))

    i_pairs = frame0["i_pairs"]
    a_pairs = frame0["a_pairs"]
    f_pairs = frame0["f_pairs"]

    if k0 == 0:
        E_pairs = frame0["E_pairs"]
        eps_occ = frame0["eps_occ"]
        eps_virt = frame0["eps_virt"]
    else:
        step_file = os.path.join(precompute_dir, f"step_{k0-1:05d}_to_{k0:05d}.npz")
        if os.path.exists(step_file):
            s = np.load(step_file)
            E_pairs = s["E_curr"]
            eps_occ = s["eps_occ_curr"]
            eps_virt = s["eps_virt_curr"]
            if "f_curr" in s:
                f_pairs = s["f_curr"]
        else:
            # Fallback to frame 0
            E_pairs = frame0["E_pairs"]
            eps_occ = frame0["eps_occ"]
            eps_virt = frame0["eps_virt"]

    qp_gap = float(eps_virt[0] - eps_occ[-1])
    return {
        "E_pairs": E_pairs,
        "f_pairs": f_pairs,
        "i_pairs": i_pairs,
        "a_pairs": a_pairs,
        "eps_occ": eps_occ,
        "eps_virt": eps_virt,
        "qp_gap": qp_gap
    }


def sample_origin_initial_states(
    precompute_dir,
    k0,
    pump_excess_ev,
    dyn_mask,
    dyn_cfg,
    n_trajectories,
    frame0=None
):
    """
    Samples initial states at frame k0 while preserving constant pump excess energy
    Delta E_excess above the instantaneous band edge E_g(t_0,k).

    Returns:
      sampled_states: 1D array of sampled pair indices
      P_init: 1D array of initial population distribution
      origin_data: dict of frame k0 energies
    """
    origin_data = load_origin_frame_data(precompute_dir, k0, frame0=frame0)
    E_pairs = origin_data["E_pairs"]
    f_pairs = origin_data["f_pairs"]
    qp_gap_k = origin_data["qp_gap"]

    # Referenced pump energy preserves excess kinetic energy
    pump_energy_k = qp_gap_k + max(0.0, pump_excess_ev)
    pulse_fwhm_ev = float(dyn_cfg.get("pulse_fwhm_ev", 0.08))
    filter_dark = bool(dyn_cfg.get("filter_dark_states", True))

    sampled_states, P_init = sample_initial_states(
        energies=E_pairs,
        f_osc=f_pairs,
        pump_energy_ev=pump_energy_k,
        pulse_fwhm_ev=pulse_fwhm_ev,
        n_trajectories=n_trajectories,
        filter_dark_states=filter_dark,
        mask=dyn_mask
    )

    origin_data["pump_energy_ev"] = pump_energy_k
    return sampled_states, P_init, origin_data


def aggregate_multi_origin_results(origin_results, times_fs):
    """
    Aggregates trajectories and observables across multiple origins,
    computing the grand ensemble mean and thermal variance confidence bands.

    Parameters:
      origin_results: list of dicts, each containing:
        - 'mean_energy': 1D array (n_frames)
        - 'mean_excess_e': 1D array (n_frames)
        - 'mean_excess_h': 1D array (n_frames)
        - 'populations': 2D array (n_frames, n_states)
        - 'trajectory_energies': optional 2D array (n_frames, n_traj)
        - 'all_energies': optional 2D array (n_frames, n_states)
      times_fs: 1D array of time points

    Returns:
      aggregated: dict containing:
        - 'mean_energy': grand mean energy vs time
        - 'std_energy': thermal standard deviation band vs time
        - 'mean_excess_e': grand mean electron excess energy vs time
        - 'std_excess_e': thermal standard deviation for electrons
        - 'mean_excess_h': grand mean hole excess energy vs time
        - 'std_excess_h': thermal standard deviation for holes
        - 'populations': grand mean populations vs time
        - 'trajectory_energies': concatenated trajectory matrix
        - 'all_energies': ensemble average orbital energies
        - 'n_origins': number of origins aggregated
    """
    n_origins = len(origin_results)
    if n_origins == 0:
        raise ValueError("origin_results list cannot be empty")

    n_frames = len(times_fs)
    energies_stack = np.zeros((n_origins, n_frames), dtype=np.float64)
    excess_e_stack = np.zeros((n_origins, n_frames), dtype=np.float64)
    excess_h_stack = np.zeros((n_origins, n_frames), dtype=np.float64)

    pop_shape = origin_results[0]["populations"].shape
    pop_stack = np.zeros((n_origins, pop_shape[0], pop_shape[1]), dtype=np.float64)

    traj_list = []
    all_e_list = []

    for idx, res in enumerate(origin_results):
        energies_stack[idx, :] = res["mean_energy"][:n_frames]
        excess_e_stack[idx, :] = res["mean_excess_e"][:n_frames]
        excess_h_stack[idx, :] = res["mean_excess_h"][:n_frames]
        pop_stack[idx, :, :] = res["populations"][:n_frames, :]

        if res.get("trajectory_energies") is not None and res["trajectory_energies"].size > 0:
            traj_list.append(res["trajectory_energies"][:n_frames, :])
        if res.get("all_energies") is not None and res["all_energies"].size > 0:
            all_e_list.append(res["all_energies"][:n_frames, :])

    mean_energy = np.mean(energies_stack, axis=0)
    std_energy = np.std(energies_stack, axis=0, ddof=1 if n_origins > 1 else 0)

    mean_excess_e = np.mean(excess_e_stack, axis=0)
    std_excess_e = np.std(excess_e_stack, axis=0, ddof=1 if n_origins > 1 else 0)

    mean_excess_h = np.mean(excess_h_stack, axis=0)
    std_excess_h = np.std(excess_h_stack, axis=0, ddof=1 if n_origins > 1 else 0)

    mean_pop = np.mean(pop_stack, axis=0)

    combined_traj = np.hstack(traj_list) if traj_list else None
    mean_all_e = np.mean(np.stack(all_e_list, axis=0), axis=0) if all_e_list else None

    return {
        "mean_energy": mean_energy,
        "std_energy": std_energy,
        "mean_excess_e": mean_excess_e,
        "std_excess_e": std_excess_e,
        "mean_excess_h": mean_excess_h,
        "std_excess_h": std_excess_h,
        "populations": mean_pop,
        "trajectory_energies": combined_traj,
        "all_energies": mean_all_e,
        "n_origins": n_origins
    }
