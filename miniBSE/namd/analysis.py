import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]


def fit_exponential_lifetime(times, values):
    """
    Fits y(t) ~ y(0) * exp(-t / tau) using linear regression on log(y).
    Returns tau in fs, or np.nan if unphysical/insufficient data.
    """
    try:
        valid = (times > 0) & (values > 1e-4)
        if np.sum(valid) >= 3:
            p = np.polyfit(times[valid], np.log(values[valid]), 1)
            if p[0] < -1e-6:
                return -1.0 / p[0]
    except Exception:
        pass
    return np.nan


def compute_nac_energy_gap_data(precompute_dir, max_steps=20, max_sample_pairs=5000):
    """
    Extracts non-adiabatic coupling magnitudes and energy differences:
      |d^virt_ab| vs |eps_b - eps_a|
      |d^occ_ij| vs |eps_j - eps_i|
    sampled across precomputed step files.
    """
    step_files = sorted(glob.glob(os.path.join(precompute_dir, "step_*.npz")), key=natural_sort_key)
    if not step_files:
        return None

    step_files = step_files[:max_steps]
    dE_virt_all = []
    nac_virt_all = []
    dE_occ_all = []
    nac_occ_all = []
    rms_nac_series = []

    for sfile in step_files:
        try:
            data = np.load(sfile)
            S_occ = data["S_occ"]
            S_virt = data["S_virt"]
            t_prev = float(data["time_prev_fs"])
            t_curr = float(data["time_curr_fs"])
            dt_fs = max(t_curr - t_prev, 1e-3)

            d_occ = np.abs((S_occ - S_occ.conj().T) / (2.0 * dt_fs))
            d_virt = np.abs((S_virt - S_virt.conj().T) / (2.0 * dt_fs))

            rms_occ = np.sqrt(np.mean(d_occ**2))
            rms_virt = np.sqrt(np.mean(d_virt**2))
            rms_nac_series.append(0.5 * (rms_occ + rms_virt))

            eps_virt = data["eps_virt_curr"] if "eps_virt_curr" in data else None
            eps_occ = data["eps_occ_curr"] if "eps_occ_curr" in data else None

            # Sample virtual transitions
            n_virt = S_virt.shape[0]
            if eps_virt is not None and n_virt > 1:
                idx_a, idx_b = np.triu_indices(n_virt, k=1)
                if len(idx_a) > max_sample_pairs:
                    sel = np.random.choice(len(idx_a), size=max_sample_pairs, replace=False)
                    idx_a = idx_a[sel]
                    idx_b = idx_b[sel]
                dE_v = np.abs(eps_virt[idx_b] - eps_virt[idx_a])
                nac_v = d_virt[idx_a, idx_b]
                dE_virt_all.append(dE_v)
                nac_virt_all.append(nac_v)

            # Sample occupied transitions
            n_occ = S_occ.shape[0]
            if eps_occ is not None and n_occ > 1:
                idx_i, idx_j = np.triu_indices(n_occ, k=1)
                if len(idx_i) > max_sample_pairs:
                    sel = np.random.choice(len(idx_i), size=max_sample_pairs, replace=False)
                    idx_i = idx_i[sel]
                    idx_j = idx_j[sel]
                dE_o = np.abs(eps_occ[idx_j] - eps_occ[idx_i])
                nac_o = d_occ[idx_i, idx_j]
                dE_occ_all.append(dE_o)
                nac_occ_all.append(nac_o)

        except Exception:
            continue

    res = {
        "rms_nac_series": np.array(rms_nac_series),
        "dE_virt": np.concatenate(dE_virt_all) if dE_virt_all else np.array([]),
        "nac_virt": np.concatenate(nac_virt_all) if nac_virt_all else np.array([]),
        "dE_occ": np.concatenate(dE_occ_all) if dE_occ_all else np.array([]),
        "nac_occ": np.concatenate(nac_occ_all) if nac_occ_all else np.array([]),
    }
    return res


def compute_band_gap_dynamics_and_spectral_density(precompute_dir, use_lowest_exciton=True):
    """
    Extracts band-edge orbital energies (HOMO and LUMO) or lowest excited state along the MD trajectory
    from precompute_dir, and computes:
      1. Band-gap fluctuation delta_Eg(t) = Eg(t) - <Eg>.
      2. Time-autocorrelation function C(tau) = <delta_Eg(t) delta_Eg(t+tau)> / sigma^2.
      3. Phonon spectral density J(omega) in cm^-1 via windowed FFT.
      4. Second-order cumulant expansion pure dephasing function:
           g(t) = (sigma^2 / hbar^2) * int_0^t dt1 int_0^t1 dt2 C(t2)
           D(t) = exp(-g(t))
      5. Pure electronic dephasing time tau_dec where D(tau_dec) = 1/e.
    """
    import glob
    from miniBSE.namd.integrator import HBAR_EV_FS

    step_files = sorted(glob.glob(os.path.join(precompute_dir, "step_*.npz")))
    lowest_exc_ens = []
    if not step_files:
        # Fallback to frame_*.npz if step files are not present
        frame_files = sorted(glob.glob(os.path.join(precompute_dir, "frame_*.npz")))
        if not frame_files:
            return None
        times = []
        homo_ens = []
        lumo_ens = []
        for ff in frame_files:
            d = np.load(ff)
            times.append(float(d["time_fs"]))
            homo_ens.append(float(d["eps_occ"][-1]))
            lumo_ens.append(float(d["eps_virt"][0]))
            if "lowest_exc" in d:
                lowest_exc_ens.append(float(d["lowest_exc"]))
    else:
        d0 = np.load(step_files[0])
        times = [float(d0["time_prev_fs"])]
        homo_ens = [float(d0["eps_occ_prev"][-1])]
        lumo_ens = [float(d0["eps_virt_prev"][0])]
        if "E_curr" in d0:
            lowest_exc_ens.append(float(np.min(d0["E_curr"])))

        for sf in step_files:
            d = np.load(sf)
            times.append(float(d["time_curr_fs"]))
            homo_ens.append(float(d["eps_occ_curr"][-1]))
            lumo_ens.append(float(d["eps_virt_curr"][0]))
            if "E_curr" in d:
                lowest_exc_ens.append(float(np.min(d["E_curr"])))

    times = np.array(times)
    if use_lowest_exciton and len(lowest_exc_ens) == len(times):
        gaps = np.array(lowest_exc_ens)
    else:
        gaps = np.array(lumo_ens) - np.array(homo_ens)
    N = len(gaps)
    if N < 3:
        return None

    dt_fs = times[1] - times[0] if N > 1 else 2.0
    delta_g = gaps - np.mean(gaps)
    var_g = np.var(delta_g)
    if var_g < 1e-18:
        return None

    # 1. Autocorrelation function C(tau)
    autocorr = np.zeros(N)
    for m in range(N):
        autocorr[m] = np.mean(delta_g[:N - m] * delta_g[m:]) / var_g
    tau_lags = np.arange(N) * dt_fs

    # 2. Windowed FFT for Phonon Spectral Density J(omega)
    window = np.hanning(N)
    n_fft = max(2048, N * 16)
    fft_vals = np.fft.rfft(autocorr * window, n=n_fft)
    freqs_fs = np.fft.rfftfreq(n_fft, d=dt_fs)
    c_cm_per_fs = 2.99792458e-5
    wavenumbers_cm = freqs_fs / c_cm_per_fs
    psd = np.abs(fft_vals) ** 2
    if np.max(psd) > 0:
        psd /= np.max(psd)

    # 3. Second-Order Cumulant Expansion for Pure Dephasing D(t)
    # g(t) = (var_g / hbar^2) * int_0^t dt1 int_0^t1 dt2 C(t2)
    I1 = np.zeros(N)
    for i in range(1, N):
        I1[i] = I1[i - 1] + 0.5 * (autocorr[i - 1] + autocorr[i]) * dt_fs

    g_t = np.zeros(N)
    coeff = var_g / (HBAR_EV_FS ** 2)
    for i in range(1, N):
        g_t[i] = g_t[i - 1] + coeff * 0.5 * (I1[i - 1] + I1[i]) * dt_fs

    D_t = np.exp(-g_t)

    # Dephasing time where D(t) drops to 1/e
    idx_dec = np.where(D_t <= 1.0 / np.e)[0]
    tau_dec_fs = float(times[idx_dec[0]]) if len(idx_dec) > 0 else np.nan

    return {
        "times": times,
        "gaps": gaps,
        "delta_g": delta_g,
        "var_g": var_g,
        "tau_lags": tau_lags,
        "autocorr": autocorr,
        "wavenumbers_cm": wavenumbers_cm,
        "psd": psd,
        "g_t": g_t,
        "D_t": D_t,
        "tau_dec_fs": tau_dec_fs,
    }


def compute_spectral_density(signal, dt_fs):
    """
    Computes normalized time-autocorrelation C(t) and power spectral density J(omega)
    in wavenumbers (cm^-1) using a Hann-windowed FFT.
    """
    N = len(signal)
    if N < 2:
        return None, None, None, None

    sig_fluc = signal - np.mean(signal)
    var = np.var(signal)
    if var < 1e-18:
        return None, None, None, None

    autocorr = np.zeros(N)
    for m in range(N):
        autocorr[m] = np.mean(sig_fluc[:N - m] * sig_fluc[m:]) / var

    tau_lags = np.arange(N) * dt_fs
    window = np.hanning(N)
    n_fft = max(512, N * 16)
    fft_vals = np.fft.rfft(autocorr * window, n=n_fft)
    freqs_fs = np.fft.rfftfreq(n_fft, d=dt_fs)
    c_cm_per_fs = 2.99792458e-5
    wavenumbers_cm = freqs_fs / c_cm_per_fs

    psd = np.abs(fft_vals) ** 2
    if np.max(psd) > 0:
        psd /= np.max(psd)

    return tau_lags, autocorr, wavenumbers_cm, psd


def analyze_and_plot_namd_results(
    times_fs,
    mean_energies_ev,
    populations,
    out_cfg,
    qp_gap_ev,
    pump_energy_ev,
    mean_excess_e=None,
    mean_excess_h=None,
    trajectory_energies=None,
    all_energies=None,
    precompute_dir=None,
    recombination_info=None
):
    """
    Comprehensive NAMD analysis module producing:
      1. Carrier cooling curve & trajectory traces (Exciton Cascade).
      2. Electron vs hole cooling decomposition.
      3. State-resolved 1Se (LUMO), 1Sh (HOMO), and 1S exciton populations vs time.
      4. Non-adiabatic coupling vs transition energy gap (Energy-Gap Law).
      5. Phonon spectral density J(omega) showing Pb-Br and cation modes.
      6. Electronic dephasing function D(t) from second-order cumulant expansion.
      7. CSV and NPZ exports.
    """
    csv_file = out_cfg.get("cooling_curve_csv", "carrier_cooling.csv")
    npz_file = out_cfg.get("populations_npz", "populations.npz")
    plot_enabled = out_cfg.get("plot_cooling", True)
    plot_file = out_cfg.get("plot_file", "carrier_cooling.png")

    excess_total = mean_energies_ev - qp_gap_ev
    if mean_excess_e is None:
        mean_excess_e = 0.5 * excess_total
    if mean_excess_h is None:
        mean_excess_h = 0.5 * excess_total

    # State-resolved population analysis: 1S_e (LUMO), 1S_h (HOMO), 1S exciton
    pop_se = None
    pop_sh = None
    pop_1s = None
    pop_hot_e = None
    pop_hot_h = None

    if precompute_dir and os.path.isdir(precompute_dir):
        frame0_path = os.path.join(precompute_dir, "frame_00000.npz")
        if os.path.exists(frame0_path):
            f0 = np.load(frame0_path)
            i_p = f0["i_pairs"]
            a_p = f0["a_pairs"]
            homo_i = np.max(i_p)
            lumo_a = 0

            mask_se = (a_p == lumo_a)
            mask_sh = (i_p == homo_i)
            mask_1s = (i_p == homo_i) & (a_p == lumo_a)

            if populations.shape[1] == len(i_p):
                pop_se = np.sum(populations[:, mask_se], axis=1)
                pop_hot_e = np.maximum(1.0 - pop_se, 0.0)
                pop_sh = np.sum(populations[:, mask_sh], axis=1)
                pop_hot_h = np.maximum(1.0 - pop_sh, 0.0)
                idx_1s = np.where(mask_1s)[0]
                pop_1s = populations[:, idx_1s[0]] if len(idx_1s) > 0 else np.zeros(len(times_fs))

    # 1. Export CSV
    csv_dict = {
        "Time_fs": times_fs,
        "E_exc_eV": mean_energies_ev,
        "Excess_Energy_eV": excess_total,
        "Excess_Electron_eV": mean_excess_e,
        "Excess_Hole_eV": mean_excess_h
    }
    if pop_se is not None:
        csv_dict["Pop_1S_e"] = pop_se
        csv_dict["Pop_Hot_e"] = pop_hot_e
        csv_dict["Pop_1S_h"] = pop_sh
        csv_dict["Pop_Hot_h"] = pop_hot_h
        csv_dict["Pop_1S_exciton"] = pop_1s

    df = pd.DataFrame(csv_dict)
    df.to_csv(csv_file, index=False)
    print(f"  [Output] Carrier cooling curve exported to: {csv_file}")

    # 2. Export Populations & Dynamics NPZ
    np.savez_compressed(
        npz_file,
        times_fs=times_fs,
        mean_energies_ev=mean_energies_ev,
        excess_total=excess_total,
        mean_excess_e=mean_excess_e,
        mean_excess_h=mean_excess_h,
        populations=populations,
        pop_se=pop_se if pop_se is not None else np.array([]),
        pop_sh=pop_sh if pop_sh is not None else np.array([]),
        pop_1s=pop_1s if pop_1s is not None else np.array([]),
        trajectory_energies=trajectory_energies if trajectory_energies is not None else np.array([]),
        qp_gap_ev=qp_gap_ev,
        pump_energy_ev=pump_energy_ev
    )
    print(f"  [Output] State populations saved to: {npz_file}")

    # 3. Fit Lifetimes
    tau_total = fit_exponential_lifetime(times_fs, excess_total)
    tau_e = fit_exponential_lifetime(times_fs, mean_excess_e)
    tau_h = fit_exponential_lifetime(times_fs, mean_excess_h)

    print("\n" + "=" * 65)
    print(" NAMD Carrier Cooling Summary")
    print("=" * 65)
    print(f"  Initial Energy <E(0)>     : {mean_energies_ev[0]:.4f} eV")
    print(f"  Final Energy   <E(end)>   : {mean_energies_ev[-1]:.4f} eV")
    print(f"  Band Edge Gap (Eg)        : {qp_gap_ev:.4f} eV")
    print(f"  Total Energy Dissipated   : {mean_energies_ev[0] - mean_energies_ev[-1]:.4f} eV")
    if np.isfinite(tau_total):
        print(f"  Exciton Cooling Lifetime  : {tau_total:.1f} fs")
    if np.isfinite(tau_e):
        print(f"  Electron Cooling Lifetime : {tau_e:.1f} fs")
    if np.isfinite(tau_h):
        print(f"  Hole Cooling Lifetime     : {tau_h:.1f} fs")
    print("=" * 65 + "\n")

    if recombination_info is not None:
        mat_name = recombination_info.get("material", "N/A")
        n_refr = recombination_info.get("refractive_index", 2.0)
        tau_rad_1 = recombination_info.get("tau_rad_1_ns", np.nan)
        tau_rad_therm = recombination_info.get("tau_rad_therm_ns", np.nan)
        tau_nr = recombination_info.get("tau_nr_ns", np.nan)
        plqy = recombination_info.get("plqy_percent", np.nan)
        k_rad_therm_s = recombination_info.get("k_rad_therm_s", 0.0)
        k_nr_s = recombination_info.get("k_nr_s", 0.0)

        print("=" * 65)
        print(" Recombination & Photoluminescence Summary")
        print("=" * 65)
        print(f"  Material                    : {mat_name} (Refractive Index n = {n_refr:.2f})")
        if np.isfinite(tau_rad_1):
            print(f"  Lowest Exciton Rad Lifetime : {tau_rad_1:.2f} ns (k_rad = {1e9/max(tau_rad_1, 1e-12):.2e} s^-1)")
        if np.isfinite(tau_rad_therm):
            print(f"  Thermalized Rad Lifetime    : {tau_rad_therm:.2f} ns (k_rad = {k_rad_therm_s:.2e} s^-1)")
        if np.isfinite(tau_nr):
            if tau_nr > 1e6:
                print(f"  Non-Radiative Lifetime      : > 1 ms (intrinsic limit, k_nr = {k_nr_s:.2e} s^-1)")
            else:
                print(f"  Non-Radiative Lifetime      : {tau_nr:.2f} ns (k_nr = {k_nr_s:.2e} s^-1)")
        if np.isfinite(plqy):
            print(f"  Predicted PL Quantum Yield  : {plqy:.1f} %")
        print("=" * 65 + "\n")

    # 4. Generate Multi-Panel Visualizations (6-panel publication layout)
    if plot_enabled:
        fig, axes = plt.subplots(2, 3, figsize=(20, 11))

        # -------------------------------------------------------------
        # Panel (a): Exciton Cascade & Trajectory Traces
        # -------------------------------------------------------------
        ax1 = axes[0, 0]
        if all_energies is not None and all_energies.shape[1] > 0:
            n_states = all_energies.shape[1]
            n_sub_states = min(150, n_states)
            e0_sorted_idx = np.argsort(all_energies[0, :])
            idx_sub = e0_sorted_idx[np.linspace(0, n_states - 1, n_sub_states, dtype=int)]
            for idx in idx_sub:
                ax1.plot(times_fs, all_energies[:, idx], color="#CBD5E1", lw=0.6, alpha=0.5)
        elif precompute_dir and os.path.isdir(precompute_dir):
            import glob
            f0_path = os.path.join(precompute_dir, "frame_00000.npz")
            step_files = sorted(glob.glob(os.path.join(precompute_dir, "step_*.npz")))
            if os.path.exists(f0_path):
                f0 = np.load(f0_path)
                E0 = f0["E_pairs"]
                n_pairs = len(E0)
                sorted_idx = np.argsort(E0)
                edge_idx = sorted_idx[:30]
                bulk_idx = sorted_idx[np.linspace(30, n_pairs - 1, 120, dtype=int)]
                tracked_idx = np.concatenate([edge_idx, bulk_idx])

                if step_files:
                    n_f = len(times_fs)
                    tracked_E = np.zeros((n_f, len(tracked_idx)))
                    tracked_E[0, :] = E0[tracked_idx]
                    for k, sf in enumerate(step_files[:n_f - 1]):
                        tracked_E[k + 1, :] = np.load(sf)["E_curr"][tracked_idx]
                    for s_i in range(len(tracked_idx)):
                        c = "#64748B" if s_i < 30 else "#94A3B8"
                        w = 0.9 if s_i < 30 else 0.5
                        al = 0.5 if s_i < 30 else 0.3
                        ax1.plot(times_fs, tracked_E[:, s_i], color=c, lw=w, alpha=al)
                else:
                    for idx in tracked_idx:
                        ax1.axhline(E0[idx], color="#CBD5E1", lw=0.6, alpha=0.5)

        if trajectory_energies is not None and trajectory_energies.size > 0:
            n_traj_plot = min(40, trajectory_energies.shape[1])
            for tr in range(n_traj_plot):
                ax1.plot(times_fs, trajectory_energies[:, tr], color="#F59E0B", lw=0.8, alpha=0.25)

        ax1.plot(times_fs, mean_energies_ev, color="#1D4ED8", lw=3.0, label=r"$\langle E_{\mathrm{exc}}(t) \rangle$")
        ax1.axhline(pump_energy_ev, color="#DC2626", linestyle="--", lw=1.8, label=f"Pump ({pump_energy_ev:.2f} eV)")
        ax1.axhline(qp_gap_ev, color="#111827", linestyle=":", lw=1.8, label=f"$E_g$ ({qp_gap_ev:.2f} eV)")

        ax1.set_xlabel("Time (fs)", fontsize=11, fontweight="bold")
        ax1.set_ylabel("Exciton Energy (eV)", fontsize=11, fontweight="bold")
        ax1.set_title("(a) Exciton State Cascade & Trajectories", fontsize=12, fontweight="bold")
        ax1.legend(loc="upper right", frameon=True, fontsize=9)
        ax1.grid(True, linestyle="--", alpha=0.4)

        # -------------------------------------------------------------
        # Panel (b): Electron vs Hole Excess Energy Relaxation
        # -------------------------------------------------------------
        ax2 = axes[0, 1]
        lbl_e = r"Electron $\Delta E_e(t)$" + (f" ($\\tau_e \\approx {tau_e:.1f}$ fs)" if np.isfinite(tau_e) else "")
        lbl_h = r"Hole $\Delta E_h(t)$" + (f" ($\\tau_h \\approx {tau_h:.1f}$ fs)" if np.isfinite(tau_h) else "")
        lbl_tot = r"Total $\Delta E_{\mathrm{exc}}(t)$" + (f" ($\\tau \\approx {tau_total:.1f}$ fs)" if np.isfinite(tau_total) else "")

        ax2.plot(times_fs, mean_excess_e, color="#2563EB", lw=2.4, label=lbl_e)
        ax2.plot(times_fs, mean_excess_h, color="#EA580C", lw=2.4, label=lbl_h)
        ax2.plot(times_fs, excess_total, color="#374151", lw=2.0, linestyle="--", label=lbl_tot)

        ax2.set_xlabel("Time (fs)", fontsize=11, fontweight="bold")
        ax2.set_ylabel("Excess Energy (eV)", fontsize=11, fontweight="bold")
        ax2.set_title("(b) Electron vs Hole Cooling Decomposition", fontsize=12, fontweight="bold")
        ax2.legend(loc="upper right", frameon=True, fontsize=9)
        ax2.grid(True, linestyle="--", alpha=0.4)

        # -------------------------------------------------------------
        # Panel (c): State-Resolved Populations vs Time (1S_e, 1S_h, 1S Exciton)
        # -------------------------------------------------------------
        ax3 = axes[0, 2]
        if pop_se is not None and pop_sh is not None:
            ax3.plot(times_fs, pop_se, color="#2563EB", lw=2.4, label=r"$1S_e$ (LUMO)")
            ax3.plot(times_fs, pop_sh, color="#EA580C", lw=2.4, label=r"$1S_h$ (HOMO)")
            ax3.plot(times_fs, pop_hot_e, color="#93C5FD", lw=1.8, linestyle="--", label=r"Hot $e^-$ ($> 1S_e$)")
            ax3.plot(times_fs, pop_hot_h, color="#FDBA74", lw=1.8, linestyle="--", label=r"Hot $h^+$ ($< 1S_h$)")
            if pop_1s is not None and np.max(pop_1s) > 1e-6:
                ax3.plot(times_fs, pop_1s, color="#059669", lw=2.0, label=r"$1S$ Exciton")
            ax3.set_xlabel("Time (fs)", fontsize=11, fontweight="bold")
            ax3.set_ylabel("Population", fontsize=11, fontweight="bold")
            ax3.set_title("(c) State-Resolved Populations vs Time", fontsize=12, fontweight="bold")
            ax3.legend(loc="center right", frameon=True, fontsize=9)
            ax3.grid(True, linestyle="--", alpha=0.4)
        else:
            ax3.text(0.5, 0.5, "Population decomposition\nrequires precompute pair mapping",
                     ha="center", va="center", transform=ax3.transAxes, color="gray", fontsize=11)
            ax3.set_title("(c) State Populations", fontsize=12, fontweight="bold")

        # -------------------------------------------------------------
        # Panel (d): Non-Adiabatic Coupling vs Transition Energy Gap
        # -------------------------------------------------------------
        ax4 = axes[1, 0]
        nac_data = None
        if precompute_dir and os.path.isdir(precompute_dir):
            nac_data = compute_nac_energy_gap_data(precompute_dir)

        if nac_data is not None and len(nac_data["dE_virt"]) > 0:
            dE_v = nac_data["dE_virt"]
            nac_v = nac_data["nac_virt"]
            dE_o = nac_data["dE_occ"]
            nac_o = nac_data["nac_occ"]

            ax4.scatter(dE_v, np.maximum(nac_v, 1e-6), color="#3B82F6", s=8, alpha=0.25, label="Virtual ($e^-$)")
            if len(dE_o) > 0:
                ax4.scatter(dE_o, np.maximum(nac_o, 1e-6), color="#EF4444", s=8, alpha=0.25, label="Occupied ($h^+$)")

            bins = np.linspace(0, min(5.0, np.max(dE_v)), 25)
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            idx_bin = np.digitize(dE_v, bins) - 1
            mean_trend = [np.mean(nac_v[idx_bin == b]) if np.sum(idx_bin == b) > 5 else np.nan for b in range(len(bin_centers))]
            ax4.plot(bin_centers, mean_trend, color="#1D4ED8", lw=2.5, label=r"Trend $\langle |d_{ab}| \rangle$")

            ax4.set_yscale("log")
            ax4.set_xlabel(r"Transition Energy Gap $\Delta E$ (eV)", fontsize=11, fontweight="bold")
            ax4.set_ylabel(r"Coupling Magnitude $|d_{IJ}|$ ($\mathrm{fs}^{-1}$)", fontsize=11, fontweight="bold")
            ax4.set_title(r"(d) NAC vs $\Delta E$ (Energy-Gap Law)", fontsize=12, fontweight="bold")
            ax4.legend(loc="upper right", frameon=True, fontsize=9)
            ax4.grid(True, linestyle="--", alpha=0.4)
        else:
            ax4.text(0.5, 0.5, "NAC vs Gap data not available",
                     ha="center", va="center", transform=ax4.transAxes, color="gray", fontsize=11)
            ax4.set_title(r"(d) NAC vs $\Delta E$", fontsize=12, fontweight="bold")

        # -------------------------------------------------------------
        # Panel (e): Phonon Spectral Density J(omega)
        # -------------------------------------------------------------
        ax5 = axes[1, 1]
        bg_data = None
        if precompute_dir and os.path.isdir(precompute_dir):
            bg_data = compute_band_gap_dynamics_and_spectral_density(precompute_dir)

        if bg_data is not None:
            wn = bg_data["wavenumbers_cm"]
            psd = bg_data["psd"]
            mask_wn = (wn >= 0) & (wn <= 400)
            ax5.plot(wn[mask_wn], psd[mask_wn], color="#059669", lw=2.4, label=r"Spectral Density $J(\omega)$")

            # Mark key perovskite phonon modes
            # 1. Pb-Br stretching mode (~155 cm^-1)
            idx_lo = np.where((wn >= 130) & (wn <= 175))[0]
            if len(idx_lo) > 0:
                p_lo = idx_lo[np.argmax(psd[idx_lo])]
                ax5.annotate(
                    f"Pb-Br LO stretch\n({wn[p_lo]:.0f} cm$^{{-1}}$)",
                    xy=(wn[p_lo], psd[p_lo]),
                    xytext=(wn[p_lo] + 25, psd[p_lo] * 0.85),
                    arrowprops=dict(facecolor="#047857", shrink=0.08, width=1.5, headwidth=6),
                    fontsize=9, fontweight="bold", color="#047857"
                )

            # 2. Cs cation / bending mode (~35-50 cm^-1)
            idx_cage = np.where((wn >= 20) & (wn <= 60))[0]
            if len(idx_cage) > 0:
                p_cage = idx_cage[np.argmax(psd[idx_cage])]
                if psd[p_cage] > 0.05:
                    ax5.annotate(
                        f"Cs / Pb-Br-Pb bend\n({wn[p_cage]:.0f} cm$^{{-1}}$)",
                        xy=(wn[p_cage], psd[p_cage]),
                        xytext=(wn[p_cage] + 20, min(1.0, psd[p_cage] + 0.15)),
                        arrowprops=dict(facecolor="#D97706", shrink=0.08, width=1.5, headwidth=6),
                        fontsize=8.5, fontweight="bold", color="#D97706"
                    )

            ax5.set_xlabel(r"Wavenumber ($\mathrm{cm}^{-1}$)", fontsize=11, fontweight="bold")
            ax5.set_ylabel(r"Spectral Density $J(\omega)$ (arb. units)", fontsize=11, fontweight="bold")
            ax5.set_title(r"(e) Phonon Spectral Density $J(\omega)$", fontsize=12, fontweight="bold")
            ax5.legend(loc="upper right", frameon=True, fontsize=9)
            ax5.grid(True, linestyle="--", alpha=0.4)
        else:
            ax5.text(0.5, 0.5, "Phonon Spectral Density\n(Requires precompute data)",
                     ha="center", va="center", transform=ax5.transAxes, color="gray", fontsize=11)
            ax5.set_title(r"(e) Phonon Spectral Density", fontsize=12, fontweight="bold")

        # -------------------------------------------------------------
        # Panel (f): Electronic Dephasing Function D(t) (Cumulant Expansion)
        # -------------------------------------------------------------
        ax6 = axes[1, 2]
        if bg_data is not None:
            t_plot = bg_data["times"]
            D_t = bg_data["D_t"]
            tau_dec = bg_data["tau_dec_fs"]

            ax6.plot(t_plot, D_t, color="#7C3AED", lw=2.4, label=r"Dephasing $D(t) = \exp(-g(t))$")
            ax6.axhline(1.0 / np.e, color="#DC2626", linestyle=":", lw=1.5, label=r"$1/e$ Threshold")
            if np.isfinite(tau_dec):
                ax6.axvline(tau_dec, color="#DC2626", linestyle="--", lw=1.5,
                            label=f"$\\tau_{{\\mathrm{{dec}}}} = {tau_dec:.1f}$ fs")

            ax6.set_xlim([0, min(100.0, np.max(t_plot))])
            ax6.set_ylim([-0.05, 1.05])
            ax6.set_xlabel("Time (fs)", fontsize=11, fontweight="bold")
            ax6.set_ylabel("Coherence Amplitude $D(t)$", fontsize=11, fontweight="bold")
            ax6.set_title(r"(f) Pure Dephasing $D(t)$ (Cumulant Expansion)", fontsize=12, fontweight="bold")
            ax6.legend(loc="upper right", frameon=True, fontsize=9)
            ax6.grid(True, linestyle="--", alpha=0.4)

            # Inset: Normalized autocorrelation C(tau)
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            ax_ins = inset_axes(ax6, width="40%", height="38%", loc="center right", borderpad=1.5)
            mask_tau = bg_data["tau_lags"] <= min(150.0, np.max(bg_data["tau_lags"]))
            ax_ins.plot(bg_data["tau_lags"][mask_tau], bg_data["autocorr"][mask_tau], color="#047857", lw=1.5)
            ax_ins.axhline(0, color="gray", linestyle=":", lw=0.8)
            ax_ins.set_xlabel(r"$\tau$ (fs)", fontsize=8)
            ax_ins.set_ylabel(r"$C(\tau)$", fontsize=8)
            ax_ins.tick_params(labelsize=7)
            ax_ins.grid(True, linestyle=":", alpha=0.4)
        else:
            ax6.text(0.5, 0.5, "Cumulant expansion\nrequires precompute data",
                     ha="center", va="center", transform=ax6.transAxes, color="gray", fontsize=11)
            ax6.set_title(r"(f) Electronic Dephasing", fontsize=12, fontweight="bold")

        fig.tight_layout()
        plt.savefig(plot_file, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  [Plot] Comprehensive carrier cooling figure saved to: {plot_file}")

        # 5. Interactive Plotly Dashboard
        html_file = out_cfg.get("html_file", plot_file.rsplit(".", 1)[0] + ".html")
        plot_interactive = out_cfg.get("plot_interactive", True)
        if plot_interactive:
            try:
                generate_interactive_plotly_dashboard(
                    times_fs=times_fs,
                    mean_energies_ev=mean_energies_ev,
                    mean_excess_e=mean_excess_e,
                    mean_excess_h=mean_excess_h,
                    excess_total=excess_total,
                    qp_gap_ev=qp_gap_ev,
                    pump_energy_ev=pump_energy_ev,
                    tau_total=tau_total,
                    tau_e=tau_e,
                    tau_h=tau_h,
                    pop_se=pop_se,
                    pop_sh=pop_sh,
                    pop_hot_e=pop_hot_e,
                    pop_hot_h=pop_hot_h,
                    pop_1s=pop_1s,
                    nac_data=nac_data,
                    bg_data=bg_data,
                    precompute_dir=precompute_dir,
                    html_file=html_file
                )
                print(f"  [Plotly] Interactive 6-panel dashboard saved to: {html_file}")
            except Exception as e:
                print(f"  [Plotly:Warn] Could not generate interactive HTML: {e}")


def generate_interactive_plotly_dashboard(
    times_fs,
    mean_energies_ev,
    mean_excess_e,
    mean_excess_h,
    excess_total,
    qp_gap_ev,
    pump_energy_ev,
    tau_total,
    tau_e,
    tau_h,
    pop_se=None,
    pop_sh=None,
    pop_hot_e=None,
    pop_hot_h=None,
    pop_1s=None,
    nac_data=None,
    bg_data=None,
    precompute_dir=None,
    html_file="carrier_cooling.html"
):
    """
    Builds a self-contained, interactive Plotly dashboard featuring:
      1. Exciton energy cascade and trajectory traces.
      2. Excess energy cooling decomposition (electron vs hole).
      3. State-resolved population dynamics (1Se, 1Sh, 1S exciton, Hot carriers).
      4. Non-adiabatic coupling vs transition energy gap (Energy-Gap Law).
      5. Phonon spectral density J(omega) showing Pb-Br and cation vibrational modes.
      6. Electronic pure dephasing D(t) from second-order cumulant expansion.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=(
            "(a) Exciton State Cascade & Trajectories",
            "(b) Electron vs Hole Cooling Decomposition",
            "(c) State-Resolved Populations vs Time",
            "(d) NAC vs &Delta;E (Energy-Gap Law)",
            "(e) Phonon Spectral Density J(&omega;)",
            "(f) Pure Dephasing D(t) (Cumulant Expansion)"
        ),
        horizontal_spacing=0.08,
        vertical_spacing=0.12
    )

    # 1. Panel (a): Exciton Cascade & Fluctuating/Intersecting States
    if precompute_dir and os.path.isdir(precompute_dir):
        import glob
        f0_path = os.path.join(precompute_dir, "frame_00000.npz")
        step_files = sorted(glob.glob(os.path.join(precompute_dir, "step_*.npz")))
        if os.path.exists(f0_path):
            f0 = np.load(f0_path)
            E0 = f0["E_pairs"]
            n_pairs = len(E0)
            sorted_idx = np.argsort(E0)
            edge_idx = sorted_idx[:30]
            bulk_idx = sorted_idx[np.linspace(30, n_pairs - 1, 120, dtype=int)]
            tracked_idx = np.concatenate([edge_idx, bulk_idx])

            if step_files:
                n_f = len(times_fs)
                tracked_E = np.zeros((n_f, len(tracked_idx)))
                tracked_E[0, :] = E0[tracked_idx]
                for k, sf in enumerate(step_files[:n_f - 1]):
                    tracked_E[k + 1, :] = np.load(sf)["E_curr"][tracked_idx]

                for s_i in range(len(tracked_idx)):
                    state_orig_idx = tracked_idx[s_i]
                    is_edge = s_i < 30
                    line_color = "rgba(100, 116, 139, 0.45)" if is_edge else "rgba(148, 163, 184, 0.25)"
                    line_w = 1.2 if is_edge else 0.8
                    fig.add_trace(go.Scatter(
                        x=times_fs, y=tracked_E[:, s_i],
                        mode="lines",
                        line=dict(color=line_color, width=line_w),
                        hoverinfo="text",
                        text=[f"State #{state_orig_idx}<br>t = {t:.1f} fs<br>E = {e:.3f} eV" for t, e in zip(times_fs, tracked_E[:, s_i])],
                        showlegend=False
                    ), row=1, col=1)
            else:
                for e_val in E0[tracked_idx]:
                    fig.add_trace(go.Scatter(
                        x=[times_fs[0], times_fs[-1]], y=[e_val, e_val],
                        mode="lines", line=dict(color="#CBD5E1", width=0.8),
                        hoverinfo="skip", showlegend=False
                    ), row=1, col=1)


    fig.add_trace(go.Scatter(
        x=times_fs, y=mean_energies_ev,
        mode="lines", name="<E_exc(t)>",
        line=dict(color="#1D4ED8", width=3.5),
        hovertemplate="Time: %{x:.1f} fs<br>Energy: %{y:.4f} eV<extra></extra>"
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=[times_fs[0], times_fs[-1]], y=[pump_energy_ev, pump_energy_ev],
        mode="lines", name=f"Pump ({pump_energy_ev:.2f} eV)",
        line=dict(color="#DC2626", width=2, dash="dash"),
        hovertemplate="Pump: %{y:.3f} eV<extra></extra>"
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=[times_fs[0], times_fs[-1]], y=[qp_gap_ev, qp_gap_ev],
        mode="lines", name=f"Eg ({qp_gap_ev:.2f} eV)",
        line=dict(color="#111827", width=2, dash="dot"),
        hovertemplate="Eg: %{y:.3f} eV<extra></extra>"
    ), row=1, col=1)

    # 2. Panel (b): Electron vs Hole Cooling
    lbl_e = f"Electron &Delta;E_e (&tau;={tau_e:.1f} fs)" if np.isfinite(tau_e) else "Electron &Delta;E_e"
    lbl_h = f"Hole &Delta;E_h (&tau;={tau_h:.1f} fs)" if np.isfinite(tau_h) else "Hole &Delta;E_h"
    lbl_tot = f"Total &Delta;E_exc (&tau;={tau_total:.1f} fs)" if np.isfinite(tau_total) else "Total &Delta;E_exc"

    fig.add_trace(go.Scatter(
        x=times_fs, y=mean_excess_e,
        mode="lines", name=lbl_e,
        line=dict(color="#2563EB", width=2.8),
        hovertemplate="Time: %{x:.1f} fs<br>&Delta;E_e: %{y:.4f} eV<extra></extra>"
    ), row=1, col=2)

    fig.add_trace(go.Scatter(
        x=times_fs, y=mean_excess_h,
        mode="lines", name=lbl_h,
        line=dict(color="#EA580C", width=2.8),
        hovertemplate="Time: %{x:.1f} fs<br>&Delta;E_h: %{y:.4f} eV<extra></extra>"
    ), row=1, col=2)

    fig.add_trace(go.Scatter(
        x=times_fs, y=excess_total,
        mode="lines", name=lbl_tot,
        line=dict(color="#374151", width=2.2, dash="dash"),
        hovertemplate="Time: %{x:.1f} fs<br>&Delta;E_total: %{y:.4f} eV<extra></extra>"
    ), row=1, col=2)

    # 3. Panel (c): State-Resolved Populations
    if pop_se is not None and pop_sh is not None:
        fig.add_trace(go.Scatter(
            x=times_fs, y=pop_se,
            mode="lines", name="1S_e (LUMO)",
            line=dict(color="#2563EB", width=2.8),
            hovertemplate="Time: %{x:.1f} fs<br>P(1S_e): %{y:.5f}<extra></extra>"
        ), row=1, col=3)

        fig.add_trace(go.Scatter(
            x=times_fs, y=pop_sh,
            mode="lines", name="1S_h (HOMO)",
            line=dict(color="#EA580C", width=2.8),
            hovertemplate="Time: %{x:.1f} fs<br>P(1S_h): %{y:.5f}<extra></extra>"
        ), row=1, col=3)

        fig.add_trace(go.Scatter(
            x=times_fs, y=pop_hot_e,
            mode="lines", name="Hot e- (> 1S_e)",
            line=dict(color="#93C5FD", width=2, dash="dash"),
            hovertemplate="Time: %{x:.1f} fs<br>P(Hot e-): %{y:.5f}<extra></extra>"
        ), row=1, col=3)

        fig.add_trace(go.Scatter(
            x=times_fs, y=pop_hot_h,
            mode="lines", name="Hot h+ (< 1S_h)",
            line=dict(color="#FDBA74", width=2, dash="dash"),
            hovertemplate="Time: %{x:.1f} fs<br>P(Hot h+): %{y:.5f}<extra></extra>"
        ), row=1, col=3)

        if pop_1s is not None and np.max(pop_1s) > 1e-8:
            fig.add_trace(go.Scatter(
                x=times_fs, y=pop_1s,
                mode="lines", name="1S Exciton",
                line=dict(color="#059669", width=2.2),
                hovertemplate="Time: %{x:.1f} fs<br>P(1S): %{y:.6e}<extra></extra>"
            ), row=1, col=3)

    # 4. Panel (d): NAC vs dE (Energy-Gap Law)
    if nac_data is not None and len(nac_data["dE_virt"]) > 0:
        n_v = len(nac_data["dE_virt"])
        sel_v = np.random.choice(n_v, size=min(2500, n_v), replace=False)
        fig.add_trace(go.Scatter(
            x=nac_data["dE_virt"][sel_v], y=np.maximum(nac_data["nac_virt"][sel_v], 1e-6),
            mode="markers", name="Virtual (e-)",
            marker=dict(color="#3B82F6", size=4, opacity=0.35),
            hovertemplate="&Delta;E: %{x:.3f} eV<br>|d_ab|: %{y:.4e} fs^-1<extra></extra>"
        ), row=2, col=1)

        n_o = len(nac_data["dE_occ"])
        if n_o > 0:
            sel_o = np.random.choice(n_o, size=min(2500, n_o), replace=False)
            fig.add_trace(go.Scatter(
                x=nac_data["dE_occ"][sel_o], y=np.maximum(nac_data["nac_occ"][sel_o], 1e-6),
                mode="markers", name="Occupied (h+)",
                marker=dict(color="#EF4444", size=4, opacity=0.35),
                hovertemplate="&Delta;E: %{x:.3f} eV<br>|d_ij|: %{y:.4e} fs^-1<extra></extra>"
            ), row=2, col=1)

        bins = np.linspace(0, min(5.0, np.max(nac_data["dE_virt"])), 25)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        idx_bin = np.digitize(nac_data["dE_virt"], bins) - 1
        mean_trend = [np.mean(nac_data["nac_virt"][idx_bin == b]) if np.sum(idx_bin == b) > 5 else np.nan for b in range(len(bin_centers))]
        fig.add_trace(go.Scatter(
            x=bin_centers, y=mean_trend,
            mode="lines", name="Trend <|d_ab|>",
            line=dict(color="#1D4ED8", width=3),
            hovertemplate="&Delta;E: %{x:.2f} eV<br>Avg |d|: %{y:.4e} fs^-1<extra></extra>"
        ), row=2, col=1)

    # 5. Panel (e): Phonon Spectral Density
    if bg_data is not None:
        wn = bg_data["wavenumbers_cm"]
        psd = bg_data["psd"]
        mask_wn = (wn >= 0) & (wn <= 400)
        fig.add_trace(go.Scatter(
            x=wn[mask_wn], y=psd[mask_wn],
            mode="lines", name="Spectral Density J(&omega;)",
            line=dict(color="#059669", width=2.8),
            hovertemplate="Wavenumber: %{x:.1f} cm^-1<br>J(&omega;): %{y:.4f}<extra></extra>"
        ), row=2, col=2)

        idx_lo = np.where((wn >= 130) & (wn <= 175))[0]
        if len(idx_lo) > 0:
            p_lo = idx_lo[np.argmax(psd[idx_lo])]
            fig.add_annotation(
                x=wn[p_lo], y=psd[p_lo],
                text=f"Pb-Br LO stretch ({wn[p_lo]:.0f} cm<sup>-1</sup>)",
                showarrow=True, arrowhead=2, arrowcolor="#047857",
                ax=30, ay=-35, font=dict(color="#047857", size=11, family="sans-serif"),
                row=2, col=2
            )

    # 6. Panel (f): Pure Dephasing D(t)
    if bg_data is not None:
        t_plot = bg_data["times"]
        D_t = bg_data["D_t"]
        tau_dec = bg_data["tau_dec_fs"]
        fig.add_trace(go.Scatter(
            x=t_plot, y=D_t,
            mode="lines", name="Dephasing D(t)",
            line=dict(color="#7C3AED", width=2.8),
            hovertemplate="Time: %{x:.1f} fs<br>D(t): %{y:.4f}<extra></extra>"
        ), row=2, col=3)

        fig.add_trace(go.Scatter(
            x=[0, min(100.0, np.max(t_plot))], y=[1.0 / np.e, 1.0 / np.e],
            mode="lines", name="1/e Threshold",
            line=dict(color="#DC2626", width=1.8, dash="dot"),
            hovertemplate="1/e = 0.368<extra></extra>"
        ), row=2, col=3)

        if np.isfinite(tau_dec):
            fig.add_trace(go.Scatter(
                x=[tau_dec, tau_dec], y=[0, 1.0],
                mode="lines", name=f"&tau;_dec = {tau_dec:.1f} fs",
                line=dict(color="#DC2626", width=1.8, dash="dash"),
                hovertemplate=f"&tau;_dec = {tau_dec:.1f} fs<extra></extra>"
            ), row=2, col=3)

    # Axis formatting
    fig.update_xaxes(title_text="Time (fs)", row=1, col=1)
    fig.update_yaxes(title_text="Exciton Energy (eV)", row=1, col=1)

    fig.update_xaxes(title_text="Time (fs)", row=1, col=2)
    fig.update_yaxes(title_text="Excess Energy (eV)", row=1, col=2)

    fig.update_xaxes(title_text="Time (fs)", row=1, col=3)
    fig.update_yaxes(title_text="Population", row=1, col=3)

    fig.update_xaxes(title_text="Transition Energy Gap &Delta;E (eV)", row=2, col=1)
    fig.update_yaxes(title_text="Coupling Magnitude |d_IJ| (fs^-1)", type="log", row=2, col=1)

    fig.update_xaxes(title_text="Wavenumber (cm^-1)", range=[0, 400], row=2, col=2)
    fig.update_yaxes(title_text="Spectral Density (arb. units)", row=2, col=2)

    fig.update_xaxes(title_text="Time (fs)", range=[0, 100], row=2, col=3)
    fig.update_yaxes(title_text="Coherence Amplitude D(t)", range=[-0.05, 1.05], row=2, col=3)

    fig.update_layout(
        height=880,
        width=1500,
        title_text="<b>miniBSE — NAMD Hot Carrier Cooling & Exciton Relaxation (CsPbBr3)</b>",
        title_font=dict(size=18, family="sans-serif"),
        template="plotly_white",
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=-0.18, xanchor="center", x=0.5, font=dict(size=10))
    )

    fig.write_html(html_file, include_plotlyjs=True)


