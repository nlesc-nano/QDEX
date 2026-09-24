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


def compute_band_edge_arrival_times(times, values, tau, temp_k=300.0):
    """
    Computes analytical estimates and actual trajectory times to reach the band edge:
      - Cooling rate k_cool = 1 / tau (in ps^-1).
      - Estimated time to dissipate 95% of excess energy: ~ 3.0 * tau (-ln(0.05) * tau).
      - Estimated time to dissipate 99% of excess energy: ~ 4.6 * tau (-ln(0.01) * tau).
      - Actual trajectory time to dissipate 95% of initial excess: first t where Delta E(t) <= 0.05 * Delta E(0).
      - Actual trajectory time to dissipate 99% of initial excess: first t where Delta E(t) <= 0.01 * Delta E(0).
      - Actual trajectory time to reach thermal energy threshold: first t where Delta E(t) <= k_B * T.
    """
    from qdex.namd.integrator import KB_EV

    times_arr = np.asarray(times, dtype=np.float64)
    val_arr = np.asarray(values, dtype=np.float64)
    v0 = float(val_arr[0]) if len(val_arr) > 0 else 0.0
    t_end = float(times_arr[-1]) if len(times_arr) > 0 else 0.0

    # 1. Cooling rate k_cool in ps^-1
    k_cool_ps = (1e3 / tau) if (np.isfinite(tau) and tau > 0) else np.nan

    # 2. Analytical estimates from fitted tau
    if np.isfinite(tau) and tau > 0:
        est_95 = 2.99573227 * tau  # -ln(0.05) * tau ≈ 3.0 * tau
        est_99 = 4.60517019 * tau  # -ln(0.01) * tau ≈ 4.6 * tau
    else:
        est_95 = np.nan
        est_99 = np.nan

    # 3. Actual numerical times from simulation trajectory
    # 95% dissipated (excess <= 5% of v0)
    idx_95 = np.where(val_arr <= 0.05 * v0)[0]
    act_95 = float(times_arr[idx_95[0]]) if len(idx_95) > 0 else np.nan

    # 99% dissipated (excess <= 1% of v0)
    idx_99 = np.where(val_arr <= 0.01 * v0)[0]
    act_99 = float(times_arr[idx_99[0]]) if len(idx_99) > 0 else np.nan

    # Thermalized: excess <= k_B * T
    thermal_thresh = KB_EV * max(float(temp_k), 1.0)
    idx_therm = np.where(val_arr <= thermal_thresh)[0]
    act_therm = float(times_arr[idx_therm[0]]) if len(idx_therm) > 0 else np.nan

    return {
        "tau_fs": tau,
        "k_cool_ps": k_cool_ps,
        "est_95_fs": est_95,
        "est_99_fs": est_99,
        "act_95_fs": act_95,
        "act_99_fs": act_99,
        "act_therm_fs": act_therm,
        "t_end_fs": t_end,
    }


def format_time_fs(val, t_max=None):
    """Formats a time in femtoseconds cleanly, handling finite, NaN, and window exceedance."""
    if np.isfinite(val):
        return f"{val:.1f} fs"
    elif t_max is not None:
        return f"> {t_max:.1f} fs (exceeds window)"
    return "N/A"


def compute_nac_energy_gap_data(precompute_dir, max_steps=20, max_sample_pairs=5000):
    """
    Extracts non-adiabatic coupling magnitudes and energy differences::

      |d_ab^virt| vs |eps_b - eps_a|
      |d_ij^occ| vs |eps_j - eps_i|

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


def compute_band_gap_dynamics_and_spectral_density(precompute_dir, use_lowest_exciton=True, temp_k=300.0):
    r"""
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
    from qdex.namd.integrator import HBAR_EV_FS

    step_files = sorted(glob.glob(os.path.join(precompute_dir, "step_*.npz")), key=natural_sort_key)
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
        # E_curr belongs at time_curr. The t = 0 exciton is the frame-0 pair list.
        frame0_path = os.path.join(precompute_dir, "frame_00000.npz")
        if os.path.exists(frame0_path):
            f0 = np.load(frame0_path)
            if "E_pairs" in f0.files and len(f0["E_pairs"]) > 0:
                lowest_exc_ens.append(float(np.min(f0["E_pairs"])))
            else:
                lowest_exc_ens.append(float(lumo_ens[0] - homo_ens[0]))
        else:
            lowest_exc_ens.append(float(lumo_ens[0] - homo_ens[0]))

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
    # J(omega) is the cosine transform of C(t), not the power of that transform.
    psd = np.maximum(np.real(fft_vals), 0.0)
    if np.max(psd) > 0:
        psd /= np.max(psd)
    rayleigh_cm = (1.0 / max((N - 1) * dt_fs, 1e-6)) / c_cm_per_fs

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

    # Extract dominant optical phonon frequency from PSD (excluding DC component < 30 cm^-1)
    mask_lo = (wavenumbers_cm >= 30.0) & (wavenumbers_cm <= 600.0)
    if np.any(mask_lo):
        sub_psd = psd[mask_lo]
        sub_wn = wavenumbers_cm[mask_lo]
        peak_idx = np.argmax(sub_psd)
        dominant_freq_cm1 = float(sub_wn[peak_idx])
    else:
        dominant_freq_cm1 = 150.0

    from qdex.hardness import extract_recombination_parameters_from_namd
    recomb_params = extract_recombination_parameters_from_namd(
        var_E_gap_ev2=var_g,
        dominant_freq_cm1=dominant_freq_cm1,
        temp_k=temp_k
    )

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
        "dominant_freq_cm1": dominant_freq_cm1,
        "rayleigh_cm": float(rayleigh_cm),
        "recomb_params": recomb_params,
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

    psd = np.maximum(np.real(fft_vals), 0.0)
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
    recombination_info=None,
    std_energies_ev=None,
    std_excess_e=None,
    std_excess_h=None,
    n_origins=1
):
    """
    Comprehensive NAMD analysis module producing:
      1. Carrier cooling curve & trajectory traces (Exciton Cascade).
      2. Electron vs hole cooling decomposition with thermal confidence intervals.
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
    if std_energies_ev is not None:
        csv_dict["E_exc_std_eV"] = std_energies_ev
    if std_excess_e is not None:
        csv_dict["Excess_Electron_std_eV"] = std_excess_e
    if std_excess_h is not None:
        csv_dict["Excess_Hole_std_eV"] = std_excess_h

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
        std_energies_ev=std_energies_ev if std_energies_ev is not None else np.array([]),
        std_excess_e=std_excess_e if std_excess_e is not None else np.array([]),
        std_excess_h=std_excess_h if std_excess_h is not None else np.array([]),
        n_origins=int(n_origins),
        populations=populations,
        pop_se=pop_se if pop_se is not None else np.array([]),
        pop_sh=pop_sh if pop_sh is not None else np.array([]),
        pop_1s=pop_1s if pop_1s is not None else np.array([]),
        trajectory_energies=trajectory_energies if trajectory_energies is not None else np.array([]),
        qp_gap_ev=qp_gap_ev,
        pump_energy_ev=pump_energy_ev
    )
    print(f"  [Output] State populations saved to: {npz_file}")

    # 3. Fit Lifetimes & Band Edge Arrival Times
    tau_total = fit_exponential_lifetime(times_fs, excess_total)
    tau_e = fit_exponential_lifetime(times_fs, mean_excess_e)
    tau_h = fit_exponential_lifetime(times_fs, mean_excess_h)

    arr_tot = compute_band_edge_arrival_times(times_fs, excess_total, tau_total)
    arr_e = compute_band_edge_arrival_times(times_fs, mean_excess_e, tau_e)
    arr_h = compute_band_edge_arrival_times(times_fs, mean_excess_h, tau_h)

    print("\n" + "=" * 68)
    print(" NAMD Carrier Cooling & Relaxation Summary")
    print("=" * 68)
    print(f"  Initial Energy <E(0)>        : {mean_energies_ev[0]:.4f} eV")
    print(f"  Final Energy   <E(end)>      : {mean_energies_ev[-1]:.4f} eV")
    print(f"  Band Edge Gap (Eg)           : {qp_gap_ev:.4f} eV")
    print(f"  Total Energy Dissipated      : {mean_energies_ev[0] - mean_energies_ev[-1]:.4f} eV")
    print()
    print("  --- Exponential Cooling Lifetimes & Rates (k_cool = 1/tau) ---")
    if np.isfinite(tau_total):
        print(f"  Exciton Lifetime (tau)       : {tau_total:.1f} fs (k_cool = {arr_tot['k_cool_ps']:.2f} ps^-1)")
    if np.isfinite(tau_e):
        print(f"  Electron Lifetime (tau_e)    : {tau_e:.1f} fs (k_cool = {arr_e['k_cool_ps']:.2f} ps^-1)")
    if np.isfinite(tau_h):
        print(f"  Hole Lifetime (tau_h)        : {tau_h:.1f} fs (k_cool = {arr_h['k_cool_ps']:.2f} ps^-1)")
    print()
    print("  --- Band Edge Arrival Time (95% excess dissipated, ~3.0*tau) ---")
    print(f"  Exciton                      : Estimated = {format_time_fs(arr_tot['est_95_fs']):<10} | Actual = {format_time_fs(arr_tot['act_95_fs'], times_fs[-1])}")
    print(f"  Electron                     : Estimated = {format_time_fs(arr_e['est_95_fs']):<10} | Actual = {format_time_fs(arr_e['act_95_fs'], times_fs[-1])}")
    print(f"  Hole                         : Estimated = {format_time_fs(arr_h['est_95_fs']):<10} | Actual = {format_time_fs(arr_h['act_95_fs'], times_fs[-1])}")
    print()
    print("  --- Complete Thermalization Time (99% excess dissipated, ~4.6*tau) ---")
    print(f"  Exciton                      : Estimated = {format_time_fs(arr_tot['est_99_fs']):<10} | Actual = {format_time_fs(arr_tot['act_99_fs'], times_fs[-1])}")
    print(f"  Electron                     : Estimated = {format_time_fs(arr_e['est_99_fs']):<10} | Actual = {format_time_fs(arr_e['act_99_fs'], times_fs[-1])}")
    print(f"  Hole                         : Estimated = {format_time_fs(arr_h['est_99_fs']):<10} | Actual = {format_time_fs(arr_h['act_99_fs'], times_fs[-1])}")
    print()
    print("  --- Lattice Thermalization Window (Excess <= k_B*T ≈ 25.8 meV) ---")
    print(f"  Exciton                      : Actual = {format_time_fs(arr_tot['act_therm_fs'], times_fs[-1])}")
    print(f"  Electron                     : Actual = {format_time_fs(arr_e['act_therm_fs'], times_fs[-1])}")
    print(f"  Hole                         : Actual = {format_time_fs(arr_h['act_therm_fs'], times_fs[-1])}")
    if n_origins > 1 and std_energies_ev is not None and len(std_energies_ev) > 0:
        print()
        print(f"  --- Multi-Origin Thermal Ensemble ({n_origins} AIMD Initial Conditions) ---")
        mean_spd = np.mean(std_energies_ev) * 1e3
        max_spd = np.max(std_energies_ev) * 1e3
        print(f"  Mean Thermal Spread (sigma)  : {mean_spd:.2f} meV (max = {max_spd:.2f} meV)")
    print("=" * 68 + "\n")

    bg_data = None
    nac_data = None
    mean_nac_fs = None
    if precompute_dir and os.path.isdir(precompute_dir):
        bg_data = compute_band_gap_dynamics_and_spectral_density(precompute_dir)
        nac_data = compute_nac_energy_gap_data(precompute_dir)
        if nac_data is not None and len(nac_data.get("rms_nac_series", [])) > 0:
            mean_nac_fs = float(np.mean(nac_data["rms_nac_series"]))

    if recombination_info is not None:
        mat_name = recombination_info.get("material", "N/A")
        n_refr = recombination_info.get("refractive_index", 2.0)
        tau_rad_1 = recombination_info.get("tau_rad_1_ns", np.nan)
        tau_rad_therm = recombination_info.get("tau_rad_therm_ns", np.nan)
        tau_nr = recombination_info.get("tau_nr_ns", np.nan)
        plqy = recombination_info.get("plqy_percent", np.nan)
        k_rad_therm_s = recombination_info.get("k_rad_therm_s", 0.0)
        k_nr_s = recombination_info.get("k_nr_s", 0.0)

        print("=" * 68)
        print(" Recombination & Photoluminescence Summary")
        print("=" * 68)
        print(f"  Material                     : {mat_name} (Refractive Index n = {n_refr:.2f})")
        if np.isfinite(tau_rad_1):
            print(f"  Lowest Exciton Rad Lifetime  : {tau_rad_1:.2f} ns (k_rad = {1e9/max(tau_rad_1, 1e-12):.2e} s^-1)")
        if np.isfinite(tau_rad_therm):
            print(f"  Thermalized Rad Lifetime     : {tau_rad_therm:.2f} ns (k_rad = {k_rad_therm_s:.2e} s^-1)")
        if np.isfinite(tau_nr):
            if tau_nr > 1e6:
                print(f"  Effective Non-Rad Lifetime   : > 1 ms (intrinsic limit, k_nr = {k_nr_s:.2e} s^-1)")
            else:
                print(f"  Effective Non-Rad Lifetime   : {tau_nr:.2f} ns (k_nr = {k_nr_s:.2e} s^-1)")
        if np.isfinite(plqy):
            print(f"  Predicted PL Quantum Yield   : {plqy:.1f} %")
        if "tau_auger_ps" in recombination_info and np.isfinite(recombination_info["tau_auger_ps"]):
            tau_aug_ps = recombination_info["tau_auger_ps"]
            tau_aug_ns = tau_aug_ps * 1.0e-3
            k_aug_s = recombination_info.get("k_auger_s", 1e12 / max(tau_aug_ps, 1e-12))
            print(f"  Biexciton Auger Lifetime     : {tau_aug_ns:.4f} ns ({tau_aug_ps:.2f} ps, k_Auger = {k_aug_s:.2e} s^-1)")

        if bg_data is not None and "recomb_params" in bg_data:
            rp = bg_data["recomb_params"]
            from qdex.hardness import compute_energy_gap_law_rate, compute_fcwd_rate
            # Intraband RMS NACs are not the exciton-to-ground coupling.
            k_jort_s, _ = compute_energy_gap_law_rate(qp_gap_ev, E_LO_ev=rp["E_LO_ev"], S_hr=rp["S_hr"])
            print()
            print("  --- Trajectory-Derived Parameters (Non-Empirical from NAMD) ---")
            rayleigh = bg_data.get("rayleigh_cm", np.nan)
            rayleigh_txt = f", resolution {rayleigh:.0f} cm^-1" if np.isfinite(rayleigh) else ""
            print(f"  Dominant Optical Phonon      : {rp['dominant_freq_cm1']:.1f} cm^-1 (hbar*omega_LO = {rp['E_LO_ev']*1e3:.1f} meV{rayleigh_txt})")
            print(f"  Nuclear Reorganization (lam) : {rp['lambda_ev']*1e3:.1f} meV")
            print(f"  Huang-Rhys Factor (S)        : {rp['S_hr']:.3f}")
            print(f"  Thermal Gap Fluctuation (sig): {rp['sigma_ev']*1e3:.1f} meV")
            if rp.get("V_el_ev") is not None:
                print(f"  Electronic Coupling (V_el)   : {rp['V_el_ev']*1e3:.2f} meV (<|d_10|> = {mean_nac_fs:.4f} fs^-1)")
                k_fcwd_s, _ = compute_fcwd_rate(qp_gap_ev, V_el_ev=rp["V_el_ev"], lambda_ev=rp["lambda_ev"], sigma_ev=rp["sigma_ev"])
                if k_fcwd_s > 0:
                    tau_fcwd = 1e9 / k_fcwd_s
                    fcwd_str = f"{tau_fcwd:.2f} ns" if tau_fcwd < 1e6 else "> 1 ms (intrinsic wide-gap limit)"
                    print(f"  FCWD Multi-Phonon Rate       : {k_fcwd_s:.2e} s^-1 (tau_nr = {fcwd_str})")
            if k_jort_s > 0:
                tau_jort = 1e9 / k_jort_s
                jort_str = f"{tau_jort:.2f} ns" if tau_jort < 1e6 else "> 1 ms (intrinsic wide-gap limit)"
                print(f"  Jortner Energy Gap Law Rate  : {k_jort_s:.2e} s^-1 (tau_nr = {jort_str})")
        print("=" * 68 + "\n")

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

        if std_energies_ev is not None and np.any(std_energies_ev > 1e-6):
            orig_lbl = f" (±1σ thermal, {n_origins} origins)" if n_origins > 1 else " (±1σ thermal)"
            ax1.fill_between(
                times_fs,
                mean_energies_ev - std_energies_ev,
                mean_energies_ev + std_energies_ev,
                color="#1D4ED8",
                alpha=0.22,
                label=r"Thermal Spread" + orig_lbl
            )

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

        if std_excess_e is not None and np.any(std_excess_e > 1e-6):
            ax2.fill_between(times_fs, mean_excess_e - std_excess_e, mean_excess_e + std_excess_e, color="#2563EB", alpha=0.18)
        if std_excess_h is not None and np.any(std_excess_h > 1e-6):
            ax2.fill_between(times_fs, mean_excess_h - std_excess_h, mean_excess_h + std_excess_h, color="#EA580C", alpha=0.18)
        if std_energies_ev is not None and np.any(std_energies_ev > 1e-6):
            ax2.fill_between(times_fs, excess_total - std_energies_ev, excess_total + std_energies_ev, color="#374151", alpha=0.15)

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
        if bg_data is None and precompute_dir and os.path.isdir(precompute_dir):
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
        title_text="<b>QDEX — NAMD Hot Carrier Cooling & Exciton Relaxation (CsPbBr3)</b>",
        title_font=dict(size=18, family="sans-serif"),
        template="plotly_white",
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=-0.18, xanchor="center", x=0.5, font=dict(size=10))
    )

    fig.write_html(html_file, include_plotlyjs=True)


# =============================================================================
# Time-Resolved Vibrational Action Spectrum / Dynamical Phonon Spectrogram
# =============================================================================

def load_trajectory_orbital_energies(precompute_dir):
    """
    Extracts orbital energy trajectories eps_occ(t) and eps_virt(t) across all frames.
    Caches results in precompute_dir/orbital_energies.npz for instantaneous subsequent loading.
    """
    cache_path = os.path.join(precompute_dir, "orbital_energies.npz")
    if os.path.exists(cache_path):
        try:
            data = np.load(cache_path)
            return data["times_fs"], data["eps_occ_traj"], data["eps_virt_traj"]
        except Exception:
            pass

    step_files = sorted(glob.glob(os.path.join(precompute_dir, "step_*.npz")), key=natural_sort_key)
    if step_files:
        s0 = np.load(step_files[0])
        n_steps = len(step_files)
        n_occ = len(s0["eps_occ_curr"])
        n_virt = len(s0["eps_virt_curr"])
        times_fs = np.zeros(n_steps + 1, dtype=np.float64)
        eps_occ_traj = np.zeros((n_steps + 1, n_occ), dtype=np.float32)
        eps_virt_traj = np.zeros((n_steps + 1, n_virt), dtype=np.float32)

        times_fs[0] = float(s0["time_prev_fs"])
        eps_occ_traj[0] = s0["eps_occ_prev"]
        eps_virt_traj[0] = s0["eps_virt_prev"]

        for k, sf in enumerate(step_files):
            d = np.load(sf)
            times_fs[k + 1] = float(d["time_curr_fs"])
            eps_occ_traj[k + 1] = d["eps_occ_curr"]
            eps_virt_traj[k + 1] = d["eps_virt_curr"]
    else:
        frame_files = sorted(glob.glob(os.path.join(precompute_dir, "frame_*.npz")), key=natural_sort_key)
        if not frame_files:
            raise FileNotFoundError(f"No step or frame npz files found in '{precompute_dir}'")
        f0 = np.load(frame_files[0])
        n_frames = len(frame_files)
        n_occ = len(f0["eps_occ"])
        n_virt = len(f0["eps_virt"])
        times_fs = np.zeros(n_frames, dtype=np.float64)
        eps_occ_traj = np.zeros((n_frames, n_occ), dtype=np.float32)
        eps_virt_traj = np.zeros((n_frames, n_virt), dtype=np.float32)
        for k, ff in enumerate(frame_files):
            d = np.load(ff)
            times_fs[k] = float(d["time_fs"])
            eps_occ_traj[k] = d["eps_occ"]
            eps_virt_traj[k] = d["eps_virt"]

    try:
        np.savez_compressed(
            cache_path,
            times_fs=times_fs,
            eps_occ_traj=eps_occ_traj,
            eps_virt_traj=eps_virt_traj
        )
    except Exception:
        pass

    return times_fs, eps_occ_traj, eps_virt_traj


def compute_pair_spectral_density(eps_1, eps_2, dt_fs, n_fft=None, w_max_cm=400.0):
    r"""
    Computes the phonon spectral density J_{uv}(omega) in cm^-1 for a pair of electronic states
    with time-dependent orbital energy trajectories eps_1(t) and eps_2(t):

      Delta eps_{uv}(t) = eps_1(t) - eps_2(t)
      delta Delta eps_{uv}(t) = Delta eps_{uv}(t) - <Delta eps_{uv}>
      sigma_{uv}^2 = Var(Delta eps_{uv})
      C_{uv}(tau) = <delta Delta eps_{uv}(0) delta Delta eps_{uv}(tau)> / sigma_{uv}^2
      J_{uv}(omega) = sigma_{uv}^2 * max(Re FFT[C_{uv}(tau) * W(tau)], 0)

    Returns:
      wavenumbers_cm: 1D array of frequencies in cm^-1 (0 to w_max_cm).
      psd: 1D array of spectral density values.
      var_d: variance of the energy difference fluctuation in eV^2.
    """
    eps_1 = np.asarray(eps_1, dtype=np.float64)
    eps_2 = np.asarray(eps_2, dtype=np.float64)
    N = len(eps_1)
    diff = eps_1 - eps_2
    delta_diff = diff - np.mean(diff)
    var_d = float(np.var(delta_diff))

    if n_fft is None:
        n_fft = max(2048, N * 16)

    freqs_fs = np.fft.rfftfreq(n_fft, d=dt_fs)
    c_cm_per_fs = 2.99792458e-5
    wavenumbers_cm = freqs_fs / c_cm_per_fs
    mask = (wavenumbers_cm >= 0) & (wavenumbers_cm <= w_max_cm)
    wn_sub = wavenumbers_cm[mask]

    if var_d < 1e-16 or N < 4:
        return wn_sub, np.zeros(len(wn_sub), dtype=np.float64), var_d

    raw_corr = np.correlate(delta_diff, delta_diff, mode="full")[N - 1:]
    norm_factors = np.arange(N, 0, -1) * var_d
    autocorr = raw_corr / np.maximum(norm_factors, 1e-20)

    window = np.hanning(N)
    fft_vals = np.fft.rfft(autocorr * window, n=n_fft)
    psd_full = np.maximum(np.real(fft_vals), 0.0)
    psd = var_d * psd_full[mask]

    return wn_sub, psd, var_d


def compute_time_resolved_spectral_density(
    times_fs,
    eps_occ_traj,
    eps_virt_traj,
    hop_records=None,
    pme_flux_records=None,
    n_trajectories=1000,
    sigma_t_fs=15.0,
    w_max_cm=400.0,
    bg_data=None,
    method_name="PME"
):
    r"""
    Synthesizes the Time-Resolved Vibrational Action Spectrum / Dynamical Phonon Spectrogram
    J(omega, t) during hot carrier cooling.

    Parameters:
      times_fs: 1D array of simulation time points (fs).
      eps_occ_traj: 2D array of occupied orbital energies (n_frames, n_occ).
      eps_virt_traj: 2D array of virtual orbital energies (n_frames, n_virt).
      hop_records: List of discrete hop dicts from stochastic surface hopping (FSSH/DISH/GDC).
      pme_flux_records: List of probability flux dicts from Pauli Master Equation (PME).
      n_trajectories: Number of trajectories (for FSSH/DISH normalization).
      sigma_t_fs: Temporal Gaussian broadening in fs.
      w_max_cm: Maximum phonon frequency in cm^-1 to include in the spectrogram.
      bg_data: Optional band-gap dynamics dict from compute_band_gap_dynamics_and_spectral_density.
      method_name: String label for the simulation method.

    Returns:
      Dictionary containing 2D arrays J_total, J_electron, J_hole, 1D frequency and time axes,
      time-integrated spectra, and metadata.
    """
    times_fs = np.asarray(times_fs, dtype=np.float64)
    N_t = len(times_fs)
    dt_fs = times_fs[1] - times_fs[0] if N_t > 1 else 2.0

    # Determine frequency grid from first virtual orbital
    sample_wn, _, _ = compute_pair_spectral_density(
        eps_virt_traj[:, 0], eps_virt_traj[:, 0], dt_fs, w_max_cm=w_max_cm
    )
    N_w = len(sample_wn)

    J_elec = np.zeros((N_w, N_t), dtype=np.float64)
    J_hole = np.zeros((N_w, N_t), dtype=np.float64)

    pair_cache = {}
    visited_pairs = {"electron": set(), "hole": set()}

    def _get_pair_psd(channel, u, v):
        key = (channel, min(u, v), max(u, v))
        if key not in pair_cache:
            if channel == "electron":
                e1 = eps_virt_traj[:, u]
                e2 = eps_virt_traj[:, v]
            else:
                e1 = eps_occ_traj[:, u]
                e2 = eps_occ_traj[:, v]
            _, psd, _ = compute_pair_spectral_density(e1, e2, dt_fs, w_max_cm=w_max_cm)
            pair_cache[key] = psd
        return pair_cache[key]

    sigma_t = max(float(sigma_t_fs), 1.0)
    # 1D Gaussian kernel for smoothing discrete transitions/fluxes along the time axis
    half_width = max(int(np.ceil(4.0 * sigma_t / max(dt_fs, 0.1))), 2)
    t_win = np.arange(-half_width, half_width + 1) * dt_fs
    gauss_kernel = np.exp(-0.5 * (t_win / sigma_t) ** 2)
    gauss_kernel /= np.maximum(np.sum(gauss_kernel) * dt_fs, 1e-30)

    pair_weights_e = {}
    pair_weights_h = {}

    # --- Mode A: Stochastic Surface Hopping (CPA-FSSH, DISH, CPA-FSSH-GDC) ---
    if hop_records is not None and len(hop_records) > 0:
        n_traj = max(int(n_trajectories), 1)
        for h in hop_records:
            t_hop = float(h["time_fs"])
            ch = str(h["channel"]).lower()
            u = int(h["from"])
            v = int(h["to"])
            if u == v:
                continue
            pair_key = (min(u, v), max(u, v))
            k_idx = int(np.clip(round((t_hop - times_fs[0]) / dt_fs), 0, N_t - 1))
            weight = 1.0 / n_traj

            if ch == "electron":
                if pair_key not in pair_weights_e:
                    pair_weights_e[pair_key] = np.zeros(N_t, dtype=np.float64)
                pair_weights_e[pair_key][k_idx] += weight
                visited_pairs["electron"].add(pair_key)
            else:
                if pair_key not in pair_weights_h:
                    pair_weights_h[pair_key] = np.zeros(N_t, dtype=np.float64)
                pair_weights_h[pair_key][k_idx] += weight
                visited_pairs["hole"].add(pair_key)

    # --- Mode B: Deterministic Master Equation (PME) ---
    elif pme_flux_records is not None and len(pme_flux_records) > 0:
        for rec in pme_flux_records:
            t_rec = float(rec["time_fs"])
            k_idx = int(np.clip(round((t_rec - times_fs[0]) / dt_fs), 0, N_t - 1))
            flux_v = rec.get("flux_virt", None)
            flux_o = rec.get("flux_occ", None)

            if flux_v is not None:
                if isinstance(flux_v, list):
                    for a, b, phi in flux_v:
                        pair_key = (min(a, b), max(a, b))
                        if pair_key not in pair_weights_e:
                            pair_weights_e[pair_key] = np.zeros(N_t, dtype=np.float64)
                        pair_weights_e[pair_key][k_idx] += phi
                        visited_pairs["electron"].add(pair_key)
                elif isinstance(flux_v, np.ndarray):
                    n_v = flux_v.shape[0]
                    np.fill_diagonal(flux_v, 0.0)
                    tot_v = np.sum(flux_v)
                    if tot_v > 1e-12:
                        thresh = 1e-4 * tot_v
                        r, c = np.where(flux_v >= thresh)
                        v_vals = flux_v[r, c]
                        if len(v_vals) > 50:
                            top_k = np.argpartition(v_vals, -50)[-50:]
                            r, c, v_vals = r[top_k], c[top_k], v_vals[top_k]
                        for a, b, phi in zip(r, c, v_vals):
                            pair_key = (min(a, b), max(a, b))
                            if pair_key not in pair_weights_e:
                                pair_weights_e[pair_key] = np.zeros(N_t, dtype=np.float64)
                            pair_weights_e[pair_key][k_idx] += phi
                            visited_pairs["electron"].add(pair_key)

            if flux_o is not None:
                if isinstance(flux_o, list):
                    for i, j, phi in flux_o:
                        pair_key = (min(i, j), max(i, j))
                        if pair_key not in pair_weights_h:
                            pair_weights_h[pair_key] = np.zeros(N_t, dtype=np.float64)
                        pair_weights_h[pair_key][k_idx] += phi
                        visited_pairs["hole"].add(pair_key)
                elif isinstance(flux_o, np.ndarray):
                    n_o = flux_o.shape[0]
                    np.fill_diagonal(flux_o, 0.0)
                    tot_o = np.sum(flux_o)
                    if tot_o > 1e-12:
                        thresh = 1e-4 * tot_o
                        r, c = np.where(flux_o >= thresh)
                        v_vals = flux_o[r, c]
                        if len(v_vals) > 50:
                            top_k = np.argpartition(v_vals, -50)[-50:]
                            r, c, v_vals = r[top_k], c[top_k], v_vals[top_k]
                        for i, j, phi in zip(r, c, v_vals):
                            pair_key = (min(i, j), max(i, j))
                            if pair_key not in pair_weights_h:
                                pair_weights_h[pair_key] = np.zeros(N_t, dtype=np.float64)
                            pair_weights_h[pair_key][k_idx] += phi
                            visited_pairs["hole"].add(pair_key)

    # Convolve 1D weights along time and synthesize 2D action spectrograms
    for (a, b), w_series in pair_weights_e.items():
        w_smooth = np.convolve(w_series, gauss_kernel, mode="same")
        psd = _get_pair_psd("electron", a, b)
        J_elec += np.outer(psd, w_smooth)

    for (i, j), w_series in pair_weights_h.items():
        w_smooth = np.convolve(w_series, gauss_kernel, mode="same")
        psd = _get_pair_psd("hole", i, j)
        J_hole += np.outer(psd, w_smooth)

    J_total = J_elec + J_hole
    max_amp = float(np.max(J_total)) if np.max(J_total) > 0 else 1.0

    # Normalized 2D spectrograms [0, 1]
    J_total_norm = J_total / max_amp
    J_elec_norm = J_elec / max_amp
    J_hole_norm = J_hole / max_amp

    # Time-integrated vibrational action spectra
    _trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
    J_int_total = _trapz(J_total, x=times_fs, axis=1)
    J_int_elec = _trapz(J_elec, x=times_fs, axis=1)
    J_int_hole = _trapz(J_hole, x=times_fs, axis=1)

    peak_int = float(np.max(J_int_total)) if np.max(J_int_total) > 0 else 1.0
    J_int_total_norm = J_int_total / peak_int
    J_int_elec_norm = J_int_elec / peak_int
    J_int_hole_norm = J_int_hole / peak_int

    # Band-gap excitation spectral density comparison
    J_bg_interp = None
    if bg_data is not None and "wavenumbers_cm" in bg_data and "psd" in bg_data:
        bg_wn = bg_data["wavenumbers_cm"]
        bg_psd = bg_data["psd"]
        J_bg_interp = np.interp(sample_wn, bg_wn, bg_psd, left=0.0, right=0.0)
        max_bg = float(np.max(J_bg_interp)) if np.max(J_bg_interp) > 0 else 1.0
        J_bg_interp /= max_bg

    return {
        "times_fs": times_fs,
        "wavenumbers_cm": sample_wn,
        "J_total": J_total_norm,
        "J_electron": J_elec_norm,
        "J_hole": J_hole_norm,
        "J_total_raw": J_total,
        "J_int_total": J_int_total_norm,
        "J_int_elec": J_int_elec_norm,
        "J_int_hole": J_int_hole_norm,
        "J_bandgap": J_bg_interp,
        "n_unique_pairs": len(visited_pairs["electron"]) + len(visited_pairs["hole"]),
        "visited_pairs": visited_pairs,
        "method": method_name,
    }


def plot_time_resolved_spectral_density(
    tr_sd_data,
    plot_file="time_resolved_spectral_density.png",
    html_file="time_resolved_spectral_density.html",
    material_name="CsPbBr3",
    method_name="PME"
):
    r"""
    Generates a 4-panel publication-quality Matplotlib figure (PNG) and an interactive Plotly HTML
    visualization of the Time-Resolved Vibrational Action Spectrum J(omega, t).
    """
    times = tr_sd_data["times_fs"]
    wn = tr_sd_data["wavenumbers_cm"]
    J_tot = tr_sd_data["J_total"]
    J_e = tr_sd_data["J_electron"]
    J_h = tr_sd_data["J_hole"]
    J_int_tot = tr_sd_data["J_int_total"]
    J_int_e = tr_sd_data["J_int_elec"]
    J_int_h = tr_sd_data["J_int_hole"]
    J_bg = tr_sd_data["J_bandgap"]
    method_lbl = tr_sd_data.get("method", method_name).upper()

    # -------------------------------------------------------------
    # 1. Matplotlib Figure (PNG)
    # -------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), constrained_layout=True)

    # Panel (a): Total Phonon Spectrogram
    im0 = axes[0, 0].pcolormesh(times, wn, J_tot, cmap="inferno", shading="auto", vmin=0, vmax=1.0)
    axes[0, 0].set_title(f"(a) Total Vibrational Spectrogram $J_{{\\mathrm{{tot}}}}(\\omega, t)$ [{method_lbl}]", fontsize=13, fontweight="bold")
    axes[0, 0].set_xlabel("Delay Time $t$ (fs)", fontsize=11)
    axes[0, 0].set_ylabel("Phonon Frequency $\\omega$ (cm$^{-1}$)", fontsize=11)
    cbar0 = fig.colorbar(im0, ax=axes[0, 0], fraction=0.046, pad=0.04)
    cbar0.set_label("Action Density (norm.)", fontsize=10)

    # Panel (b): Time-Integrated Footprint vs Band-Edge Exciton PSD
    axes[0, 1].plot(wn, J_int_tot, color="#111827", lw=2.4, label="Cooling Total $J_{\\mathrm{int}}(\\omega)$")
    axes[0, 1].plot(wn, J_int_e, color="#2563EB", lw=1.8, ls="--", label="Electron Cooling $J_e$")
    axes[0, 1].plot(wn, J_int_h, color="#DC2626", lw=1.8, ls="--", label="Hole Cooling $J_h$")
    if J_bg is not None:
        axes[0, 1].plot(wn, J_bg, color="#059669", lw=2.2, ls=":", label="Band-Gap Exciton $J_{\\mathrm{gap}}(\\omega)$")
        axes[0, 1].fill_between(wn, J_bg, color="#10B981", alpha=0.15)

    # Identify dominant peak in cooling spectrum (excluding DC < 30 cm^-1)
    mask_opt = wn >= 30.0
    if np.any(mask_opt):
        p_opt = np.where(mask_opt)[0][np.argmax(J_int_tot[mask_opt])]
        axes[0, 1].annotate(
            f"Active Mode ({wn[p_opt]:.0f} cm$^{{-1}}$)",
            xy=(wn[p_opt], J_int_tot[p_opt]),
            xytext=(wn[p_opt] + 25, J_int_tot[p_opt] * 0.85),
            arrowprops=dict(facecolor="#111827", shrink=0.08, width=1.5, headwidth=6),
            fontsize=10, fontweight="semibold"
        )

    axes[0, 1].set_title("(b) Dynamical Action Footprint vs. Band-Gap PSD", fontsize=13, fontweight="bold")
    axes[0, 1].set_xlabel("Phonon Frequency $\\omega$ (cm$^{-1}$)", fontsize=11)
    axes[0, 1].set_ylabel("Spectral Density (arb. units)", fontsize=11)
    axes[0, 1].set_ylim(-0.02, 1.08)
    axes[0, 1].legend(loc="upper right", frameon=True, fontsize=10)
    axes[0, 1].grid(True, ls=":", alpha=0.6)

    # Panel (c): Electron Cooling Channel
    im1 = axes[1, 0].pcolormesh(times, wn, J_e, cmap="viridis", shading="auto", vmin=0, vmax=max(0.01, np.max(J_e)))
    axes[1, 0].set_title("(c) Electron Cooling Channel $J_e(\\omega, t)$", fontsize=13, fontweight="bold")
    axes[1, 0].set_xlabel("Delay Time $t$ (fs)", fontsize=11)
    axes[1, 0].set_ylabel("Phonon Frequency $\\omega$ (cm$^{-1}$)", fontsize=11)
    cbar1 = fig.colorbar(im1, ax=axes[1, 0], fraction=0.046, pad=0.04)
    cbar1.set_label("Action Density", fontsize=10)

    # Panel (d): Hole Cooling Channel
    im2 = axes[1, 1].pcolormesh(times, wn, J_h, cmap="plasma", shading="auto", vmin=0, vmax=max(0.01, np.max(J_h)))
    axes[1, 1].set_title("(d) Hole Cooling Channel $J_h(\\omega, t)$", fontsize=13, fontweight="bold")
    axes[1, 1].set_xlabel("Delay Time $t$ (fs)", fontsize=11)
    axes[1, 1].set_ylabel("Phonon Frequency $\\omega$ (cm$^{-1}$)", fontsize=11)
    cbar2 = fig.colorbar(im2, ax=axes[1, 1], fraction=0.046, pad=0.04)
    cbar2.set_label("Action Density", fontsize=10)

    fig.suptitle(
        f"QDEX — Time-Resolved Vibrational Action Spectrum $J(\\omega, t)$ ({material_name} | {method_lbl})",
        fontsize=16, fontweight="bold", y=1.02
    )

    fig.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [NAMD:Spectrogram] Saved publication plot to: {plot_file}")

    # -------------------------------------------------------------
    # 2. Interactive Plotly Figure (HTML)
    # -------------------------------------------------------------
    try:
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig_plotly = make_subplots(
            rows=2, cols=2,
            subplot_titles=[
                f"<b>(a) Total Vibrational Spectrogram J(w, t) [{method_lbl}]</b>",
                "<b>(b) Time-Integrated Action Footprint vs Bandgap PSD</b>",
                "<b>(c) Electron Cooling Channel J_e(w, t)</b>",
                "<b>(d) Hole Cooling Channel J_h(w, t)</b>"
            ],
            horizontal_spacing=0.09,
            vertical_spacing=0.12
        )

        fig_plotly.add_trace(
            go.Heatmap(
                x=times, y=wn, z=J_tot,
                colorscale="Inferno", zmin=0.0, zmax=1.0,
                colorbar=dict(title="Action Density", x=0.45, len=0.45, y=0.8),
                hovertemplate="Time: %{x:.1f} fs<br>Freq: %{y:.1f} cm^-1<br>J_tot: %{z:.4f}<extra></extra>"
            ), row=1, col=1
        )

        fig_plotly.add_trace(
            go.Scatter(
                x=wn, y=J_int_tot,
                mode="lines", name="Total Cooling Footprint",
                line=dict(color="#111827", width=2.6),
                hovertemplate="Freq: %{x:.1f} cm^-1<br>J_int: %{y:.4f}<extra></extra>"
            ), row=1, col=2
        )
        fig_plotly.add_trace(
            go.Scatter(
                x=wn, y=J_int_e,
                mode="lines", name="Electron Cooling Footprint",
                line=dict(color="#2563EB", width=2.0, dash="dash"),
                hovertemplate="Freq: %{x:.1f} cm^-1<br>J_e: %{y:.4f}<extra></extra>"
            ), row=1, col=2
        )
        fig_plotly.add_trace(
            go.Scatter(
                x=wn, y=J_int_h,
                mode="lines", name="Hole Cooling Footprint",
                line=dict(color="#DC2626", width=2.0, dash="dash"),
                hovertemplate="Freq: %{x:.1f} cm^-1<br>J_h: %{y:.4f}<extra></extra>"
            ), row=1, col=2
        )
        if J_bg is not None:
            fig_plotly.add_trace(
                go.Scatter(
                    x=wn, y=J_bg,
                    mode="lines", name="Band-Gap Exciton PSD",
                    line=dict(color="#059669", width=2.2, dash="dot"),
                    fill="tozeroy", fillcolor="rgba(16, 185, 129, 0.15)",
                    hovertemplate="Freq: %{x:.1f} cm^-1<br>J_gap: %{y:.4f}<extra></extra>"
                ), row=1, col=2
            )

        fig_plotly.add_trace(
            go.Heatmap(
                x=times, y=wn, z=J_e,
                colorscale="Viridis",
                colorbar=dict(title="J_e", x=0.45, len=0.45, y=0.2),
                hovertemplate="Time: %{x:.1f} fs<br>Freq: %{y:.1f} cm^-1<br>J_e: %{z:.4f}<extra></extra>"
            ), row=2, col=1
        )

        fig_plotly.add_trace(
            go.Heatmap(
                x=times, y=wn, z=J_h,
                colorscale="Plasma",
                colorbar=dict(title="J_h", x=1.02, len=0.45, y=0.2),
                hovertemplate="Time: %{x:.1f} fs<br>Freq: %{y:.1f} cm^-1<br>J_h: %{z:.4f}<extra></extra>"
            ), row=2, col=2
        )

        fig_plotly.update_xaxes(title_text="Delay Time t (fs)", row=1, col=1)
        fig_plotly.update_yaxes(title_text="Phonon Frequency (cm^-1)", row=1, col=1)

        fig_plotly.update_xaxes(title_text="Phonon Frequency (cm^-1)", row=1, col=2)
        fig_plotly.update_yaxes(title_text="Spectral Density (arb. units)", row=1, col=2)

        fig_plotly.update_xaxes(title_text="Delay Time t (fs)", row=2, col=1)
        fig_plotly.update_yaxes(title_text="Phonon Frequency (cm^-1)", row=2, col=1)

        fig_plotly.update_xaxes(title_text="Delay Time t (fs)", row=2, col=2)
        fig_plotly.update_yaxes(title_text="Phonon Frequency (cm^-1)", row=2, col=2)

        fig_plotly.update_layout(
            height=900,
            width=1500,
            title_text=f"<b>QDEX — Time-Resolved Vibrational Action Spectrum J(w, t) ({material_name} | {method_lbl})</b>",
            title_font=dict(size=18, family="sans-serif"),
            template="plotly_white",
            legend=dict(orientation="h", yanchor="bottom", y=-0.12, xanchor="center", x=0.5)
        )

        fig_plotly.write_html(html_file, include_plotlyjs=True)
        print(f"  [NAMD:Spectrogram] Saved interactive HTML to: {html_file}")
    except Exception as e:
        print(f"  [NAMD:Warn] Failed to create Plotly HTML widget: {e}")


def export_time_resolved_spectral_density(
    tr_sd_data,
    output_npz="time_resolved_spectral_density.npz",
    output_csv="time_resolved_spectral_density.csv"
):
    """
    Exports 2D spectrogram matrices and 1D time-integrated spectra to compressed NPZ and CSV.
    """
    np.savez_compressed(
        output_npz,
        times_fs=tr_sd_data["times_fs"],
        wavenumbers_cm=tr_sd_data["wavenumbers_cm"],
        J_total=tr_sd_data["J_total"],
        J_electron=tr_sd_data["J_electron"],
        J_hole=tr_sd_data["J_hole"],
        J_int_total=tr_sd_data["J_int_total"],
        J_int_elec=tr_sd_data["J_int_elec"],
        J_int_hole=tr_sd_data["J_int_hole"],
        J_bandgap=tr_sd_data["J_bandgap"] if tr_sd_data["J_bandgap"] is not None else np.zeros_like(tr_sd_data["wavenumbers_cm"]),
        method=str(tr_sd_data.get("method", "NAMD"))
    )
    print(f"  [NAMD:Spectrogram] Exported 2D spectrogram arrays to: {output_npz}")

    # Export 1D time-integrated spectrum to CSV
    cols = [
        tr_sd_data["wavenumbers_cm"],
        tr_sd_data["J_int_total"],
        tr_sd_data["J_int_elec"],
        tr_sd_data["J_int_hole"]
    ]
    hdr = "wavenumber_cm,J_int_total,J_int_electron,J_int_hole"
    if tr_sd_data["J_bandgap"] is not None:
        cols.append(tr_sd_data["J_bandgap"])
        hdr += ",J_bandgap"

    np.savetxt(
        output_csv,
        np.column_stack(cols),
        header=hdr,
        delimiter=",",
        comments=""
    )
    print(f"  [NAMD:Spectrogram] Exported 1D integrated spectra to: {output_csv}")


def compute_2d_vibronic_action_map(
    tr_sd_data,
    n_fft=2048,
    detrend_mode="linear",
    wmax_acc=250.0,
    wmax_prom=350.0
):
    r"""
    Decomposes the time-resolved vibrational action spectrum J(omega_acc, t)
    along the delay time axis t into coherent promoting / driving phonon frequencies:
    S(omega_acc, Omega_prom) = abs(FFT_t[ J_osc(omega_acc, t) * W(t) ])

    Parameters
    ----------
    tr_sd_data : dict
        Output dictionary from compute_time_resolved_spectral_density.
    n_fft : int
        Zero-padded FFT points for high promoting mode frequency resolution.
    detrend_mode : str
        'linear' for standard linear detrending or 'savgol' for Savitzky-Golay baseline subtraction.
    wmax_acc : float
        Maximum frequency cut-off (cm^-1) for accepting modes.
    wmax_prom : float
        Maximum frequency cut-off (cm^-1) for promoting modes.

    Returns
    -------
    dict
        Contains:
          - 'wavenumbers_acc_cm': 1D array of accepting frequencies
          - 'wavenumbers_prom_cm': 1D array of promoting frequencies
          - 'S_2d_total': 2D array [n_acc, n_prom]
          - 'S_2d_elec': 2D array [n_acc, n_prom]
          - 'S_2d_hole': 2D array [n_acc, n_prom]
          - 'proj_acc': 1D marginal accepting action
          - 'proj_prom': 1D marginal promoting action
          - 'dominant_promoting_modes': list of peak tuples (freq_cm, amplitude)
          - 'J_bandgap': 1D array of optical exciton spectral density
          - 'method': method name string
    """
    from scipy.signal import detrend, savgol_filter, find_peaks

    times_fs = tr_sd_data["times_fs"]
    wn_acc_full = tr_sd_data["wavenumbers_cm"]
    J_tot = tr_sd_data["J_total"]
    J_e = tr_sd_data["J_electron"]
    J_h = tr_sd_data["J_hole"]
    J_bg_full = tr_sd_data.get("J_bandgap")

    dt_fs = float(times_fs[1] - times_fs[0]) if len(times_fs) > 1 else 1.0
    n_times = len(times_fs)

    # Filter accepting frequency range
    mask_acc = (wn_acc_full >= 0.0) & (wn_acc_full <= wmax_acc)
    wn_acc = wn_acc_full[mask_acc]
    J_tot_sub = J_tot[mask_acc, :]
    J_e_sub = J_e[mask_acc, :]
    J_h_sub = J_h[mask_acc, :]

    # Detrend along the time axis (axis=1) to isolate coherent oscillations from macroscopic cooling decay
    if detrend_mode == "savgol" and n_times > 15:
        win_len = min(n_times if n_times % 2 != 0 else n_times - 1, max(7, int(150.0 / dt_fs) | 1))
        J_tot_osc = J_tot_sub - savgol_filter(J_tot_sub, window_length=win_len, polyorder=2, axis=1)
        J_e_osc = J_e_sub - savgol_filter(J_e_sub, window_length=win_len, polyorder=2, axis=1)
        J_h_osc = J_h_sub - savgol_filter(J_h_sub, window_length=win_len, polyorder=2, axis=1)
    else:
        J_tot_osc = detrend(J_tot_sub, axis=1)
        J_e_osc = detrend(J_e_sub, axis=1)
        J_h_osc = detrend(J_h_sub, axis=1)

    # Windowing along time axis
    w_hann = np.hanning(n_times)[None, :]
    S_2d_tot_full = np.abs(np.fft.rfft(J_tot_osc * w_hann, n=n_fft, axis=1))
    S_2d_e_full = np.abs(np.fft.rfft(J_e_osc * w_hann, n=n_fft, axis=1))
    S_2d_h_full = np.abs(np.fft.rfft(J_h_osc * w_hann, n=n_fft, axis=1))

    # Frequency axis for promoting modes
    freqs_thz = np.fft.rfftfreq(n_fft, d=dt_fs * 1e-15) * 1e-12
    wn_prom_full = freqs_thz * 33.35641

    mask_prom = (wn_prom_full >= 5.0) & (wn_prom_full <= wmax_prom)
    wn_prom = wn_prom_full[mask_prom]

    S_2d_tot = S_2d_tot_full[:, mask_prom]
    S_2d_e = S_2d_e_full[:, mask_prom]
    S_2d_h = S_2d_h_full[:, mask_prom]

    # Marginal 1D projections
    proj_acc = np.sum(S_2d_tot, axis=1)
    proj_prom = np.sum(S_2d_tot, axis=0)

    # Detect dominant promoting peaks
    pks, _ = find_peaks(proj_prom, height=np.max(proj_prom) * 0.15, distance=max(1, int(15.0 / (wn_prom[1] - wn_prom[0]))))
    dom_modes = [(float(wn_prom[p]), float(proj_prom[p])) for p in pks]

    # Interpolate bandgap PSD onto accepting frequencies
    J_bg_acc = np.interp(wn_acc, wn_acc_full, J_bg_full) if J_bg_full is not None else None

    return {
        "wavenumbers_acc_cm": wn_acc,
        "wavenumbers_prom_cm": wn_prom,
        "S_2d_total": S_2d_tot,
        "S_2d_elec": S_2d_e,
        "S_2d_hole": S_2d_h,
        "proj_acc": proj_acc,
        "proj_prom": proj_prom,
        "dominant_promoting_modes": dom_modes,
        "J_bandgap": J_bg_acc,
        "J_bandgap_full": J_bg_full,
        "wavenumbers_full": wn_acc_full,
        "method": tr_sd_data.get("method", "NAMD")
    }


def plot_2d_vibronic_action_map(
    vib2d_data,
    output_png="2d_vibronic_action_map.png",
    output_html="2d_vibronic_action_map.html",
    material_name="CsPbBr3",
    method_name=None
):
    """
    Renders publication-grade 2D Vibronic Action Map with marginal projections and Plotly HTML widget.
    """
    wn_acc = vib2d_data["wavenumbers_acc_cm"]
    wn_prom = vib2d_data["wavenumbers_prom_cm"]
    S_tot = vib2d_data["S_2d_total"]
    proj_acc = vib2d_data["proj_acc"]
    proj_prom = vib2d_data["proj_prom"]
    J_bg_full = vib2d_data.get("J_bandgap_full")
    wn_full = vib2d_data.get("wavenumbers_full", wn_acc)
    method_lbl = method_name or vib2d_data.get("method", "NAMD").upper()
    dom_modes = vib2d_data.get("dominant_promoting_modes", [])

    fig = plt.figure(figsize=(10, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[4, 1.25], height_ratios=[1.25, 4],
                           hspace=0.08, wspace=0.08)

    ax_main = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)

    norm_max = np.max(S_tot) if np.max(S_tot) > 0 else 1.0
    c = ax_main.pcolormesh(
        wn_acc, wn_prom, S_tot.T / norm_max,
        shading="gouraud", cmap="magma", vmin=0, vmax=1
    )
    ax_main.set_xlabel(r"Accepting Mode Frequency $\omega_{\mathrm{acc}}$ (cm$^{-1}$)", fontsize=12, fontweight="bold")
    ax_main.set_ylabel(r"Promoting Mode Frequency $\Omega_{\mathrm{prom}}$ (cm$^{-1}$)", fontsize=12, fontweight="bold")
    ax_main.set_xlim(wn_acc[0], wn_acc[-1])
    ax_main.set_ylim(wn_prom[0], wn_prom[-1])
    ax_main.grid(True, ls=":", alpha=0.5)

    # Highlight dominant dual-phonon cross peak if available
    if len(dom_modes) > 0:
        opt_modes = [m for m in dom_modes if 120.0 <= m[0] <= 180.0]
        ref_mode = opt_modes[0] if len(opt_modes) > 0 else dom_modes[0]
        pk_acc = wn_acc[np.argmax(proj_acc)]
        pk_prom = ref_mode[0]
        ax_main.scatter([pk_acc], [pk_prom], color="cyan", s=90, marker="+", lw=2, zorder=5)
        ax_main.annotate(
            f"Dual-Phonon Cross-Peak\n({pk_acc:.0f} cm$^{{-1}}$ acc, {pk_prom:.0f} cm$^{{-1}}$ prom)",
            xy=(pk_acc, pk_prom), xytext=(pk_acc + 30, pk_prom + 20),
            arrowprops=dict(facecolor="cyan", shrink=0.05, width=1.5, headwidth=6),
            fontweight="bold", fontsize=9, color="white",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.8, edgecolor="cyan")
        )

    # Top Panel: Accepting Mode Projection
    p_acc_norm = proj_acc / np.max(proj_acc) if np.max(proj_acc) > 0 else proj_acc
    ax_top.plot(wn_acc, p_acc_norm, color="#DC2626", lw=2.2)
    ax_top.fill_between(wn_acc, 0, p_acc_norm, color="#DC2626", alpha=0.2)
    ax_top.set_ylabel("Accepting\nAction", fontsize=10, fontweight="bold")
    ax_top.set_title(
        rf"QDEX — 2D Non-Adiabatic Vibronic Action Map ({material_name} | {method_lbl})",
        fontsize=13, fontweight="bold", pad=12
    )
    ax_top.grid(True, ls=":", alpha=0.5)
    plt.setp(ax_top.get_xticklabels(), visible=False)

    # Right Panel: Promoting Mode Projection overlaid with Optical Band-Gap Spectrum
    p_prom_norm = proj_prom / np.max(proj_prom) if np.max(proj_prom) > 0 else proj_prom
    ax_right.plot(p_prom_norm, wn_prom, color="#2563EB", lw=2.2, label=r"Promoting $\Omega_{\mathrm{prom}}$")

    if J_bg_full is not None:
        j_bg_interp = np.interp(wn_prom, wn_full, J_bg_full)
        j_bg_norm = j_bg_interp / np.max(j_bg_interp) if np.max(j_bg_interp) > 0 else j_bg_interp
        ax_right.plot(j_bg_norm, wn_prom, color="#059669", lw=1.8, ls="--", label=r"Optical $J_{\mathrm{gap}}$")
        ax_right.fill_betweenx(wn_prom, 0, j_bg_norm, color="#059669", alpha=0.15)

    ax_right.set_xlabel("Promoting\nAction", fontsize=10, fontweight="bold")
    ax_right.grid(True, ls=":", alpha=0.5)
    ax_right.legend(loc="upper right", fontsize=8)
    plt.setp(ax_right.get_yticklabels(), visible=False)

    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [NAMD:2D-Map] Saved publication plot to: {output_png}")

    # Plotly interactive widget
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        fig_p = make_subplots(
            rows=2, cols=2,
            column_widths=[0.75, 0.25],
            row_heights=[0.25, 0.75],
            vertical_spacing=0.04,
            horizontal_spacing=0.04,
            shared_xaxes=True,
            shared_yaxes=True
        )

        fig_p.add_trace(
            go.Scatter(
                x=wn_acc, y=p_acc_norm,
                mode="lines", fill="tozeroy",
                line=dict(color="#DC2626", width=2.5),
                name="Accepting Modes (Cooling)",
                hovertemplate="w_acc: %{x:.1f} cm^-1<br>Action: %{y:.3f}<extra></extra>"
            ), row=1, col=1
        )

        fig_p.add_trace(
            go.Heatmap(
                x=wn_acc, y=wn_prom, z=(S_tot.T / norm_max),
                colorscale="Magma",
                colorbar=dict(title="Action", x=1.02, len=0.75, y=0.38),
                hovertemplate="w_acc: %{x:.1f} cm^-1<br>Omega_prom: %{y:.1f} cm^-1<br>Coupling: %{z:.3f}<extra></extra>"
            ), row=2, col=1
        )

        fig_p.add_trace(
            go.Scatter(
                x=p_prom_norm, y=wn_prom,
                mode="lines",
                line=dict(color="#2563EB", width=2.5),
                name="Promoting Modes (NAC Driving)",
                hovertemplate="Omega_prom: %{y:.1f} cm^-1<br>Promoting Action: %{x:.3f}<extra></extra>"
            ), row=2, col=2
        )

        if J_bg_full is not None:
            fig_p.add_trace(
                go.Scatter(
                    x=j_bg_norm, y=wn_prom,
                    mode="lines", fill="tozerox",
                    line=dict(color="#059669", width=2, dash="dash"),
                    name="Bandgap Optical Phonon PSD",
                    hovertemplate="Omega: %{y:.1f} cm^-1<br>J_gap: %{x:.3f}<extra></extra>"
                ), row=2, col=2
            )

        fig_p.update_xaxes(title_text="Accepting Mode Frequency (cm^-1)", row=2, col=1)
        fig_p.update_yaxes(title_text="Promoting Mode Frequency (cm^-1)", row=2, col=1)
        fig_p.update_yaxes(title_text="Accepting Action", row=1, col=1)
        fig_p.update_xaxes(title_text="Promoting Action", row=2, col=2)

        fig_p.update_layout(
            width=1000, height=900,
            title_text=f"<b>QDEX — 2D Non-Adiabatic Vibronic Action Map ({material_name} | {method_lbl})</b>",
            template="plotly_white",
            legend=dict(orientation="h", yanchor="bottom", y=-0.12, xanchor="center", x=0.5)
        )
        fig_p.write_html(output_html, include_plotlyjs=True)
        print(f"  [NAMD:2D-Map] Saved interactive HTML to: {output_html}")
    except Exception as e:
        print(f"  [NAMD:Warn] Failed to create Plotly HTML 2D map: {e}")


def export_2d_vibronic_action_map(
    vib2d_data,
    output_npz="2d_vibronic_action_map.npz",
    output_csv="2d_vibronic_action_projections.csv"
):
    """
    Exports 2D vibronic action matrices and 1D marginal projections to NPZ and CSV.
    """
    np.savez_compressed(
        output_npz,
        wavenumbers_acc_cm=vib2d_data["wavenumbers_acc_cm"],
        wavenumbers_prom_cm=vib2d_data["wavenumbers_prom_cm"],
        S_2d_total=vib2d_data["S_2d_total"],
        S_2d_elec=vib2d_data["S_2d_elec"],
        S_2d_hole=vib2d_data["S_2d_hole"],
        proj_acc=vib2d_data["proj_acc"],
        proj_prom=vib2d_data["proj_prom"],
        method=str(vib2d_data.get("method", "NAMD"))
    )
    print(f"  [NAMD:2D-Map] Exported 2D vibronic arrays to: {output_npz}")

    # Export 1D marginal projections to CSV
    # Pad to equal length if needed
    len_acc = len(vib2d_data["wavenumbers_acc_cm"])
    len_prom = len(vib2d_data["wavenumbers_prom_cm"])
    max_len = max(len_acc, len_prom)

    w_acc_col = np.pad(vib2d_data["wavenumbers_acc_cm"], (0, max_len - len_acc), constant_values=np.nan)
    p_acc_col = np.pad(vib2d_data["proj_acc"], (0, max_len - len_acc), constant_values=np.nan)
    w_prom_col = np.pad(vib2d_data["wavenumbers_prom_cm"], (0, max_len - len_prom), constant_values=np.nan)
    p_prom_col = np.pad(vib2d_data["proj_prom"], (0, max_len - len_prom), constant_values=np.nan)

    np.savetxt(
        output_csv,
        np.column_stack([w_acc_col, p_acc_col, w_prom_col, p_prom_col]),
        header="wavenumber_acc_cm,action_acc,wavenumber_prom_cm,action_prom",
        delimiter=",",
        comments=""
    )
    print(f"  [NAMD:2D-Map] Exported 1D marginal projections to: {output_csv}")




