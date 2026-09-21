r"""
QDEX Ultrafast Pump-Probe Transient Absorption (TA) Spectroscopy Module.

Calculates time-resolved differential absorption spectra :math:`\Delta A(E, t)` from
NAMD non-adiabatic dynamics trajectories, modeling:

#. Ground-State Bleach (GSB) via dynamic state-filling (Pauli blocking) of conduction
   and valence orbitals.
#. Stimulated Emission (SE) from populated excited states.
#. Band-edge 1S bleach kinetic rise profiling to extract carrier cooling rates :math:`k_C`.
#. 2D false-color pump-probe maps and spectral slices at specified delay times.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Dict, Any, Tuple, List
from scipy.optimize import curve_fit


def compute_transient_absorption(
    times_fs: np.ndarray,
    populations: np.ndarray,
    E_pairs: np.ndarray,
    f_pairs: np.ndarray,
    i_pairs: np.ndarray,
    a_pairs: np.ndarray,
    sigma_ev: float = 0.03,
    e_range: Optional[Tuple[float, float]] = None,
    n_e_points: int = 300,
    include_se: bool = True,
    all_energies: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Computes time-resolved differential absorption spectra Delta A(E, t).

    Parameters:
      times_fs: (n_times,) array of trajectory delay times in fs.
      populations: (n_times, n_pairs) array of exciton/pair populations over time.
      E_pairs: (n_pairs,) initial pair transition energies in eV, or (n_times, n_pairs).
      f_pairs: (n_pairs,) ground-state oscillator strengths.
      i_pairs: (n_pairs,) occupied MO indices for each pair.
      a_pairs: (n_pairs,) virtual MO indices for each pair.
      sigma_ev: Gaussian spectral line broadening in eV (default: 0.03 eV).
      e_range: (E_min, E_max) probe photon energy window in eV.
      n_e_points: Number of probe energy grid points.
      include_se: Whether to include stimulated emission (default: True).
      all_energies: Optional (n_times, n_pairs) instantaneous pair energies.

    Returns:
      dict with keys:
        'probe_energies_ev': (n_e_points,) array of probe photon energies.
        'times_fs': (n_times,) array of delay times.
        'delta_A': (n_times, n_e_points) 2D differential absorption matrix.
        'delta_A_1s': (n_times,) band-edge 1S bleach kinetic trace.
        'e_1s_ev': Energy of the 1S band-edge transition in eV.
        'fit_results': Dictionary with 1S rise time tau_rise_fs and rate k_c_ps.
        'p_1s': (n_times,) population of the lowest band-edge state.
    """
    times_arr = np.asarray(times_fs, dtype=np.float64)
    pop_arr = np.asarray(populations, dtype=np.float64)
    n_times, n_pairs = pop_arr.shape

    # Determine instantaneous energies E_mat (n_times, n_pairs)
    if all_energies is not None and all_energies.shape == (n_times, n_pairs):
        E_mat = np.asarray(all_energies, dtype=np.float64)
    elif E_pairs.ndim == 2 and E_pairs.shape == (n_times, n_pairs):
        E_mat = np.asarray(E_pairs, dtype=np.float64)
    else:
        E_base = np.asarray(E_pairs, dtype=np.float64)
        if E_base.ndim == 1:
            E_mat = np.tile(E_base, (n_times, 1))
        else:
            E_mat = E_base

    f_base = np.asarray(f_pairs, dtype=np.float64)
    if f_base.ndim == 2:
        f_base = f_base[0]

    i_arr = np.asarray(i_pairs, dtype=int)
    a_arr = np.asarray(a_pairs, dtype=int)

    # Unique active orbitals
    n_occ_max = int(np.max(i_arr)) + 1
    n_virt_max = int(np.max(a_arr)) + 1

    # 1. Probe Energy Grid
    e_min_data = float(np.min(E_mat))
    e_max_data = float(np.max(E_mat))
    if e_range is not None:
        e_grid_min = float(e_range[0])
        e_grid_max = float(e_range[1])
    else:
        e_grid_min = max(0.5, e_min_data - 0.4)
        e_grid_max = e_max_data + 0.4

    probe_energies = np.linspace(e_grid_min, e_grid_max, n_e_points)

    # 2. Identify 1S Band-Edge State (lowest energy pair)
    e0_pairs = E_mat[0]
    idx_1s = int(np.argmin(e0_pairs))
    e_1s_ev = float(e0_pairs[idx_1s])

    # 3. Compute Dynamic Differential Oscillator Strengths Delta f_ia(t)
    delta_A = np.zeros((n_times, n_e_points), dtype=np.float64)
    p_1s = pop_arr[:, idx_1s]

    inv_sqrt2pi_sigma = 1.0 / (np.sqrt(2.0 * np.pi) * max(sigma_ev, 1e-4))
    two_sigma_sq = 2.0 * (max(sigma_ev, 1e-4) ** 2)

    for k in range(n_times):
        p_k = pop_arr[k]  # (n_pairs,)
        E_k = E_mat[k]    # (n_pairs,)

        # Compute single-particle occupations from pair populations
        # n_a = sum_i P_ia (fractional occupation in virtual orbital a)
        # p_i = sum_a P_ia (fractional hole occupation in occupied orbital i)
        n_a = np.bincount(a_arr, weights=p_k, minlength=n_virt_max)
        p_i = np.bincount(i_arr, weights=p_k, minlength=n_occ_max)

        # Ground-State Bleach (GSB): Pauli blocking reduces absorption
        # Delta f_ia^GSB = - f_0 * (n_a + p_i)
        blocking_factor = n_a[a_arr] + p_i[i_arr]

        if include_se:
            # Stimulated Emission (SE): populated pairs emit coherently
            # Delta f_ia^SE = - f_0 * P_ia
            delta_f = - f_base * (blocking_factor + p_k)
        else:
            delta_f = - f_base * blocking_factor

        # Gaussian spectral convolution
        # delta_A[k, :] = sum_ia delta_f_ia * G(E - E_ia)
        # Vectorized over probe grid, evaluated only on pairs with non-negligible |delta_f|
        active_mask = np.abs(delta_f) > 1e-12
        if np.any(active_mask):
            E_active = E_k[active_mask]
            df_active = delta_f[active_mask]
            diff_sq = (probe_energies[np.newaxis, :] - E_active[:, np.newaxis]) ** 2
            gauss_profiles = inv_sqrt2pi_sigma * np.exp(- diff_sq / two_sigma_sq)
            delta_A[k, :] = np.dot(df_active, gauss_profiles)
        else:
            delta_A[k, :] = 0.0

    # 4. Extract 1S Bleach Kinetic Profile
    # Find probe grid point closest to the 1S transition
    idx_probe_1s = int(np.argmin(np.abs(probe_energies - e_1s_ev)))
    # Or track the local minimum near e_1s
    delta_A_1s = delta_A[:, idx_probe_1s]

    # 5. Fit 1S Bleach Rise Kinetics to Extract Cooling Rate k_C
    fit_results = fit_bleach_rise_kinetics(times_arr, delta_A_1s)

    return {
        "probe_energies_ev": probe_energies,
        "times_fs": times_arr,
        "delta_A": delta_A,
        "delta_A_1s": delta_A_1s,
        "e_1s_ev": e_1s_ev,
        "idx_probe_1s": idx_probe_1s,
        "fit_results": fit_results,
        "p_1s": p_1s,
        "sigma_ev": sigma_ev,
    }


def fit_bleach_rise_kinetics(times_fs: np.ndarray, delta_A_1s: np.ndarray) -> Dict[str, Any]:
    """
    Fits the negative 1S bleach rise profile -Delta A_1S(t) to an exponential rise::

        S(t) = A_0 * (1 - exp(- t / tau_rise)) + offset

    extracting the hot-carrier cooling / 1S arrival lifetime tau_rise (in fs/ps)
    and cooling rate k_C = 1 / tau_rise (in ps^-1), matching Figure 2a,c of the manuscript.
    """
    times = np.asarray(times_fs, dtype=np.float64)
    # Bleach is negative; invert so S(t) is positive rising signal
    S = - np.asarray(delta_A_1s, dtype=np.float64)

    if len(times) < 4 or np.max(S) - np.min(S) < 1e-12:
        return {
            "tau_rise_fs": np.nan,
            "tau_rise_ps": np.nan,
            "k_cool_ps": np.nan,
            "A0": 0.0,
            "fit_curve": np.zeros_like(S),
            "success": False
        }

    # Normalize for robust fitting
    S_max = float(np.max(S))
    S_min = float(S[0])
    delta_S = max(S_max - S_min, 1e-10)

    # Initial guess for tau: time where S reaches 63% of maximum
    target = S_min + 0.632 * delta_S
    idx_guess = np.where(S >= target)[0]
    tau_guess = float(times[idx_guess[0]]) if len(idx_guess) > 0 else (0.3 * times[-1])
    tau_guess = max(10.0, min(tau_guess, times[-1]))

    def rise_model(t, a0, tau, c0):
        return a0 * (1.0 - np.exp(- np.maximum(t, 0.0) / max(tau, 1e-3))) + c0

    p0 = [delta_S, tau_guess, S_min]
    bounds = (
        [0.0, 1.0, -np.inf],
        [10.0 * delta_S, 10.0 * times[-1], np.inf]
    )

    try:
        popt, _ = curve_fit(rise_model, times, S, p0=p0, bounds=bounds, maxfev=5000)
        a0_fit, tau_fit, c0_fit = popt
        fit_curve = - rise_model(times, a0_fit, tau_fit, c0_fit)  # Return in original Delta A sign
        tau_rise_fs = float(tau_fit)
        tau_rise_ps = tau_rise_fs * 1.0e-3
        k_cool_ps = 1.0 / tau_rise_ps if tau_rise_ps > 0 else np.nan
        success = True
    except Exception:
        # Fallback to analytical 1/e estimation
        tau_rise_fs = float(tau_guess)
        tau_rise_ps = tau_rise_fs * 1.0e-3
        k_cool_ps = 1.0 / tau_rise_ps if tau_rise_ps > 0 else np.nan
        fit_curve = delta_A_1s
        success = False

    return {
        "tau_rise_fs": tau_rise_fs,
        "tau_rise_ps": tau_rise_ps,
        "k_cool_ps": k_cool_ps,
        "fit_curve": fit_curve,
        "success": success
    }


def plot_transient_absorption(
    ta_data: Dict[str, Any],
    plot_file: str = "transient_absorption_map.png",
    delay_slices_fs: Optional[List[float]] = None,
    material_name: str = "CSPBBR3",
) -> None:
    r"""
    Generates a publication-quality 3-panel Transient Absorption dashboard:

    - **Panel (a)**: 2D False-Color TA Map :math:`\Delta A(E, t)` (Probe Energy vs Delay Time).
    - **Panel (b)**: 1S Bleach Kinetic Rise Trace with exponential fit and cooling rate :math:`k_C`.
    - **Panel (c)**: Differential Absorption Spectra :math:`\Delta A(E)` at selected delay times.
    """
    probe_e = ta_data["probe_energies_ev"]
    times = ta_data["times_fs"]
    delta_A = ta_data["delta_A"]
    delta_A_1s = ta_data["delta_A_1s"]
    e_1s = ta_data["e_1s_ev"]
    fit = ta_data["fit_results"]

    fig, axes = plt.subplots(1, 3, figsize=(20, 6), gridspec_kw={"width_ratios": [1.2, 1.0, 1.0]})

    # -----------------------------------------------------------------
    # Panel (a): 2D False-Color Transient Absorption Map
    # -----------------------------------------------------------------
    ax1 = axes[0]
    # Symmetrize colormap limits around 0
    max_abs = float(np.max(np.abs(delta_A))) if np.max(np.abs(delta_A)) > 0 else 1.0
    vmin = - max_abs
    vmax = max_abs

    # Use RdBu_r: Blue = Negative Bleach (GSB/SE), Red = Positive (ESA)
    cax = ax1.pcolormesh(
        probe_e, times, delta_A,
        cmap="RdBu_r", vmin=vmin, vmax=vmax, shading="auto"
    )
    cb = fig.colorbar(cax, ax=ax1, orientation="vertical", pad=0.02)
    cb.set_label(r"Differential Absorbance $\Delta A$ (arb. units)", fontsize=10, fontweight="bold")

    ax1.axvline(e_1s, color="#111827", linestyle="--", lw=1.5, alpha=0.8, label=f"1S Band Edge ({e_1s:.2f} eV)")
    ax1.set_xlabel("Probe Photon Energy (eV)", fontsize=11, fontweight="bold")
    ax1.set_ylabel("Pump–Probe Delay Time (fs)", fontsize=11, fontweight="bold")
    ax1.set_title(r"(a) Ultrafast Transient Absorption 2D Map $\Delta A(E, t)$", fontsize=12, fontweight="bold")
    ax1.legend(loc="upper right", frameon=True, fontsize=9)

    # -----------------------------------------------------------------
    # Panel (b): 1S Bleach Kinetic Rise Profile (Figure 2a/c Analog)
    # -----------------------------------------------------------------
    ax2 = axes[1]
    # Plot -Delta A_1S so bleach appears as a rising positive signal
    bleach_signal = - delta_A_1s
    fit_signal = - fit["fit_curve"]

    ax2.scatter(times, bleach_signal, color="#2563EB", s=25, alpha=0.8, edgecolors="none", label=r"1S Bleach $-\Delta A_{1S}(t)$")
    if fit.get("success", False):
        lbl_fit = f"Rise Fit ($\\tau = {fit['tau_rise_fs']:.1f}$ fs, $k_C = {fit['k_cool_ps']:.2f}$ ps$^{{-1}}$)"
        ax2.plot(times, fit_signal, color="#DC2626", lw=2.4, label=lbl_fit)

    ax2.set_xlabel("Pump–Probe Delay Time (fs)", fontsize=11, fontweight="bold")
    ax2.set_ylabel(r"Bleach Amplitude $-\Delta A_{1S}(t)$ (arb. units)", fontsize=11, fontweight="bold")
    ax2.set_title(r"(b) 1S Bleach Rise Dynamics (Carrier Cooling)", fontsize=12, fontweight="bold")
    ax2.legend(loc="lower right", frameon=True, fontsize=9)
    ax2.grid(True, linestyle="--", alpha=0.4)

    # Annotation box with cooling metrics
    text_info = (
        f"1S Energy: {e_1s:.3f} eV\n"
        f"Rise Time $\\tau_C$: {fit['tau_rise_fs']:.1f} fs ({fit['tau_rise_ps']:.3f} ps)\n"
        f"Cooling Rate $k_C$: {fit['k_cool_ps']:.2f} ps$^{{-1}}$"
    )
    ax2.text(
        0.05, 0.95, text_info,
        transform=ax2.transAxes, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8FAFC", edgecolor="#CBD5E1", alpha=0.9),
        fontsize=9.5, fontweight="bold"
    )

    # -----------------------------------------------------------------
    # Panel (c): Spectral Slices at Selected Delay Times
    # -----------------------------------------------------------------
    ax3 = axes[2]
    if delay_slices_fs is None:
        # Pick 5 representative delays along the trajectory
        t_max = float(times[-1])
        delay_slices_fs = [0.0, 0.1 * t_max, 0.25 * t_max, 0.5 * t_max, t_max]

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(delay_slices_fs)))
    for d_t, col in zip(delay_slices_fs, colors):
        idx_t = int(np.argmin(np.abs(times - d_t)))
        actual_t = times[idx_t]
        ax3.plot(probe_e, delta_A[idx_t, :], color=col, lw=2.0, label=f"t = {actual_t:.1f} fs")

    ax3.axhline(0.0, color="#6B7280", linestyle=":", lw=1.2)
    ax3.axvline(e_1s, color="#111827", linestyle="--", lw=1.2, alpha=0.7)
    ax3.set_xlabel("Probe Photon Energy (eV)", fontsize=11, fontweight="bold")
    ax3.set_ylabel(r"$\Delta A(E)$ (arb. units)", fontsize=11, fontweight="bold")
    ax3.set_title(r"(c) Transient Spectra at Selected Delays", fontsize=12, fontweight="bold")
    ax3.legend(loc="lower right", frameon=True, fontsize=9)
    ax3.grid(True, linestyle="--", alpha=0.4)

    plt.suptitle(
        f"QDEX Ultrafast Pump–Probe Transient Absorption Spectroscopy: {material_name.upper()}",
        fontsize=14, fontweight="bold", y=0.98
    )
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()
    print(f"  [Transient Absorption] Publication dashboard saved to: {plot_file}")


def export_transient_absorption_data(
    ta_data: Dict[str, Any],
    kinetics_csv: str = "ta_bleach_kinetics.csv",
    map_npz: Optional[str] = "ta_2d_map.npz",
) -> None:
    """Exports kinetic bleach traces to CSV and full 2D map to compressed NPZ."""
    times = ta_data["times_fs"]
    delta_A_1s = ta_data["delta_A_1s"]
    fit = ta_data["fit_results"]

    df = pd.DataFrame({
        "time_fs": times,
        "time_ps": times * 1.0e-3,
        "delta_A_1s": delta_A_1s,
        "bleach_amplitude": - delta_A_1s,
        "fit_curve": fit.get("fit_curve", np.zeros_like(delta_A_1s)),
        "pop_1s": ta_data.get("p_1s", np.zeros_like(delta_A_1s))
    })
    df.to_csv(kinetics_csv, index=False)
    print(f"  [Transient Absorption] 1S Bleach kinetics exported to: {kinetics_csv}")

    if map_npz:
        np.savez_compressed(
            map_npz,
            probe_energies_ev=ta_data["probe_energies_ev"],
            times_fs=ta_data["times_fs"],
            delta_A=ta_data["delta_A"],
            e_1s_ev=ta_data["e_1s_ev"],
            tau_rise_fs=fit.get("tau_rise_fs", np.nan),
            k_cool_ps=fit.get("k_cool_ps", np.nan)
        )
        print(f"  [Transient Absorption] 2D Map array exported to: {map_npz}")
