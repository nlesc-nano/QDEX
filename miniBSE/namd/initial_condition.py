import numpy as np

def sample_initial_states(
    energies,
    f_osc,
    pump_energy_ev,
    pulse_fwhm_ev=0.08,
    n_trajectories=1000,
    filter_dark_states=True,
    min_f_thresh=1e-5,
    random_seed=None
):
    """
    Samples initial exciton states based on laser pulse excitation:
      P_I(0) ~ f_I * exp(-(E_I - E_pump)^2 / (2 * sigma^2))
    
    Parameters:
      energies: 1D array of exciton energies (eV)
      f_osc: 1D array of oscillator strengths
      pump_energy_ev: center frequency of the laser pulse (e.g. 2 * Eg)
      pulse_fwhm_ev: full-width at half-maximum of the laser pulse (eV)
      n_trajectories: number of trajectories to sample
      filter_dark_states: whether to weight by oscillator strength
      min_f_thresh: threshold below which oscillator strength is treated as dark
      random_seed: optional integer for reproducible sampling
    
    Returns:
      sampled_states: 1D array of length n_trajectories containing state indices
      P_norm: normalized probability distribution across all states
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    energies = np.asarray(energies, dtype=np.float64)
    f_osc = np.asarray(f_osc, dtype=np.float64)

    sigma = pulse_fwhm_ev / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    spectral_profile = np.exp(-0.5 * ((energies - pump_energy_ev) / sigma) ** 2)

    if filter_dark_states:
        weights = np.maximum(f_osc, 0.0) * spectral_profile
        weights[weights < min_f_thresh * np.max(spectral_profile)] = 0.0
    else:
        weights = spectral_profile

    total_weight = np.sum(weights)
    if total_weight <= 1e-15:
        # Fallback: choose the closest state in energy to pump_energy_ev
        idx_closest = np.argmin(np.abs(energies - pump_energy_ev))
        P_norm = np.zeros_like(energies)
        P_norm[idx_closest] = 1.0
        sampled_states = np.full(n_trajectories, idx_closest, dtype=int)
        return sampled_states, P_norm

    P_norm = weights / total_weight
    state_indices = np.arange(len(energies))
    sampled_states = np.random.choice(state_indices, size=n_trajectories, p=P_norm)

    return sampled_states, P_norm
