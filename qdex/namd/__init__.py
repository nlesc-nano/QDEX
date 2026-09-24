"""
QDEX Non-Adiabatic Molecular Dynamics (NAMD) Module
"""

from qdex.namd.precompute import precompute_namd_data, compact_precomputed_data
from qdex.namd.surface_hopping import run_namd_dynamics
from qdex.namd.transient_absorption import (
    compute_transient_absorption,
    fit_bleach_rise_kinetics,
    plot_transient_absorption,
    export_transient_absorption_data,
)
from qdex.namd.analysis import (
    load_trajectory_orbital_energies,
    compute_pair_spectral_density,
    compute_time_resolved_spectral_density,
    plot_time_resolved_spectral_density,
    export_time_resolved_spectral_density,
    compute_2d_vibronic_action_map,
    plot_2d_vibronic_action_map,
    export_2d_vibronic_action_map,
)
from qdex.namd.ensemble import (
    compute_energy_autocorrelation_time,
    auto_calibrate_ensemble_origins,
    sample_origin_initial_states,
    aggregate_multi_origin_results,
)

__all__ = [
    "precompute_namd_data",
    "compact_precomputed_data",
    "run_namd_dynamics",
    "compute_transient_absorption",
    "fit_bleach_rise_kinetics",
    "plot_transient_absorption",
    "export_transient_absorption_data",
    "load_trajectory_orbital_energies",
    "compute_pair_spectral_density",
    "compute_time_resolved_spectral_density",
    "plot_time_resolved_spectral_density",
    "export_time_resolved_spectral_density",
    "compute_2d_vibronic_action_map",
    "plot_2d_vibronic_action_map",
    "export_2d_vibronic_action_map",
    "compute_energy_autocorrelation_time",
    "auto_calibrate_ensemble_origins",
    "sample_origin_initial_states",
    "aggregate_multi_origin_results",
]


