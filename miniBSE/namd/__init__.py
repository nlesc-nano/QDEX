"""
miniBSE Non-Adiabatic Molecular Dynamics (NAMD) Module
"""

from miniBSE.namd.precompute import precompute_namd_data, compact_precomputed_data
from miniBSE.namd.surface_hopping import run_namd_dynamics

__all__ = ["precompute_namd_data", "compact_precomputed_data", "run_namd_dynamics"]

