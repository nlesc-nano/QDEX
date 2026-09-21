"""
QDEX Non-Adiabatic Molecular Dynamics (NAMD) Module
"""

from qdex.namd.precompute import precompute_namd_data, compact_precomputed_data
from qdex.namd.surface_hopping import run_namd_dynamics

__all__ = ["precompute_namd_data", "compact_precomputed_data", "run_namd_dynamics"]

