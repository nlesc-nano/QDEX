import numpy as np
import pytest
from qdex.orbital_analysis import compute_orbital_centroids, print_orbital_summary


def test_compute_orbital_centroids_core_vs_trap():
    # Model system: 1 center atom (Cd, index 0 at [0, 0, 0]),
    # 4 core atoms (Se, indices 1..4 at tetrahedral positions at radius R=5.0 A),
    # 4 surface ligand atoms (Cl, indices 5..8 at radius R=8.0 A).
    coords = np.array([
        [0.0, 0.0, 0.0],       # 0: Cd (center)
        [2.88, 2.88, 2.88],    # 1: Se
        [-2.88, -2.88, 2.88],  # 2: Se
        [-2.88, 2.88, -2.88],  # 3: Se
        [2.88, -2.88, -2.88],  # 4: Se
        [8.0, 0.0, 0.0],       # 5: Cl (surface +x)
        [-8.0, 0.0, 0.0],      # 6: Cl (surface -x)
        [0.0, 8.0, 0.0],       # 7: Cl (surface +y)
        [0.0, -8.0, 0.0],      # 8: Cl (surface -y)
    ])
    syms = ["Cd", "Se", "Se", "Se", "Se", "Cl", "Cl", "Cl", "Cl"]

    # One s-orbital per atom (n_ao = 9)
    shells = [{"sym": syms[i], "atom_idx": i, "l": 0} for i in range(len(syms))]

    # State 0: S-like envelope delocalized symmetrically across core atoms (0..4)
    # State 1: Localized surface trap on atom 5 (Cl at [8.0, 0, 0])
    # State 2: Mixed state (partially on core, partially on surface)
    pops = np.zeros((9, 3))
    # State 0: 20% on each of 5 core atoms
    pops[0:5, 0] = 0.20
    # State 1: 95% on atom 5 (Cl), 5% elsewhere
    pops[5, 1] = 0.95
    pops[0, 1] = 0.05
    # State 2: 50% on atom 1 (Se), 50% on atom 5 (Cl)
    pops[1, 2] = 0.50
    pops[5, 2] = 0.50

    info = compute_orbital_centroids(pops, shells, coords, syms)

    # Core COM should be at [0, 0, 0]
    assert np.allclose(info["com_core"], [0.0, 0.0, 0.0], atol=1e-3)
    assert np.isclose(info["r_core"], np.linalg.norm(coords[1]), atol=1e-2)

    # State 0: centered, low xi, high f_core -> [Core]
    assert info["xi"][0] < 0.10
    assert info["f_core"][0] == 1.0
    assert info["labels"][0] == "[Core]"

    # State 1: located at +x surface, high xi, low f_core -> [Surf/Trap]
    assert info["xi"][1] > 1.0
    assert info["f_core"][1] < 0.10
    assert info["labels"][1] == "[Surf/Trap]"

    # State 2: intermediate -> [Mixed] or [Surf/Trap]
    assert info["labels"][2] in ("[Mixed]", "[Surf/Trap]")


def test_print_orbital_summary_with_coords(caplog):
    import logging
    coords = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [-2.0, 0.0, 0.0],
    ])
    syms = ["Cd", "Se", "Cl"]
    shells = [{"sym": syms[i], "atom_idx": i, "l": 0} for i in range(len(syms))]
    pops = np.eye(3)
    energies_eV = np.array([-2.0, -1.0, 0.5])
    occ = np.array([2.0, 2.0, 0.0])

    with caplog.at_level(logging.INFO):
        info = print_orbital_summary(
            energies_eV, occ, homo_idx=1, pops=pops, syms=syms, shells=shells,
            coords_ang=coords, print_range=2
        )
    assert info is not None
    assert len(info["labels"]) == 3
