import numpy as np
import pytest
from qdex.exciton_hamiltonian import ExcitonHamiltonian
from qdex.orbital_analysis import compute_spin_character


def test_trap_filter_spin_free():
    # 4 MOs: 2 occupied (0, 1), 2 virtual (2, 3)
    # homo_index = 1
    n_ao = 4
    C = np.eye(n_ao)
    eps = np.array([-4.0, -2.0, 1.0, 3.0])  # gaps: (1-(-2)=3), (3-(-2)=5), (1-(-4)=5), (3-(-4)=7)
    overlap = np.eye(n_ao)
    atom_ao_ranges = [(0, 1), (1, 2), (2, 3), (3, 4)]
    
    # Flag HOMO (index 1 in global, which is index 1 in occ_idx [0, 1]) as trap
    hole_trap_mask = np.array([False, True])  # hole 0 is core, hole 1 is trap
    elec_trap_mask = np.array([False, False])

    # Unfiltered Hamiltonian
    ham_full = ExcitonHamiltonian(
        C=C, eps=eps, overlap=overlap, atom_ao_ranges=atom_ao_ranges,
        homo_index=1, n_occ=2, n_virt=2, scissor_ev=0.0,
        gamma_qp=np.ones((4, 4)), gamma_bse=np.ones((4, 4)),
        excitation_mode="diagonal_bse",
    )
    assert ham_full.dim == 4  # 2 occ x 2 virt = 4 transitions

    # Filtered Hamiltonian: transitions from hole 1 are excluded!
    ham_filtered = ExcitonHamiltonian(
        C=C, eps=eps, overlap=overlap, atom_ao_ranges=atom_ao_ranges,
        homo_index=1, n_occ=2, n_virt=2, scissor_ev=0.0,
        gamma_qp=np.ones((4, 4)), gamma_bse=np.ones((4, 4)),
        excitation_mode="diagonal_bse",
        hole_trap_mask=hole_trap_mask, elec_trap_mask=elec_trap_mask,
    )
    assert ham_filtered.dim == 2  # Only hole 0 -> virt 0, virt 1 remain
    assert all(ham_filtered.valid_i == 0)


def test_trap_filter_soc_and_spin_character():
    # 2 spatial MOs: 1 occ, 1 virt -> 2 occupied spinors, 2 virtual spinors
    n_ao = 2
    C = np.eye(n_ao)
    eps = np.array([-2.0, 2.0])
    overlap = np.eye(n_ao)
    atom_ao_ranges = [(0, 1), (1, 2)]

    # 4 spinors: 2 occ (indices 0, 1), 2 virt (indices 2, 3)
    soc_E = np.array([-2.1, -1.9, 1.9, 2.1])
    soc_U = np.eye(4, dtype=complex)  # 2 alpha + 2 beta basis

    # Flag spinor 1 (second occupied spinor) as trap
    hole_trap_mask_soc = np.array([False, True])
    elec_trap_mask_soc = np.array([False, False])

    ham = ExcitonHamiltonian(
        C=C, eps=eps, overlap=overlap, atom_ao_ranges=atom_ao_ranges,
        homo_index=0, n_occ=1, n_virt=1, scissor_ev=0.0,
        gamma_qp=np.ones((2, 2)), gamma_bse=np.ones((2, 2)),
        include_exchange=False, include_direct_eh=False,
        soc_U=soc_U, soc_E=soc_E,
        hole_trap_mask_soc=hole_trap_mask_soc, elec_trap_mask_soc=elec_trap_mask_soc,
    )
    # Total spinor transitions: 2 occ x 2 virt = 4.
    # Excluding hole spinor 1 removes 2 transitions -> exactly 2 transitions remain
    assert ham.dim == 2
    assert ham.valid_spinor_mask.sum() == 2

    # Verify compute_spin_character with valid_spinor_mask works without dimension mismatch
    dummy_vec = np.array([1.0, 0.0], dtype=complex)
    s_pct, t_pct = compute_spin_character(dummy_vec, soc_U, ham.n_occ_spinor, ham.n_virt_spinor, valid_mask=ham.valid_spinor_mask)
    assert np.isclose(s_pct + t_pct, 100.0)
