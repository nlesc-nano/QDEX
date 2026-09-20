import unittest

import numpy as np

from miniBSE.exciton_hamiltonian import ExcitonHamiltonian
from miniBSE.solver import ExcitonSolver


def make_two_transition_hamiltonian():
    ham = object.__new__(ExcitonHamiltonian)
    ham.dim = 2
    ham.soc_flag = False
    ham.spin = "singlet"
    ham.D_dft = np.array([2.0, 3.0])
    ham.D = np.array([2.8, 3.8])
    ham.D_spatial = ham.D
    ham.gamma = np.array([[4.0, 0.5], [0.5, 3.0]])
    ham.q_flat = np.array([[0.20, -0.20], [0.10, -0.10]])
    ham.include_direct_eh = True
    ham.n_occ_act = 1
    ham.n_virt_act = 2
    ham.valid_i = np.array([0, 0])
    ham.valid_a = np.array([0, 1])
    ham.q_occ = np.array([[[0.6, 0.4]]])
    ham.W_virt = np.zeros((2, 2, 2))
    ham.W_virt[0, 0] = [1.1, 0.7]
    ham.W_virt[1, 1] = [0.8, 0.9]
    ham.W_virt[0, 1] = ham.W_virt[1, 0] = [0.2, 0.1]
    return ham


class IndependentTransitionModeTests(unittest.TestCase):
    def test_kernel_separation_and_diagonal_approximation(self):
        ham = make_two_transition_hamiltonian()
        basis = np.eye(ham.dim)
        dense = np.column_stack([ham.matvec(basis[:, i]) for i in range(ham.dim)])
        energies, kx, kd = ham.independent_transition_energies("diagonal_bse")
        np.testing.assert_allclose(energies, np.diag(dense))
        np.testing.assert_allclose(energies, ham.D + kx - kd)

    def test_dft_and_qp_modes_do_not_diagonalize(self):
        ham = make_two_transition_hamiltonian()
        solver = object.__new__(ExcitonSolver)
        solver.ham = ham
        solver.soc_flag = False

        e_dft, x_dft = solver.solve(nroots=2, excitation_mode="independent_dft")
        e_qp, x_qp = solver.solve(nroots=2, excitation_mode="independent_qp")
        np.testing.assert_allclose(e_dft, [2.0, 3.0])
        np.testing.assert_allclose(e_qp, [2.8, 3.8])
        np.testing.assert_allclose(x_dft, np.eye(2))
        np.testing.assert_allclose(x_qp, np.eye(2))

    def test_dense_and_matrix_free_are_equivalent(self):
        ham = make_two_transition_hamiltonian()
        basis = np.eye(ham.dim)
        dense = np.column_stack([ham.matvec(basis[:, i]) for i in range(ham.dim)])
        probe = np.array([0.3, -0.7])
        np.testing.assert_allclose(ham.matvec(probe), dense @ probe)
        np.testing.assert_allclose(dense, dense.T, atol=1.0e-12)

    def test_diagonal_bse_vectors_and_oscillator_strengths(self):
        from miniBSE.solver import DiagonalBSEVectors
        from miniBSE.oscillator import compute_oscillator_strengths

        dim = 1000
        order = np.arange(dim)[::-1]
        vecs = DiagonalBSEVectors(dim, order, dtype=complex)

        self.assertEqual(vecs.shape, (dim, dim))
        # 1D column slicing
        v0 = vecs[:, 0]
        self.assertEqual(v0.shape, (dim,))
        self.assertEqual(v0[order[0]], 1.0)
        self.assertEqual(np.sum(np.abs(v0)), 1.0)

        # 2D submatrix slicing
        v_sub = vecs[:, :10]
        self.assertEqual(v_sub.shape, (dim, 10))
        np.testing.assert_allclose(v_sub[order[:10], np.arange(10)], 1.0)

        # Matmul via transpose
        mu_ia = np.random.randn(dim, 3)
        mu_trans = vecs.T @ mu_ia
        np.testing.assert_allclose(mu_trans, mu_ia[order])

        # Oscillator strengths
        energies_ev = np.linspace(1.0, 3.0, dim)
        f_diag = compute_oscillator_strengths(energies_ev, vecs, mu_ia, is_spinor=False)
        self.assertEqual(len(f_diag), dim)
        # Check first state manually: prefactor=4/3, HA_TO_EV
        from miniBSE.constants import HA_TO_EV
        expected_f0 = (4.0 / 3.0) * (energies_ev[0] / HA_TO_EV) * np.sum(mu_ia[order[0]] ** 2)
        self.assertAlmostEqual(f_diag[0], expected_f0, places=10)

    def test_diagonal_bse_expectation_components(self):
        ham = make_two_transition_hamiltonian()
        solver = object.__new__(ExcitonSolver)
        solver.ham = ham
        solver.soc_flag = False

        e, vecs = solver.solve(nroots=2, excitation_mode="diagonal_bse")
        d_qp, kx, minus_kd = solver.expectation_components(vecs[:, 0])
        self.assertAlmostEqual(d_qp + kx + minus_kd, e[0], places=10)

    def test_diagonal_bse_low_memory_path(self):
        np.random.seed(42)
        n_ao, n_occ, n_virt, n_atoms = 15, 3, 3, 2
        atom_ranges = [(0, 7), (7, 15)]

        S = np.random.randn(n_ao, n_ao)
        S = S.T @ S + np.eye(n_ao)
        C = np.random.randn(n_ao, 10)
        C, _ = np.linalg.qr(C)
        eps = np.linspace(-5.0, 5.0, 10)

        gamma = np.random.randn(n_atoms, n_atoms)
        gamma = gamma.T @ gamma

        k = n_occ + n_virt
        soc_U = np.random.randn(2*k, 2*k) + 1j * np.random.randn(2*k, 2*k)
        soc_U, _ = np.linalg.qr(soc_U)
        soc_E = np.linspace(-4.0, 4.0, 2*k)

        # Spatial test
        ham_full = ExcitonHamiltonian(
            C=C, eps=eps, overlap=S, atom_ao_ranges=atom_ranges,
            homo_index=4, n_occ=n_occ, n_virt=n_virt, scissor_ev=1.0,
            gamma_qp=gamma, gamma_bse=gamma, gamma_bare=gamma,
            include_exchange=True, excitation_mode='bse'
        )
        ham_diag = ExcitonHamiltonian(
            C=C, eps=eps, overlap=S, atom_ao_ranges=atom_ranges,
            homo_index=4, n_occ=n_occ, n_virt=n_virt, scissor_ev=1.0,
            gamma_qp=gamma, gamma_bse=gamma, gamma_bare=gamma,
            include_exchange=True, excitation_mode='diagonal_bse'
        )

        self.assertFalse(hasattr(ham_diag, 'q_occ'))
        self.assertFalse(hasattr(ham_diag, 'q_virt'))
        self.assertFalse(hasattr(ham_diag, 'q_flat'))

        e_full, kx_full, kd_full = ham_full.independent_transition_energies('diagonal_bse')
        e_diag, kx_diag, kd_diag = ham_diag.independent_transition_energies('diagonal_bse')

        np.testing.assert_allclose(e_diag, e_full, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(kx_diag, kx_full, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(kd_diag, kd_full, rtol=1e-10, atol=1e-10)

        # SOC test
        ham_full_soc = ExcitonHamiltonian(
            C=C, eps=eps, overlap=S, atom_ao_ranges=atom_ranges,
            homo_index=4, n_occ=n_occ, n_virt=n_virt, scissor_ev=1.0,
            gamma_qp=gamma, gamma_bse=gamma, gamma_bare=gamma,
            include_exchange=True, soc_U=soc_U, soc_E=soc_E,
            excitation_mode='bse'
        )
        ham_diag_soc = ExcitonHamiltonian(
            C=C, eps=eps, overlap=S, atom_ao_ranges=atom_ranges,
            homo_index=4, n_occ=n_occ, n_virt=n_virt, scissor_ev=1.0,
            gamma_qp=gamma, gamma_bse=gamma, gamma_bare=gamma,
            include_exchange=True, soc_U=soc_U, soc_E=soc_E,
            excitation_mode='diagonal_bse'
        )

        self.assertFalse(hasattr(ham_diag_soc, 'q_occ'))
        self.assertFalse(hasattr(ham_diag_soc, 'q_virt'))
        self.assertFalse(hasattr(ham_diag_soc, 'q_hole_spinor'))
        self.assertFalse(hasattr(ham_diag_soc, 'q_elec_spinor'))
        self.assertFalse(hasattr(ham_diag_soc, 'q_spinor'))

        e_full_soc, kx_full_soc, kd_full_soc = ham_full_soc.independent_transition_energies('diagonal_bse')
        e_diag_soc, kx_diag_soc, kd_diag_soc = ham_diag_soc.independent_transition_energies('diagonal_bse')

        np.testing.assert_allclose(e_diag_soc, e_full_soc, rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(kx_diag_soc, kx_full_soc, rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(kd_diag_soc, kd_full_soc, rtol=1e-8, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
