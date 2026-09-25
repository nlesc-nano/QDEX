import unittest
import numpy as np

import libint_cpp
from qdex.integrals import compute_two_electron_ao
from qdex.hardness import build_xs_kernel, build_sbse_kernel, estimate_sgw_qp_gap
from qdex.exciton_hamiltonian import ExcitonHamiltonian
from qdex.solver import ExcitonSolver
from qdex.constants import HA_TO_EV, BOHR_PER_ANG


def create_test_h2_shells():
    """Create a minimal H2 system with s and p shells for unit testing."""
    # Two H atoms at (0, 0, 0) and (0, 0, 1.4) Bohr
    shells = [
        # Atom 0: s shell (1 AO)
        {
            'center': np.array([0.0, 0.0, 0.0], dtype=float),
            'l': 0,
            'exps': np.array([0.5, 0.2], dtype=float),
            'coefs': np.array([0.6, 0.4], dtype=float),
            'pure': True
        },
        # Atom 0: p shell (3 AOs)
        {
            'center': np.array([0.0, 0.0, 0.0], dtype=float),
            'l': 1,
            'exps': np.array([0.8, 0.3], dtype=float),
            'coefs': np.array([0.5, 0.5], dtype=float),
            'pure': True
        },
        # Atom 1: s shell (1 AO)
        {
            'center': np.array([0.0, 0.0, 1.4], dtype=float),
            'l': 0,
            'exps': np.array([0.5, 0.2], dtype=float),
            'coefs': np.array([0.6, 0.4], dtype=float),
            'pure': True
        },
        # Atom 1: p shell (3 AOs)
        {
            'center': np.array([0.0, 0.0, 1.4], dtype=float),
            'l': 1,
            'exps': np.array([0.8, 0.3], dtype=float),
            'coefs': np.array([0.5, 0.5], dtype=float),
            'pure': True
        }
    ]
    atom_symbols = ["H", "H"]
    coords_ang = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.4 / BOHR_PER_ANG]])
    atom_ao_ranges = [(0, 4), (4, 8)]
    return shells, atom_symbols, coords_ang, atom_ao_ranges


class TestXsKernel(unittest.TestCase):
    def setUp(self):
        self.shells, self.symbols, self.coords, self.ao_ranges = create_test_h2_shells()
        self.n_ao = 8

    def test_compute_two_electron_ao_properties(self):
        """Test analytical (mu mu | nu nu) integrals computed via Libint2."""
        gamma_au = compute_two_electron_ao(self.shells, nthreads=1)
        self.assertEqual(gamma_au.shape, (self.n_ao, self.n_ao))

        # 1. Symmetry: (mu mu | nu nu) = (nu nu | mu mu)
        np.testing.assert_allclose(gamma_au, gamma_au.T, atol=1e-14,
                                   err_msg="AO Coulomb matrix Gamma must be strictly symmetric")

        # 2. Strict Positivity of diagonal (self-repulsion)
        diag = np.diag(gamma_au)
        self.assertTrue(np.all(diag > 0.0), "Self-repulsion integrals (mu mu | mu mu) must be strictly positive")

        # 3. Positive semi-definiteness (Coulomb metric)
        eigvals = np.linalg.eigvalsh(gamma_au)
        self.assertTrue(np.all(eigvals >= -1e-12), f"Coulomb matrix eigenvalues must be non-negative: min={eigvals.min()}")

    def test_build_xs_kernel_uniform(self):
        """Test Xs-QDEX kernel builder with uniform screening."""
        alpha = 0.65
        w_bse, w_qp, gamma_bare = build_xs_kernel(
            shells=self.shells,
            atom_symbols=self.symbols,
            coords=self.coords,
            atom_ao_ranges=self.ao_ranges,
            material_name="DEFAULT",
            kernel_mode="uniform",
            alpha=alpha,
            nthreads=1
        )

        self.assertEqual(w_bse.shape, (self.n_ao, self.n_ao))
        self.assertEqual(gamma_bare.shape, (self.n_ao, self.n_ao))

        # w_bse must equal alpha * gamma_bare
        np.testing.assert_allclose(w_bse, alpha * gamma_bare, atol=1e-12)
        np.testing.assert_allclose(w_qp, w_bse, atol=1e-12)

    def test_build_xs_kernel_resta(self):
        """Test Xs-QDEX kernel builder with Resta microscopic screening."""
        w_bse, w_qp, gamma_bare = build_xs_kernel(
            shells=self.shells,
            atom_symbols=self.symbols,
            coords=self.coords,
            atom_ao_ranges=self.ao_ranges,
            material_name="CDSE",
            kernel_mode="resta",
            nthreads=1
        )

        self.assertEqual(w_bse.shape, (self.n_ao, self.n_ao))

        # Check on-site block (atom 0 with atom 0: AOs 0..3)
        # On-site Resta screening factor is 1.0 (unscreened atomic core)
        for a0, a1 in self.ao_ranges:
            np.testing.assert_allclose(
                w_bse[a0:a1, a0:a1], gamma_bare[a0:a1, a0:a1], atol=1e-12,
                err_msg="On-site AO repulsion must remain unscreened under Resta model"
            )

        # Check cross-atom block (atom 0 with atom 1: AOs 0..3 vs 4..7)
        # Screened repulsion must be strictly less than bare repulsion
        cross_w = w_bse[0:4, 4:8]
        cross_bare = gamma_bare[0:4, 4:8]
        self.assertTrue(np.all(cross_w < cross_bare), "Inter-atomic repulsion must be screened: W < Gamma_bare")

    def test_build_xs_kernel_dim(self):
        """Test Xs-QDEX kernel builder with Polarizable Dipole Interaction Model (DIM)."""
        w_bse, w_qp, gamma_bare, eps_info = build_xs_kernel(
            shells=self.shells,
            atom_symbols=self.symbols,
            coords=self.coords,
            atom_ao_ranges=self.ao_ranges,
            material_name="CDSE",
            kernel_mode="dim",
            eps_out=2.4,
            return_eps_info=True,
            nthreads=1
        )

        self.assertEqual(w_bse.shape, (self.n_ao, self.n_ao))
        self.assertEqual(eps_info["kernel_mode"], "dim")
        self.assertGreater(eps_info["eps_eff_exciton"], 1.0)

        # On-site blocks must remain unscreened (core invariance)
        for a0, a1 in self.ao_ranges:
            np.testing.assert_allclose(
                w_bse[a0:a1, a0:a1], gamma_bare[a0:a1, a0:a1], atol=1e-12,
                err_msg="On-site AO repulsion must remain unscreened under DIM model"
            )

        # Cross-atom block must be screened
        cross_w = w_bse[0:4, 4:8]
        cross_bare = gamma_bare[0:4, 4:8]
        self.assertTrue(np.all(cross_w < cross_bare), "Inter-atomic repulsion must be screened: W < Gamma_bare")

    def test_build_dim_mnok(self):
        """Test MNOK kernel builder with Polarizable Dipole Interaction Model (DIM)."""
        from qdex.hardness import build_dim_mnok
        w_dim, w_qp, gamma_bare, eps_info = build_dim_mnok(
            atom_symbols=self.symbols,
            coords=self.coords,
            material_name="CDSE",
            eps_out=2.4
        )

        self.assertEqual(w_dim.shape, (len(self.symbols), len(self.symbols)))
        self.assertEqual(eps_info["kernel_mode"], "dim")
        # On-site must remain bare
        np.testing.assert_allclose(np.diag(w_dim), np.diag(gamma_bare), atol=1e-12)
        # Off-diagonal must be screened
        self.assertTrue(w_dim[0, 1] < gamma_bare[0, 1])

    def test_exciton_hamiltonian_xs_mode(self):
        """Test ExcitonHamiltonian initialization and matvec in Xs-QDEX mode."""
        S = libint_cpp.overlap(self.shells, 1)
        # Build orthogonal MOs using Löwdin
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals)) # C.T @ S @ C = I

        eps = np.linspace(-10.0, 5.0, self.n_ao)
        homo_index = 3 # 4 occ, 4 virt
        n_occ, n_virt = 2, 2

        w_bse, _, gamma_bare = build_xs_kernel(
            shells=self.shells,
            atom_symbols=self.symbols,
            coords=self.coords,
            atom_ao_ranges=self.ao_ranges,
            material_name="CDSE",
            kernel_mode="resta",
            nthreads=1
        )

        # Singlet Hamiltonian
        ham_singlet = ExcitonHamiltonian(
            C=C, eps=eps, overlap=S, atom_ao_ranges=self.ao_ranges,
            homo_index=homo_index, n_occ=n_occ, n_virt=n_virt,
            scissor_ev=1.0, gamma_qp=w_bse, gamma_bse=w_bse,
            gamma_bare=gamma_bare, include_exchange=True,
            spin='singlet', kernel_type="xs"
        )

        self.assertTrue(ham_singlet.is_xs)
        self.assertEqual(ham_singlet.n_features, self.n_ao)
        self.assertEqual(ham_singlet.q_flat.shape, (ham_singlet.dim, self.n_ao))

        # Test matvec hermiticity: <x|H|y> = <y|H|x>^*
        x = np.random.randn(ham_singlet.dim)
        y = np.random.randn(ham_singlet.dim)
        Hx = ham_singlet.matvec(x)
        Hy = ham_singlet.matvec(y)
        np.testing.assert_allclose(np.dot(y, Hx), np.dot(x, Hy), atol=1e-10)

        # Triplet Hamiltonian (no Kx bare exchange)
        ham_triplet = ExcitonHamiltonian(
            C=C, eps=eps, overlap=S, atom_ao_ranges=self.ao_ranges,
            homo_index=homo_index, n_occ=n_occ, n_virt=n_virt,
            scissor_ev=1.0, gamma_qp=w_bse, gamma_bse=w_bse,
            gamma_bare=gamma_bare, include_exchange=True,
            spin='triplet', kernel_type="xs"
        )
        kx_trip, kd_trip = ham_triplet.kernel_actions(x)
        np.testing.assert_allclose(kx_trip, 0.0, atol=1e-14, err_msg="Triplet must have zero exchange Kx")

    def test_exciton_solver_xs_qdex(self):
        """Test full ExcitonSolver workflow with Xs-QDEX and Davidson diagonalization."""
        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        eps = np.array([-15.0, -10.0, -8.0, -5.0, 1.0, 3.0, 6.0, 10.0])
        occ = np.array([2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0])
        homo_index = 3
        n_occ, n_virt = 2, 2

        # 1. Xs-QDEX with Resta screening
        solver_resta = ExcitonSolver(
            C=C, eps=eps, occ=occ, overlap=S,
            atom_symbols=self.symbols, atom_coords=self.coords,
            atom_ao_ranges=self.ao_ranges, homo_index=homo_index,
            n_occ=n_occ, n_virt=n_virt, scissor_ev=1.5,
            kernel="xs-resta", alpha=1.0, material="CDSE",
            include_exchange=True, spin='singlet',
            kernel_type="xs", shells=self.shells
        )

        energies_resta, vecs_resta = solver_resta.solve(nroots=3, tol=1e-5)
        self.assertEqual(len(energies_resta), 3)
        self.assertTrue(np.all(np.diff(energies_resta) >= -1e-8), "Eigenvalues must be sorted")

        # 2. Legacy MNOK comparison (user choice)
        solver_mnok = ExcitonSolver(
            C=C, eps=eps, occ=occ, overlap=S,
            atom_symbols=self.symbols, atom_coords=self.coords,
            atom_ao_ranges=self.ao_ranges, homo_index=homo_index,
            n_occ=n_occ, n_virt=n_virt, scissor_ev=1.5,
            kernel="resta", alpha=1.0, material="CDSE",
            include_exchange=True, spin='singlet',
            kernel_type="mnok"
        )
        energies_mnok, vecs_mnok = solver_mnok.solve(nroots=3, tol=1e-5)
        self.assertEqual(len(energies_mnok), 3)
        self.assertFalse(solver_mnok.ham.is_xs)
        self.assertEqual(solver_mnok.ham.n_features, len(self.symbols))

    def test_xs_uks_spin_preserving(self):
        """Test UKS spin-preserving mode with Xs-QDEX kernel."""
        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C_alpha = U @ np.diag(1.0 / np.sqrt(eigvals))
        C_beta = C_alpha.copy()

        eps_alpha = np.array([-15.0, -10.0, -8.0, -5.0, 1.0, 3.0, 6.0, 10.0])
        eps_beta = eps_alpha + 0.2
        occ_alpha = np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        occ_beta = occ_alpha.copy()
        homo_index = 3
        n_occ, n_virt = 2, 2

        solver_uks = ExcitonSolver(
            C=C_alpha, eps=eps_alpha, occ=occ_alpha, overlap=S,
            atom_symbols=self.symbols, atom_coords=self.coords,
            atom_ao_ranges=self.ao_ranges, homo_index=homo_index,
            n_occ=n_occ, n_virt=n_virt, scissor_ev=1.0,
            kernel="xs", alpha=1.0, material="DEFAULT",
            include_exchange=True, spin='uks_spin_preserving',
            C_beta=C_beta, eps_beta=eps_beta, homo_index_beta=homo_index,
            n_occ_beta=n_occ, n_virt_beta=n_virt,
            kernel_type="xs", shells=self.shells
        )

        energies, vecs = solver_uks.solve(nroots=4, tol=1e-5)
        self.assertEqual(len(energies), 4)
        self.assertTrue(solver_uks.ham.is_xs)
        self.assertEqual(solver_uks.ham.n_features, self.n_ao)

    def test_xs_soc_spinor(self):
        """Test SOC spinor basis transformation with Xs-QDEX kernel."""
        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        eps = np.array([-15.0, -10.0, -8.0, -5.0, 1.0, 3.0, 6.0, 10.0])
        occ = np.array([2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0])
        homo_index = 3
        n_occ, n_virt = 2, 2

        k = n_occ + n_virt
        soc_U = np.eye(2 * k, dtype=complex)
        soc_E = np.linspace(-6.0, 6.0, 2 * k)

        solver_soc = ExcitonSolver(
            C=C, eps=eps, occ=occ, overlap=S,
            atom_symbols=self.symbols, atom_coords=self.coords,
            atom_ao_ranges=self.ao_ranges, homo_index=homo_index,
            n_occ=n_occ, n_virt=n_virt, scissor_ev=1.0,
            kernel="xs-resta", alpha=1.0, material="CDSE",
            include_exchange=True, soc_U=soc_U, soc_E=soc_E,
            kernel_type="xs", shells=self.shells
        )

        energies, vecs = solver_soc.solve(nroots=4, tol=1e-5)
        self.assertEqual(len(energies), 4)
        self.assertTrue(solver_soc.ham.is_xs)
        self.assertEqual(solver_soc.ham.n_features, self.n_ao)
        self.assertEqual(solver_soc.ham.dim, (2 * n_occ) * (2 * n_virt))

    def test_xs_rpa_kernel_and_eps_logging(self):
        """Test parameter-free microscopic ZDO-RPA screening and eps_eff logging."""
        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        S_half = U @ np.diag(np.sqrt(eigvals)) @ U.T

        eps = np.array([-15.0, -10.0, -8.0, -5.0, 1.0, 3.0, 6.0, 10.0])
        homo_index = 3
        n_occ, n_virt = 2, 2

        occ_idx = np.arange(homo_index - n_occ + 1, homo_index + 1)
        virt_idx = np.arange(homo_index + 1, homo_index + 1 + n_virt)
        C_occ_low = S_half @ C[:, occ_idx]
        C_virt_low = S_half @ C[:, virt_idx]

        w_rpa, _, gamma_bare, eps_info = build_xs_kernel(
            shells=self.shells,
            atom_symbols=self.symbols,
            coords=self.coords,
            atom_ao_ranges=self.ao_ranges,
            material_name="CDSE",
            kernel_mode="rpa",
            alpha=1.0,
            nthreads=1,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            eps_occ=eps[occ_idx],
            eps_virt=eps[virt_idx],
            return_eps_info=True
        )

        # 1. Symmetry
        np.testing.assert_allclose(w_rpa, w_rpa.T, atol=1e-12, err_msg="RPA W must be symmetric")

        # 2. Dielectric Screening: W must be screened relative to bare Gamma
        self.assertTrue(np.all(w_rpa <= gamma_bare + 1e-12), "RPA screening must reduce Coulomb repulsion (W <= Gamma)")

        # 3. Check eps_info diagnostics
        self.assertIn("eps_eff_exciton", eps_info)
        self.assertIn("eps_eff_hole", eps_info)
        self.assertIn("eps_eff_elec", eps_info)
        self.assertGreater(eps_info["eps_eff_exciton"], 1.0, "Effective exciton dielectric constant must be > 1")
        self.assertGreaterEqual(eps_info["eps_eff_hole"], 1.0 - 1e-6)
        self.assertGreaterEqual(eps_info["eps_eff_elec"], 1.0 - 1e-6)

    def test_xs_rpa_solver(self):
        """Test full ExcitonSolver workflow with kernel='xs-rpa'."""
        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        eps = np.array([-15.0, -10.0, -8.0, -5.0, 1.0, 3.0, 6.0, 10.0])
        occ = np.array([2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0])
        homo_index = 3
        n_occ, n_virt = 2, 2

        solver_rpa = ExcitonSolver(
            C=C, eps=eps, occ=occ, overlap=S,
            atom_symbols=self.symbols, atom_coords=self.coords,
            atom_ao_ranges=self.ao_ranges, homo_index=homo_index,
            n_occ=n_occ, n_virt=n_virt, scissor_ev=1.5,
            kernel="xs-rpa", alpha=1.0, material="CDSE",
            include_exchange=True, spin='singlet',
            kernel_type="xs", shells=self.shells
        )

        energies_rpa, vecs_rpa = solver_rpa.solve(nroots=3, tol=1e-5)
        self.assertEqual(len(energies_rpa), 3)
        self.assertTrue(np.all(np.diff(energies_rpa) >= -1e-8), "Eigenvalues must be sorted")
        self.assertIsNotNone(solver_rpa.eps_info)
        self.assertGreater(solver_rpa.eps_info["eps_eff_exciton"], 1.0)

    def test_build_sbse_kernel_atom(self):
        """Test atom-resolved sBSE kernel builder (Cho et al. 2022)."""
        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        S_half = U @ np.diag(np.sqrt(eigvals)) @ U.T

        eps = np.array([-15.0, -10.0, -8.0, -5.0, 1.0, 3.0, 6.0, 10.0])
        homo_index = 3
        n_occ, n_virt = 2, 2

        occ_idx = np.arange(homo_index - n_occ + 1, homo_index + 1)
        virt_idx = np.arange(homo_index + 1, homo_index + 1 + n_virt)
        C_occ_low = S_half @ C[:, occ_idx]
        C_virt_low = S_half @ C[:, virt_idx]

        # 1. Vacuum screening (eps_out = 1.0)
        w_vac, _, j_bare, info_vac = build_sbse_kernel(
            atom_symbols=self.symbols,
            coords=self.coords,
            atom_ao_ranges=self.ao_ranges,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            eps_occ=eps[occ_idx],
            eps_virt=eps[virt_idx],
            mode="atom",
            eps_out=1.0,
            material_name="CDSE",
            return_eps_info=True
        )

        n_atoms = len(self.symbols)
        self.assertEqual(w_vac.shape, (n_atoms, n_atoms))
        np.testing.assert_allclose(w_vac, w_vac.T, atol=1e-12, err_msg="Atom-resolved W must be symmetric")
        self.assertTrue(np.all(np.diag(w_vac) <= np.diag(j_bare) + 1e-12), "Diagonal W must be <= bare J")
        diff_eigs = np.linalg.eigvalsh(j_bare - w_vac)
        self.assertTrue(np.all(diff_eigs >= -1e-10), "Screening must be positive semi-definite: J - W >= 0")
        self.assertEqual(info_vac["sbse_mode"], "atom")
        self.assertGreaterEqual(info_vac["eps_eff_exciton"], 1.0 - 1e-8)

        # 2. Solvent screening (eps_out = 2.4): solvent response should increase screening (reduce W)
        w_solv, _, _, info_solv = build_sbse_kernel(
            atom_symbols=self.symbols,
            coords=self.coords,
            atom_ao_ranges=self.ao_ranges,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            eps_occ=eps[occ_idx],
            eps_virt=eps[virt_idx],
            mode="atom",
            eps_out=2.4,
            material_name="CDSE",
            return_eps_info=True
        )
        self.assertTrue(np.all(np.diag(w_solv) <= np.diag(w_vac) + 1e-12), "Solvent screening must reduce diagonal W relative to vacuum")
        self.assertGreater(info_solv["eps_eff_exciton"], info_vac["eps_eff_exciton"])

    def test_build_sbse_kernel_ao(self):
        """Test AO-resolved sBSE kernel builder (Cho et al. 2022)."""
        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        S_half = U @ np.diag(np.sqrt(eigvals)) @ U.T

        eps = np.array([-15.0, -10.0, -8.0, -5.0, 1.0, 3.0, 6.0, 10.0])
        homo_index = 3
        n_occ, n_virt = 2, 2

        occ_idx = np.arange(homo_index - n_occ + 1, homo_index + 1)
        virt_idx = np.arange(homo_index + 1, homo_index + 1 + n_virt)
        C_occ_low = S_half @ C[:, occ_idx]
        C_virt_low = S_half @ C[:, virt_idx]

        w_ao, _, j_bare_ao, info_ao = build_sbse_kernel(
            atom_symbols=self.symbols,
            coords=self.coords,
            atom_ao_ranges=self.ao_ranges,
            shells=self.shells,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            eps_occ=eps[occ_idx],
            eps_virt=eps[virt_idx],
            mode="ao",
            eps_out=2.4,
            material_name="CDSE",
            return_eps_info=True
        )

        self.assertEqual(w_ao.shape, (self.n_ao, self.n_ao))
        np.testing.assert_allclose(w_ao, w_ao.T, atol=1e-12, err_msg="AO-resolved W must be symmetric")
        self.assertTrue(np.all(np.diag(w_ao) <= np.diag(j_bare_ao) + 1e-12), "Diagonal AO W must be <= bare J")
        diff_eigs_ao = np.linalg.eigvalsh(j_bare_ao - w_ao)
        self.assertTrue(np.all(diff_eigs_ao >= -1e-10), "Screening must be positive semi-definite: J - W >= 0")
        self.assertEqual(info_ao["sbse_mode"], "ao")
        self.assertGreater(info_ao["eps_eff_exciton"], 1.0)

    def test_estimate_sgw_qp_gap(self):
        """Test microscopic sGW Delta-W QP gap shift."""
        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        S_half = U @ np.diag(np.sqrt(eigvals)) @ U.T

        eps = np.array([-15.0, -10.0, -8.0, -5.0, 1.0, 3.0, 6.0, 10.0])
        homo_index = 3
        n_occ, n_virt = 2, 2

        occ_idx = np.arange(homo_index - n_occ + 1, homo_index + 1)
        virt_idx = np.arange(homo_index + 1, homo_index + 1 + n_virt)
        C_occ_low = S_half @ C[:, occ_idx]
        C_virt_low = S_half @ C[:, virt_idx]

        # Test atom mode
        scissor_atom, prov_atom = estimate_sgw_qp_gap(
            coords=self.coords,
            atom_symbols=self.symbols,
            material_name="CDSE",
            eps_out=2.4,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            eps_occ=eps[occ_idx],
            eps_virt=eps[virt_idx],
            atom_ao_ranges=self.ao_ranges,
            mode="atom",
            return_details=True
        )
        self.assertGreater(scissor_atom, 0.0)
        self.assertIn("bulk_shift_ev", prov_atom)
        self.assertIn("confinement_shift_ev", prov_atom)
        self.assertIn("total_scissor_ev", prov_atom)
        self.assertEqual(prov_atom["qp_model"], "sgw_dw")

        # Test AO mode
        scissor_ao, prov_ao = estimate_sgw_qp_gap(
            coords=self.coords,
            atom_symbols=self.symbols,
            material_name="CDSE",
            eps_out=2.4,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            eps_occ=eps[occ_idx],
            eps_virt=eps[virt_idx],
            atom_ao_ranges=self.ao_ranges,
            shells=self.shells,
            mode="ao",
            return_details=True
        )
        self.assertGreater(scissor_ao, 0.0)
        self.assertEqual(prov_ao["qp_model"], "sgw_dw")

    def test_exciton_solver_sbse(self):
        """Test full ExcitonSolver workflow with kernel='sbse' (atom) and 'sbse-ao' (AO)."""
        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        eps = np.array([-15.0, -10.0, -8.0, -5.0, 1.0, 3.0, 6.0, 10.0])
        occ = np.array([2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0])
        homo_index = 3
        n_occ, n_virt = 2, 2

        # 1. Atom-resolved sBSE
        solver_atom = ExcitonSolver(
            C=C, eps=eps, occ=occ, overlap=S,
            atom_symbols=self.symbols, atom_coords=self.coords,
            atom_ao_ranges=self.ao_ranges, homo_index=homo_index,
            n_occ=n_occ, n_virt=n_virt, scissor_ev=1.5,
            kernel="sbse-atom", alpha=1.0, material="CDSE",
            include_exchange=True, spin='singlet',
            kernel_type="mnok", shells=self.shells
        )
        energies_atom, vecs_atom = solver_atom.solve(nroots=3, tol=1e-5)
        self.assertEqual(len(energies_atom), 3)
        self.assertTrue(np.all(np.diff(energies_atom) >= -1e-8))
        self.assertIsNotNone(solver_atom.eps_info)
        self.assertEqual(solver_atom.eps_info["sbse_mode"], "atom")

        # 2. AO-resolved sBSE
        solver_ao = ExcitonSolver(
            C=C, eps=eps, occ=occ, overlap=S,
            atom_symbols=self.symbols, atom_coords=self.coords,
            atom_ao_ranges=self.ao_ranges, homo_index=homo_index,
            n_occ=n_occ, n_virt=n_virt, scissor_ev=1.5,
            kernel="sbse-ao", alpha=1.0, material="CDSE",
            include_exchange=True, spin='singlet',
            kernel_type="xs", shells=self.shells
        )
        energies_ao, vecs_ao = solver_ao.solve(nroots=3, tol=1e-5)
        self.assertEqual(len(energies_ao), 3)
        self.assertTrue(np.all(np.diff(energies_ao) >= -1e-8))
        self.assertIsNotNone(solver_ao.eps_info)
        self.assertEqual(solver_ao.eps_info["sbse_mode"], "ao")

    def test_estimate_sgw_dim_qp_gap(self):
        """Test sGW-DIM Delta-W quasiparticle scissor estimation."""
        from qdex.hardness import estimate_sgw_dim_qp_gap
        from qdex.lowdin import lowdin_sqrt

        S = libint_cpp.overlap(self.shells, 1)
        S_half = lowdin_sqrt(S)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        homo_index = 3
        occ_idx = np.arange(homo_index - 1, homo_index + 1)
        virt_idx = np.arange(homo_index + 1, homo_index + 3)
        C_occ_low = S_half @ C[:, occ_idx]
        C_virt_low = S_half @ C[:, virt_idx]

        # 1. Solvent eps_out = 2.4
        scissor_solv, prov_solv = estimate_sgw_dim_qp_gap(
            coords=self.coords,
            atom_symbols=self.symbols,
            material_name="CDSE",
            eps_out=2.4,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            atom_ao_ranges=self.ao_ranges,
            return_details=True
        )
        self.assertGreater(scissor_solv, 0.0)
        self.assertEqual(prov_solv["qp_model"], "sgw_dim")
        self.assertIn("confinement_shift_internal_ev", prov_solv)
        self.assertIn("confinement_shift_solvent_ev", prov_solv)
        self.assertGreater(prov_solv["confinement_shift_internal_ev"], 0.0)
        self.assertGreater(prov_solv["confinement_shift_solvent_ev"], 0.0)

        # 2. Vacuum eps_out = 1.0 (must have larger scissor than solvent)
        scissor_vac, prov_vac = estimate_sgw_dim_qp_gap(
            coords=self.coords,
            atom_symbols=self.symbols,
            material_name="CDSE",
            eps_out=1.0,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            atom_ao_ranges=self.ao_ranges,
            return_details=True
        )
        self.assertGreater(scissor_vac, scissor_solv)

    def test_estimate_sgw_resta_qp_gap(self):
        """Test sGW-Resta Delta-W quasiparticle scissor estimation (Penn-scaled vs Pure boundary)."""
        from qdex.hardness import estimate_sgw_resta_qp_gap
        from qdex.lowdin import lowdin_sqrt

        S = libint_cpp.overlap(self.shells, 1)
        S_half = lowdin_sqrt(S)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        homo_index = 3
        occ_idx = np.arange(homo_index - 1, homo_index + 1)
        virt_idx = np.arange(homo_index + 1, homo_index + 3)
        C_occ_low = S_half @ C[:, occ_idx]
        C_virt_low = S_half @ C[:, virt_idx]

        # 1. Penn-scaled Resta
        scissor_penn, prov_penn = estimate_sgw_resta_qp_gap(
            coords=self.coords,
            atom_symbols=self.symbols,
            material_name="CDSE",
            eps_out=2.4,
            dft_gap=2.5,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            atom_ao_ranges=self.ao_ranges,
            penn_scaling=True,
            return_details=True
        )
        self.assertGreater(scissor_penn, 0.0)
        self.assertEqual(prov_penn["qp_model"], "sgw_resta")
        self.assertTrue(prov_penn["penn_scaling"])
        self.assertGreater(prov_penn["confinement_shift_internal_ev"], 0.0)

        # 2. Pure boundary Resta (penn_scaling=False)
        scissor_pure, prov_pure = estimate_sgw_resta_qp_gap(
            coords=self.coords,
            atom_symbols=self.symbols,
            material_name="CDSE",
            eps_out=2.4,
            dft_gap=2.5,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            atom_ao_ranges=self.ao_ranges,
            penn_scaling=False,
            return_details=True
        )
        self.assertGreater(scissor_pure, 0.0)
        self.assertFalse(prov_pure["penn_scaling"])
        self.assertEqual(prov_pure["confinement_shift_internal_ev"], 0.0)

        # Penn scaling must produce larger scissor than pure boundary
        self.assertGreater(scissor_penn, scissor_pure)

    def test_compute_dynamic_z(self):
        """Test physical behavior of dynamic Z calculation via plasmon-pole model."""
        from qdex.hardness import compute_dynamic_z

        # 1. Zero self-energy shift gives Z = 1
        z_zero = compute_dynamic_z(0.0, gap_ev=2.5, eps_eff=5.0)
        self.assertAlmostEqual(z_zero, 1.0, places=12)

        # One plasmon pole: Z = 1 / (1 + dSigma / omega_tilde), omega_tilde = omega_p / sqrt(1 - 1/eps)
        omega = 15.0 / np.sqrt(1.0 - 1.0 / 5.0)
        self.assertAlmostEqual(compute_dynamic_z(1.0, eps_eff=5.0, omega_p_ev=15.0), 1.0 / (1.0 + 1.0 / omega), places=12)

        # 2. Typical self-energy shifts (0.3 - 1.5 eV) yield physical Z in [0.75, 0.98]
        z_small = compute_dynamic_z(0.3, gap_ev=2.5, eps_eff=5.0)
        z_med = compute_dynamic_z(0.8, gap_ev=2.5, eps_eff=5.0)
        z_large = compute_dynamic_z(1.5, gap_ev=2.5, eps_eff=5.0)

        self.assertTrue(0.70 < z_small < 0.99)
        self.assertTrue(0.70 < z_med < 0.99)
        self.assertTrue(0.60 < z_large < 0.99)

        # 3. Monotonicity: larger self-energy leads to stronger renormalization (lower Z)
        self.assertGreater(z_small, z_med)
        self.assertGreater(z_med, z_large)

    def test_evgw_dim_qp_gap(self):
        """Test Eigenvalue Self-Consistent evGW with Polarizable Dipole Interaction Model (DIM)."""
        from qdex.hardness import estimate_evgw_dim_qp_gap
        from qdex.lowdin import lowdin_sqrt

        S = libint_cpp.overlap(self.shells, 1)
        S_half = lowdin_sqrt(S)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        homo_index = 3
        occ_idx = np.arange(homo_index - 1, homo_index + 1)
        virt_idx = np.arange(homo_index + 1, homo_index + 3)
        C_occ_low = S_half @ C[:, occ_idx]
        C_virt_low = S_half @ C[:, virt_idx]

        scissor, prov = estimate_evgw_dim_qp_gap(
            coords=self.coords,
            atom_symbols=self.symbols,
            material_name="CDSE",
            eps_out=2.4,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            atom_ao_ranges=self.ao_ranges,
            return_details=True
        )

        self.assertGreater(scissor, 0.0)
        self.assertEqual(prov["qp_model"], "evgw_dim")
        self.assertTrue(prov["self_consistent"])
        self.assertTrue(prov["dynamic_z"])
        self.assertTrue(prov["evgw_converged"])
        self.assertLessEqual(prov["evgw_iterations"], 25)
        self.assertTrue(0.5 <= prov["z_homo"] <= 0.99)
        self.assertTrue(0.5 <= prov["z_lumo"] <= 0.99)

    def test_evgw_resta_qp_gap(self):
        """Test Eigenvalue Self-Consistent evGW with Resta-Penn Model."""
        from qdex.hardness import estimate_evgw_resta_qp_gap
        from qdex.lowdin import lowdin_sqrt

        S = libint_cpp.overlap(self.shells, 1)
        S_half = lowdin_sqrt(S)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        homo_index = 3
        occ_idx = np.arange(homo_index - 1, homo_index + 1)
        virt_idx = np.arange(homo_index + 1, homo_index + 3)
        C_occ_low = S_half @ C[:, occ_idx]
        C_virt_low = S_half @ C[:, virt_idx]

        scissor, prov = estimate_evgw_resta_qp_gap(
            coords=self.coords,
            atom_symbols=self.symbols,
            material_name="CDSE",
            eps_out=2.4,
            dft_gap=2.5,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            atom_ao_ranges=self.ao_ranges,
            return_details=True
        )

        self.assertGreater(scissor, 0.0)
        self.assertEqual(prov["qp_model"], "evgw_resta")
        self.assertTrue(prov["self_consistent"])
        self.assertTrue(prov["dynamic_z"])
        self.assertTrue(prov["evgw_converged"])
        self.assertLessEqual(prov["evgw_iterations"], 25)
        self.assertTrue(0.5 <= prov["z_homo"] <= 0.99)
        self.assertTrue(0.5 <= prov["z_lumo"] <= 0.99)

    def test_qsgw_dim_qp_gap(self):
        """Test Full AO-basis Quasiparticle Self-Consistent qsGW with DIM."""
        from qdex.hardness import estimate_qsgw_dim_qp_gap

        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        homo_index = 3
        eps = np.array([-15.0, -12.0, -9.0, -6.0, 1.0, 3.0, 5.0, 7.0])

        scissor, prov, C_qp, eps_qp = estimate_qsgw_dim_qp_gap(
            coords=self.coords,
            atom_symbols=self.symbols,
            C=C,
            eps=eps,
            S=S,
            atom_ao_ranges=self.ao_ranges,
            homo_index=homo_index,
            material_name="CDSE",
            eps_out=2.4,
            return_details=True
        )

        self.assertGreater(scissor, 0.0)
        self.assertEqual(prov["qp_model"], "qsgw_dim")
        self.assertTrue(prov["orbital_update"])
        self.assertTrue(prov["qsgw_converged"])
        self.assertLessEqual(prov["qsgw_iterations"], 25)
        self.assertGreater(prov["homo_fidelity"], 0.90)
        self.assertGreater(prov["lumo_fidelity"], 0.90)

        # Orthonormality check of updated C_qp
        err_ortho = np.max(np.abs(C_qp.T @ S @ C_qp - np.eye(self.n_ao)))
        self.assertLess(err_ortho, 1e-10)

    def test_qsgw_resta_qp_gap(self):
        """Test Full AO-basis Quasiparticle Self-Consistent qsGW with Resta."""
        from qdex.hardness import estimate_qsgw_resta_qp_gap

        S = libint_cpp.overlap(self.shells, 1)
        eigvals, U = np.linalg.eigh(S)
        C = U @ np.diag(1.0 / np.sqrt(eigvals))
        homo_index = 3
        eps = np.array([-15.0, -12.0, -9.0, -6.0, 1.0, 3.0, 5.0, 7.0])

        scissor, prov, C_qp, eps_qp = estimate_qsgw_resta_qp_gap(
            coords=self.coords,
            atom_symbols=self.symbols,
            C=C,
            eps=eps,
            S=S,
            atom_ao_ranges=self.ao_ranges,
            homo_index=homo_index,
            material_name="CDSE",
            eps_out=2.4,
            return_details=True
        )

        self.assertGreater(scissor, 0.0)
        self.assertEqual(prov["qp_model"], "qsgw_resta")
        self.assertTrue(prov["orbital_update"])
        self.assertTrue(prov["qsgw_converged"])
        self.assertLessEqual(prov["qsgw_iterations"], 25)
        self.assertGreater(prov["homo_fidelity"], 0.90)
        self.assertGreater(prov["lumo_fidelity"], 0.90)

        # Orthonormality check of updated C_qp
        err_ortho = np.max(np.abs(C_qp.T @ S @ C_qp - np.eye(self.n_ao)))
        self.assertLess(err_ortho, 1e-10)

    def test_cli_sgw_anchor_and_2e_integrals(self):
        """Test CLI and config parsing for sgw-anchor and 2e-integrals / two-electron-integrals."""
        from qdex.cli import _apply_config
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--qp_gap", type=str, default="brus")
        parser.add_argument("--two-electron-integrals", "--two_electron_integrals", "--2e-integrals", "--2e_integrals", "--kernel_type", "--kernel-type",
                            dest="kernel_type", choices=["mnok", "xs", "xs-qdex"], default="mnok")
        
        # Test 1: YAML with 2e-integrals: xs and qp_gap: sgw-anchor
        args1 = parser.parse_args([])
        cfg1 = {
            "physics": {
                "qp_gap": "sgw-anchor",
                "2e-integrals": "xs",
            }
        }
        _apply_config(args1, cfg1, explicit_cli_args=set())
        self.assertEqual(args1.qp_gap, "sgw-anchor")
        self.assertEqual(args1.kernel_type, "xs")

        # Test 2: YAML with two_electron_integrals: xs
        args2 = parser.parse_args([])
        cfg2 = {
            "physics": {
                "two_electron_integrals": "xs",
            }
        }
        _apply_config(args2, cfg2, explicit_cli_args=set())
        self.assertEqual(args2.kernel_type, "xs")

        # Test 3: CLI explicit flag overrides YAML
        args3 = parser.parse_args(["--2e-integrals", "xs"])
        explicit = {"2e_integrals"}
        cfg3 = {
            "physics": {
                "two_electron_integrals": "mnok",
            }
        }
        _apply_config(args3, cfg3, explicit_cli_args=explicit)
        self.assertEqual(args3.kernel_type, "xs")


if __name__ == '__main__':
    unittest.main()




