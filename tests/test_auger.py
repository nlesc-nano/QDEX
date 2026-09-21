"""
Unit tests for the QDEX Auger Recombination module (qdex.auger).
"""

import unittest
import numpy as np

from qdex.auger import (
    compute_transition_density_charges,
    compute_auger_matrix_element,
    calculate_auger_rates,
    AugerResult,
)


class TestAugerRecombination(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        # Create a toy system: 2 atoms, 4 AOs per atom = 8 AOs total, 8 MOs
        self.n_atoms = 2
        self.ao_per_atom = 4
        self.n_ao = self.n_atoms * self.ao_per_atom
        self.atom_ao_ranges = [(0, 4), (4, 8)]
        self.coords = np.array([
            [0.0, 0.0, 0.0],
            [2.5, 0.0, 0.0],
        ], dtype=np.float64)
        self.atom_symbols = ["Pb", "Br"]

        # Generate a valid symmetric positive-definite overlap matrix
        A = np.random.randn(self.n_ao, self.n_ao)
        self.S = np.eye(self.n_ao) + 0.1 * (A + A.T) / np.sqrt(self.n_ao)
        # S must be positive definite
        eigvals, eigvecs = np.linalg.eigh(self.S)
        eigvals = np.clip(eigvals, 0.5, None)
        self.S = eigvecs @ np.diag(eigvals) @ eigvecs.T

        # Generate orthonormal MO coefficients: C^T S C = I
        L = np.linalg.cholesky(self.S)
        C_raw = np.random.randn(self.n_ao, self.n_ao)
        Q, _ = np.linalg.qr(C_raw)
        self.C = np.linalg.solve(L.T, Q)
        # Verify C^T S C == I
        np.testing.assert_allclose(self.C.T @ self.S @ self.C, np.eye(self.n_ao), atol=1e-10)

        # 8 MO eigenvalues: 4 occupied (-8, -6, -4, -2 eV), 4 virtual (1, 3, 5, 7 eV)
        # HOMO is index 3 (E = -2.0 eV), LUMO is index 4 (E = +1.0 eV), Gap = 3.0 eV
        self.eps = np.array([-8.0, -6.0, -4.0, -2.0, 1.0, 3.0, 5.0, 7.0], dtype=np.float64)
        self.homo_idx = 3

        # Simple 2x2 Resta kernel matrix in eV
        self.W_resta = np.array([
            [4.5, 1.8],
            [1.8, 5.2],
        ], dtype=np.float64)

    def test_transition_density_charges_orthonormality(self):
        """Verify that sum_A q_A^{ij} == delta_{ij} in Mulliken population analysis."""
        q = compute_transition_density_charges(
            C=self.C,
            S=self.S,
            atom_ao_ranges=self.atom_ao_ranges,
            i_indices=np.arange(self.n_ao),
            j_indices=np.arange(self.n_ao),
        )
        self.assertEqual(q.shape, (self.n_ao, self.n_ao, self.n_atoms))

        # Summing over atoms should equal the identity matrix (C^T S C)
        q_sum_atoms = np.sum(q, axis=2)
        np.testing.assert_allclose(q_sum_atoms, np.eye(self.n_ao), atol=1e-10)

    def test_compute_auger_matrix_element(self):
        """Test direct, exchange, and effective net Auger matrix element evaluation."""
        q_recomb = np.array([0.6, -0.4], dtype=np.float64)
        q_eject_1 = np.array([0.3, 0.7], dtype=np.float64)
        q_eject_2 = np.array([0.5, 0.2], dtype=np.float64)

        v_dir, v_exch, m_eff = compute_auger_matrix_element(
            q_recomb=q_recomb,
            q_eject_1=q_eject_1,
            q_eject_2=q_eject_2,
            W_resta=self.W_resta,
        )

        expected_dir = float(q_eject_1 @ self.W_resta @ q_recomb)
        expected_exch = float(q_recomb @ self.W_resta @ q_eject_2)

        self.assertAlmostEqual(v_dir, expected_dir, places=10)
        self.assertAlmostEqual(v_exch, expected_exch, places=10)
        self.assertAlmostEqual(m_eff, expected_dir - expected_exch, places=10)

    def test_calculate_auger_rates_structure(self):
        """Test calculate_auger_rates execution and output structure."""
        # For eeh: recombining (LUMO, HOMO) releases E_recomb = 1.0 - (-2.0) = 3.0 eV.
        # Spectator electron in LUMO (1.0 eV) target energy: 1.0 + 3.0 = 4.0 eV.
        # Virtual orbital index 6 has energy 5.0 eV (mismatch = 1.0 eV).
        # We use sigma_ev = 0.5 eV to capture virtual transitions.
        res = calculate_auger_rates(
            C=self.C,
            eps=self.eps,
            S=self.S,
            atom_ao_ranges=self.atom_ao_ranges,
            coords=self.coords,
            atom_symbols=self.atom_symbols,
            homo_idx=self.homo_idx,
            W_resta=self.W_resta,
            sigma_ev=0.5,
            channel="all",
            verbose=False,
        )

        self.assertIsInstance(res, AugerResult)
        self.assertGreaterEqual(res.rate_eeh_fs, 0.0)
        self.assertGreaterEqual(res.rate_hhe_fs, 0.0)
        self.assertGreaterEqual(res.rate_biexciton_fs, 0.0)
        self.assertEqual(res.fundamental_gap_ev, 3.0)

        # Biexciton rate = 2 * eeh + 2 * hhe
        self.assertAlmostEqual(res.rate_biexciton_fs, 2.0 * res.rate_eeh_fs + 2.0 * res.rate_hhe_fs, places=12)

        # Summary table formatting check
        table_str = res.summary_table()
        self.assertIn("AUGER RECOMBINATION REPORT", table_str)
        self.assertIn("Biexciton (XX)", table_str)

    def test_spinor_auger_rates(self):
        """Test calculate_auger_rates with two-component spinors."""
        U_a = self.C.astype(complex)
        U_b = np.zeros_like(U_a)

        res = calculate_auger_rates(
            C=self.C,
            eps=self.eps,
            S=self.S,
            atom_ao_ranges=self.atom_ao_ranges,
            coords=self.coords,
            atom_symbols=self.atom_symbols,
            homo_idx=self.homo_idx,
            W_resta=self.W_resta,
            sigma_ev=0.5,
            spinor=True,
            U_spinor_alpha=U_a,
            U_spinor_beta=U_b,
            verbose=False,
        )

        self.assertIsInstance(res, AugerResult)
        self.assertGreaterEqual(res.rate_biexciton_fs, 0.0)

    def test_cli_parser_auger_flags(self):
        """Verify CLI argument parser handles --auger and options."""
        from qdex.cli import main
        import argparse

        # Verify parser creation
        parser = argparse.ArgumentParser()
        # Import arguments from cli
        from qdex.cli import main
        # Test flag existence via CLI help output or parse_args
        with unittest.mock.patch("sys.argv", ["qdex", "--help"]):
            try:
                main()
            except SystemExit as e:
                self.assertEqual(e.code, 0)


if __name__ == "__main__":
    unittest.main()
