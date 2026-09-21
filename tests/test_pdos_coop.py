import unittest
import numpy as np
import os
import shutil
import tempfile
from qdex.pdos_coop import compute_pdos_and_coop, export_pdos_coop_data
from qdex.fuzzy_bands import _matmul_real_matrix

class TestPdosCoop(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.test_dir)

    def tearDown(self):
        os.chdir(self.orig_dir)
        shutil.rmtree(self.test_dir)

    def test_matmul_real_matrix(self):
        rng = np.random.default_rng(42)
        A = rng.standard_normal((20, 30))
        B_real = rng.standard_normal((30, 10))
        B_cplx = rng.standard_normal((30, 10)) + 1j * rng.standard_normal((30, 10))

        res_real = _matmul_real_matrix(A, B_real)
        np.testing.assert_allclose(res_real, A @ B_real, rtol=1e-12)

        res_cplx = _matmul_real_matrix(A, B_cplx)
        np.testing.assert_allclose(res_cplx, A @ B_cplx, rtol=1e-12)

    def test_compute_pdos_and_coop_basic(self):
        rng = np.random.default_rng(123)
        n_ao = 6
        n_mo = 4
        # Two atoms: Pb (index 0, l=1 -> 3 AOs), I (index 1, l=1 -> 3 AOs)
        shells = [
            {"sym": "Pb", "atom_idx": 0, "l": 1, "O": [0.0, 0.0, 0.0]},
            {"sym": "I",  "atom_idx": 1, "l": 1, "O": [2.0, 0.0, 0.0]},
        ]

        S = np.eye(n_ao)
        C = rng.standard_normal((n_ao, n_mo))
        eps = np.array([-2.0, -1.0, 1.0, 2.0])
        ewin = [-3.0, 3.0]

        analysis = compute_pdos_and_coop(
            C, S, eps, shells,
            pdos_atoms=["Pb", "I"],
            coop_pairs=["Pb-I"],
            ewin=ewin,
            sigma=0.1,
            is_soc=False,
            prefix="test_sf"
        )
        self.assertIn("P_weights", analysis)
        self.assertIn("coop_results", analysis)
        self.assertTrue(os.path.exists("pdos_data_test_sf.csv"))
        self.assertTrue(os.path.exists("coop_data_test_sf.csv"))

    def test_compute_pdos_and_coop_soc(self):
        rng = np.random.default_rng(456)
        n_ao = 6
        n_spinor = 4
        shells = [
            {"sym": "Pb", "atom_idx": 0, "l": 1, "O": [0.0, 0.0, 0.0]},
            {"sym": "I",  "atom_idx": 1, "l": 1, "O": [2.0, 0.0, 0.0]},
        ]

        S = np.eye(n_ao)
        C_spinor = rng.standard_normal((2 * n_ao, n_spinor)) + 1j * rng.standard_normal((2 * n_ao, n_spinor))
        eps_soc = np.array([-1.5, -0.5, 0.5, 1.5])
        ewin = [-2.0, 2.0]

        analysis = compute_pdos_and_coop(
            C_spinor, S, eps_soc, shells,
            pdos_atoms=["Pb", "I"],
            coop_pairs=["Pb-I"],
            ewin=ewin,
            sigma=0.1,
            is_soc=True,
            prefix="test_soc"
        )
        self.assertIn("P_weights", analysis)
        self.assertIn("coop_results", analysis)
        self.assertTrue(os.path.exists("pdos_data_test_soc.csv"))
        self.assertTrue(os.path.exists("coop_data_test_soc.csv"))

if __name__ == "__main__":
    unittest.main()
