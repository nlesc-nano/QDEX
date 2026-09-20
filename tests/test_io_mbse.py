import unittest
import os
import tempfile
import numpy as np

from miniBSE.io_utils import read_mos_mbse, write_mos_mbse, read_mos_auto


class TestMBSEIO(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_write_read_roundtrip(self):
        n_ao = 120
        n_mo = 50
        n_occ = 25
        n_elec = 50

        np.random.seed(42)
        C_orig = np.random.randn(n_ao, n_mo)
        eps_orig = np.linspace(-1.0, 0.5, n_mo)
        occ_orig = np.zeros(n_mo)
        occ_orig[:n_occ] = 2.0

        file_path = os.path.join(self.temp_dir.name, "test_orbitals.mbse")
        write_mos_mbse(file_path, C_orig, eps_orig, occ_orig, n_occ=n_occ, n_elec=n_elec)

        # Read back using read_mos_mbse
        C_read, eps_read, occ_read = read_mos_mbse(file_path, n_ao, verbose=False)

        self.assertEqual(C_read.shape, (n_ao, n_mo))
        self.assertEqual(eps_read.shape, (n_mo,))
        self.assertEqual(occ_read.shape, (n_mo,))

        np.testing.assert_allclose(C_read, C_orig, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(eps_read, eps_orig, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(occ_read, occ_orig, rtol=1e-14, atol=1e-14)

    def test_read_mos_auto_detection(self):
        n_ao = 80
        n_mo = 30
        n_occ = 15

        np.random.seed(123)
        C_orig = np.random.randn(n_ao, n_mo)
        eps_orig = np.linspace(-0.8, 0.3, n_mo)
        occ_orig = np.zeros(n_mo)
        occ_orig[:n_occ] = 2.0

        # Test 1: extension .mbse
        file_path1 = os.path.join(self.temp_dir.name, "mos.mbse")
        write_mos_mbse(file_path1, C_orig, eps_orig, occ_orig)
        C1, eps1, occ1 = read_mos_auto(file_path1, n_ao)
        np.testing.assert_allclose(C1, C_orig)

        # Test 2: arbitrary extension (magic bytes detection)
        file_path2 = os.path.join(self.temp_dir.name, "mos_custom.dat")
        write_mos_mbse(file_path2, C_orig, eps_orig, occ_orig)
        C2, eps2, occ2 = read_mos_auto(file_path2, n_ao)
        np.testing.assert_allclose(C2, C_orig)

    def test_ao_mismatch_error(self):
        n_ao = 50
        n_mo = 20
        C = np.zeros((n_ao, n_mo))
        eps = np.zeros(n_mo)
        occ = np.zeros(n_mo)

        file_path = os.path.join(self.temp_dir.name, "mismatch.mbse")
        write_mos_mbse(file_path, C, eps, occ)

        with self.assertRaises(ValueError) as ctx:
            read_mos_mbse(file_path, n_ao_total=60)
        self.assertIn("AO count mismatch", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
