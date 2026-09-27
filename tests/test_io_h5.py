import unittest
import os
import tempfile
import numpy as np

from qdex.io_utils import (
    read_mos_h5,
    write_mos_h5,
    read_mos_auto,
    read_mos_uks,
    read_geometry_h5,
    read_xyz,
    is_h5_file,
    build_trexio_to_cp2k_ao_map,
)

CDSE_H5_PATH = "/Users/ivaninfante/Documents/University/Programs/miniBSE/CdSe/2.4nm/orbitals.h5"
CDSE_BASIS_PATH = "/Users/ivaninfante/Documents/University/Programs/miniBSE/CdSe/2.4nm/BASIS_MOLOPT_UZH"


class TestH5IO(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_is_h5_file(self):
        self.assertFalse(is_h5_file("nonexistent_file.xyz"))
        dummy_txt = os.path.join(self.temp_dir.name, "test.txt")
        with open(dummy_txt, "w") as f:
            f.write("hello world\n")
        self.assertFalse(is_h5_file(dummy_txt))

        if os.path.exists(CDSE_H5_PATH):
            self.assertTrue(is_h5_file(CDSE_H5_PATH))

    def test_build_trexio_to_cp2k_ao_map(self):
        # Shells: s (l=0, 1 AO), p (l=1, 3 AOs), d (l=2, 5 AOs)
        # Total AOs = 1 + 3 + 5 = 9
        shells = [0, 1, 2]
        mapping = build_trexio_to_cp2k_ao_map(shells)
        self.assertEqual(len(mapping), 9)
        # Every index 0..8 should appear exactly once
        self.assertEqual(sorted(mapping.tolist()), list(range(9)))

        # l=0: [0]
        # l=1: [2, 3, 1] relative to offset 1
        # l=2: [6, 7, 5, 8, 4] relative to offset 4
        expected = [0, 2, 3, 1, 6, 7, 5, 8, 4]
        self.assertEqual(mapping.tolist(), expected)

    def test_synthetic_roundtrip_trexio(self):
        # 1 s-shell, 2 p-shells, 1 d-shell, 1 f-shell
        # 1 + 3 + 3 + 5 + 7 = 19 AOs
        shell_ang_mom = [0, 1, 1, 2, 3]
        n_ao = 19
        n_mo = 8
        n_occ = 4

        np.random.seed(42)
        C_orig = np.random.randn(n_ao, n_mo)
        eps_orig = np.linspace(-0.9, -0.1, n_mo)
        occ_orig = np.zeros(n_mo)
        occ_orig[:n_occ] = 2.0

        labels = ["Cd", "Se"]
        coords_ang = np.array([[0.0, 0.0, 0.0], [1.5, 1.5, 1.5]])

        out_path = os.path.join(self.temp_dir.name, "synthetic.h5")
        write_mos_h5(
            out_path,
            C=C_orig,
            eps=eps_orig,
            occ=occ_orig,
            shell_ang_mom=shell_ang_mom,
            labels=labels,
            coords_ang=coords_ang,
        )

        self.assertTrue(is_h5_file(out_path))

        # Read MOs back
        C_read, eps_read, occ_read = read_mos_h5(out_path, n_ao_total=n_ao)
        self.assertEqual(C_read.shape, (n_ao, n_mo))
        self.assertEqual(eps_read.shape, (n_mo,))
        self.assertEqual(occ_read.shape, (n_mo,))

        np.testing.assert_allclose(C_read, C_orig, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(eps_read, eps_orig, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(occ_read, occ_orig, rtol=1e-14, atol=1e-14)

        # Read geometry back
        syms_read, coords_read = read_geometry_h5(out_path)
        self.assertEqual(syms_read, labels)
        np.testing.assert_allclose(coords_read, coords_ang, rtol=1e-10, atol=1e-10)

        # Verify read_xyz route
        syms_xyz, coords_xyz = read_xyz(out_path)
        self.assertEqual(syms_xyz, labels)
        np.testing.assert_allclose(coords_xyz, coords_ang, rtol=1e-10, atol=1e-10)

    def test_synthetic_uks_spin_channels(self):
        shell_ang_mom = [0, 1]  # 4 AOs
        n_ao = 4
        n_mo = 6  # 3 alpha, 3 beta
        C_orig = np.random.randn(n_ao, n_mo)
        eps_orig = np.array([-0.5, -0.3, 0.1, -0.45, -0.25, 0.15])
        occ_orig = np.array([1.0, 1.0, 0.0, 1.0, 1.0, 0.0])
        mo_spin = np.array([0, 0, 0, 1, 1, 1])

        out_path = os.path.join(self.temp_dir.name, "synthetic_uks.h5")
        write_mos_h5(
            out_path,
            C=C_orig,
            eps=eps_orig,
            occ=occ_orig,
            shell_ang_mom=shell_ang_mom,
            mo_spin=mo_spin,
        )

        # Read alpha
        Ca, ea, oa = read_mos_h5(out_path, n_ao, spin="alpha")
        self.assertEqual(Ca.shape, (n_ao, 3))
        np.testing.assert_allclose(Ca, C_orig[:, :3])
        np.testing.assert_allclose(ea, eps_orig[:3])
        np.testing.assert_allclose(oa, occ_orig[:3])

        # Read beta
        Cb, eb, ob = read_mos_h5(out_path, n_ao, spin="beta")
        self.assertEqual(Cb.shape, (n_ao, 3))
        np.testing.assert_allclose(Cb, C_orig[:, 3:])
        np.testing.assert_allclose(eb, eps_orig[3:])
        np.testing.assert_allclose(ob, occ_orig[3:])

        # Test read_mos_uks with single file
        C_a, e_a, o_a, C_b, e_b, o_b = read_mos_uks(out_path, out_path, n_ao)
        np.testing.assert_allclose(C_a, Ca)
        np.testing.assert_allclose(C_b, Cb)

    def test_auto_detection_and_ao_mismatch(self):
        n_ao = 10
        n_mo = 4
        C_orig = np.random.randn(n_ao, n_mo)
        eps_orig = np.linspace(-0.5, 0.1, n_mo)
        occ_orig = np.array([2.0, 2.0, 0.0, 0.0])

        # Save with arbitrary extension
        custom_path = os.path.join(self.temp_dir.name, "orbitals.custom_ext")
        write_mos_h5(custom_path, C_orig, eps_orig, occ_orig, is_cartesian=True)

        C_read, eps_read, occ_read = read_mos_auto(custom_path, n_ao_total=n_ao)
        np.testing.assert_allclose(C_read, C_orig)

        # Test AO count mismatch error
        with self.assertRaises(ValueError) as ctx:
            read_mos_auto(custom_path, n_ao_total=15)
        self.assertIn("AO count mismatch", str(ctx.exception))

    @unittest.skipUnless(os.path.exists(CDSE_H5_PATH), "CdSe 2.4nm orbitals.h5 not available")
    def test_cdse_2_4nm_file(self):
        # 1. Read geometry
        syms, coords_ang = read_geometry_h5(CDSE_H5_PATH)
        self.assertEqual(len(syms), 381)
        self.assertEqual(coords_ang.shape, (381, 3))
        self.assertEqual(syms.count("Cd"), 176)
        self.assertEqual(syms.count("Se"), 147)
        self.assertEqual(syms.count("Cl"), 58)

        # 2. Read MOs via read_mos_h5 and read_mos_auto
        C, eps, occ = read_mos_auto(CDSE_H5_PATH, n_ao_total=7065, verbose=False)
        self.assertEqual(C.shape, (7065, 1700))
        self.assertEqual(eps.shape, (1700,))
        self.assertEqual(occ.shape, (1700,))

        # All 1700 MOs are occupied
        self.assertTrue(np.all(occ == 2.0))
        self.assertAlmostEqual(eps[0], -0.76308301, places=5)
        self.assertAlmostEqual(eps[-1], -0.20612560, places=5)

        # 3. Test MO orthonormality with libint overlap matrix if basis exists
        if os.path.exists(CDSE_BASIS_PATH):
            import libint_cpp
            from qdex.io_utils import parse_basis, build_shell_dicts

            basis_dict = parse_basis(CDSE_BASIS_PATH, "DZVP-MOLOPT-PBE-GTH", required_elements=set(syms))
            shells = build_shell_dicts(syms, coords_ang, basis_dict)
            shells = [{**sh, "pure": True} for sh in shells]
            S = libint_cpp.overlap(shells, 8)

            # Test highest 10 occupied MOs
            C_sub = C[:, -10:]
            overlap_mo = C_sub.T @ S @ C_sub
            np.testing.assert_allclose(np.diag(overlap_mo), np.ones(10), atol=1e-6)
            off_diag = overlap_mo - np.diag(np.diag(overlap_mo))
            self.assertLess(np.max(np.abs(off_diag)), 1e-6)


if __name__ == "__main__":
    unittest.main()
