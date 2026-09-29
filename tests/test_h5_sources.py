"""HDF5 MO files as geometry source and in the readers used by NAMD and Auger."""
import os
import tempfile
import unittest

import numpy as np

try:
    import h5py  # noqa: F401
    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

from qdex.io_utils import geometry_source, read_mos_dense, read_mos_h5, write_mos_h5, write_mos_mbse, read_xyz


@unittest.skipUnless(HAVE_H5PY, "h5py not installed")
class H5SourceTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="qdex_h5src_")
        rng = np.random.default_rng(0)
        self.l = [0, 1, 0, 2]                       # 1 + 3 + 1 + 5 = 10 AOs
        self.n_ao = 10
        self.C = rng.standard_normal((self.n_ao, self.n_ao))
        self.eps = np.sort(rng.standard_normal(self.n_ao))
        self.occ = np.array([2.0] * 4 + [0.0] * 6)
        self.syms = ["Cd", "Se"]
        self.xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.63]])
        self.h5 = os.path.join(self.tmp, "orbitals.h5")
        write_mos_h5(self.h5, self.C, self.eps, self.occ, shell_ang_mom=self.l, labels=self.syms, coords_ang=self.xyz)

    def test_geometry_source(self):
        self.assertEqual(geometry_source("geom.xyz", self.h5), "geom.xyz")
        self.assertEqual(geometry_source(None, self.h5), self.h5)
        self.assertEqual(geometry_source("", self.h5), self.h5)
        self.assertIsNone(geometry_source(None, "MOs.mbse"))
        self.assertIsNone(geometry_source(None, None))
        syms, xyz = read_xyz(geometry_source(None, self.h5))
        self.assertEqual(syms, self.syms)
        np.testing.assert_allclose(xyz, self.xyz, atol=1e-12)

    def test_read_mos_dense_h5_and_mbse_agree(self):
        mbse = os.path.join(self.tmp, "MOs.mbse")
        write_mos_mbse(mbse, self.C, self.eps, self.occ)
        for path in (self.h5, mbse):
            C, eps, occ = read_mos_dense(path, self.n_ao)
            self.assertIsInstance(C, np.ndarray)
            np.testing.assert_allclose(C, self.C, atol=1e-14)
            np.testing.assert_allclose(eps, self.eps)
            np.testing.assert_allclose(occ, self.occ)

    def test_unrestricted_file_needs_a_spin_channel(self):
        path = os.path.join(self.tmp, "uks.h5")
        C2 = np.hstack([self.C, self.C])
        write_mos_h5(path, C2, np.tile(self.eps, 2), np.tile(self.occ / 2, 2), shell_ang_mom=self.l,
                     mo_spin=np.array([0] * self.n_ao + [1] * self.n_ao))
        with self.assertRaisesRegex(ValueError, "unrestricted"):
            read_mos_h5(path, self.n_ao)
        Ca, _, _ = read_mos_h5(path, self.n_ao, spin="alpha")
        np.testing.assert_allclose(Ca, self.C, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
