import os
import tempfile
import unittest

import numpy as np

import libint_cpp
from qdex.io_utils import read_mos_dense, read_xyz, count_ao_from_shells, build_atom_ao_ranges
from qdex.xtb.molden import build_frame_shells, convert_molden, parse_molden, molden_to_spherical

DATA = os.path.join(os.path.dirname(__file__), "data", "xtb")


class XtbMoldenTests(unittest.TestCase):
    """xtb molden -> per-atom basis + MOs.mbse (CdSeCl2, Cartesian d shells on Se and Cl)."""

    def _roundtrip(self, name):
        with tempfile.TemporaryDirectory() as tmp:
            info = convert_molden(os.path.join(DATA, name), tmp)
            syms, X = read_xyz(os.path.join(tmp, "frame.xyz"))
            shells = build_frame_shells(os.path.join(tmp, "BASIS_GXTB"), syms, X)
            C, eps, occ = read_mos_dense(os.path.join(tmp, "MOs.mbse"), count_ao_from_shells(shells))
            S = libint_cpp.cross_overlap_geometries(shells, shells, 1)
        return info, syms, shells, C, occ, S

    def test_gxtb_conventions_and_orthonormality(self):
        # The g-xTB (tblite) writer: raw primitives, axis-normalized Cartesian components
        info, syms, shells, C, occ, S = self._roundtrip("cdsecl2_gxtb.molden")
        self.assertEqual(info["primitive_convention"], "raw")
        self.assertEqual(info["cartesian_convention"], "axis")
        self.assertEqual(C.shape, (count_ao_from_shells(shells), C.shape[1]))
        np.testing.assert_allclose(C.T @ S @ C, np.eye(C.shape[1]), atol=1e-7)

    def test_gfn2_conventions_and_orthonormality(self):
        # The classic xtb (GFN2) writer: raw primitives, individually normalized components
        info, _syms, _shells, C, _occ, S = self._roundtrip("cdsecl2_gfn2.molden")
        self.assertEqual(info["primitive_convention"], "raw")
        self.assertEqual(info["cartesian_convention"], "component")
        np.testing.assert_allclose(C.T @ S @ C, np.eye(C.shape[1]), atol=1e-6)

    def test_mulliken_charges_match_xtb(self):
        info, syms, shells, C, occ, S = self._roundtrip("cdsecl2_gxtb.molden")
        P = (C * occ) @ C.T
        PS = P @ S
        pop = np.array([np.trace(PS[a:b, a:b]) for a, b in build_atom_ao_ranges(shells)])
        zval = {"Cd": 2, "Se": 6, "Cl": 7}
        q = np.array([zval[s] for s in syms]) - pop
        np.testing.assert_allclose(q, np.loadtxt(os.path.join(DATA, "cdsecl2_gxtb.charges")), atol=1e-4)

    def test_per_atom_basis_is_kept(self):
        # q-vSZP contractions depend on the atomic charge: the two Cl atoms must keep their own
        mol = parse_molden(os.path.join(DATA, "cdsecl2_gxtb.molden"))
        atom_shells, _, _, _ = molden_to_spherical(mol)
        with tempfile.TemporaryDirectory() as tmp:
            convert_molden(os.path.join(DATA, "cdsecl2_gxtb.molden"), tmp)
            syms, X = read_xyz(os.path.join(tmp, "frame.xyz"))
            shells = build_frame_shells(os.path.join(tmp, "BASIS_GXTB"), syms, X)
        cl = [sh for sh in shells if sh["sym"] == "Cl" and sh["l"] == 0]
        self.assertEqual(len(cl), 2)
        for sh, ref in zip(cl, [atom_shells[2][0], atom_shells[3][0]]):
            np.testing.assert_allclose(sh["coefs"], ref[2], rtol=1e-12)


if __name__ == "__main__":
    unittest.main()
