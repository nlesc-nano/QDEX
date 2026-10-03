import os
import tempfile
import unittest

import numpy as np

import libint_cpp
from qdex.config_schema import flatten_config
from qdex.hardness import build_stda_gammas, stda_parameters
from qdex.io_utils import build_atom_ao_ranges, count_ao_from_shells, read_mos_dense, read_xyz
from qdex.namd.precompute import compute_frame_diagonal_bse
from qdex.xtb.molden import build_frame_shells, convert_molden

DATA = os.path.join(os.path.dirname(__file__), "data", "xtb")
HA_TO_EV = 27.211386245988


class StdaParameterTests(unittest.TestCase):
    def test_range_separated_sets(self):
        # stda main.f: -wB97MV / -wB97XD3 (a_x 0.51, alpha 4.51, beta 8.0), -CAMB3LYP (0.38, 0.90, 1.86)
        for name in ("gxtb", "wB97M-V", "wb97x-d3"):
            ax, a, b, _ = stda_parameters(name)
            self.assertEqual((ax, a, b), (0.51, 4.51, 8.0))
        self.assertEqual(stda_parameters("cam-b3lyp")[:3], (0.38, 0.90, 1.86))

    def test_xtb_set(self):
        # stda main.f, -xtb: a_x 0.50, alpha 2.0, beta 4.0; 'gfn2' borrows it (no GFN2 set exists)
        for name in ("stda-xtb", "gfn2", "GFN2-xTB"):
            ax, a, b, src = stda_parameters(name)
            self.assertEqual((ax, a, b), (0.50, 2.0, 4.0))
            self.assertIn("sTDA-xTB set", src)

    def test_global_hybrid_formulas(self):
        ax, a, b, _ = stda_parameters("pbe0")
        self.assertAlmostEqual(ax, 0.25)
        self.assertAlmostEqual(a, 1.42 + 0.48 * 0.25)
        self.assertAlmostEqual(b, 0.20 + 1.83 * 0.25)

    def test_explicit_values_override(self):
        ax, a, b, src = stda_parameters("gxtb", ax=1.0, alpha=2.0, beta=3.0)
        self.assertEqual((ax, a, b), (1.0, 2.0, 3.0))
        self.assertIn("explicit", src)
        ax, a, b, _ = stda_parameters(None, ax=1.0)
        self.assertAlmostEqual(a, 1.90)
        self.assertAlmostEqual(b, 2.03)

    def test_yaml_keys(self):
        flat = flatten_config({"excitations": {"functional": "gxtb", "stda_alpha": 4.0, "stda_beta": 7.0}})
        self.assertEqual((flat["stda_functional"], flat["stda_alpha"], flat["stda_beta"]), ("gxtb", 4.0, 7.0))

    def test_gammas_use_alpha_beta(self):
        syms, X = ["Cd", "Se"], np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.6]])
        gj, gk, info = build_stda_gammas(syms, X, 0.51, alpha=4.51, beta=8.0)
        self.assertEqual((info["alpha_K"], info["beta_J"]), (4.51, 8.0))
        gj2, gk2, _ = build_stda_gammas(syms, X, 0.51)
        self.assertFalse(np.allclose(gj, gj2))


class DiagonalStdaPrecomputeTests(unittest.TestCase):
    """NAMD diagonal_stda: E_ia = de + 2 (ia|ia)_K - (ii|aa)_J with Loewdin charges and Grimme gammas."""

    def test_matches_explicit_formula(self):
        ax, a, b, _ = stda_parameters("gxtb")
        with tempfile.TemporaryDirectory() as tmp:
            convert_molden(os.path.join(DATA, "cdsecl2_gxtb.molden"), tmp)
            xyz, basis, mo = (os.path.join(tmp, f) for f in ("frame.xyz", "BASIS_GXTB", "MOs.mbse"))
            syms, X = read_xyz(xyz)
            shells = build_frame_shells(basis, syms, X)
            n_ao, ranges = count_ao_from_shells(shells), build_atom_ao_ranges(shells)
            res = compute_frame_diagonal_bse(
                xyz_path=xyz, mo_path=mo, basis_dict=None, n_ao=n_ao, atom_ao_ranges=ranges, scissor=0.0,
                nhomos=3, nlumos=3, excitation_mode="diagonal_stda", material="CDSE", compute_dipoles=False,
                frame_basis_path=basis, stda=dict(ax=ax, alpha=a, beta=b))
            C, eps, occ = read_mos_dense(mo, n_ao)

        # Reference: Loewdin charges from S^1/2 C, sTDA gammas
        S = libint_cpp.cross_overlap_geometries(shells, shells, 1)
        w, V = np.linalg.eigh(S)
        Xl = (V * np.sqrt(w)) @ V.T @ C
        n_occ = int((occ > 0.5).sum())
        o, v = np.arange(n_occ - 3, n_occ), np.arange(n_occ, n_occ + 3)
        q = lambda p, r: np.array([np.sum(Xl[a0:a1, p] * Xl[a0:a1, r]) for a0, a1 in ranges])
        gj, gk, _ = build_stda_gammas(syms, X, ax, alpha=a, beta=b)
        ref = []
        for i in o:
            for c in v:
                de = (eps[c] - eps[i]) * HA_TO_EV
                ref.append(de + 2.0 * q(i, c) @ gk @ q(i, c) - q(i, i) @ gj @ q(c, c))
        np.testing.assert_allclose(res["E_pairs"], np.array(ref), atol=1e-8)


if __name__ == "__main__":
    unittest.main()
