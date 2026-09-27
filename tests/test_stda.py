import unittest

import numpy as np

from qdex.hardness import HARDNESS_DICT, HA_TO_EV, MATERIAL_DB, build_stda_gammas, stda_ax


class STDATests(unittest.TestCase):
    def test_parameters_and_limits(self):
        coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]])
        gj, gk, info = build_stda_gammas(["C", "C"], coords, 0.25)
        self.assertAlmostEqual(info["beta_J"], 0.20 + 1.83 * 0.25)
        self.assertAlmostEqual(info["alpha_K"], 1.42 + 0.48 * 0.25)
        eta = 2.0 * HARDNESS_DICT["c"]                    # Grimme: (ii|ii) = IP - EA = 2 eta
        self.assertAlmostEqual(gk[0, 0], eta, places=10)   # R = 0: gamma^K = eta
        self.assertAlmostEqual(gj[0, 0], 0.25 * eta, places=10)
        far = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 500.0]])
        gj, gk, _ = build_stda_gammas(["C", "C"], far, 0.25)
        r_bohr = 500.0 / 0.529177210903
        self.assertAlmostEqual(gk[0, 1] * r_bohr / HA_TO_EV, 1.0, places=2)   # 1/R tail
        # gamma^J reaches 1/R slowly (beta = 0.66): 7 % low at 500 A, 0.4 % at 50,000 A.
        gj, _, _ = build_stda_gammas(["C", "C"], 100.0 * far, 0.25)
        self.assertAlmostEqual(gj[0, 1] * 100.0 * r_bohr / HA_TO_EV, 1.0, places=2)

    def test_no_attraction_for_pure_functionals(self):
        gj, gk, _ = build_stda_gammas(["C", "O"], np.array([[0, 0, 0], [0, 0, 1.2]]), 0.0)
        self.assertTrue(np.all(gj == 0.0))
        self.assertTrue(np.all(gk > 0.0))

    def test_ax_resolution(self):
        self.assertEqual(stda_ax("PBE0")[0], 0.25)
        self.assertEqual(stda_ax("pbe")[0], 0.0)
        self.assertEqual(stda_ax(None, "0.3")[0], 0.3)
        self.assertAlmostEqual(stda_ax(None, "dielectric", "CDSE")[0], 1.0 / MATERIAL_DB["CDSE"][0])
        with self.assertRaises(ValueError):
            stda_ax(None, None)
        with self.assertRaises(ValueError):
            stda_ax("mystery")


if __name__ == "__main__":
    unittest.main()
