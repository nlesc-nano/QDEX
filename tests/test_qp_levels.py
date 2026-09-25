import unittest

import numpy as np

from qdex.hardness import expand_atom_to_ao
from qdex.qp_levels import atom_delta_w, orbital_qp_energies, orbital_sigma, xs_shared_w

RANGES = [(0, 2), (2, 3)]


def _parts():
    gamma = np.array([[5.0, 3.0], [3.0, 6.0]])
    w_bulk = gamma * np.array([[1.0, 0.3], [0.3, 1.0]])
    w_qd = gamma * np.array([[1.0, 0.5], [0.5, 1.0]])
    w_add = np.array([[0.2, 0.1], [0.1, 0.2]])
    return {"w_qd": w_qd, "w_bulk": w_bulk, "w_add": w_add, "gamma": gamma, "eps_z": 4.0, "bulk_shift": 1.0}


class QPLevelTests(unittest.TestCase):
    def test_xs_kernel_reduces_to_mnok_blocks(self):
        # With AO integrals equal to the expanded MNOK gamma, xs W is the expanded atom W.
        p = _parts()
        gamma_ao = expand_atom_to_ao(p["gamma"], RANGES, 3)
        W_ao, _, _ = xs_shared_w(p, 1.0, gamma_ao, RANGES)
        np.testing.assert_allclose(W_ao, expand_atom_to_ao(p["w_qd"] + p["w_add"], RANGES, 3), atol=1e-12)

    def test_orbital_shifts_follow_the_formula(self):
        p = _parts()
        dW = atom_delta_w(p)
        C = np.eye(3)
        eps = np.array([-6.0, -5.0, -2.0])
        eps_qp, info = orbital_qp_energies(eps, C[:, :2], C[:, 2:], np.array([0, 1]), np.array([2]), dW, RANGES,
                                           "atom", 1.0, "fixed", 1.0, None, "CDSE")
        sig = orbital_sigma(C, dW, RANGES, "atom")
        np.testing.assert_allclose(eps_qp, eps + np.array([-0.5 - sig[0], -0.5 - sig[1], 0.5 + sig[2]]), atol=1e-12)
        self.assertEqual(info["qp_levels"], "orbital")

    def test_edge_pinning_keeps_frontier_shifts(self):
        p = _parts()
        C = np.eye(3)
        eps_qp, _ = orbital_qp_energies(np.array([-6.0, -5.0, -2.0]), C[:, :2], C[:, 2:], np.array([0, 1]),
                                        np.array([2]), atom_delta_w(p), RANGES, "atom", 1.0, "fixed", 1.0,
                                        None, "CDSE", edge_shifts=(1.4, 2.0))
        self.assertAlmostEqual(eps_qp[1], -6.4, places=12)
        self.assertAlmostEqual(eps_qp[2], 0.0, places=12)


if __name__ == "__main__":
    unittest.main()
