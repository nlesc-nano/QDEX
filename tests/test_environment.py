import unittest

import numpy as np

from qdex.environment import (
    COULOMB_EV_ANG,
    environment_qp_energies,
    polarization_self_energies,
    reaction_field_matrix,
    sphere_cavity,
)


class ReactionFieldTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(4)
        v = rng.normal(size=(40, 3))
        self.coords = 6.0 * v / np.linalg.norm(v, axis=1)[:, None] * rng.uniform(0, 1, (40, 1)) ** (1 / 3)
        self.R = 8.0

    def test_vanishes_for_matched_media(self):
        W = reaction_field_matrix(self.coords, 6.2, 6.2, self.R, center=np.zeros(3))
        self.assertEqual(np.abs(W).max(), 0.0)

    def test_born_term_at_center(self):
        W = reaction_field_matrix(np.zeros((1, 3)), 6.2, 2.0, self.R, center=np.zeros(3))
        self.assertAlmostEqual(W[0, 0], COULOMB_EV_ANG * (1 / 2.0 - 1 / 6.2) / self.R, places=10)

    def test_kelvin_image_for_conducting_exterior(self):
        # eps_in = 1, eps_out -> infinity: image charge -R/r' at R^2 r'/|r'|^2.
        a = np.array([[1.0, 2.0, -0.5]])
        b = np.array([[-2.0, 0.5, 3.0]])
        W = reaction_field_matrix(np.vstack([a, b]), 1.0, 1.0e9, self.R, center=np.zeros(3), tol_ev=1e-12)
        rb = np.linalg.norm(b)
        image = (self.R ** 2 / rb ** 2) * b[0]
        exact = -COULOMB_EV_ANG * (self.R / rb) / np.linalg.norm(a[0] - image)
        self.assertAlmostEqual(W[0, 1], exact, places=6)

    def test_symmetric_and_repulsive_for_low_exterior_permittivity(self):
        W = reaction_field_matrix(self.coords, 6.2, 1.0, self.R, center=np.zeros(3))
        np.testing.assert_allclose(W, W.T, atol=1e-12)
        q = np.full(len(self.coords), 1.0 / len(self.coords))
        self.assertGreater(polarization_self_energies(q[None, :], W)[0], 0.0)

    def test_neutral_pair_cancellation_identity(self):
        W = reaction_field_matrix(self.coords, 6.2, 1.0, self.R, center=np.zeros(3))
        rng = np.random.default_rng(0)
        qe = rng.random(len(self.coords)); qe /= qe.sum()
        qh = rng.random(len(self.coords)); qh /= qh.sum()
        # self-energies (QP gap) minus e-h reaction attraction (BSE) = 1/2 (qe-qh)^T W (qe-qh)
        lhs = 0.5 * qe @ W @ qe + 0.5 * qh @ W @ qh - qh @ W @ qe
        self.assertAlmostEqual(lhs, 0.5 * (qe - qh) @ W @ (qe - qh), places=12)
        self.assertLess(abs(lhs), 0.1 * (0.5 * qe @ W @ qe + 0.5 * qh @ W @ qh))

    def test_atoms_outside_cavity_are_rejected(self):
        with self.assertRaises(ValueError):
            reaction_field_matrix(np.array([[9.0, 0, 0]]), 6.2, 1.0, self.R, center=np.zeros(3))
        c, R = sphere_cavity(self.coords, buffer_ang=1.0)
        self.assertGreater(R, np.max(np.linalg.norm(self.coords - c, axis=1)))


class EnvironmentQPTests(unittest.TestCase):
    def test_qp_levels_open_by_bulk_shift_plus_polarization(self):
        # Two orthonormal "orbitals" on two atoms, one AO per atom.
        coords = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
        C = np.eye(2)
        S = np.eye(2)
        ranges = [(0, 1), (1, 2)]
        W = reaction_field_matrix(coords, 6.2, 1.0, 5.0, center=np.zeros(3))
        eps = np.array([-1.0, 1.0])
        eps_qp, det = environment_qp_energies(eps, C, S, ranges, 0, W, 1.27, n_window=5)
        sig = 0.5 * W[0, 0]
        self.assertAlmostEqual(eps_qp[0], -1.0 - sig, places=12)
        self.assertAlmostEqual(eps_qp[1], 1.0 + 1.27 + 0.5 * W[1, 1], places=12)
        self.assertAlmostEqual(det["sigma_pol_homo_ev"], sig, places=12)


if __name__ == "__main__":
    unittest.main()
