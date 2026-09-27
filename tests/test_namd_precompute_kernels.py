import unittest

import numpy as np

from qdex.namd.precompute import _diag_exchange


class DiagExchangeTests(unittest.TestCase):
    def test_matches_explicit_mulliken_contraction(self):
        rng = np.random.default_rng(0)
        n_ao, n_o, n_v = 7, 3, 4
        ranges = [(0, 3), (3, 5), (5, 7)]
        A = rng.normal(size=(n_ao, n_ao))
        S = A @ A.T / n_ao + np.eye(n_ao)
        C = rng.normal(size=(n_ao, n_o + n_v))
        Co, Cv = C[:, :n_o], C[:, n_o:]
        g = rng.normal(size=(3, 3)); g = g @ g.T + 3 * np.eye(3)
        kx = _diag_exchange([Co], [Cv], [S @ Co], [S @ Cv], ranges, g, chunk_elems=5)  # force chunking
        q = np.zeros((n_o, n_v, 3))
        for A_, (a0, a1) in enumerate(ranges):
            q[:, :, A_] = 0.5 * (Co[a0:a1].T @ (S @ Cv)[a0:a1] + (S @ Co)[a0:a1].T @ Cv[a0:a1])
        ref = np.einsum("iam,mn,ian->ia", q, g, q)
        np.testing.assert_allclose(kx, ref, atol=1e-12)

    def test_spin_components_add_before_contraction(self):
        rng = np.random.default_rng(1)
        n_ao, ranges = 4, [(0, 2), (2, 4)]
        S = np.eye(n_ao)
        Co = rng.normal(size=(n_ao, 2)); Cv = rng.normal(size=(n_ao, 2))
        g = np.eye(2)
        one = _diag_exchange([Co], [Cv], [Co], [Cv], ranges, g)
        two = _diag_exchange([Co, Co], [Cv, Cv], [Co, Co], [Cv, Cv], ranges, g)
        np.testing.assert_allclose(two, 4.0 * one, atol=1e-12)   # doubled charge -> 4x


if __name__ == "__main__":
    unittest.main()
