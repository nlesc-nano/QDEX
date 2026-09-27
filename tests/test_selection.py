import unittest

import numpy as np

from qdex.selection import HARTREE_EV, perturbative_selection


def _model(seed=0, n_o=4, n_v=5, m=3):
    rng = np.random.default_rng(seed)
    D = 2.0 + np.add.outer(np.arange(n_o)[::-1] * 0.4, np.arange(n_v) * 0.5)
    q_ov = 0.3 * rng.normal(size=(n_o, n_v, m))
    q_occ = 0.3 * rng.normal(size=(n_o, n_o, m)); q_occ = 0.5 * (q_occ + q_occ.transpose(1, 0, 2))
    q_virt = 0.3 * rng.normal(size=(n_v, n_v, m)); q_virt = 0.5 * (q_virt + q_virt.transpose(1, 0, 2))
    g = rng.normal(size=(m, m)); g = g @ g.T + m * np.eye(m)
    W = 0.3 * g
    W_virt = q_virt @ W.T
    return D, q_ov, q_occ, W_virt, g


def _A(D, q_ov, q_occ, W_virt, g, kx):
    n_o, n_v = D.shape
    A = np.zeros((n_o * n_v, n_o * n_v))
    for i in range(n_o):
        for a in range(n_v):
            for j in range(n_o):
                for b in range(n_v):
                    A[i * n_v + a, j * n_v + b] = (kx * q_ov[i, a] @ g @ q_ov[j, b]
                                                   - q_occ[i, j] @ W_virt[a, b]
                                                   + (D[i, a] if (i, a) == (j, b) else 0.0))
    return A


class SelectionTests(unittest.TestCase):
    def test_matches_brute_force(self):
        D, q_ov, q_occ, W_virt, g = _model()
        A = _A(D, q_ov, q_occ, W_virt, g, 2.0)
        diag = np.diag(A)
        e_thr = np.sort(diag)[6]
        t = 1e-4
        mask, shift, info = perturbative_selection(D, q_ov, q_occ, W_virt, g, 2.0, True, e_thr, t)
        prim = diag <= e_thr
        cand = ~prim
        pt = A[np.ix_(cand, prim)] ** 2 / (diag[cand][:, None] - diag[prim][None, :] + 1e-10 * HARTREE_EV)
        keep = pt.sum(1) > t * HARTREE_EV
        ref = prim.copy()
        ref[np.flatnonzero(cand)[keep]] = True
        np.testing.assert_array_equal(mask.ravel(), ref)
        np.testing.assert_allclose(shift.ravel()[prim], -pt[~keep].sum(0), atol=1e-12)
        self.assertEqual(info["n_primary"], int(prim.sum()))

    def test_large_threshold_keeps_everything(self):
        D, q_ov, q_occ, W_virt, g = _model(1)
        mask, shift, info = perturbative_selection(D, q_ov, q_occ, W_virt, g, 2.0, True, 1e3)
        self.assertTrue(mask.all())
        self.assertTrue(np.all(shift == 0.0))

    def test_empty_primary_space_raises(self):
        D, q_ov, q_occ, W_virt, g = _model(2)
        with self.assertRaises(ValueError):
            perturbative_selection(D, q_ov, q_occ, W_virt, g, 2.0, True, -100.0)


if __name__ == "__main__":
    unittest.main()
