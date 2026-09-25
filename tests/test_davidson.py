import unittest

import numpy as np

from qdex.davidson import davidson


class DavidsonTests(unittest.TestCase):
    def _problem(self, n=1500, seed=1):
        rng = np.random.default_rng(seed)
        diag = np.sort(rng.uniform(3.0, 6.0, n))
        off = rng.normal(0.0, 0.02, (n, n))
        A = np.diag(diag) + 0.5 * (off + off.T)
        return A

    def test_matches_dense_at_tight_tolerance(self):
        # Regression: an absolute 1e-5 linear-independence threshold on
        # unnormalized corrections stalled the solver near tol ~ 1e-5.
        A = self._problem()
        ref = np.linalg.eigvalsh(A)[:6]
        vals, vecs = davidson(lambda x: A @ x, np.diag(A), 6, tol=1e-7, max_iter=300)
        np.testing.assert_allclose(np.real(vals), ref, atol=1e-9)
        res = A @ vecs - vecs * np.real(vals)
        self.assertLess(np.max(np.linalg.norm(res, axis=0)), 1e-7)

    def test_degenerate_diagonal(self):
        n = 400
        diag = np.repeat(np.linspace(1.0, 4.0, n // 4), 4)
        A = np.diag(diag) + 1e-3 * np.eye(n, k=1) + 1e-3 * np.eye(n, k=-1)
        ref = np.linalg.eigvalsh(A)[:8]
        vals, _ = davidson(lambda x: A @ x, diag, 8, tol=1e-8, max_iter=300)
        np.testing.assert_allclose(np.real(vals), ref, atol=1e-10)


if __name__ == "__main__":
    unittest.main()
