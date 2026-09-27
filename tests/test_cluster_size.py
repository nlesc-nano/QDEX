import unittest

import numpy as np

from qdex.cluster_size import cluster_size


def _sphere(n, radius, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.uniform(-radius, radius, (6 * n, 3))
    return x[np.linalg.norm(x, axis=1) < radius][:n]


class ClusterSizeTests(unittest.TestCase):
    def test_uniform_sphere(self):
        x = _sphere(3000, 20.0)                       # 4 nm diameter
        d = cluster_size(x, ["Cd", "Se"] * 1500, "CDSE")
        self.assertAlmostEqual(d["d_saxs_nm"], 4.0, delta=0.1)
        self.assertAlmostEqual(d["d_guinier_nm"], 4.0, delta=0.1)

    def test_organic_shell_is_invisible(self):
        core = _sphere(1000, 10.0, 1)
        shell = _sphere(4000, 18.0, 2)
        shell = shell[np.linalg.norm(shell, axis=1) > 12.0]
        syms = ["Cd", "Se"] * 500 + ["C"] * len(shell)
        d = cluster_size(np.vstack([core, shell]), syms, "CDSE")
        self.assertAlmostEqual(d["d_saxs_nm"], 2.0, delta=0.1)
        self.assertEqual(d["n_inorganic"], 1000)
        # the user can make the ligands visible
        d2 = cluster_size(np.vstack([core, shell]), syms, "CDSE", inorganic_elements=["Cd", "Se", "C"])
        self.assertGreater(d2["d_saxs_nm"], 2.5)


if __name__ == "__main__":
    unittest.main()
