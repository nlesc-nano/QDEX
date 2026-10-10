import unittest

import numpy as np

from qdex.hardness import expand_atom_to_ao
from qdex.qp_levels import (atom_delta_w, cohsex_diagonal, orbital_populations, orbital_qp_energies,
                             orbital_sigma, xs_shared_w)
from qdex.hardness import MATERIAL_DB, anchor_residual_scale

RANGES = [(0, 2), (2, 3)]


def _parts():
    gamma = np.array([[5.0, 3.0], [3.0, 6.0]])
    w_bulk = gamma * np.array([[1.0, 0.3], [0.3, 1.0]])
    w_qd = gamma * np.array([[1.0, 0.5], [0.5, 1.0]])
    w_add = np.array([[0.2, 0.1], [0.1, 0.2]])
    return {"w_qd": w_qd, "w_bulk": w_bulk, "w_add": w_add, "gamma": gamma, "eps_z": 4.0, "bulk_shift": 1.0}


class QPLevelTests(unittest.TestCase):
    def test_mulliken_and_lowdin_populations_sum_to_one(self):
        rng = np.random.default_rng(0)
        A = rng.normal(size=(3, 3))
        S = A @ A.T / 3 + np.eye(3)
        w, V = np.linalg.eigh(S)
        C = V / np.sqrt(w)            # S-orthonormal orbitals
        for mode in ("mulliken", "lowdin"):
            q = orbital_populations(C, S, RANGES, mode, "atom")
            np.testing.assert_allclose(q.sum(axis=0), 1.0, atol=1e-12)

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
        q = orbital_populations(C, np.eye(3), RANGES, "mulliken", "atom")
        eps = np.array([-6.0, -5.0, -2.0])
        eps_qp, info = orbital_qp_energies(eps, q[:, :2], q[:, 2:], np.array([0, 1]), np.array([2]), dW,
                                           1.0, "fixed", 1.0, None, "CDSE")
        sig = orbital_sigma(q, dW)
        np.testing.assert_allclose(eps_qp, eps + np.array([-0.5 - sig[0], -0.5 - sig[1], 0.5 + sig[2]]), atol=1e-12)
        self.assertEqual(info["qp_levels"], "orbital")

    def test_edge_pinning_keeps_frontier_shifts(self):
        p = _parts()
        q = orbital_populations(np.eye(3), np.eye(3), RANGES, "lowdin", "atom")
        eps_qp, _ = orbital_qp_energies(np.array([-6.0, -5.0, -2.0]), q[:, :2], q[:, 2:], np.array([0, 1]),
                                        np.array([2]), atom_delta_w(p), 1.0, "fixed", 1.0,
                                        None, "CDSE", edge_shifts=(1.4, 2.0))
        self.assertAlmostEqual(eps_qp[1], -6.4, places=12)
        self.assertAlmostEqual(eps_qp[2], 0.0, places=12)

    def test_cohsex_constant_dw_is_the_classical_limit(self):
        # For a constant Delta W = c: COH = c/2 for every orbital, SEX = -c (occupied) and 0 (virtual),
        # so the gap opens by c, the classical result.
        rng = np.random.default_rng(1)
        A = rng.normal(size=(3, 3))
        S = A @ A.T / 3 + np.eye(3)
        w, V = np.linalg.eigh(S)
        C = V / np.sqrt(w)
        c = 0.7
        coh, sex = cohsex_diagonal(C, S, 0, RANGES, dW_atom=np.full((2, 2), c))
        np.testing.assert_allclose(coh, 0.5 * c, atol=1e-12)
        np.testing.assert_allclose(sex, [-c, 0.0, 0.0], atol=1e-12)

    def test_cohsex_xs_matches_atom_blocks(self):
        rng = np.random.default_rng(2)
        A = rng.normal(size=(3, 3))
        S = A @ A.T / 3 + np.eye(3)
        w, V = np.linalg.eigh(S)
        C = V / np.sqrt(w)
        dW = np.array([[0.4, 0.2], [0.2, 0.6]])
        a = cohsex_diagonal(C, S, 1, RANGES, dW_atom=dW)
        b = cohsex_diagonal(C, S, 1, RANGES, dW_ao=expand_atom_to_ao(dW, RANGES, 3))
        np.testing.assert_allclose(a, b, atol=1e-12)

    def test_cohsex_active_window_exact_match(self):
        rng = np.random.default_rng(42)
        n_ao = 12
        homo = 4
        A = rng.normal(size=(n_ao, n_ao))
        S = A @ A.T / n_ao + np.eye(n_ao)
        w, V = np.linalg.eigh(S)
        C = V / np.sqrt(w)
        # Atom ranges
        ranges = [(0, 4), (4, 8), (8, 12)]
        dW = rng.uniform(0.5, 2.0, (3, 3))
        dW = 0.5 * (dW + dW.T)

        coh_full, sex_full = cohsex_diagonal(C, S, homo, ranges, dW_atom=dW)
        eval_idx = np.array([3, 4, 5, 6])
        coh_act, sex_act = cohsex_diagonal(C, S, homo, ranges, dW_atom=dW, eval_indices=eval_idx)

        # On eval_idx, must be numerically identical to machine precision
        np.testing.assert_allclose(coh_act[eval_idx], coh_full[eval_idx], atol=1e-14)
        np.testing.assert_allclose(sex_act[eval_idx], sex_full[eval_idx], atol=1e-14)
        # Outside window, should be clamped to edge values
        self.assertEqual(coh_act[0], coh_act[3])
        self.assertEqual(sex_act[0], sex_act[3])
        self.assertEqual(coh_act[-1], coh_act[6])
        self.assertEqual(sex_act[-1], sex_act[6])

    def test_cohsex_occ_density_window(self):
        rng = np.random.default_rng(7)
        n_ao, homo = 12, 4
        A = rng.normal(size=(n_ao, n_ao))
        S = A @ A.T / n_ao + np.eye(n_ao)
        w, V = np.linalg.eigh(S)
        C = V / np.sqrt(w)
        ranges = [(0, 4), (4, 8), (8, 12)]
        dW = rng.uniform(0.5, 2.0, (3, 3))
        dW = 0.5 * (dW + dW.T)
        eval_idx = np.array([3, 4, 5, 6])

        coh_ref, sex_ref = cohsex_diagonal(C, S, homo, ranges, dW_atom=dW, eval_indices=eval_idx)
        coh_all, sex_all = cohsex_diagonal(C, S, homo, ranges, dW_atom=dW, eval_indices=eval_idx,
                                           occ_density=np.arange(homo + 1))
        np.testing.assert_allclose(coh_all[eval_idx], coh_ref[eval_idx], atol=1e-14)
        np.testing.assert_allclose(sex_all[eval_idx], sex_ref[eval_idx], atol=1e-14)

        # window {2, 3, 4}: COH unchanged, SEX_n = -sum_{i in window} q^{ni} dW q^{ni} (Löwdin transition charges)
        occ_d = np.array([2, 3, 4])
        coh_w, sex_w = cohsex_diagonal(C, S, homo, ranges, dW_atom=dW, eval_indices=eval_idx, occ_density=occ_d)
        np.testing.assert_allclose(coh_w[eval_idx], coh_ref[eval_idx], atol=1e-14)
        c = (V * np.sqrt(w)) @ V.T @ C
        for n in eval_idx:
            q = np.array([[np.sum(c[a0:a1, n] * c[a0:a1, i]) for a0, a1 in ranges] for i in occ_d])
            self.assertAlmostEqual(sex_w[n], -np.einsum("ia,ab,ib->", q, dW, q), places=12)

    def test_residual_scale_econf(self):
        e = MATERIAL_DB["CDSE"]
        gap0 = e[11] - e[10]
        self.assertAlmostEqual(anchor_residual_scale("CDSE", 20.0, gap0)[0], 1.0, places=12)
        self.assertAlmostEqual(anchor_residual_scale("CDSE", 20.0, e[7])[0], 0.0, places=12)
        mid = e[7] + 0.5 * (gap0 - e[7])
        self.assertAlmostEqual(anchor_residual_scale("CDSE", 20.0, mid)[0], 0.5, places=12)
        self.assertEqual(anchor_residual_scale("CDSE", 20.0, None)[1], "power")


if __name__ == "__main__":
    unittest.main()
