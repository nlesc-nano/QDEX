import unittest
from unittest.mock import patch

import numpy as np

from miniBSE.hardness import (
    MATERIAL_DB,
    build_resta_mnok,
    estimate_gw_qp_gap,
    get_cluster_size_metrics,
    compute_radiative_rates,
    compute_energy_gap_law_rate,
    compute_fcwd_rate,
    extract_recombination_parameters_from_namd,
)


class GeometryAndQPModelTests(unittest.TestCase):
    def setUp(self):
        self.coords = np.array(
            [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.5]]
        )
        self.symbols = ["Cd", "Se", "Cd", "Se"]

    def test_radius_is_rigid_motion_and_permutation_invariant(self):
        angle = 0.731
        rotation = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0.0],
                [np.sin(angle), np.cos(angle), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        permutation = np.array([2, 0, 3, 1])
        reference = get_cluster_size_metrics(self.coords, self.symbols, "CDSE")
        moved = self.coords[permutation] @ rotation.T + np.array([7.0, -3.0, 2.0])
        symbols = [self.symbols[i] for i in permutation]
        transformed = get_cluster_size_metrics(moved, symbols, "CDSE")
        self.assertAlmostEqual(reference["R_eff_hull"], transformed["R_eff_hull"], places=12)
        np.testing.assert_allclose(
            sorted(reference["principal_extents_ang"]),
            sorted(transformed["principal_extents_ang"]),
            atol=1.0e-12,
        )

    def test_vacuum_anchor_is_reproduced(self):
        r0 = MATERIAL_DB["CDSE"][9]
        # Scale a tetrahedron so its equivalent-volume radius plus the documented
        # 1.25 A surface offset equals the database anchor radius.
        base = get_cluster_size_metrics(self.coords, self.symbols, "CDSE")
        geometric_radius = base["R_eff_hull"] - base["surface_offset_ang"]
        scale = (r0 - base["surface_offset_ang"]) / geometric_radius
        coords = self.coords * scale
        shift, details = estimate_gw_qp_gap(
            coords, self.symbols, "CDSE", eps_out=1.0, return_details=True
        )
        expected = details["monomer_gw_gap_ev"] - details["monomer_pbe_gap_ev"]
        self.assertAlmostEqual(shift, expected, places=10)

    def test_resta_uses_eps_in_not_eps_out(self):
        _, w_vac = build_resta_mnok(self.symbols, self.coords, 1.0, "CDSE", eps_out=1.0)
        _, w_solvent = build_resta_mnok(self.symbols, self.coords, 1.0, "CDSE", eps_out=8.0)
        np.testing.assert_allclose(w_vac, w_solvent, atol=0.0, rtol=0.0)

    def _qp_at_radius(self, radius, eps_out=1.0):
        metrics = {
            "R_eff_hull": float(radius),
            "anisotropy": 1.0,
            "anisotropy_ratio": 1.0,
            "principal_extents_ang": [2.0 * radius] * 3,
            "selected_atom_indices": [0, 1],
            "radius_definition_version": "test",
            "surface_offset_ang": 0.0,
        }
        with patch("miniBSE.hardness.get_cluster_size_metrics", return_value=metrics):
            return estimate_gw_qp_gap(
                self.coords, self.symbols, "CDSE", eps_out=eps_out, return_details=True
            )

    def test_qp_model_has_bulk_limit(self):
        shift, details = self._qp_at_radius(1.0e9, eps_out=1.0)
        self.assertAlmostEqual(shift, details["bulk_gw_shift_ev"], places=7)

    def test_dielectric_matching_removes_only_asymptotic_term(self):
        eps_inf = MATERIAL_DB["CDSE"][0]
        shift, details = self._qp_at_radius(2.0 * MATERIAL_DB["CDSE"][9], eps_out=eps_inf)
        expected = details["bulk_gw_shift_ev"] + details["anchor_residual_ev"] / 4.0
        self.assertAlmostEqual(shift, expected, places=12)
        self.assertAlmostEqual(details["kappa_solvent_ev_ang"], 0.0, places=12)

    def test_qp_model_is_continuous_at_anchor_from_valid_domain(self):
        r0 = MATERIAL_DB["CDSE"][9]
        at_anchor, _ = self._qp_at_radius(r0, eps_out=1.0)
        just_above, _ = self._qp_at_radius(r0 + 1.0e-7, eps_out=1.0)
        self.assertLess(abs(just_above - at_anchor), 1.0e-6)


class RecombinationModelTests(unittest.TestCase):
    def test_radiative_rates_scaling(self):
        k_s, k_fs = compute_radiative_rates(2.0, 1.0, refractive_index=2.0)
        self.assertGreater(k_s, 0.0)
        self.assertAlmostEqual(k_fs, k_s * 1e-15, places=20)
        # Double energy -> 4x rate
        k_s2, _ = compute_radiative_rates(4.0, 1.0, refractive_index=2.0)
        self.assertAlmostEqual(k_s2 / k_s, 4.0, places=5)

    def test_energy_gap_law_decay(self):
        k1_s, _ = compute_energy_gap_law_rate(0.2, E_LO_ev=0.018, S_hr=1.0)
        k2_s, _ = compute_energy_gap_law_rate(0.4, E_LO_ev=0.018, S_hr=1.0)
        self.assertGreater(k1_s, 0.0)
        self.assertGreater(k2_s, 0.0)
        self.assertGreater(k1_s, k2_s)  # larger gap -> slower non-radiative decay

    def test_extract_recombination_parameters(self):
        var_g = 0.05 ** 2  # 50 meV fluctuation
        params = extract_recombination_parameters_from_namd(
            var_E_gap_ev2=var_g,
            dominant_freq_cm1=150.0,
            temp_k=300.0,
            mean_nac_fs=0.001
        )
        self.assertAlmostEqual(params["dominant_freq_cm1"], 150.0)
        self.assertAlmostEqual(params["sigma_ev"], 0.05, places=5)
        self.assertGreater(params["E_LO_ev"], 0.015)
        self.assertGreater(params["lambda_ev"], 0.0)
        self.assertGreater(params["S_hr"], 0.0)
        self.assertIsNotNone(params["V_el_ev"])
        self.assertGreater(params["V_el_ev"], 0.0)

    def test_fcwd_rate_scaling(self):
        k_s, k_fs = compute_fcwd_rate(E_gap_ev=0.5, V_el_ev=0.001, lambda_ev=0.05, sigma_ev=0.05)
        self.assertGreater(k_s, 0.0)
        self.assertAlmostEqual(k_fs, k_s * 1e-15, places=20)
        # Double V_el -> 4x rate
        k_s2, _ = compute_fcwd_rate(E_gap_ev=0.5, V_el_ev=0.002, lambda_ev=0.05, sigma_ev=0.05)
        self.assertAlmostEqual(k_s2 / k_s, 4.0, places=5)

    def test_band_edge_arrival_times(self):
        from miniBSE.namd.analysis import compute_band_edge_arrival_times
        times = np.linspace(0, 1000, 1001)
        tau = 100.0
        # Ideal exponential decay: excess(t) = 1.0 * exp(-t / 100)
        excess = 1.0 * np.exp(-times / tau)
        res = compute_band_edge_arrival_times(times, excess, tau, temp_k=300.0)
        self.assertAlmostEqual(res["est_95_fs"], 300.0, delta=2.0)
        self.assertAlmostEqual(res["est_99_fs"], 460.5, delta=2.0)
        self.assertAlmostEqual(res["act_95_fs"], 300.0, delta=2.0)
        self.assertAlmostEqual(res["act_99_fs"], 461.0, delta=2.0)
        self.assertAlmostEqual(res["k_cool_ps"], 10.0, places=3)


if __name__ == "__main__":
    unittest.main()
