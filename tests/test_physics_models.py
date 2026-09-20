import unittest
from unittest.mock import patch

import numpy as np

from miniBSE.hardness import (
    MATERIAL_DB,
    build_resta_mnok,
    estimate_gw_qp_gap,
    get_cluster_size_metrics,
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


if __name__ == "__main__":
    unittest.main()
