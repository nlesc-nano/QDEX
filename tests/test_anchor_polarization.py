"""Two-anchor QP curve with sphere polarization, Penn eps(R), plasmon-pole Z."""
import os
import unittest
from unittest.mock import patch

import numpy as np

from qdex.hardness import (
    COULOMB_EV_ANG,
    MATERIAL_DB,
    anchor_edge_curves,
    build_sphere_reaction_field,
    compute_dynamic_z,
    estimate_gw_qp_gap,
    get_cluster_size_metrics,
    penn_eps_eff,
    penn_gap_ev,
    scale_w_difference,
    sphere_polarization_factor,
    valence_plasmon_ev,
)

HERE = os.path.dirname(__file__)


def _metrics(radius):
    return {
        "R_eff_hull": float(radius), "anisotropy": 1.0, "anisotropy_ratio": 1.0,
        "principal_extents_ang": [2.0 * radius] * 3, "selected_atom_indices": [0, 1],
        "radius_definition_version": "test", "surface_offset_ang": 0.0,
    }


class SpherePolarizationTests(unittest.TestCase):
    def test_matched_dielectric_has_no_polarization(self):
        self.assertEqual(sphere_polarization_factor(6.2, 6.2), 0.0)

    def test_multipoles_add_to_born_term(self):
        born = 1.0 - 1.0 / 6.2
        F = sphere_polarization_factor(6.2, 1.0)
        self.assertGreater(F, born)
        self.assertLess(F, 1.5 * born)

    def test_solvent_reduces_polarization(self):
        self.assertLess(sphere_polarization_factor(6.2, 2.24), sphere_polarization_factor(6.2, 1.0))

    def test_reaction_field_center_is_born(self):
        coords = np.array([[0.0, 0.0, 0.0], [2.6, 0.0, 0.0], [-2.6, 0.0, 0.0], [0.0, 2.6, 0.0],
                           [0.0, -2.6, 0.0], [0.0, 0.0, 2.6], [0.0, 0.0, -2.6]])
        syms = ["Cd"] + ["Se"] * 6
        G = build_sphere_reaction_field(coords, syms, "CDSE", eps_out=1.0)
        R = get_cluster_size_metrics(coords, syms, "CDSE")["R_eff_hull"]
        eps = MATERIAL_DB["CDSE"][0]
        self.assertAlmostEqual(G[0, 0], COULOMB_EV_ANG / R * (1.0 - 1.0 / eps), places=10)
        np.testing.assert_allclose(G, G.T, atol=1e-12)
        self.assertTrue(np.all(np.diag(G) > 0.0))


class TwoAnchorCurveTests(unittest.TestCase):
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    syms = ["Cd", "Se"]

    def _run(self, radius, eps_out=1.0, **kw):
        with patch("qdex.hardness.get_cluster_size_metrics", return_value=_metrics(radius)):
            return estimate_gw_qp_gap(self.coords, self.syms, "CDSE", eps_out=eps_out, return_details=True, **kw)

    def test_anchor_gap_and_edges_are_exact(self):
        e = MATERIAL_DB["CDSE"]
        shift, det = self._run(e[9])
        self.assertAlmostEqual(shift, (e[13] - e[12]) - (e[11] - e[10]), places=10)
        d_h, d_l = -(e[12] - e[10]), e[13] - e[11]
        self.assertAlmostEqual(det["f_homo"], d_h / (d_h + d_l), places=10)

    def test_bulk_limit(self):
        shift, det = self._run(1.0e9)
        self.assertAlmostEqual(shift, det["bulk_gw_shift_ev"], places=6)
        edges = anchor_edge_curves("CDSE", 1.0e9, 1.0)
        self.assertAlmostEqual(edges["f_homo"], edges["bulk_homo_fraction"], places=6)

    def test_polarization_is_one_over_r_at_large_size(self):
        eps = MATERIAL_DB["CDSE"][0]
        R = 200.0
        shift, det = self._run(R)
        pol = shift - det["bulk_gw_shift_ev"] - det["anchor_residual_ev"] * (det["monomer_radius_ang"] / R) ** 2
        self.assertAlmostEqual(pol, sphere_polarization_factor(eps, 1.0) * COULOMB_EV_ANG / R, places=10)

    def test_legacy_curve_still_available(self):
        e = MATERIAL_DB["CDSE"]
        shift, det = self._run(e[9], polarization_model="legacy")
        self.assertAlmostEqual(shift, (e[13] - e[12]) - (e[11] - e[10]), places=10)
        self.assertEqual(det["polarization_model"], "legacy")

    def test_anchor_radius_matches_monomer_geometry(self):
        path = os.path.join(HERE, "CdSe", "1.2nm", "geom.xyz")
        lines = open(path).read().split("\n")
        n = int(lines[0])
        rows = [ln.split() for ln in lines[2:2 + n]]
        syms = [r[0] for r in rows]
        xyz = np.array([[float(v) for v in r[1:4]] for r in rows])
        R = get_cluster_size_metrics(xyz, syms, "CDSE")["R_eff_hull"]
        self.assertAlmostEqual(R, MATERIAL_DB["CDSE"][9], places=4)


class PennAndZTests(unittest.TestCase):
    def test_plasmon_and_penn_gap(self):
        self.assertAlmostEqual(valence_plasmon_ev("CDSE"), 14.12, places=1)
        E_P = penn_gap_ev("CDSE")
        self.assertAlmostEqual(E_P, valence_plasmon_ev("CDSE") / np.sqrt(MATERIAL_DB["CDSE"][0] - 1.0), places=12)

    def test_penn_eps(self):
        self.assertAlmostEqual(penn_eps_eff(6.2, 0.0, 6.0), 6.2, places=12)
        self.assertAlmostEqual(penn_eps_eff(6.2, 6.0, 6.0), 1.0 + 5.2 / 4.0, places=12)
        self.assertGreater(penn_eps_eff(6.2, 1.0, 6.0), penn_eps_eff(6.2, 2.0, 6.0))

    def test_z_uses_material_plasmon(self):
        wp = valence_plasmon_ev("CDSE")
        omega = wp / np.sqrt(1.0 - 1.0 / 4.0)
        self.assertAlmostEqual(compute_dynamic_z(1.0, eps_eff=4.0, material_name="CDSE"),
                               1.0 / (1.0 + 1.0 / omega), places=12)

    def test_kernel_scaling(self):
        W_qd, W_b = np.array([[3.0, 1.0], [1.0, 3.0]]), np.array([[2.0, 0.5], [0.5, 2.0]])
        np.testing.assert_allclose(scale_w_difference(W_qd, W_b, 1.0), W_qd)
        np.testing.assert_allclose(scale_w_difference(W_qd, W_b, 0.5), 0.5 * (W_qd + W_b))


if __name__ == "__main__":
    unittest.main()
