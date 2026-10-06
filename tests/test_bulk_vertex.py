"""Vertex correction of the bulk QSGW shift (quasiparticles.bulk_vertex)."""
import unittest

import qdex.hardness as h

PBE, QSGW = h.MATERIAL_DB["CDSE"][7], h.MATERIAL_DB["CDSE"][8]     # bulk spin-free PBE and QSGW gaps
D = QSGW - PBE


class TestBulkVertex(unittest.TestCase):
    def tearDown(self):
        h.set_bulk_vertex("none", 0.8, verbose=False)

    def test_none_is_pure_qsgw(self):
        h.set_bulk_vertex("none", verbose=False)
        shift, info = h.bulk_qp_shift("CDSE", 1.457)
        self.assertAlmostEqual(shift, D)
        self.assertEqual(info["bulk_vertex_fraction"], 0.0)

    def test_full_applies_the_bulk_factor(self):
        h.set_bulk_vertex("full", 0.8, verbose=False)
        shift, info = h.bulk_qp_shift("CDSE", 2.639)
        self.assertAlmostEqual(shift, 0.8 * D)
        self.assertAlmostEqual(info["bulk_vertex_correction_ev"], -0.2 * D)

    def test_scaled_follows_the_penn_screening_fraction(self):
        h.set_bulk_vertex("scaled", 0.8, verbose=False)
        d = D
        bulk, _ = h.bulk_qp_shift("CDSE", PBE)              # no confinement: full bulk correction
        self.assertAlmostEqual(bulk, 0.8 * d)
        shifts = [h.bulk_qp_shift("CDSE", g)[0] for g in (PBE, 1.0, 1.457, 2.639, 6.0)]
        self.assertTrue(all(a <= b for a, b in zip(shifts, shifts[1:])))   # smaller dots keep more of the QSGW shift
        s12, info = h.bulk_qp_shift("CDSE", 2.639)
        f = (info["bulk_vertex_eps_eff"] - 1.0) / (6.2 - 1.0)
        self.assertAlmostEqual(info["bulk_vertex_fraction"], f)
        self.assertAlmostEqual(s12, d - 0.2 * d * f)

    def test_invalid_settings(self):
        with self.assertRaises(ValueError):
            h.set_bulk_vertex("partial", verbose=False)
        with self.assertRaises(ValueError):
            h.set_bulk_vertex("full", 1.5, verbose=False)


if __name__ == "__main__":
    unittest.main()
