import unittest

from qdex.config_schema import flatten_config
from qdex.hardness import GXTB_BULK, MATERIAL_DB, gxtb_bulk_shift


class GxtbBulkShiftTests(unittest.TestCase):
    def test_formula(self):
        # Delta = (gap_exp + Delta_so/3) - (gap_pbe_bulk + <E_gxtb - E_pbe>_dots)
        d, e = GXTB_BULK["CDSE"], MATERIAL_DB["CDSE"]
        expected = (e[3] + d["delta_so"] / 3.0) - (e[7] + d["delta"])
        shift, info = gxtb_bulk_shift("CdSe")
        self.assertAlmostEqual(shift, expected, places=12)
        self.assertLess(shift, 0.0)                       # g-xTB overestimates the gap
        self.assertAlmostEqual(info["gap_ref_spin_free_ev"], 1.74 + 0.14, places=12)

    def test_periodic_value_takes_precedence(self):
        GXTB_BULK["CDSE"]["gap_gxtb_bulk"] = 4.0
        try:
            shift, info = gxtb_bulk_shift("CDSE")
            self.assertAlmostEqual(shift, 1.88 - 4.0, places=12)
            self.assertEqual(info["gxtb_bulk_source"], "periodic g-xTB")
        finally:
            del GXTB_BULK["CDSE"]["gap_gxtb_bulk"]

    def test_unknown_material(self):
        with self.assertRaises(ValueError):
            gxtb_bulk_shift("CSPBBR3")

    def test_yaml_key(self):
        flat = flatten_config({"quasiparticles": {"model": "bulk", "reference": "gxtb"}})
        self.assertEqual(flat["qp_reference"], "gxtb")


if __name__ == "__main__":
    unittest.main()
