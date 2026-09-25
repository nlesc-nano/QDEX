import types
import unittest

import numpy as np

from qdex.cli import resolve_bse_kernel, resolve_qp_z


def _args(**kw):
    base = dict(qp_gap="sgw-resta", kernel=None, kernel_type="mnok", allow_inconsistent_kernel=False,
                qp_z=None, dynamic_z=False)
    base.update(kw)
    return types.SimpleNamespace(**base)


W = (np.eye(2), "test W")


class KernelConsistencyTests(unittest.TestCase):
    def test_w_based_model_defaults_to_shared_w(self):
        self.assertEqual(resolve_bse_kernel(_args(), W), "qp")

    def test_w_based_model_rejects_a_different_kernel(self):
        for k in ["resta", "dim", "bse", "xs-resta", "sbse"]:
            with self.assertRaises(ValueError):
                resolve_bse_kernel(_args(kernel=k), W)

    def test_legacy_escape_hatch(self):
        self.assertEqual(resolve_bse_kernel(_args(kernel="resta", allow_inconsistent_kernel=True), W), "resta")

    def test_sbse_kernel_is_the_same_w_for_sgw(self):
        self.assertEqual(resolve_bse_kernel(_args(qp_gap="sgw", kernel="sbse"), W), "qp")

    def test_shared_w_requires_atom_representation(self):
        with self.assertRaises(ValueError):
            resolve_bse_kernel(_args(kernel_type="xs"), W)

    def test_gap_only_models(self):
        self.assertEqual(resolve_bse_kernel(_args(qp_gap="gw", kernel="resta"), None), "resta")
        self.assertEqual(resolve_bse_kernel(_args(qp_gap="pbe"), None), "bse")
        with self.assertRaises(ValueError):
            resolve_bse_kernel(_args(qp_gap="gw", kernel="qp"), None)


class ZOptionTests(unittest.TestCase):
    def test_model_default(self):
        self.assertEqual(resolve_qp_z(_args(), False), (False, 0.8))
        self.assertEqual(resolve_qp_z(_args(), True), (True, 0.8))

    def test_fixed_and_derived(self):
        self.assertEqual(resolve_qp_z(_args(qp_z="1.0"), True), (False, 1.0))
        self.assertEqual(resolve_qp_z(_args(qp_z="derived"), False)[0], True)
        self.assertEqual(resolve_qp_z(_args(dynamic_z=True), False)[0], True)

    def test_invalid_values(self):
        for bad in ["0", "1.2", "-0.5", "abc"]:
            with self.assertRaises(ValueError):
                resolve_qp_z(_args(qp_z=bad), False)


if __name__ == "__main__":
    unittest.main()
