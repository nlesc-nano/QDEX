import os
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from miniBSE.cli import transform_ao_operator
from miniBSE.exciton_analysis import coherent_mulliken_ao_populations
from miniBSE.io_utils import parse_basis, read_mos_auto
from miniBSE.oscillator import compute_oscillator_strengths
from miniBSE.soc_utils import _ortho_active_coeffs


class NumericalSafetyTests(unittest.TestCase):
    def test_mps_dipole_transform_uses_float64_reference_path(self):
        mu = np.array([[2.0, 0.5], [0.5, 1.0]])
        left = np.array([[1.0], [0.2]])
        right = np.array([[0.3], [0.8]])
        result = transform_ao_operator(mu, left, right, device="mps")
        self.assertEqual(result.dtype, np.float64)
        np.testing.assert_allclose(result, left.T @ mu @ right)

    def test_complex_cholesky_orthogonalization(self):
        c = np.array([[1.0, 0.2j], [0.3 + 0.1j, 1.1], [0.2, -0.4j]])
        s = np.array([[1.4, 0.1, 0.0], [0.1, 1.2, 0.05], [0.0, 0.05, 1.1]])
        c_ortho = _ortho_active_coeffs(c, s)
        np.testing.assert_allclose(c_ortho.conj().T @ s @ c_ortho, np.eye(2), atol=1.0e-12)

    def test_explicit_uks_uses_two_thirds_prefactor(self):
        energy = np.array([1.0])
        vectors = np.array([[1.0]])
        dipole = np.array([[1.0, 0.0, 0.0]])
        spatial = compute_oscillator_strengths(energy, vectors, dipole)
        uks = compute_oscillator_strengths(energy, vectors, dipole, spin_resolved=True)
        np.testing.assert_allclose(uks, 0.5 * spatial)

    def test_active_space_mulliken_contraction_matches_dense_reference(self):
        rng = np.random.default_rng(7)
        coeff = rng.normal(size=(7, 3)) + 1j * rng.normal(size=(7, 3))
        raw_s = rng.normal(size=(7, 7)) + 1j * rng.normal(size=(7, 7))
        overlap = raw_s.conj().T @ raw_s
        raw_d = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        density = raw_d @ raw_d.conj().T
        overlap_coeff = overlap @ coeff
        compact = coherent_mulliken_ao_populations(coeff, overlap_coeff, density)
        p_ao = coeff @ density @ coeff.conj().T
        dense = np.real(np.sum(p_ao * overlap.T, axis=1))
        np.testing.assert_allclose(compact, dense, atol=1.0e-11)

    def test_validated_mo_cache_avoids_reparsing_text(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "mos.txt")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("placeholder")
            expected = (
                np.arange(6.0).reshape(3, 2),
                np.array([-0.2, 0.3]),
                np.array([2.0, 0.0]),
            )
            with patch("miniBSE.io_utils.read_mos_txt_cc", return_value=expected) as parser:
                first = read_mos_auto(path, 3, cache=True)
                self.assertEqual(parser.call_count, 1)
            with patch("miniBSE.io_utils.read_mos_txt_cc", side_effect=AssertionError("reparsed")):
                second = read_mos_auto(path, 3, cache=True)
            for actual, reference in zip(first, second):
                np.testing.assert_array_equal(actual, reference)

    def test_ambiguous_basis_prefix_is_rejected(self):
        content = """Cd DZVP-q2\n1\n1 0 0 1 1\n1.0 1.0\nCd DZVP-q4\n1\n1 0 0 1 1\n1.0 1.0\n"""
        with tempfile.NamedTemporaryFile("w", delete=False) as handle:
            handle.write(content)
            path = handle.name
        try:
            with self.assertRaisesRegex(ValueError, "Ambiguous basis selection"):
                parse_basis(path, "DZVP", required_elements={"Cd"})
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
