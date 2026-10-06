from math import gamma, pi, sqrt

import numpy as np
import pytest

libint_cpp = pytest.importorskip("libint_cpp")

from qdex.io_utils import parse_gth_soc_potentials
from qdex.soc_utils import _build_soc_projectors, compute_spinor_subspace, get_angular_momentum_matrices

# Se GTH-PBE-q6 and a dummy GTH-BLYP-q6 entry with different SOC constants
GTH_TEXT = """\
Se GTH-BLYP-q6 GTH-BLYP
    2    4
     0.51000000    0
    3  SOC
     0.43367514    3     6.43369922    -0.22063995    -1.17545323
                                       -1.62702987     3.03500718
                                                      -2.40896227
     0.47248310    2     2.23970327     0.40191857
                                       -0.47555647
                         9.99000000     0.00000000
                                        0.00000000
     0.60911821    1     0.49494679
                         0.00563899
#
Se GTH-PBE-q6 GTH-PBE
    2    4
     0.51000000    0
    3  SOC
     0.43246005    3     6.51810990    -0.22271639    -1.19612899
                                       -1.65797833     3.08839178
                                                      -2.45133498
     0.47049162    2     2.28126223     0.36533529
                                       -0.43227055
                         0.15386757    -0.09508784
                                        0.11250945
     0.62560034    1     0.43979948
                         0.00489231
#
"""


@pytest.fixture
def gth_file(tmp_path):
    path = tmp_path / "GTH_SOC_POTENTIALS"
    path.write_text(GTH_TEXT)
    return str(path)


def _shell(l, exp, center=(0.0, 0.0, 0.0)):
    return dict(l=l, exps=np.array([exp]), coefs=np.array([1.0]), center=np.asarray(center, float), pure=True)


def _radial(shell, r):
    """Radial part R(r) of a libint pure shell, from its m = 0 component on the z axis."""
    l = shell["l"]
    A = libint_cpp.evaluate_basis_on_grid([shell], np.c_[0 * r, 0 * r, r], 1)
    return A[:, l] / sqrt((2 * l + 1) / (4 * pi))


def _hgh_radial(l, i, r_l, r):
    k = l + (4 * i - 1) / 2
    return sqrt(2) * r ** (l + 2 * (i - 1)) * np.exp(-r ** 2 / (2 * r_l ** 2)) / (r_l ** k * sqrt(gamma(k)))


def test_gth_parser_selects_the_functional(gth_file):
    pbe = parse_gth_soc_potentials(gth_file, {"Se": 6})["Se"]
    blyp = parse_gth_soc_potentials(gth_file, {"Se": 6}, functional="BLYP")["Se"]
    assert pbe["name"] == "GTH-PBE-q6" and blyp["name"] == "GTH-BLYP-q6"
    assert pbe["so"][1]["k_coeffs"] == pytest.approx([0.15386757, -0.09508784, 0.11250945])
    assert blyp["so"][1]["k_coeffs"][0] == pytest.approx(9.99)


@pytest.mark.parametrize("l", [1, 2, 3])
def test_angular_momentum_matrices_act_on_libint_real_harmonics(l):
    """<R_m|L_k|R_m'> from -i r x grad applied to libint's own pure functions."""
    sh = [_shell(l, 0.7)]
    rng = np.random.default_rng(1)
    P = rng.normal(size=(300, 3))
    h = 1e-5
    f = lambda Q: libint_cpp.evaluate_basis_on_grid(sh, Q, 1)
    F = f(P)
    d = [(f(P + h * e) - f(P - h * e)) / (2 * h) for e in np.eye(3)]
    x, y, z = P.T[:, :, None]
    act = [-1j * (y * d[2] - z * d[1]), -1j * (z * d[0] - x * d[2]), -1j * (x * d[1] - y * d[0])]
    for L, a in zip(get_angular_momentum_matrices(l), act):
        assert np.allclose(np.linalg.lstsq(F, a, rcond=None)[0], L, atol=1e-6)


def test_projector_overlaps_match_radial_integrals(gth_file):
    tbl = parse_gth_soc_potentials(gth_file, {"Se": 6})
    shells = [_shell(1, 0.9), _shell(1, 0.25), _shell(2, 0.6)]
    projectors, _ = _build_soc_projectors(["Se"], [[0.0, 0.0, 0.0]], tbl)
    B = libint_cpp.compute_hgh_overlaps(shells, projectors, 1)
    r = np.linspace(1e-6, 16.0, 160001)
    ao0 = {0: 0, 1: 3, 2: 6}                       # first AO of each shell
    col = 0
    for p in projectors:
        l = p["l"]
        for s_idx, sh in enumerate(shells):
            if sh["l"] != l:
                continue
            ref = np.trapezoid(_radial(sh, r) * _hgh_radial(l, p["i"], p["r_l"], r) * r ** 2, r)
            block = B[ao0[s_idx]:ao0[s_idx] + 2 * l + 1, col:col + 2 * l + 1]
            assert np.allclose(block, ref * np.eye(2 * l + 1), rtol=1e-6, atol=1e-9)
        col += 2 * l + 1


def test_single_atom_p_shell_splits_into_j_multiplets(gth_file):
    """V = sum_ij |p_i> k_ij L.S <p_j|: a p shell gives j = 1/2 at -lam and j = 3/2 at +lam/2."""
    shells = [_shell(1, 0.3)]
    S = libint_cpp.overlap(shells, 1)
    E, U, _ = compute_spinor_subspace(["Se"], [[0.0, 0.0, 0.0]], shells, np.eye(3), np.zeros(3), S,
                                      np.arange(3), gth_file, verbose=False)
    so = parse_gth_soc_potentials(gth_file, {"Se": 6})["Se"]["so"][1]
    r = np.linspace(1e-6, 16.0, 160001)
    b = [np.trapezoid(_radial(shells[0], r) * _hgh_radial(1, i, so["r"], r) * r ** 2, r) for i in (1, 2)]
    k11, k12, k22 = so["k_coeffs"]
    lam = k11 * b[0] ** 2 + 2 * k12 * b[0] * b[1] + k22 * b[1] ** 2
    assert np.allclose(E, [-lam, -lam, lam / 2, lam / 2, lam / 2, lam / 2], atol=1e-8)


def test_spinor_levels_come_in_kramers_pairs(gth_file):
    shells = [_shell(0, 0.5), _shell(1, 0.4), _shell(2, 0.5),
              _shell(1, 0.3, (2.1, -1.3, 3.0)), _shell(2, 0.7, (2.1, -1.3, 3.0))]
    for sh, a in zip(shells, [0, 0, 0, 1, 1]):
        sh["atom_idx"] = a
    coords_ang = np.array([[0.0, 0.0, 0.0], [2.1, -1.3, 3.0]]) / 1.8897259886
    S = libint_cpp.overlap(shells, 1)
    w, V = np.linalg.eigh(S)
    C = V @ np.diag(w ** -0.5) @ V.T                 # Löwdin-orthonormal AOs as MOs
    eps = np.random.default_rng(3).uniform(-0.3, 0.3, C.shape[1])
    E, _, _ = compute_spinor_subspace(["Se", "Se"], coords_ang, shells, C, eps, S, np.arange(C.shape[1]),
                                      gth_file, verbose=False)
    assert np.allclose(E[0::2], E[1::2], atol=1e-10)
    assert np.abs(E - np.sort(np.repeat(eps, 2))).max() > 1e-3   # SOC is switched on


def test_spinor_fuzzy_weights_add_the_two_spin_components():
    from qdex.fuzzy_bands import folded_plane_wave_weights
    shells = [_shell(0, 0.4), _shell(1, 0.6, (1.0, 0.5, -0.7))]
    k = np.random.default_rng(0).normal(size=(20, 3))
    C = np.eye(4)
    U = np.zeros((8, 1), complex)                   # (phi_0 alpha + i phi_2 beta) / sqrt2
    U[0, 0] = 1 / np.sqrt(2)
    U[4 + 2, 0] = 1j / np.sqrt(2)
    G = np.zeros((1, 3))
    W, W_spin = folded_plane_wave_weights(shells, k, G, 1, [C], [(0, range(4), U[:4]), (0, range(4), U[4:])])
    assert np.allclose(W_spin[0], 0.5 * (W[0][0] + W[0][2]))
