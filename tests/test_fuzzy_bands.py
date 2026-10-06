import numpy as np
import pytest

libint_cpp = pytest.importorskip("libint_cpp")


def _shells(l_values, centers):
    return [dict(l=l, exps=np.array([0.8]), coefs=np.array([1.0]), center=np.asarray(c, float), pure=True)
            for l, c in zip(l_values, centers)]


def test_ao_fourier_transform_reproduces_libint_overlap():
    """Parseval: S_mn = int F_m(k)* F_n(k) d3k / (2 pi)^3 for s, p, d and f solid harmonics."""
    shells = _shells([0, 1, 2, 3, 2, 3], [[0, 0, 0]] * 4 + [[0.7, -0.4, 0.9]] * 2)
    S = libint_cpp.overlap(shells, 1)
    h = 0.25
    g = np.arange(-10.0, 10.0 + h / 2, h)
    K = np.stack(np.meshgrid(g, g, g, indexing="ij"), axis=-1).reshape(-1, 3)
    F = libint_cpp.ao_ft_complex(shells, K, 4)
    S_k = np.real(F.conj() @ F.T) * h ** 3 / (2 * np.pi) ** 3
    assert np.allclose(S_k, S, atol=1e-8)


def test_ao_grid_values_reproduce_libint_overlap():
    shells = _shells([2, 3, 2, 3], [[0, 0, 0], [0, 0, 0], [0.7, -0.4, 0.9], [0.7, -0.4, 0.9]])
    S = libint_cpp.overlap(shells, 1)
    h = 0.12
    g = np.arange(-7.0, 8.0 + h / 2, h)
    P = np.stack(np.meshgrid(g, g, g, indexing="ij"), axis=-1).reshape(-1, 3)
    A = libint_cpp.evaluate_basis_on_grid(shells, P, 4)
    assert np.allclose(A.T @ A * h ** 3, S, atol=1e-5)


def _zincblende():
    from pymatgen.core import Lattice, Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    conv = Structure.from_spacegroup("F-43m", Lattice.cubic(6.08), ["Cd", "Se"], [[0, 0, 0], [0.25, 0.25, 0.25]])
    return conv, SpacegroupAnalyzer(conv).get_primitive_standard_structure()


def test_lattice_orientation_recovers_a_rotated_cluster():
    from scipy.spatial.transform import Rotation
    from qdex.fuzzy_bands import fit_lattice_orientation

    conv, prim = _zincblende()
    sites = conv.get_sites_in_sphere(conv[0].coords, 9.0)
    coords = np.array([s.coords for s in sites])
    syms = [s.specie.symbol for s in sites]
    Q = Rotation.from_euler("zyx", [37.0, -21.0, 64.0], degrees=True).as_matrix()
    info = fit_lattice_orientation(prim, 1.02 * coords @ Q.T, syms)

    assert info["match"] > 0.999
    assert info["scale"] == pytest.approx(1.02, abs=1e-6)
    # the fitted frame differs from Q at most by a cubic operation (signed permutation)
    M = Q.T @ info["rotation"]
    assert np.allclose(np.sort(np.abs(M), axis=1), [[0, 0, 1]] * 3, atol=1e-6)


def test_bulk_bands_land_on_the_fuzzy_path():
    from pymatgen.symmetry.bandstructure import HighSymmKpath
    from qdex.bulk_bands import find_bulk_bs, parse_cp2k_bs, map_bulk_to_path

    _, prim = _zincblende()
    kfrac, labels = HighSymmKpath(prim).get_kpoints(line_density=50, coords_are_cartesian=False)
    kfrac = np.asarray(kfrac)
    bs = parse_cp2k_bs(find_bulk_bs("CdSe"))
    mapped = map_bulk_to_path(bs["segments"], kfrac)

    assert len(mapped) == len(bs["segments"])
    x0, _ = mapped[0]                     # Gamma -> X starts at Gamma and ends at the first X
    assert x0[0] == pytest.approx(0.0)
    assert labels[int(round(x0[-1]))] == "X"
    assert bs["n_k"] == sum(len(s["kfrac"]) for s in bs["segments"]) - (len(bs["segments"]) - 2)
