import os

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


def test_bulk_segments_skip_path_breaks():
    """CsPbBr3 path G-X-M-G-R-X|M-R: the X->M segment is the stretch 27..55, not the X|M jump."""
    from pymatgen.core import Lattice, Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.symmetry.bandstructure import HighSymmKpath
    from qdex.bulk_bands import find_bulk_bs, parse_cp2k_bs, map_bulk_to_path
    st = Structure.from_spacegroup("Pm-3m", Lattice.cubic(5.95), ["Cs", "Pb", "Br"],
                                   [[0.5, 0.5, 0.5], [0, 0, 0], [0.5, 0, 0]])
    prim = SpacegroupAnalyzer(st).get_primitive_standard_structure()
    kf, labels = HighSymmKpath(prim).get_kpoints(line_density=50, coords_are_cartesian=False)
    bs = parse_cp2k_bs(find_bulk_bs("CsPbBr3"))
    mapped = map_bulk_to_path(bs["segments"], np.asarray(kf))
    assert len(mapped) == len(bs["segments"])
    assert all(x[-1] - x[0] > 2 for x, _ in mapped)


def test_soc_bulk_bands_of_cdse():
    """PBE+SOC bulk CdSe (qdex.bulk_soc): Gamma8/Gamma7 split by Delta_so, VBM up by Delta/3, same Cd 4d mean."""
    from qdex.bulk_bands import find_bulk_bs, parse_cp2k_bs, bulk_semicore_level, load_bulk_meta
    sf_path = find_bulk_bs("CdSe")
    assert os.path.basename(sf_path).startswith("CdSe_zb")          # zinc blende preferred
    sf = parse_cp2k_bs(sf_path)
    so = parse_cp2k_bs(find_bulk_bs("CdSe", soc=True))
    meta = load_bulk_meta(sf_path)
    G = so["bands"][0]
    n_occ = 9
    assert np.allclose(G[0::2], G[1::2], atol=1e-6)                    # Kramers pairs at Gamma
    d_so = G[2 * n_occ - 1] - G[2 * n_occ - 6]
    assert 0.33 < d_so < 0.40
    assert meta["gamma_vb_triplet"]["delta_so"] == pytest.approx(d_so, abs=0.002)
    assert so["vbm"] - sf["vbm"] == pytest.approx(d_so / 3, abs=0.005)
    assert so["cbm"] == pytest.approx(sf["cbm"], abs=0.005)
    assert meta["semicore"]["label"] == "Cd-d" and meta["semicore"]["bands_sf"] == [1, 6]
    assert bulk_semicore_level(meta, so, spinor=True) == pytest.approx(bulk_semicore_level(meta, sf), abs=0.01)


def test_bulk_band_files_by_cif_and_material():
    from qdex.bulk_bands import find_bulk_bs, load_bulk_meta
    assert os.path.basename(find_bulk_bs("CdSe", cif="/any/where/CdSe_wz.cif")).startswith("CdSe_wz")
    assert os.path.basename(find_bulk_bs("CdSe", cif="CdSe_wz.cif", soc=True)).startswith("CdSe_wz_soc")
    assert os.path.basename(find_bulk_bs("PbSe")).startswith("PbSe_rs")
    assert os.path.basename(find_bulk_bs("CsPbBr3")).startswith("CsPbBr3_cubic")
    assert find_bulk_bs("Unobtainium") is None
    # every material has spin-free and SOC bands and a semicore anchor
    from qdex.bulk_bands import DEFAULT_BULK_DIR
    for f in os.listdir(DEFAULT_BULK_DIR):
        if f.endswith(".json"):
            name = f[:-5]
            meta = load_bulk_meta(find_bulk_bs(None, cif=name + ".cif"))
            assert meta["name"] == name and meta["semicore"] is not None
            assert find_bulk_bs(None, cif=name + ".cif", soc=True) is not None


@pytest.mark.parametrize("name", ["HgS_zb", "HgSe_zb", "HgTe_zb", "InAs_zb"])
def test_inverted_zincblende_band_order(name):
    """PBE puts Gamma1 (s) below Gamma15 (p) in these semimetals: negative s-p order, Gamma15 partly filled."""
    from qdex.bulk_bands import find_bulk_bs, load_bulk_meta
    meta = load_bulk_meta(find_bulk_bs(None, cif=name + ".cif"))
    assert meta["gamma_band_order"]["inverted"] and meta["gamma_band_order"]["sf"] < -0.2
    assert not meta["gamma_vb_triplet"]["filled"]
    assert meta["gamma_vb_triplet"]["delta_so"] > 0.0


def test_write_bs_round_trip(tmp_path):
    from qdex.bulk_soc import write_bs
    from qdex.bulk_bands import parse_cp2k_bs
    kf = [np.linspace([0, 0, 0], [0.5, 0, 0.5], 5), np.linspace([0.5, 0, 0.5], [0.5, 0.25, 0.75], 3)]
    bands = [np.sort(np.random.default_rng(i).normal(size=(len(k), 6)), axis=1) for i, k in enumerate(kf)]
    write_bs(tmp_path / "x.bs", kf, bands, n_occupied=2)
    d = parse_cp2k_bs(tmp_path / "x.bs")
    assert d["n_k"] == 7 and d["n_bands"] == 6 and d["vbm_band"] == 1
    assert np.allclose(d["segments"][1]["bands"], bands[1]) and np.allclose(d["segments"][0]["kfrac"], kf[0])


def test_material_db_bulk_pbe_gaps():
    """Index 7 is the PBE gap at the experimental lattice; MATERIAL_DB_BULK_PBE keeps both lattices."""
    from qdex.hardness import MATERIAL_DB, MATERIAL_DB_BULK_PBE
    for key, v in MATERIAL_DB_BULK_PBE.items():
        if v["a_exp"] is None:
            continue
        assert MATERIAL_DB[key][7] == pytest.approx(v["gap_exp_lattice"], abs=1e-6)
        assert MATERIAL_DB[key][2] == pytest.approx(v["a_exp"], abs=1e-3)
        assert v["a_pbe"] > 0.99 * v["a_exp"]                       # PBE lattices are not smaller
    assert MATERIAL_DB_BULK_PBE["GAAS"]["gap_exp_lattice"] > MATERIAL_DB_BULK_PBE["GAAS"]["gap_pbe_lattice"] + 0.3
