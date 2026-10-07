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


def test_perovskite_scale_from_the_b_sublattice():
    """Tilted octahedra: the k scale follows Pb-Pb (the lattice), not the Pb-Br bond."""
    from pymatgen.core import Lattice, Structure
    from qdex.fuzzy_bands import fit_lattice_orientation

    a = 5.95
    prim = Structure(Lattice.cubic(a), ["Cs", "Pb", "Br", "Br", "Br"],
                     [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    sc = prim * (5, 5, 5)
    coords = np.array([s.coords for s in sc]) * 0.99
    syms = [s.specie.symbol for s in sc]
    # buckle the Br off the Pb-Pb lines (alternating, as a tilt): bonds lengthen, Pb-Pb is unchanged
    for k, (s_, x) in enumerate(zip(syms, coords)):
        if s_ == "Br":
            frac = np.round(x / (0.99 * a) * 2).astype(int)
            axis = int(np.argmax(frac % 2 == 1))   # perpendicular to the Pb-Br-Pb line
            coords[k, axis] += 0.35 * (-1) ** int(frac.sum() // 2)
    info = fit_lattice_orientation(prim, coords, syms)

    assert info["sublattice_element"] == "Pb"
    assert info["scale"] == pytest.approx(0.99, abs=1e-3)
    assert info["qd_bond_ang"] > 0.99 * a / 2 * 1.005


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
        if f.endswith(".json") and not f.endswith("_unfolded.json"):   # unfolded overlays: test_bulk_unfold
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
    # Delta_so = E(Gamma8) - E(Gamma7): negative in HgS (doublet above the quartet; LDA -0.11 eV, Svane 2011)
    assert (meta["gamma_vb_triplet"]["delta_so"] < 0.0) == (name == "HgS_zb")


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
        if not key.endswith("_WZ"):                  # wurtzite: index 2 is the zinc-blende-equivalent a
            assert MATERIAL_DB[key][2] == pytest.approx(v["a_exp"], abs=1e-3)
        assert v["a_pbe"] > 0.99 * v["a_exp"]                       # PBE lattices are not smaller
    assert MATERIAL_DB_BULK_PBE["GAAS"]["gap_exp_lattice"] > MATERIAL_DB_BULK_PBE["GAAS"]["gap_pbe_lattice"] + 0.3


def test_state_weights_and_k_peaks():
    from qdex.fuzzy_bands import state_weights, fuzzy_state_peaks
    k = np.arange(60)
    I = np.vstack([5.0 * np.exp(-0.5 * ((k - 20) / 3.0) ** 2) + 0.01,                  # one peak at k = 20
                   np.exp(-0.5 * ((k - 10) / 2.0) ** 2) + np.exp(-0.5 * ((k - 45) / 2.0) ** 2),
                   np.ones(60)])                                                      # flat: no peak
    P = state_weights(I)
    assert np.allclose(P.mean(axis=1), 1.0)
    kk, ee, ww, nn = fuzzy_state_peaks(P, np.array([-1.0, 0.5, 2.0]), (-5, 5))
    assert sorted(zip(nn.tolist(), kk.tolist())) == [(0, 20), (1, 10), (1, 45)]
    assert np.all(ww > 1.5) and np.allclose(ee[nn == 1], 0.5)
    # a jump in the path (K|U) splits the smoothing: a peak at the last point of a segment stays there
    frac = np.zeros((60, 3))
    frac[:, 0] = np.r_[np.linspace(0, 0.3, 30), np.linspace(0.6, 0.9, 30)]
    edge = np.where(k < 30, np.exp(-0.5 * ((k - 29) / 3.0) ** 2), np.exp(-0.5 * ((k - 30) / 3.0) ** 2))[None, :] + 0.01
    kk, *_ = fuzzy_state_peaks(state_weights(edge), np.array([0.0]), (-5, 5), kpts_frac=frac)
    assert sorted(kk.tolist()) == [29, 30]                       # one peak per side of the jump
    kk, *_ = fuzzy_state_peaks(state_weights(edge), np.array([0.0]), (-5, 5))
    assert len(kk) == 1                                          # merged when the path is taken as continuous


def test_spinor_soc_energy_is_the_soc_expectation_value():
    from qdex.fuzzy_bands import spinor_soc_energy
    rng = np.random.default_rng(3)
    n = 12
    eps = np.sort(rng.normal(size=n))
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    V = 0.1 * (A + A.conj().T)
    V -= np.trace(V).real / n * np.eye(n)          # V_SOC is traceless
    E, U = np.linalg.eigh(np.diag(eps) + V)
    v_ref = np.einsum("in,ij,jn->n", U.conj(), V, U).real
    assert np.allclose(spinor_soc_energy(E, U, eps), v_ref)
    assert np.allclose(spinor_soc_energy(E - 0.7, U, eps), v_ref)   # frames shifted by a constant


def test_fuzzy_export_and_state_norm_display(tmp_path, monkeypatch):
    from qdex.fuzzy_bands import smear_and_export_fuzzy
    from qdex.plot_fuzzy import load_fuzzy, prepare_fuzzy_display
    monkeypatch.chdir(tmp_path)
    k = np.arange(40)
    E = np.array([-1.0, -0.5, 0.5, 1.0])
    I = np.vstack([np.exp(-0.5 * ((k - c) / 3.0) ** 2) * s + 1e-3 for c, s in zip((0, 10, 20, 30), (1.0, 50.0, 1.0, 2.0))])
    labels = ["Γ"] + [""] * 38 + ["X"]
    v = np.array([np.nan, 0.05, -0.1, np.nan])
    smear_and_export_fuzzy(I, E, labels, (-2.0, 2.0), 0.02, prefix="soc", kpts_frac=np.c_[k / 40.0, k * 0, k * 0],
                           soc_energy=v)
    fz = load_fuzzy("fuzzy_data_soc.npz")
    assert fz["Z_norm"] is not None and fz["peaks"]["k"].size >= 4
    assert np.isnan(fz["peaks"]["soc_energy"][fz["peaks"]["state"] == 0]).all()
    assert np.isclose(np.nanmax(fz["soc_energy_map"]), 0.05, atol=1e-3)
    Zd, zmin, zmax, *_ = prepare_fuzzy_display(fz["Z"], "state_norm", Z_norm=fz["Z_norm"])
    assert zmin == 0.0 and zmax == 1.0 and 0.0 <= Zd.min() and Zd.max() <= 1.0
    # per state the 50x heavier state no longer dominates: both peaks reach comparable brightness
    rows = [np.argmin(abs(fz["centres"] - e)) for e in (-0.5, 0.5)]
    assert Zd[rows[0]].max() / Zd[rows[1]].max() < 1.5
