import numpy as np

from qdex.fuzzy_bands import kspace_descriptors


def test_flat_and_gamma_profiles():
    k = np.zeros((10, 3))
    k[:, 0] = np.linspace(0.0, 0.9, 10)          # first point at Gamma, spacing 0.1 1/A
    flat = np.ones((1, 10))
    peak = np.zeros((1, 10)); peak[0, 0] = 2.0
    g, s, uniform = kspace_descriptors(np.vstack([flat, peak]), k, gamma_radius=0.12)
    assert np.isclose(uniform, 0.2)              # 0.0 and 0.1 lie within 0.12 1/A
    assert np.isclose(g[0], 0.2) and np.isclose(s[0], 1.0)
    assert np.isclose(g[1], 1.0) and np.isclose(s[1], 0.0)


def test_k_participation_and_band_edges():
    from qdex.fuzzy_bands import classify_band_edges, k_participation
    flat, spot = np.ones((1, 8)), np.eye(1, 8)
    assert np.allclose(k_participation(np.vstack([flat, spot])), [1.0, 1.0 / 8])
    # occupied: two localized states above a band state; empty: an off-band state in the bulk gap, then a band state
    E = np.array([-0.5, -0.6, -1.0, 0.1, 0.6])
    occ = np.array([True, True, True, False, False])
    kpart = np.array([0.9, 0.8, 0.1, 0.1, 0.05])
    onband = np.array([0.1, 0.1, 0.7, 0.05, 0.3])
    bulk = {"vbm": -0.9, "cbm": 0.4}
    r = classify_band_edges(E, occ, kpart, onband, bulk)
    assert r["homo"] == 0 and r["deloc_homo"] == 2 and r["lumo"] == 3 and r["deloc_lumo"] == 4
    assert list(r["flag"]) == [1, 1, 2, 1, 3]


def test_deloc_edge_mos_skips_traps():
    import numpy as np
    from qdex.fuzzy_bands import deloc_edge_mos
    E = np.array([-1.2, -1.0, -0.8, -0.5, 0.5, 0.7, 0.9])          # MOs 10..16, HOMO = 13
    mo = np.arange(10, 17)
    localized = np.array([False, True, False, True, False, True, False])
    cls = {"deloc_homo": 2, "deloc_lumo": 4, "localized": localized}
    out = deloc_edge_mos(cls, E, mo, 2, 2)
    assert out == {12: "dHOMO", 10: "dHOMO-1", 14: "dLUMO", 16: "dLUMO+1"}
