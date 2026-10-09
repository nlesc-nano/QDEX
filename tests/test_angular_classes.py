"""Angular coverage Omega (qdex.angular), the band-edge classes and the delocalized-edge cubes."""
import numpy as np

from qdex.angular import harmonics, omega_of
from qdex.fuzzy_bands import classify_edges, deloc_edge_mos


def _sphere(n=2000):
    i = np.arange(n) + 0.5
    pol, az = np.arccos(1 - 2 * i / n), np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(az) * np.sin(pol), np.sin(az) * np.sin(pol), np.cos(pol)], axis=1)


def test_omega_reference_values():
    u = _sphere()
    Y, ls = harmonics(u, 4)
    flat = np.full(len(u), 1.0 / len(u))
    pz = u[:, 2] ** 2 / np.sum(u[:, 2] ** 2)                    # |Y_10|^2: a 1P envelope, no density at the centre
    point = np.zeros(len(u)); point[0] = 1.0
    om = omega_of(np.stack([flat, pz, point], axis=1), Y, ls)
    assert np.allclose(om, [1.0, 5.0 / 9.0, 1.0 / 15.0], atol=1e-3)   # even l only: 1/(1+5+9)


def test_edges_from_classes():
    E = np.array([-1.2, -1.0, -0.8, -0.5, 0.5, 0.7, 0.9])
    occ = E < 0
    labels = ["S-like", "D-like", "Facet", "Localized", "S-like", "Facet", "P-like"]
    r = classify_edges(E, occ, labels)
    assert (r["homo"], r["deloc_homo"], r["lumo"], r["deloc_lumo"]) == (3, 1, 4, 4)
    assert list(r["flag"]) == [0, 2, 4, 1, 3, 0, 0]               # 1 Localized, 4 Facet between the edges
    assert list(r["localized"]) == [False, False, True, True, False, True, False]


def test_deloc_edge_mos_skips_traps():
    E = np.array([-1.2, -1.0, -0.8, -0.5, 0.5, 0.7, 0.9])          # MOs 10..16, HOMO = 13
    mo = np.arange(10, 17)
    localized = np.array([False, True, False, True, False, True, False])
    cls = {"deloc_homo": 2, "deloc_lumo": 4, "localized": localized}
    out = deloc_edge_mos(cls, E, mo, 2, 2)
    assert out == {12: "dHOMO", 10: "dHOMO-1", 14: "dLUMO", 16: "dLUMO+1"}


def test_omega_trap_sets_band_boundary():
    import qdex.angular as qa
    u = _sphere(400) * 10.0
    syms = ["Cd"] * len(u)
    try:
        qa.set_omega_trap(0.40)
        assert abs(qa.AngularCoverage(u, syms).thresholds["D"] - 0.40) < 1e-12
        qa.set_omega_trap("auto")
        a = qa.AngularCoverage(u, syms)
        assert abs(a.thresholds["D"] - 0.5 * (a.references["D"] + a.references["facet"])) < 1e-12
    finally:
        qa.set_omega_trap(0.40)
