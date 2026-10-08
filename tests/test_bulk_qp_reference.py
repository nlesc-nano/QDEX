"""Bulk QP reference: Delta_Sigma from literature QSGW+SOC and PBE+SOC at a_exp, and the geometry correction."""
import math

import numpy as np
import pytest

import qdex.hardness as h


def test_index8_is_index7_plus_delta_sigma():
    for key, v in h.MATERIAL_DB_BULK_PBE.items():
        if "delta_sigma" not in v:
            continue
        entry = h.MATERIAL_DB[key]
        assert entry[8] - entry[7] == pytest.approx(v["delta_sigma"], abs=1.5e-3)
        assert v["qsgw_soc_exp_lattice"] - v["gap_soc_exp_lattice"] == pytest.approx(v["delta_sigma"], abs=1.5e-3)


def test_qsgw_reference_is_the_literature_value_at_a_exp():
    for key, (gap, a_lit, _, _, dlnv) in h.QSGW_SOC_LITERATURE.items():
        v = h.MATERIAL_DB_BULK_PBE[key]
        if dlnv is not None:                         # carried with the source's deformation potential
            assert v["qsgw_soc_exp_lattice"] == pytest.approx(gap + dlnv * 3 * math.log(v["a_exp"] / a_lit), abs=2e-3)
        else:                                        # a_lit within 0.6 % of a_exp, carried with the PBE dEg/dlnV
            assert abs(a_lit / v["a_exp"] - 1) < 6e-3
            assert v["qsgw_soc_exp_lattice"] == pytest.approx(gap, abs=0.05)
    # CdSe: Deguchi et al. QSGW+SO 2.16 eV; the spin-free equivalent is close to their spin-free QSGW (2.28)
    assert h.MATERIAL_DB["CDSE"][8] == pytest.approx(2.28, abs=0.02)


def _zb_cluster(a, n=3, cation="Cd", anion="Se"):
    """Bulk-terminated zinc-blende block with lattice constant a."""
    fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    syms, X = [], []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for b in fcc:
                    for s, off in ((cation, 0.0), (anion, 0.25)):
                        syms.append(s)
                        X.append((np.array([i, j, k]) + b + off) * a)
    return syms, np.array(X)


def test_dot_lattice_strain_fraction():
    v = h.MATERIAL_DB_BULK_PBE["CDSE"]
    for a, s_ref in ((v["a_exp"], 0.0), (v["a_pbe"], 1.0), (0.5 * (v["a_exp"] + v["a_pbe"]), 0.5)):
        s, info = h.dot_lattice_strain("CDSE", *_zb_cluster(a))
        assert s == pytest.approx(s_ref, abs=1e-3)
        assert info["bond"] == "Cd-Se" and info["n_bonds"] >= 8


def test_bulk_shift_is_delta_sigma_plus_delta_geom():
    v = h.MATERIAL_DB_BULK_PBE["CDSE"]
    d_sigma = h.MATERIAL_DB["CDSE"][8] - h.MATERIAL_DB["CDSE"][7]
    try:
        h.set_bulk_vertex("none", verbose=False)
        h.set_bulk_geometry("none", verbose=False)
        assert h.bulk_qp_shift("CDSE")[0] == pytest.approx(d_sigma)
        h.set_bulk_geometry("full", verbose=False)
        shift, info = h.bulk_qp_shift("CDSE")
        assert shift == pytest.approx(d_sigma + v["gap_exp_lattice"] - v["gap_pbe_lattice"])
        assert h.bulk_pbe_gap_dot("CDSE") == pytest.approx(v["gap_pbe_lattice"])
        h.set_bulk_geometry("strain", verbose=False)
        h.set_dot_strain("CDSE", *_zb_cluster(0.5 * (v["a_exp"] + v["a_pbe"])))
        shift, info = h.bulk_qp_shift("CDSE")
        assert info["bulk_strain_fraction"] == pytest.approx(0.5, abs=1e-3)
        assert info["bulk_geometry_shift_ev"] == pytest.approx(0.5 * (v["gap_exp_lattice"] - v["gap_pbe_lattice"]), abs=1e-3)
        # the vertex correction scales Delta_Sigma only
        h.set_bulk_vertex("full", 0.8, verbose=False)
        shift_v, info_v = h.bulk_qp_shift("CDSE")
        assert shift_v == pytest.approx(0.8 * d_sigma + info["bulk_geometry_shift_ev"])
    finally:
        h.set_bulk_vertex("none", verbose=False)
        h.set_bulk_geometry("none", verbose=False)


def test_material_vertex_factor_puts_the_bulk_on_experiment():
    try:
        h.set_bulk_geometry("none", verbose=False)
        h.set_bulk_vertex("full", "material", verbose=False)
        for key in ("CDSE", "GAAS", "INP", "PBS", "HGTE"):
            b = h.MATERIAL_DB_BULK_PBE[key]
            shift, info = h.bulk_qp_shift(key)
            assert b["gap_soc_exp_lattice"] + shift == pytest.approx(h.experimental_bulk_gap(key), abs=1e-9)
        # CsPbBr3: the partner of the measured gap is the PBE+SOC gap of the measured (Pnma) structure
        a, info = h.material_vertex_factor("CSPBBR3")
        assert info["bulk_gap_pbe_soc_ev"] == pytest.approx(h.MATERIAL_DB_BULK_PBE["CSPBBR3"]["gap_soc_exp_structure"])
        assert a == pytest.approx((2.30 - 0.797) / h.MATERIAL_DB_BULK_PBE["CSPBBR3"]["delta_sigma"], abs=1e-3)
        assert h.material_vertex_factor("CDSE")[0] == pytest.approx(0.702, abs=0.005)     # zinc blende, E_exp 1.675 eV
        assert h.material_vertex_factor("CDSE_WZ")[0] == pytest.approx(0.729, abs=0.005)  # wurtzite, E_exp 1.751 eV
    finally:
        h.set_bulk_vertex("none", verbose=False)
    with pytest.raises(ValueError):
        h.set_bulk_vertex("scaled", 1.5, verbose=False)


def _cubic_perovskite(a, n=4):
    syms, X = [], []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                o = np.array([i, j, k]) * a
                syms += ["Pb", "Br", "Br", "Br", "Cs"]
                X += [o, o + [a / 2, 0, 0], o + [0, a / 2, 0], o + [0, 0, a / 2], o + [a / 2, a / 2, a / 2]]
    return syms, np.array(X)


def test_perovskite_strain_from_pb_pb_distances():
    v = h.MATERIAL_DB_BULK_PBE["CSPBBR3"]
    for a, s_ref in ((v["a_exp"], 0.0), (v["a_pbe"], 1.0)):
        s, info = h.dot_lattice_strain("CSPBBR3", *_cubic_perovskite(a))
        assert s == pytest.approx(s_ref, abs=1e-3) and info["bond"] == "Pb-Pb"
    try:                                       # thermal expansion beyond a_PBE is not corrected: s capped at 1
        h.set_bulk_geometry("strain", verbose=False)
        assert h.set_dot_strain("CSPBBR3", *_cubic_perovskite(1.01 * v["a_pbe"]))[0] == 1.0
    finally:
        h.set_bulk_geometry("none", verbose=False)


def test_split_model_vertex_scales_and_residual_does_not():
    b = h.MATERIAL_DB_BULK_PBE["CDSE"]
    ds = b["delta_sigma"]
    try:
        h.set_bulk_geometry("none", verbose=False)
        h.set_bulk_vertex("scaled", 0.8, verbose=False)
        h.set_bulk_residual("experimental", verbose=False)
        bulk, _ = h.bulk_qp_shift("CDSE")                         # no DFT gap: bulk limit (f = 1)
        assert b["gap_soc_exp_lattice"] + bulk == pytest.approx(h.experimental_bulk_gap("CDSE"), abs=1e-9)
        res = h.bulk_residual_shift("CDSE", 0.8)[0]
        assert res == pytest.approx(1.675 - (b["gap_soc_exp_lattice"] + 0.8 * ds), abs=1e-9)
        for gap in (1.2, 2.0, 3.0):                              # the residual is the same at every size
            sh, info = h.bulk_qp_shift("CDSE", gap)
            assert info["bulk_residual_shift_ev"] == pytest.approx(res)
            f = info["bulk_vertex_fraction"]
            assert sh == pytest.approx(ds * (1 - 0.2 * f) + res)
        # with the material factor there is nothing left for the residual
        h.set_bulk_vertex("full", "material", verbose=False)
        assert h.bulk_qp_shift("CDSE")[1]["bulk_residual_shift_ev"] == pytest.approx(0.0, abs=1e-9)
    finally:
        h.set_bulk_vertex("none", verbose=False)
        h.set_bulk_residual("none", verbose=False)
    with pytest.raises(ValueError):
        h.set_bulk_residual("bogus", verbose=False)


def test_perovskite_confinement_is_measured_from_the_tilted_bulk():
    b = h.MATERIAL_DB_BULK_PBE["CSPBBR3"]
    t_exp = b["gap_exp_structure"] - b["gap_exp_lattice"]
    t_pbe = b["gap_pbe_structure"] - b["gap_pbe_lattice"]
    try:
        h.set_bulk_geometry("none", verbose=False)               # dot at a_exp: the tilt opening at a_exp
        assert h.bulk_structure_opening("CSPBBR3")[0] == pytest.approx(t_exp)
        assert h.bulk_pbe_confinement_ref("CSPBBR3") == pytest.approx(b["gap_exp_structure"])
        h.set_bulk_geometry("full", verbose=False)               # dot at the PBE lattice
        assert h.bulk_pbe_confinement_ref("CSPBBR3") == pytest.approx(b["gap_pbe_structure"])
        assert h.bulk_structure_opening("CDSE")[0] == 0.0        # no tilts: cubic reference unchanged
        assert h.bulk_pbe_confinement_ref("CDSE") == pytest.approx(h.bulk_pbe_gap_dot("CDSE"))
        # a dot exactly at the tilted bulk gap has no confinement: full bulk screening, f = 1
        h.set_bulk_vertex("scaled", 0.8, verbose=False)
        info = h.bulk_qp_shift("CSPBBR3", b["gap_pbe_structure"])[1]
        assert info["bulk_vertex_fraction"] == pytest.approx(1.0)
        assert info["bulk_pbe_confinement_ref_ev"] == pytest.approx(b["gap_pbe_structure"])
        # ... and the self-energy shift keeps the cubic reference
        assert info["bulk_pbe_gap_dot_lattice_ev"] == pytest.approx(b["gap_pbe_lattice"])
        assert t_pbe > 0
    finally:
        h.set_bulk_vertex("none", verbose=False)
        h.set_bulk_geometry("none", verbose=False)


def test_bulk_edge_split_modes():
    import qdex.hardness as hh
    try:
        hh.set_bulk_edge_split("symmetric", verbose=False)
        assert hh.bulk_homo_fraction("CdSe")[0] == 0.5
        hh.set_bulk_edge_split("bulk", verbose=False)
        f, src = hh.bulk_homo_fraction("CdSe")
        assert f == 0.95 and "Grueneis" in src
        assert hh.bulk_homo_fraction("PBS")[0] == 0.5            # no bulk value: default 1/2
        hh.set_bulk_edge_split("cluster", verbose=False)
        assert hh.bulk_homo_fraction("CdSe")[0] == 0.07
        assert hh.bulk_homo_fraction("CSPBBR3")[0] == 0.60
        hh.set_bulk_edge_split(0.3, verbose=False)
        assert hh.bulk_homo_fraction("GaAs")[0] == 0.3
        with pytest.raises(ValueError):
            hh.set_bulk_edge_split("vacuum", verbose=False)
        for table in (hh.BULK_EDGE_SPLIT, hh.CLUSTER_EDGE_SPLIT):
            assert all(m in hh.MATERIAL_DB for m in table)
    finally:
        hh.set_bulk_edge_split("symmetric", verbose=False)
