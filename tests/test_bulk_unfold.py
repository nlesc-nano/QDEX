"""Unfolded supercell bands (qdex.bulk_unfold) and their overlay on the cubic path."""
import json
import os

import numpy as np
import pytest

from qdex.bulk_bands import DEFAULT_BULK_DIR, find_unfolded, get_aligned_unfolded_bands, parse_cp2k_bs
from qdex.bulk_unfold import pseudocubic_cell


def _orthorhombic_cubic_supercell(a=6.0):
    """Untilted cubic CsPbBr3 written in the sqrt2 x 2 x sqrt2 (Pnma-like) supercell."""
    prim = np.eye(3) * a
    frac = {"Pb": [0, 0, 0], "Br": [[.5, 0, 0], [0, .5, 0], [0, 0, .5]], "Cs": [.5, .5, .5]}
    M = np.array([[1, 1, 0], [0, 0, 2], [1, -1, 0]])
    cell = M @ prim
    syms, X = [], []
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                o = np.array([i, j, k]) @ prim
                for s, fs in frac.items():
                    for f in (fs if isinstance(fs[0], list) else [fs]):
                        r = o + np.array(f) @ prim
                        u = r @ np.linalg.inv(cell)
                        if np.all(u >= -1e-9) and np.all(u < 1 - 1e-9):
                            syms.append(s)
                            X.append(r)
    return cell, syms, np.array(X), M


def test_pseudocubic_cell_of_a_supercell():
    cell, syms, X, M = _orthorhombic_cubic_supercell()
    assert len(syms) == 20
    A_p, M_found = pseudocubic_cell(cell, syms, X, "Pb")
    assert abs(round(np.linalg.det(M_found))) == 4
    assert np.allclose(np.linalg.norm(A_p, axis=1), 6.0)
    assert np.allclose(M_found @ A_p, cell)                    # exact supercell of the primitive cell


@pytest.mark.parametrize("name", ["CsPbBr3_cubic", "CsPbI3_cubic", "CsPbCl3_cubic"])
def test_unfolded_overlay_files(name):
    npz = find_unfolded(cif=name + ".cif")
    assert npz and npz.endswith(f"{name}_unfolded.npz")
    meta = json.load(open(os.path.splitext(npz)[0] + ".json"))
    assert meta["semicore"]["label"] == "Cs-p"
    d = np.load(npz)
    for tag in ("sf", "soc"):
        P = d[f"P_{tag}"]
        assert P.min() >= -1e-9 and P.max() <= 1 + 1e-6
    # on the cubic path of the primitive-cell file, anchored on the Cs 5p level
    path = np.concatenate([s["kfrac"] for s in parse_cp2k_bs(os.path.join(DEFAULT_BULK_DIR, name + ".bs.gz"))["segments"]])
    level = meta["semicore"]["level_soc"]
    res = get_aligned_unfolded_bands(npz, qd_semicore_rel=level - 5.0, path_frac=path, soc=True, qd_semicore_label="Cs-p")
    assert res["unfolded"] and res["alignment_mode"] == "core_level"
    assert res["shift"] == pytest.approx(-5.0)
    assert res["vbm_aligned"] == pytest.approx(meta["soc"]["vbm"] - 5.0)
    x, b, w = res["segments"][0]
    assert b.shape == w.shape and len(x) == len(b)
