"""excitations.active_space: delocalized -- a BSE active space of selected, non-contiguous MOs.

The solver with explicit occ / virt MO indices must give exactly what the contiguous window gives when the same MOs
are moved next to the gap; the selection takes the first S/P/D-like MOs counted outward from the gap.
"""
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from qdex.solver import ExcitonSolver


def _system(seed=5, n_occ_mo=7, n_virt_mo=7):
    """4 atoms x 4 AOs, S-orthonormal random MOs, sorted energies."""
    rng = np.random.default_rng(seed)
    syms = ["Cd", "Se", "Cd", "Se"]
    coords = np.array([[0.0, 0.0, 0.0], [2.6, 0.0, 0.0], [0.0, 2.6, 0.0], [2.6, 2.6, 0.4]])
    n_ao = 16
    ranges = [(4 * a, 4 * a + 4) for a in range(4)]
    A = rng.normal(size=(n_ao, n_ao))
    S = np.eye(n_ao) + 0.04 * (A + A.T)
    w, V = np.linalg.eigh(S)
    Q, _ = np.linalg.qr(rng.normal(size=(n_ao, n_ao)))
    C = (V @ np.diag(w ** -0.5) @ V.T) @ Q
    eps = np.sort(np.concatenate([rng.uniform(-9, -5, n_occ_mo), rng.uniform(-3, 1, n_virt_mo)]))
    occ = np.array([2.0] * n_occ_mo + [0.0] * n_virt_mo)
    return dict(C=C, eps=eps, occ=occ, overlap=S, atom_symbols=syms, atom_coords=coords, atom_ao_ranges=ranges,
                homo_index=n_occ_mo - 1, scissor_ev=0.5, kernel="resta", alpha=1.0, material="CDSE",
                include_exchange=True, spin="singlet")


OCC_SEL = np.array([0, 2, 3, 6])          # skips MOs 1, 4, 5 ("traps")
VIRT_SEL = np.array([7, 9, 10, 13])       # skips MOs 8, 11, 12


def _contiguous_copy(kw):
    """The same MOs, reordered so that the selected ones form the window next to the gap."""
    homo = kw["homo_index"]
    rest_o = [i for i in range(homo + 1) if i not in OCC_SEL]
    rest_v = [a for a in range(homo + 1, len(kw["eps"])) if a not in VIRT_SEL]
    perm = np.array(rest_o + list(OCC_SEL) + list(VIRT_SEL) + rest_v)
    out = dict(kw)
    out.update(C=kw["C"][:, perm], eps=kw["eps"][perm], occ=kw["occ"][perm])
    return out


def _soc(rng_seed=9):
    rng = np.random.default_rng(rng_seed)
    k = len(OCC_SEL) + len(VIRT_SEL)
    M = rng.normal(size=(2 * k, 2 * k)) + 1j * rng.normal(size=(2 * k, 2 * k))
    U, _ = np.linalg.qr(M)
    E = np.sort(np.concatenate([rng.uniform(-9, -5, k), rng.uniform(-3, 1, k)]))
    return U, E


@pytest.mark.parametrize("charge_type", ["mulliken", "lowdin"])
@pytest.mark.parametrize("mode", ["bse", "diagonal_bse"])
@pytest.mark.parametrize("soc", [False, True])
def test_explicit_indices_match_the_contiguous_window(charge_type, mode, soc):
    kw = _system()
    extra = dict(charge_type=charge_type, excitation_mode=mode, n_occ=len(OCC_SEL), n_virt=len(VIRT_SEL))
    if soc:
        U, E = _soc()
        extra.update(soc_U=U, soc_E=E)
    ref = ExcitonSolver(**_contiguous_copy(kw), **extra)
    sel = ExcitonSolver(**kw, **extra, occ_indices=OCC_SEL, virt_indices=VIRT_SEL)
    e_ref, _ = ref.solve(nroots=ref.ham.dim, full_diag=True, excitation_mode=mode)
    e_sel, _ = sel.solve(nroots=sel.ham.dim, full_diag=True, excitation_mode=mode)
    np.testing.assert_allclose(e_sel, e_ref, atol=1e-10)
    np.testing.assert_array_equal(sel.ham.occ_idx, OCC_SEL)
    np.testing.assert_array_equal(sel.ham.virt_idx, VIRT_SEL)
    assert (sel.n_occ, sel.n_virt) == (len(OCC_SEL), len(VIRT_SEL))


def test_explicit_indices_must_be_occupied_and_virtual():
    kw = _system()
    with pytest.raises(ValueError, match="ascending occupied / virtual"):
        ExcitonSolver(**kw, n_occ=2, n_virt=2, occ_indices=[2, 7], virt_indices=[8, 9])


def test_selection_skips_traps_and_stays_in_the_window():
    from qdex import cli
    n_mo, homo = 40, 19
    eps = np.concatenate([np.linspace(-12.0, -6.0, homo + 1), np.linspace(-4.0, 3.0, n_mo - homo - 1)])
    trap = {19, 18, 17, 15, 20, 23}                       # HOMO, HOMO-1, HOMO-2, HOMO-4; LUMO, LUMO+3
    window_ev = 5.0
    idx = np.where((eps >= eps[homo] - window_ev) & (eps <= eps[homo + 1] + window_ev))[0]
    cls = ["Facet" if i in trap else "S-like" for i in idx]
    args = SimpleNamespace(centroid_core_elements=None)
    with mock.patch("qdex.angular.states_omega", lambda *a, **k: (None, cls, None, None)):
        occ, virt, info = cli._delocalized_active_mos(np.eye(n_mo), np.eye(n_mo), eps, homo, 4, 3, [], None, [],
                                                       args, window_ev=window_ev)
    # 4 holes asked: 16, 14, 13, 12 counted down from the HOMO (19, 18, 17, 15 are traps)
    np.testing.assert_array_equal(occ, [12, 13, 14, 16])
    np.testing.assert_array_equal(virt, [21, 22, 24])
    assert info["n_trap_occ"] == 4 and info["n_trap_virt"] == 2
    assert info["n_classified"] == len(idx) and not info["short_occ"] and not info["short_virt"]
    # only the traps between the delocalized edges (HOMO 16, LUMO 21) leave the SOC window; 15 and 23 lie inside
    # the bands and stay in it
    np.testing.assert_array_equal(info["gap_trap_idx"], [17, 18, 19, 20])


def test_selection_reports_a_short_window():
    from qdex import cli
    n_mo, homo = 20, 9
    eps = np.concatenate([np.linspace(-8.0, -6.0, homo + 1), np.linspace(-4.0, -2.0, n_mo - homo - 1)])
    idx = np.where((eps >= eps[homo] - 0.5) & (eps <= eps[homo + 1] + 0.5))[0]
    cls = ["P-like"] * len(idx)
    args = SimpleNamespace(centroid_core_elements=None)
    with mock.patch("qdex.angular.states_omega", lambda *a, **k: (None, cls, None, None)):
        occ, virt, info = cli._delocalized_active_mos(np.eye(n_mo), np.eye(n_mo), eps, homo, 25, 25, [], None, [],
                                                       args, window_ev=0.5)
    assert info["short_occ"] and info["short_virt"]
    assert set(occ) | set(virt) == set(idx)
