"""One-centre extension of xs (qdex.onecentre): matrix-free products against an explicit dense build."""
import os

import numpy as np
import pytest

libint_cpp = pytest.importorskip("libint_cpp")
if not hasattr(libint_cpp, "compute_onecentre_coulomb"):
    pytest.skip("libint_cpp built without compute_onecentre_coulomb", allow_module_level=True)

from qdex.io_utils import parse_basis, build_shell_dicts, build_atom_ao_ranges, count_ao_from_shells
from qdex.onecentre import OneCentreKernel

BASIS = os.path.join(os.path.dirname(__file__), "CdSe", "2.0nm", "BASIS_MOLOPT_UZH")


def _system(cutoff, tail="dipole"):
    syms = ["Cd", "Se", "Cl", "Se"]
    xyz = np.array([[0.0, 0.0, 0.0], [2.62, 0.0, 0.0], [-1.0, 2.3, 0.4], [0.3, -0.5, 3.9]])
    bd = parse_basis(BASIS, "DZVP-MOLOPT-GGA-GTH", required_elements=set(syms))
    shells = build_shell_dicts(syms, xyz, bd)
    ranges = build_atom_ao_ranges(shells)
    S = np.asarray(libint_cpp.overlap(shells, 1))
    G = np.asarray(libint_cpp.compute_aabb_coulomb(shells, 1)) * 27.211386245988
    return (OneCentreKernel(shells, ranges, xyz, S, G, cutoff_ang=cutoff, nthreads=2, tail=tail),
            count_ao_from_shells(shells))


def _dense(oc, screened=False):
    """Delta (+ dipole tail) over ordered same-atom AO pairs (all atoms), and the pair -> (mu, nu) map.
    screened: the direct-term version, every atom-pair block scaled by S_AB."""
    offs = np.concatenate([[0], np.cumsum(oc.size ** 2)])
    D = np.zeros((offs[-1], offs[-1]))
    Sab = oc.screen_ratio if screened else np.ones((oc.n_atoms, oc.n_atoms))
    if getattr(oc, "tail", "none") == "dipole":
        # moment vectors of the pair basis: m[(A, k), (mu nu)] = m_k^A[mu, nu]
        Mvec = np.zeros((4 * oc.n_atoms, offs[-1]))
        for nA, idx in oc.kinds.items():
            atoms = np.where(oc.size == nA)[0]
            for j, A in enumerate(atoms):
                Mvec[4 * A:4 * A + 4, offs[A]:offs[A + 1]] = oc.kmom[nA][j].reshape(4, -1)
        N = oc.n_atoms
        D += Mvec.T @ (oc.M.reshape(N, 4, N, 4) * Sab[:, None, :, None]).reshape(4 * N, 4 * N) @ Mvec
    for g in oc.groups:
        for k, (A, B) in enumerate(zip(g["A"], g["B"])):
            blk = g["dx"][k] * Sab[A, B]
            D[offs[A]:offs[A + 1], offs[B]:offs[B + 1]] += blk
            if A != B:
                D[offs[B]:offs[B + 1], offs[A]:offs[A + 1]] += blk.T
    mu = np.concatenate([np.repeat(np.arange(a0, a1), a1 - a0) for (a0, a1) in oc.ranges])
    nu = np.concatenate([np.tile(np.arange(a0, a1), a1 - a0) for (a0, a1) in oc.ranges])
    return D, mu, nu


def _pairs(comps, mu, nu, left, right):
    """F[(mu nu), (p q)] = sum_s conj(L_s[mu, p]) R_s[nu, q]."""
    F = 0.0
    for c in comps:
        L, R = c[left], c[right]
        F = F + np.conj(L[mu])[:, :, None] * R[nu][:, None, :]
    return F.reshape(len(mu), -1)


@pytest.mark.parametrize("cutoff,tail", [(0.0, "none"), (3.0, "none"), (0.0, "dipole"), (3.0, "dipole")])
@pytest.mark.parametrize("ncomp,cplx", [(1, False), (2, True)])
def test_products_match_dense(cutoff, tail, ncomp, cplx):
    oc, nao = _system(cutoff, tail)
    D, mu, nu = _dense(oc)
    rng = np.random.default_rng(3)
    no, nv = 4, 5

    def rand(n, m):
        a = rng.normal(size=(n, m))
        return a + 1j * rng.normal(size=(n, m)) if cplx else a

    comps = [(rand(nao, no), rand(nao, nv)) for _ in range(ncomp)]
    X = rand(no, nv)
    # exchange: K_ia,jb = conj(F_ov)^T D F_ov
    Fov = _pairs(comps, mu, nu, 0, 1)
    Kx = np.conj(Fov).T @ D @ Fov
    assert np.allclose(oc.exchange_action([(comps, X)])[0].reshape(-1), Kx @ X.reshape(-1))
    # direct: K_ia,jb = (ij|D|ab) = conj(F_oo[ij]) D F_vv[ab]
    Foo = _pairs([(T, T) for T, _ in comps], mu, nu, 0, 1).reshape(len(mu), no, no)
    Fvv = _pairs([(V, V) for _, V in comps], mu, nu, 0, 1).reshape(len(mu), nv, nv)
    Kd = np.einsum("fij,fg,gab->iajb", np.conj(Foo), D, Fvv).reshape(no * nv, no * nv)
    # the spinor/UKS convention: sum over spin components on both sides (comps paired independently)
    Foo_s = [_pairs([(T, T)], mu, nu, 0, 1).reshape(len(mu), no, no) for T, _ in comps]
    Fvv_s = [_pairs([(V, V)], mu, nu, 0, 1).reshape(len(mu), nv, nv) for _, V in comps]
    Kd = sum(np.einsum("fij,fg,gab->iajb", np.conj(a), D, b) for a in Foo_s for b in Fvv_s).reshape(no * nv, -1)
    assert np.allclose(oc.direct_action(comps, X).reshape(-1), Kd @ X.reshape(-1))
    # diagonals
    vi, va = np.divmod(np.arange(no * nv), nv)
    assert np.allclose(oc.exchange_diag(comps, vi, va), np.diag(Kx).real)
    assert np.allclose(oc.direct_diag(comps).reshape(-1), np.diag(Kd).real)


def test_exchange_couples_channels():
    oc, nao = _system(3.0)
    D, mu, nu = _dense(oc)
    rng = np.random.default_rng(5)
    ca = [(rng.normal(size=(nao, 3)), rng.normal(size=(nao, 4)))]
    cb = [(rng.normal(size=(nao, 2)), rng.normal(size=(nao, 4)))]
    Xa, Xb = rng.normal(size=(3, 4)), rng.normal(size=(2, 4))
    Fa, Fb = _pairs(ca, mu, nu, 0, 1), _pairs(cb, mu, nu, 0, 1)
    F = np.hstack([Fa, Fb])
    ref = F.T @ D @ F @ np.concatenate([Xa.ravel(), Xb.ravel()])
    ka, kb = oc.exchange_action([(ca, Xa), (cb, Xb)])
    assert np.allclose(np.concatenate([ka.ravel(), kb.ravel()]), ref)


def test_full_cutoff_restores_exact_one_centre_on_diagonal_products():
    """With every atom pair inside R_c the diagonal products carry (lambda lambda|lambda lambda), not xs."""
    oc, _ = _system(50.0)
    for g in oc.groups:
        assert g["dx"].shape[0] > 0
    n_pairs = sum(len(g["A"]) for g in oc.groups)
    assert n_pairs == oc.n_atoms + oc.n_atoms * (oc.n_atoms - 1) // 2


def test_onecentre_layout_matches_xs_integrals():
    """raw (mu mu|ka ka) of the libint one-centre blocks = compute_aabb_coulomb."""
    syms = ["Cd", "Se"]
    xyz = np.array([[0.0, 0.0, 0.0], [2.62, 0.3, -0.2]])
    bd = parse_basis(BASIS, "DZVP-MOLOPT-GGA-GTH", required_elements=set(syms))
    shells = build_shell_dicts(syms, xyz, bd)
    ranges = build_atom_ao_ranges(shells)
    start = [0, sum(1 for sh in shells if sh["atom_idx"] == 0), len(shells)]
    blocks = libint_cpp.compute_onecentre_coulomb(shells, start, [(0, 0), (0, 1), (1, 1)], 1)
    G = np.asarray(libint_cpp.compute_aabb_coulomb(shells, 1))
    for (A, B), blk in zip([(0, 0), (0, 1), (1, 1)], blocks):
        (a0, a1), (b0, b1) = ranges[A], ranges[B]
        nA, nB = a1 - a0, b1 - b0
        dA, dB = np.arange(nA) * (nA + 1), np.arange(nB) * (nB + 1)
        assert np.allclose(np.asarray(blk)[np.ix_(dA, dB)], G[a0:a1, b0:b1], atol=1e-10)
        # (mu nu|ka la) = (nu mu|ka la) = (mu nu|la ka)
        b4 = np.asarray(blk).reshape(nA, nA, nB, nB)
        assert np.allclose(b4, b4.transpose(1, 0, 2, 3)) and np.allclose(b4, b4.transpose(0, 1, 3, 2))


def test_dipole_tail_is_symmetric_and_excludes_exact_pairs():
    oc, _ = _system(3.0, "dipole")
    assert np.allclose(oc.M, oc.M.T)
    M4 = oc.M.reshape(oc.n_atoms, 4, oc.n_atoms, 4)
    for A in range(oc.n_atoms):
        assert np.all(M4[A, :, A, :] == 0.0)
    assert np.all(M4[:, 0, :, 0] == 0.0)       # charge-charge stays in xs


def test_direct_term_is_screened_by_the_atom_pair_ratio():
    oc, nao = _system(3.0, "dipole")
    rng = np.random.default_rng(11)
    r = rng.uniform(0.2, 1.0, size=(oc.n_atoms, oc.n_atoms))
    r = 0.5 * (r + r.T)
    G = np.zeros((nao, nao))
    for A, (a0, a1) in enumerate(oc.ranges):
        for B, (b0, b1) in enumerate(oc.ranges):
            G[a0:a1, b0:b1] = r[A, B]
    W = G * oc._G_atom.repeat(oc.size, 0).repeat(oc.size, 1) / np.repeat(np.repeat(
        np.array([[oc.size[A] * oc.size[B] for B in range(oc.n_atoms)] for A in range(oc.n_atoms)]), oc.size, 0), oc.size, 1)
    oc.set_screening(W)
    assert np.allclose(oc.screen_ratio, r)
    Dd, mu, nu = _dense(oc, screened=True)
    Dx, _, _ = _dense(oc, screened=False)
    no, nv = 3, 4
    T, V = rng.normal(size=(nao, no)), rng.normal(size=(nao, nv))
    X = rng.normal(size=(no, nv))
    Foo = _pairs([(T, T)], mu, nu, 0, 1).reshape(len(mu), no, no)
    Fvv = _pairs([(V, V)], mu, nu, 0, 1).reshape(len(mu), nv, nv)
    Kd = np.einsum("fij,fg,gab->iajb", Foo, Dd, Fvv).reshape(no * nv, -1)
    assert np.allclose(oc.direct_action([(T, V)], X).reshape(-1), Kd @ X.reshape(-1))
    assert np.allclose(oc.direct_diag([(T, V)]).reshape(-1), np.diag(Kd))
    Fov = _pairs([(T, V)], mu, nu, 0, 1)
    assert np.allclose(oc.exchange_action([([(T, V)], X)])[0].reshape(-1), Fov.T @ Dx @ Fov @ X.reshape(-1))
