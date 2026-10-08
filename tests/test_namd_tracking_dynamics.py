"""NAMD fixes: K^x in tracked frames, Kramers-pair alignment, DISH targets, PME line widths."""
import numpy as np
import pytest

from qdex.namd.precompute import (_pair_energies, _degenerate_groups, _group_average_pairs,
                                  align_degenerate_blocks, align_spinor_phases_and_crossings)
from qdex.namd.integrator import step_dish_batch
from qdex.namd.master_equation import propagate_pme_tensor


def test_tracked_pair_energies_keep_the_exchange():
    rng = np.random.default_rng(1)
    eps_o, eps_v = np.sort(rng.normal(-1, 0.3, 4)), np.sort(rng.normal(1, 0.3, 3))
    Kx, Kd = rng.uniform(0, 0.02, (4, 3)), rng.uniform(0.1, 0.3, (4, 3))
    i, a = np.repeat(np.arange(4), 3), np.tile(np.arange(3), 4)
    d = dict(eps_occ=eps_o, eps_virt=eps_v, Kx_mat=Kx, Kd_mat=Kd, kx_factor=2.0)
    E0 = eps_v[a] - eps_o[i] + 2 * Kx[i, a] - Kd[i, a]                    # frame 0 (singlets)
    assert np.allclose(_pair_energies(d, i, a), E0)
    # a relabelling of the states permutes K^x with them: same pair energies, permuted
    po, pv = np.array([1, 0, 3, 2]), np.array([2, 0, 1])
    d2 = dict(eps_occ=eps_o[po], eps_virt=eps_v[pv], Kx_mat=Kx[np.ix_(po, pv)], Kd_mat=Kd[np.ix_(po, pv)],
              kx_factor=2.0)
    E_perm = _pair_energies(d2, i, a).reshape(4, 3)
    assert np.allclose(E_perm, E0.reshape(4, 3)[np.ix_(po, pv)])


def test_kramers_pairs_are_parallel_transported():
    rng = np.random.default_rng(0)
    U_prev = np.linalg.qr(rng.normal(size=(8, 4)) + 1j * rng.normal(size=(8, 4)))[0]

    def su2():
        a, b = rng.normal(size=2) + 1j * rng.normal(size=2)
        n = np.sqrt(abs(a) ** 2 + abs(b) ** 2)
        return np.array([[a, -b.conj()], [b, a.conj()]]) / n
    R = np.zeros((4, 4), complex)
    R[:2, :2], R[2:, 2:] = su2(), su2()
    U_next = U_prev @ R                                   # same states, new SU(2) frame per pair
    S = U_prev.conj().T @ U_next
    S1, (U1,), _ = align_spinor_phases_and_crossings(S, [U_next.copy()])
    assert np.abs(S1 - np.diag(np.diag(S1))).max() > 0.1  # U(1) alone leaves a spurious coupling
    S2, (U2,), n = align_degenerate_blocks(S1, [U1], np.array([0.0, 0.0, 1.0, 1.0]))
    assert n == 2
    assert np.allclose(S2, np.eye(4), atol=1e-12)
    assert np.allclose(U_prev.conj().T @ U2, S2)          # the coefficients carry the same rotation


def test_group_average_is_the_basis_independent_part():
    M = np.arange(12.0).reshape(4, 3)
    g = _degenerate_groups(np.array([0.0, 0.0, 1.0, 2.0]), 1e-5)
    A = _group_average_pairs(M, g, [])
    assert np.allclose(A[:2], M[:2].mean(axis=0)) and np.allclose(A[2:], M[2:])
    assert A[:2].sum() == pytest.approx(M[:2].sum())


def test_dish_never_lands_on_a_pair_that_is_not_stored():
    np.random.seed(0)
    n_states, n_traj = 3, 2000
    C = np.zeros((n_states, n_traj), complex)
    C[0], C[1], C[2] = np.sqrt(0.4), np.sqrt(0.3), np.sqrt(0.3)
    allowed = np.ones((n_states, n_traj), bool)
    allowed[2] = False
    tau = np.full((n_states, n_states), 1.0)
    _, act, hopped = step_dish_batch(C, np.zeros(n_traj, int), np.zeros((n_states, n_traj)), 5.0,
                                     tau_mat=tau, beta=None, detailed_balance=False, allowed=allowed)
    assert hopped.any() and not np.any(act == 2)


def test_pme_pair_line_widths_and_detailed_balance():
    d = np.array([[0.0, 0.01], [-0.01, 0.0]])
    eps_v = np.array([1.0, 1.05])
    P = np.array([[0.0, 1.0]])                            # one hole, electron in the upper state
    kw = dict(E_mat=np.zeros((1, 2)), d_occ=np.zeros((1, 1)), d_virt=d, dt_fs=1.0, temp_k=300.0,
              eps_occ=np.array([0.0]), eps_virt=eps_v)
    tau = np.full((2, 2), 10.0)
    P_pair = propagate_pme_tensor(P, tau_virt_mat=tau, tau_occ_mat=np.full((1, 1), 10.0), **kw)
    P_unif = propagate_pme_tensor(P, tau_dec_fs=10.0, **kw)
    assert np.allclose(P_pair, P_unif)                     # a uniform matrix equals the single tau
    for _ in range(20000):                                 # long-time limit: Boltzmann
        P = propagate_pme_tensor(P, tau_virt_mat=tau, tau_occ_mat=np.full((1, 1), 10.0), **kw)
    assert P.sum() == pytest.approx(1.0)
    assert P[0, 1] / P[0, 0] == pytest.approx(np.exp(-0.05 / (8.617333e-5 * 300.0)), rel=1e-3)


@pytest.mark.parametrize("cplx", [False, True])
def test_logm_nac_is_the_exact_generator_of_the_step(cplx):
    from scipy.linalg import expm, logm
    from qdex.namd.nac import nac_logm, nac_hst
    rng = np.random.default_rng(3)
    n, dt = 40, 2.0
    A = rng.normal(size=(n, n)) * 0.15
    if cplx:
        A = A + 1j * rng.normal(size=(n, n)) * 0.15
    G = A - A.conj().T                                    # anti-Hermitian generator
    G *= 2.5 / np.abs(np.linalg.eigvals(G)).max()         # largest rotation 2.5 rad (< pi): exercises the
    U = expm(G)                                           # scipy fallback for angles beyond pi/2
    d = nac_logm(U, dt)
    assert np.allclose(d * dt, G, atol=1e-10)
    assert np.allclose(d * dt, 0.5 * (logm(U) - logm(U).conj().T), atol=1e-10)
    assert np.iscomplexobj(d) == cplx
    # propagating with exp(-d dt) maps the coefficients exactly: c(t+dt) = U^+ c(t)
    c = rng.normal(size=n) + 0j
    assert np.allclose(expm(-d * dt) @ c, U.conj().T @ c)
    # HST is first order: smaller coupling for large rotations
    assert np.linalg.norm(nac_hst(U, dt)) < np.linalg.norm(d)
    # norm leaking out of the window (scaled overlap) does not change the generator
    assert np.allclose(nac_logm(0.9 * U, dt), d, atol=1e-10)


def test_stored_nacs_are_read_back(tmp_path):
    from qdex.namd.nac import compute_nac_files, step_nacs, nac_logm
    rng = np.random.default_rng(5)
    from scipy.linalg import expm
    for k in range(2):
        So = expm(np.triu(rng.normal(size=(6, 6)) * 0.2, 1) - np.triu(rng.normal(size=(6, 6)) * 0.2, 1).T)
        Sv = np.eye(4)
        np.savez(tmp_path / f"step_{k:05d}_to_{k + 1:05d}.npz", S_occ=So, S_virt=Sv,
                 time_prev_fs=2.0 * k, time_curr_fs=2.0 * (k + 1))
    assert compute_nac_files(str(tmp_path), workers=1, threads_per_worker=1) == 2
    z = np.load(tmp_path / "step_00001_to_00002.npz")
    d_o, d_v = step_nacs(str(tmp_path), 1, z["S_occ"], z["S_virt"], 2.0, "logm")
    assert np.allclose(d_o, nac_logm(z["S_occ"], 2.0)) and np.allclose(d_v, 0.0)
    assert compute_nac_files(str(tmp_path), workers=1, threads_per_worker=1) == 0     # restartable


def test_pair_dephasing_cumulant_limits():
    from qdex.namd.precompute import _pair_dephasing_times
    hbar, dt, N = 0.6582119569, 2.0, 4000
    rng = np.random.default_rng(7)
    t = np.arange(N) * dt
    # state 0 fixed; state 1: slow, large fluctuation (period 2 ps, std 0.05 eV) -> short-time Gaussian
    # state 2: fast fluctuation (white noise, std 0.03 eV) -> motional narrowing, slower dephasing
    E = np.zeros((N, 3))
    E[:, 1] = 0.05 * np.sqrt(2) * np.sin(2 * np.pi * t / 2000.0 + 0.3)
    E[:, 2] = 0.03 * rng.normal(size=N)
    t1e, sig = _pair_dephasing_times(E, dt, 1500.0, hbar)
    tau = t1e / np.sqrt(2)
    assert tau[0, 1] == pytest.approx(hbar / sig[0, 1], rel=0.05)            # Gaussian limit
    assert tau[0, 2] > 3 * hbar / sig[0, 2]                                 # narrowed: much slower than hbar/sigma
    assert np.allclose(tau, tau.T, equal_nan=True)
    # direct single-pair cumulant
    x = E[:, 1] - E[:, 0]; x = x - x.mean(); L = 300
    A = np.array([np.mean(x[:N - s] * x[s:]) for s in range(L)])
    I1 = np.concatenate([[0], np.cumsum(0.5 * (A[1:] + A[:-1]) * dt)])
    g = np.concatenate([[0], np.cumsum(0.5 * (I1[1:] + I1[:-1]) * dt)]) / hbar ** 2
    k = np.argmax(g >= 1.0)
    assert t1e[0, 1] == pytest.approx((k - 1 + (1 - g[k - 1]) / (g[k] - g[k - 1])) * dt, rel=1e-6)


def test_nac_correlation_time_of_an_exponential_coupling(tmp_path):
    from qdex.namd.nac import nac_correlation
    rng = np.random.default_rng(11)
    n, T, dt, tau_c = 12, 3000, 2.0, 6.0
    phi = np.exp(-dt / tau_c)                       # AR(1): C(s) = phi^(s/dt) = exp(-s/tau_c)
    x = np.zeros((T, n, n))
    for t in range(1, T):
        x[t] = phi * x[t - 1] + np.sqrt(1 - phi ** 2) * rng.normal(size=(n, n))
    for t in range(T):
        d = np.triu(x[t], 1); d = d - d.T            # antisymmetric couplings
        np.savez(tmp_path / f"nac_{t:05d}_to_{t + 1:05d}.npz", d_occ=d, d_virt=d, dt_fs=dt)
    out = nac_correlation(str(tmp_path), max_steps=T, max_lag=30, write=False)
    # exact integral of exp(-s/tau_c) up to the first zero crossing ~ tau_c (minus the trapezoid offset)
    assert out["tau_c_occ_fs"] == pytest.approx(tau_c, rel=0.15)
    assert out["C_virt"][1] == pytest.approx(phi, abs=0.05)


def test_pair_coupling_reduction(tmp_path):
    from qdex.namd.nac import pair_coupling_reduction
    rng = np.random.default_rng(3)
    n, T, dt = 3, 4000, 2.0
    # pair (0,1): constant coupling -> r = 1; pair (0,2): white noise (decorrelates within a step) -> r small
    d = np.zeros((T, n, n))
    d[:, 0, 1] = 0.05
    d[:, 0, 2] = 0.05 * rng.normal(size=T)
    d = d - np.transpose(d, (0, 2, 1))
    for t in range(T):
        np.savez(tmp_path / f"nac_{t:05d}_to_{t + 1:05d}.npz", d_occ=d[t], d_virt=d[t], dt_fs=dt)
    tau = np.full((n, n), 20.0)
    np.savez(tmp_path / "decoherence_times.npz", tau_occ=tau, tau_virt=tau)
    r_o, _ = pair_coupling_reduction(str(tmp_path), max_lag_fs=60.0)
    assert r_o[0, 1] == pytest.approx(1.0, abs=1e-6)
    # white noise: only the s = 0 term survives -> r = 0.5 / sum_s w_s D(s)
    s = np.arange(31) * dt; w = np.ones(31); w[0] = 0.5
    assert r_o[0, 2] == pytest.approx(0.5 / np.sum(w * np.exp(-0.5 * (s / 20.0) ** 2)), rel=0.25)
