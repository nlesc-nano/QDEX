"""
Non-adiabatic couplings between consecutive frames from their orbital overlaps.

S_ij = <i(t)|j(t+dt)> (tracked, phase-aligned). Two schemes:

* ``hst`` (Hammes-Schiffer-Tully): d = (S - S^+)/(2 dt). First order in the rotation angle: a pair of
  states rotating by theta within the step gets sin(theta)/dt instead of theta/dt.
* ``logm`` (default): the overlap is first made unitary (Loewdin, U = S (S^+ S)^{-1/2}, which also
  absorbs the norm leaking out of the active window) and d = log(U)/dt is the constant generator that
  carries the states from t to t+dt exactly: propagating c with exp(-d dt), as the Strang integrator
  does with d held constant over the step, reproduces c(t+dt) = U^+ c(t).

In the dense valence band of the CsPbX3 dots the states rotate by up to ~50 deg in a 2 fs step and
HST underestimates the coupling of the hole states by 16 % (CsPbBr3) to 37 % (CsPbCl3), the rates
(~|d|^2) by 35-120 %; at the band edges the two agree within 3-6 %.

The logm couplings cost an SVD and an eigendecomposition per step (seconds for 1300 real or 2600
complex states), so they are computed once, in parallel over the steps (``qdex --namd-nac``), and
stored as nac_<k>_to_<k+1>.npz next to the step files; the dynamics reads them.
"""
import glob
import logging
import os
import time

import numpy as np

logger = logging.getLogger(__name__)

NAC_SCHEMES = ("logm", "hst")


def nac_hst(S, dt):
    """Hammes-Schiffer-Tully finite difference (S - S^+)/(2 dt)."""
    S = np.asarray(S)
    return (S - S.conj().T) / (2.0 * dt)


def loewdin_unitary(S):
    """Closest unitary matrix to S: U = W V^+ from S = W s V^+."""
    W, _, Vh = np.linalg.svd(np.asarray(S))
    return W @ Vh


def nac_logm(S, dt, check_tol=1e-8):
    """Generator d = log(U)/dt of the Loewdin-orthonormalised overlap (anti-Hermitian; real for real S).

    U is normal, so it shares its eigenvectors with the Hermitian H = (U - U^+)/(2i), whose eigh is
    four times cheaper than a general eig. The eigenvalues of H are sin(phi): two eigenphases phi and
    pi - phi (rotations beyond 90 deg) would mix, which the residual check catches (scipy logm then)."""
    S = np.asarray(S)
    U = loewdin_unitary(S)
    H = (U - U.conj().T) / 2j
    _, V = np.linalg.eigh(H)
    UV = U @ V
    lam = np.einsum("ij,ij->j", V.conj(), UV)
    if np.max(np.abs(UV - V * lam)) > check_tol * max(1.0, np.sqrt(U.shape[0])):
        from scipy.linalg import logm
        L = logm(U)
    else:
        L = (V * (1j * np.angle(lam))) @ V.conj().T
    L = 0.5 * (L - L.conj().T)
    if not np.iscomplexobj(S):
        L = L.real
    return L / dt


def nac_from_overlap(S, dt, scheme="logm"):
    scheme = str(scheme).lower()
    if scheme == "hst":
        return nac_hst(S, dt)
    if scheme == "logm":
        return nac_logm(S, dt)
    raise ValueError(f"nac_scheme must be one of {NAC_SCHEMES}, not '{scheme}'.")


def nac_file(precompute_dir, k):
    return os.path.join(precompute_dir, f"nac_{k:05d}_to_{k + 1:05d}.npz")


def _one_step(args):
    step_file, out_file, scheme = args
    z = np.load(step_file)
    dt = float(z["time_curr_fs"] - z["time_prev_fs"]) if "time_prev_fs" in z else 2.0
    d_occ = nac_from_overlap(z["S_occ"], dt, scheme)
    d_virt = nac_from_overlap(z["S_virt"], dt, scheme)
    tmp = out_file + ".tmp.npz"
    np.savez(tmp, d_occ=d_occ, d_virt=d_virt, scheme=np.array(scheme), dt_fs=dt)
    os.replace(tmp, out_file)
    return out_file


def compute_nac_files(precompute_dir, scheme="logm", workers=None, threads_per_worker=None):
    """Write nac_<k>_to_<k+1>.npz for every step file (skips existing ones: restartable).

    workers processes run in parallel; each uses threads_per_worker BLAS threads."""
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    steps = sorted(glob.glob(os.path.join(precompute_dir, "step_*_to_*.npz")))
    todo = []
    for s in steps:
        k = int(os.path.basename(s).split("_")[1])
        out = nac_file(precompute_dir, k)
        if not os.path.exists(out):
            todo.append((s, out, scheme))
    n_cpu = os.cpu_count() or 1
    workers = int(workers or max(1, min(len(todo), n_cpu // 8)))
    threads = int(threads_per_worker or max(1, n_cpu // workers))
    logger.info(f"[NAMD NAC] {scheme} couplings for {len(todo)} of {len(steps)} steps in {precompute_dir} "
                f"({workers} workers x {threads} threads)")
    if not todo:
        return 0
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = str(threads)        # inherited by the spawned workers
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        for i, _ in enumerate(ex.map(_one_step, todo, chunksize=1)):
            if (i + 1) % max(1, len(todo) // 20) == 0:
                logger.info(f"  [NAMD NAC] {i + 1}/{len(todo)} steps ({time.time() - t0:.0f} s)")
    logger.info(f"[NAMD NAC] Done in {time.time() - t0:.1f} s")
    try:
        nac_correlation(precompute_dir)
        if os.path.exists(os.path.join(precompute_dir, "decoherence_times.npz")):
            pair_coupling_reduction(precompute_dir)
    except Exception as exc:
        logger.warning(f"[NAMD NAC] Coupling correlation not computed: {exc}")
    return len(todo)


_warned = set()


def step_nacs(precompute_dir, k, S_occ, S_virt, dt, scheme="logm"):
    """d_occ, d_virt of step k -> k+1: the stored couplings of this scheme when present, else computed."""
    scheme = str(scheme).lower()
    if scheme == "logm":
        f = nac_file(precompute_dir, k)
        if os.path.exists(f):
            z = np.load(f)
            return z["d_occ"], z["d_virt"]
        if precompute_dir not in _warned:
            _warned.add(precompute_dir)
            logger.warning(f"  [NAMD:Warn] No stored logm couplings in {precompute_dir}: computing them on the fly "
                           f"(slow for large windows; run 'qdex --namd-nac' once after the precompute).")
    return nac_from_overlap(S_occ, dt, scheme), nac_from_overlap(S_virt, dt, scheme)


def nac_correlation(precompute_dir, max_steps=500, max_lag=25, neighbours=5, n_pairs=4000, write=True):
    """Time correlation of the couplings of each channel and its correlation time.

    For the pairs that carry the flux (labels i, i+1 .. i+neighbours, the strongest n_pairs by <|d|^2>),
    C(s) = <d_ij(t) d_ij(t+s)> / <d_ij^2>, averaged over the pairs (weighted by <d_ij^2>), and
    tau_c = int_0^s0 C(s) ds up to the first zero of C. In the golden rule of a coupling that fluctuates
    (pure dephasing D(s) = e^{-s/tau}, C(s) = e^{-s/tau_c}) the Lorentzian width is 1/tau + 1/tau_c:
    k = 2 <|d|^2> tau_eff / (1 + (dE tau_eff/hbar)^2), tau_eff = (1/tau + 1/tau_c)^-1. Constant-coupling
    rates (tau_eff = tau) overestimate the transfer when tau_c << tau, as in a dense valence band whose
    states reshuffle within femtoseconds. Writes nac_correlation.npz (C_occ, C_virt, tau_c_occ_fs,
    tau_c_virt_fs, lags_fs). Returns that dict."""
    files = sorted(glob.glob(os.path.join(precompute_dir, "nac_*_to_*.npz")))[:max_steps]
    if len(files) < max_lag + 5:
        return None
    z0 = np.load(files[0])
    dt = float(z0["dt_fs"]) if "dt_fs" in z0.files else 2.0
    out = {"lags_fs": np.arange(max_lag + 1) * dt}
    for ch, key in (("occ", "d_occ"), ("virt", "d_virt")):
        n = z0[key].shape[0]
        ii = np.concatenate([np.arange(n - k) for k in range(1, neighbours + 1)])
        jj = np.concatenate([np.arange(k, n) for k in range(1, neighbours + 1)])
        X = np.array([np.load(f)[key][ii, jj] for f in files])          # (T, pairs), complex for spinors
        w = np.mean(np.abs(X) ** 2, axis=0)
        keep = np.argsort(w)[-min(n_pairs, len(w)):]
        X, w = X[:, keep], w[keep]
        X = X - X.mean(axis=0)
        T = X.shape[0]
        C = np.array([np.real(np.sum(np.mean(np.conj(X[:T - s]) * X[s:], axis=0))) / np.sum(np.mean(np.abs(X) ** 2, axis=0))
                      for s in range(max_lag + 1)])
        s0 = int(np.argmax(C <= 0)) if np.any(C <= 0) else len(C)
        trap = getattr(np, "trapezoid", None) or np.trapz
        tau_c = float(trap(C[:max(s0, 2)], dx=dt))
        out[f"C_{ch}"] = C
        out[f"tau_c_{ch}_fs"] = max(tau_c, 1e-3)
    if write:
        np.savez(os.path.join(precompute_dir, "nac_correlation.npz"), **out)
    logger.info(f"[NAMD NAC] Coupling correlation time: holes {out['tau_c_occ_fs']:.2f} fs, "
                f"electrons {out['tau_c_virt_fs']:.2f} fs")
    return out


def pair_coupling_reduction(precompute_dir, max_lag_fs=60.0, max_steps=None, write=True):
    """Pair-resolved effect of the coupling fluctuations on the golden-rule rate.

    For every pair of each channel, with C_ij(s) = Re<d_ij*(t) d_ij(t+s)>_t / <|d_ij|^2>_t (truncated at its
    first non-positive lag) and the pair dephasing D_ij(s) = exp(-s^2 / 2 tau_ij^2) of decoherence_times.npz,

        r_ij = sum_s w_s C_ij(s) D_ij(s) / sum_s w_s D_ij(s)      (trapezoid weights, s <= max_lag_fs)

    is the fraction of the constant-coupling rate that survives; the PME (pme_tau: pairs_nac) uses
    tau_eff_ij = r_ij tau_ij, which is tau_ij for a coupling that keeps its phase (r = 1). The lagged
    products are accumulated streaming over the steps with a buffer of the last max_lag_fs/dt couplings.
    Adds r_occ, r_virt (float32) to nac_correlation.npz. Returns (r_occ, r_virt)."""
    from collections import deque
    files = sorted(glob.glob(os.path.join(precompute_dir, "nac_*_to_*.npz")))
    if max_steps:
        files = files[:max_steps]
    dec = np.load(os.path.join(precompute_dir, "decoherence_times.npz"))
    z0 = np.load(files[0])
    dt = float(z0["dt_fs"]) if "dt_fs" in z0.files else 2.0
    L = int(max(1, round(max_lag_fs / dt)))
    w = np.ones(L + 1); w[0] = 0.5
    result = {}
    for ch, key, tkey in (("occ", "d_occ", "tau_occ"), ("virt", "d_virt", "tau_virt")):
        n = z0[key].shape[0]
        acc = np.zeros((L + 1, n, n))
        buf = deque(maxlen=L + 1)
        T = 0
        for f in files:
            d = np.load(f)[key]
            buf.appendleft(d)                       # buf[s] = d(t - s)
            for s in range(len(buf)):
                acc[s] += np.real(np.conj(buf[s]) * d)
            T += 1
        counts = np.array([T - s for s in range(L + 1)], dtype=float)
        C = acc / counts[:, None, None]
        C0 = np.maximum(C[0], 1e-30)
        C = C / C0
        alive = np.cumprod(C > 0, axis=0).astype(bool)        # truncate each pair at its first zero
        C = np.where(alive, C, 0.0)
        tau = np.asarray(dec[tkey], float)[:n, :n]
        s_fs = (np.arange(L + 1) * dt)[:, None, None]
        D = np.exp(-0.5 * (s_fs / np.maximum(tau, 1e-6)) ** 2)
        r = np.einsum("s,sij->ij", w, C * D) / np.einsum("s,sij->ij", w, D)
        r = np.clip(r, 0.0, 1.0)
        np.fill_diagonal(r, 1.0)
        result[f"r_{ch}"] = r.astype(np.float32)
        off = ~np.eye(n, dtype=bool)
        wgt = C0[off]
        logger.info(f"[NAMD NAC] Pair coupling reduction ({ch}): median {np.median(r[off]):.3f}, "
                    f"|d|^2-weighted mean {np.sum(r[off] * wgt) / np.sum(wgt):.3f} (max lag {L * dt:.0f} fs)")
        del acc, C, D
    if write:
        path = os.path.join(precompute_dir, "nac_correlation.npz")
        old = dict(np.load(path)) if os.path.exists(path) else {}
        old.update(result)
        np.savez(path, **old)
    return result["r_occ"], result["r_virt"]
