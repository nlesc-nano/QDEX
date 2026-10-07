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
