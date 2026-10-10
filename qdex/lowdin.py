import numpy as np
from qdex.device_utils import is_gpu, to_tensor, to_numpy

_FACTOR_CACHE = {}
_SQRT_CACHE = {}
# Device of the eigen-decomposition of S (system.lowdin_device): 'cpu', or a CUDA device ('cuda', 'cuda:1'), which
# uses cuSOLVER's 64-bit cusolverDnXsyevd when the matrix fits in GPU memory (CPU otherwise).
LOWDIN_DEVICE = "cpu"


def set_lowdin_device(device=None):
    global LOWDIN_DEVICE
    d = str(device or "cpu").strip().lower()
    LOWDIN_DEVICE = "cpu" if d in ("", "none", "numpy", "cpu") else d


def _key(S):
    """Cheap fingerprint of an overlap matrix (same S -> same key within a run)."""
    S = np.asarray(S)
    return (S.shape, float(np.trace(S)), float(S[::97, ::89].sum()), float(S[-1, ::53].sum()))


def _disk_cache_path(S, k):
    """File of the eigen-decomposition of S in $QDEX_LOWDIN_CACHE (a directory), or None. Runs on the same orbital
    file (e.g. several QP models or trap settings of one dot) then diagonalize S once."""
    import os
    d = os.environ.get("QDEX_LOWDIN_CACHE")
    if not d:
        return None
    import hashlib
    h = hashlib.sha1(repr(k).encode() + np.ascontiguousarray(S[:: max(1, S.shape[0] // 64)]).tobytes()).hexdigest()[:16]
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"lowdin_{S.shape[0]}_{h}")


def _eigh(S):
    """Eigen-decomposition of S: cuSOLVER on LOWDIN_DEVICE if it is a CUDA device and S fits, else LAPACK divide and
    conquer on the CPU. numpy's syevd needs a 1 + 6n + 2n^2 workspace, which overflows 32-bit LAPACK integers past
    n = 32767; there the ILP64 dsyevd_64 of MKL is called, and only without it scipy's slower syevr."""
    import logging, time
    log = logging.getLogger(__name__)
    n = S.shape[0]
    t0 = time.time()
    if LOWDIN_DEVICE.startswith("cuda"):
        try:
            w, V = _eigh_cusolver(S, LOWDIN_DEVICE)
            log.info(f"  [Lowdin] S ({n} AOs) diagonalized on {LOWDIN_DEVICE} (cusolverDnXsyevd) in {time.time() - t0:.1f} s")
            return w, V
        except Exception as e:  # no CUDA, too large for the GPU, cuSOLVER without 64-bit syevd
            log.info(f"  [Lowdin] GPU diagonalization not used ({e}); diagonalizing S on the CPU")
            t0 = time.time()
    if n <= 32000:
        w, V = np.linalg.eigh(S)
        how = "LAPACK syevd"
    else:
        try:
            w, V = _eigh_mkl_ilp64(S)
            how = "MKL ILP64 dsyevd_64"
        except (OSError, AttributeError, RuntimeError) as e:
            log.info(f"  [Lowdin] MKL ILP64 dsyevd_64 not available ({e}); using scipy syevr")
            from scipy.linalg import eigh
            w, V = eigh(S, driver="evr", overwrite_a=False, check_finite=False)
            how = "scipy syevr"
    log.info(f"  [Lowdin] S ({n} AOs) diagonalized on the CPU ({how}) in {time.time() - t0:.1f} s")
    return w, V


def _mkl_rt():
    import ctypes, glob, os, sys
    for d in (os.path.join(sys.prefix, "lib"), os.environ.get("MKLROOT", "") and os.path.join(os.environ["MKLROOT"], "lib")):
        libs = sorted(glob.glob(os.path.join(d, "libmkl_rt.so*"))) if d else []
        if libs:
            return ctypes.CDLL(libs[0])
    return ctypes.CDLL("libmkl_rt.so")


def _eigh_mkl_ilp64(S):
    """w, V by MKL's Fortran dsyevd_64 (64-bit integers, oneMKL >= 2023), from the libmkl_rt numpy is linked to. The
    symmetric S is passed as a column-major copy, so the eigenvectors end up in the rows of the C array."""
    import ctypes
    f = _mkl_rt().dsyevd_64
    f.restype = None
    n = S.shape[0]
    A = np.array(S, dtype=np.float64, order="C", copy=True)
    w = np.empty(n)
    i64 = ctypes.c_int64
    N, info = i64(n), i64(0)

    def call(work, lwork, iwork, liwork):
        f(ctypes.c_char_p(b"V"), ctypes.c_char_p(b"L"), ctypes.byref(N), A.ctypes.data_as(ctypes.c_void_p),
          ctypes.byref(N), w.ctypes.data_as(ctypes.c_void_p), work.ctypes.data_as(ctypes.c_void_p), ctypes.byref(lwork),
          iwork.ctypes.data_as(ctypes.c_void_p), ctypes.byref(liwork), ctypes.byref(info), ctypes.c_size_t(1),
          ctypes.c_size_t(1))

    wq, iq = np.zeros(1), np.zeros(1, dtype=np.int64)
    call(wq, i64(-1), iq, i64(-1))  # workspace query
    lwork, liwork = int(wq[0]) + 1, int(iq[0]) + 1
    call(np.empty(lwork), i64(lwork), np.empty(liwork, dtype=np.int64), i64(liwork))
    if info.value != 0:
        raise RuntimeError(f"dsyevd_64 info = {info.value}")
    return w, A.T


def _cusolver():
    """libcusolver with the 64-bit cusolverDnXsyevd: $QDEX_CUSOLVER if set, else the newest found next to torch
    (nvidia/cu13 from the nvidia-cusolver wheel before nvidia/cusolver from nvidia-cusolver-cu12)."""
    import ctypes, glob, os
    path = os.environ.get("QDEX_CUSOLVER")
    if not path:
        import torch
        base = os.path.join(os.path.dirname(torch.__file__), "..", "nvidia")
        libs = sorted(glob.glob(os.path.join(base, "cu1*", "lib", "libcusolver.so*")), reverse=True) + \
            sorted(glob.glob(os.path.join(base, "cusolver", "lib", "libcusolver.so*")))
        if not libs:
            raise OSError("no libcusolver next to torch; set QDEX_CUSOLVER")
        path = libs[0]
    return ctypes.CDLL(path)


def _eigh_cusolver(S, device):
    """w, V of S on a CUDA device by cusolverDnXsyevd (64-bit API; torch.linalg.eigh fails for n of a few 10^4).
    Needs about 8 n^2 bytes for S plus the solver workspace (~32 n^2 bytes) in free GPU memory; cuSOLVER 11.x
    rejects n > ~32k, 12.x (CUDA 13) does not. Raises if the matrix does not fit or the call fails."""
    import ctypes
    import torch
    if not torch.cuda.is_available():
        raise RuntimeError("torch has no CUDA device")
    dev = torch.device(device)
    lib = _cusolver()
    c_p, c_i64, c_sz = ctypes.c_void_p, ctypes.c_int64, ctypes.c_size_t
    R64F, VECTOR, LOWER = 1, 1, 0
    n = S.shape[0]
    h, prm = c_p(), c_p()
    with torch.cuda.device(dev):
        if lib.cusolverDnCreate(ctypes.byref(h)) or lib.cusolverDnCreateParams(ctypes.byref(prm)):
            raise RuntimeError("cusolverDnCreate failed")
        try:
            lib.cusolverDnSetStream(h, c_p(torch.cuda.current_stream(dev).cuda_stream))
            probe = torch.zeros(16, dtype=torch.float64, device=dev)
            dws, hws = c_sz(), c_sz()
            st = lib.cusolverDnXsyevd_bufferSize(h, prm, VECTOR, LOWER, c_i64(n), R64F, c_p(probe.data_ptr()), c_i64(n),
                                                 R64F, c_p(probe.data_ptr()), R64F, ctypes.byref(dws), ctypes.byref(hws))
            if st != 0:
                raise RuntimeError(f"cusolverDnXsyevd rejects n = {n} (status {st}; cuSOLVER >= 12 needed)")
            free = torch.cuda.mem_get_info(dev)[0]
            need = 8 * n * n + 8 * n + dws.value
            if need > 0.95 * free:
                raise RuntimeError(f"needs {need / 1e9:.0f} GB, {free / 1e9:.0f} GB free on {device}")
            A = torch.from_numpy(np.ascontiguousarray(S, dtype=np.float64)).to(dev)
            w = torch.empty(n, dtype=torch.float64, device=dev)
            dbuf = torch.empty(dws.value, dtype=torch.uint8, device=dev)
            hbuf = np.empty(max(hws.value, 1), dtype=np.uint8)
            info = torch.zeros(1, dtype=torch.int32, device=dev)
            st = lib.cusolverDnXsyevd(h, prm, VECTOR, LOWER, c_i64(n), R64F, c_p(A.data_ptr()), c_i64(n), R64F,
                                      c_p(w.data_ptr()), R64F, c_p(dbuf.data_ptr()), c_sz(dws.value),
                                      hbuf.ctypes.data_as(c_p), c_sz(hws.value), c_p(info.data_ptr()))
            torch.cuda.synchronize(dev)
            if st != 0 or int(info.item()) != 0:
                raise RuntimeError(f"cusolverDnXsyevd status {st}, info {int(info.item())}")
            del dbuf
            # column-major result on the row-major buffer: eigenvector k is row k of A
            return w.cpu().numpy(), A.cpu().numpy().T
        finally:
            lib.cusolverDnDestroyParams(prm)
            lib.cusolverDnDestroy(h)
            torch.cuda.empty_cache()


def lowdin_factor(S):
    """Eigen-decomposition of S (eigenvalues, eigenvectors), computed once per S and cached (in memory, and on disk
    when $QDEX_LOWDIN_CACHE names a directory)."""
    import os
    S = S.toarray() if hasattr(S, "toarray") else np.asarray(S, dtype=np.float64)
    k = _key(S)
    if k not in _FACTOR_CACHE:
        _FACTOR_CACHE.clear()
        path = _disk_cache_path(S, k)
        if path and os.path.exists(path + "_V.npy") and os.path.exists(path + "_w.npy"):
            _FACTOR_CACHE[k] = (np.load(path + "_w.npy"), np.load(path + "_V.npy"))
            return _FACTOR_CACHE[k]
        w, V = _eigh(S)
        _FACTOR_CACHE[k] = (np.clip(w, 1e-15, None), V)
        if path:
            tmp = f"{path}.{os.getpid()}"
            np.save(tmp + "_V.npy", _FACTOR_CACHE[k][1]); np.save(tmp + "_w.npy", _FACTOR_CACHE[k][0])
            os.replace(tmp + "_V.npy", path + "_V.npy"); os.replace(tmp + "_w.npy", path + "_w.npy")
    return _FACTOR_CACHE[k]


def lowdin_apply(S, X):
    """S^{1/2} X for selected columns X without forming S^{1/2}: V (sqrt(w) * (V^T X)) after one cached
    diagonalization of S."""
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    w, V = lowdin_factor(S)
    return V @ (np.sqrt(w)[:, None] * (V.T @ X))


def lowdin_sqrt(S, device="numpy"):
    """Full S^{1/2} (numpy) from the cached factor of lowdin_factor; device is accepted for the callers, the
    diagonalization itself runs on LOWDIN_DEVICE."""
    S = S.toarray() if hasattr(S, "toarray") else np.asarray(S, dtype=np.float64)
    k = _key(S)
    if k not in _SQRT_CACHE:
        _SQRT_CACHE.clear()
        w, V = lowdin_factor(S)
        # one GEMM; the old V @ diag(sqrt w) @ V.T did two n^3 products and an n x n diagonal
        _SQRT_CACHE[k] = (V * np.sqrt(w)[None, :]) @ V.T
    return _SQRT_CACHE[k]


def lowdin_columns(S, X, device="numpy"):
    """S^{1/2} X for the given columns, without forming S^{1/2} (lowdin_apply)."""
    return lowdin_apply(S, X)


def transform_mos(C, S, device="numpy"):
    return lowdin_apply(S, C)


def build_lowdin_transition_charges_flat(C_occ_act, C_virt_act, S, atom_ao_ranges, valid_i, valid_a, device="numpy"):
    """
    Computes transition charges using Löwdin symmetric orthogonalization:
    C^L = S^{1/2} * C.
    q^L_{ia, A} = sum_{mu in A} C^L_{mu i} * C^L_{mu a}.
    """
    if is_gpu(device):
        import torch
        dev = torch.device(device) if isinstance(device, str) else device
        S_half = to_tensor(lowdin_sqrt(S), dev, dtype=torch.float64)

        C_occ_t = to_tensor(C_occ_act, dev, dtype=torch.float64)
        C_virt_t = to_tensor(C_virt_act, dev, dtype=torch.float64)
        C_occ_lowdin = S_half @ C_occ_t
        C_virt_lowdin = S_half @ C_virt_t

        dim = len(valid_i)
        n_atoms = len(atom_ao_ranges)
        q_flat_t = torch.zeros((dim, n_atoms), dtype=torch.float64, device=dev)
        vi_t = torch.as_tensor(valid_i, device=dev, dtype=torch.long)
        va_t = torch.as_tensor(valid_a, device=dev, dtype=torch.long)

        for A, (a0, a1) in enumerate(atom_ao_ranges):
            Ci_A = C_occ_lowdin[a0:a1, :][:, vi_t]
            Ca_A = C_virt_lowdin[a0:a1, :][:, va_t]
            q_flat_t[:, A] = torch.sum(Ci_A * Ca_A, dim=0)

        return to_numpy(q_flat_t)

    C_occ_lowdin = lowdin_apply(S, C_occ_act)
    C_virt_lowdin = lowdin_apply(S, C_virt_act)
    
    dim = len(valid_i)
    n_atoms = len(atom_ao_ranges)
    q_flat = np.zeros((dim, n_atoms))
    
    for A, (a0, a1) in enumerate(atom_ao_ranges):
        Ci_A = C_occ_lowdin[a0:a1, :][:, valid_i]
        Ca_A = C_virt_lowdin[a0:a1, :][:, valid_a]
        q_flat[:, A] = np.sum(Ci_A * Ca_A, axis=0)
        
    return q_flat


def build_xs_transition_densities_flat(C_occ_act, C_virt_act, S, valid_i, valid_a, device="numpy"):
    r"""
    Computes AO-resolved transition densities for XsTD-DFT:
    C^L = S^{1/2} * C
    Q_{ia}(\mu) = C^L_{\mu i} * C^L_{\mu a}
    Returns shape (dim, n_ao).
    """
    if is_gpu(device):
        import torch
        dev = torch.device(device) if isinstance(device, str) else device
        S_half = to_tensor(lowdin_sqrt(S), dev, dtype=torch.float64)

        C_occ_t = to_tensor(C_occ_act, dev, dtype=torch.float64)
        C_virt_t = to_tensor(C_virt_act, dev, dtype=torch.float64)
        C_occ_lowdin = S_half @ C_occ_t
        C_virt_lowdin = S_half @ C_virt_t

        vi_t = torch.as_tensor(valid_i, device=dev, dtype=torch.long)
        va_t = torch.as_tensor(valid_a, device=dev, dtype=torch.long)
        Q_t = (C_occ_lowdin[:, vi_t] * C_virt_lowdin[:, va_t]).T
        return to_numpy(Q_t)

    C_occ_lowdin = lowdin_apply(S, C_occ_act)
    C_virt_lowdin = lowdin_apply(S, C_virt_act)
    Q = (C_occ_lowdin[:, valid_i] * C_virt_lowdin[:, valid_a]).T
    return Q


def build_xs_state_densities(C_occ_act, C_virt_act, S, device="numpy"):
    """
    Computes AO-resolved pair densities for occupied and virtual MOs for XsTD-DFT:
    Q_occ[i, j, mu] = C^L_{mu i} * C^L_{mu j}
    Q_virt[a, b, mu] = C^L_{mu a} * C^L_{mu b}
    Q_ov[i, a, mu] = C^L_{mu i} * C^L_{mu a}
    """
    C_occ_lowdin = lowdin_apply(S, C_occ_act)
    C_virt_lowdin = lowdin_apply(S, C_virt_act)

    Q_occ = np.einsum("mi,mj->ijm", C_occ_lowdin, C_occ_lowdin, optimize=True)
    Q_virt = np.einsum("ma,mb->abm", C_virt_lowdin, C_virt_lowdin, optimize=True)
    Q_ov = np.einsum("mi,ma->iam", C_occ_lowdin, C_virt_lowdin, optimize=True)
    return Q_occ, Q_virt, Q_ov
