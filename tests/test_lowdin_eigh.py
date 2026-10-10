import numpy as np
import pytest

import qdex.lowdin as lw


def _spd(n, seed=0):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return A @ A.T / n + np.eye(n)


def _check(S, w, V):
    np.testing.assert_allclose(w, np.linalg.eigvalsh(S), atol=1e-12)
    np.testing.assert_allclose(S @ V, V * w, atol=1e-12)
    np.testing.assert_allclose(V.T @ V, np.eye(len(w)), atol=1e-12)


def test_mkl_ilp64_dsyevd_matches_numpy():
    S = _spd(200)
    try:
        w, V = lw._eigh_mkl_ilp64(S)
    except (OSError, AttributeError) as e:
        pytest.skip(f"MKL with dsyevd_64 not available: {e}")
    _check(S, w, V)


def test_cusolver_xsyevd_matches_numpy():
    torch = pytest.importorskip("torch")
    if not torch.cuda.is_available():
        pytest.skip("no CUDA device")
    S = _spd(300, seed=1)
    w, V = lw._eigh_cusolver(S, "cuda")
    _check(S, w, V)


def test_gpu_request_without_gpu_falls_back_to_cpu(monkeypatch):
    def no_gpu(S, device):
        raise RuntimeError("torch has no CUDA device")
    monkeypatch.setattr(lw, "_eigh_cusolver", no_gpu)
    lw.set_lowdin_device("cuda")
    try:
        S = _spd(120, seed=2)
        w, V = lw._eigh(S)
        _check(S, w, V)
    finally:
        lw.set_lowdin_device("cpu")


def test_lowdin_columns_equal_full_sqrt():
    S = _spd(150, seed=3)
    X = np.random.default_rng(4).standard_normal((150, 7))
    lw._FACTOR_CACHE.clear(); lw._SQRT_CACHE.clear()
    np.testing.assert_allclose(lw.lowdin_columns(S, X), lw.lowdin_sqrt(S) @ X, atol=1e-12)
