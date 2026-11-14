"""The torch (GPU) solver path with SOC against the NumPy path.

torch's matmul and einsum do not promote real operands to complex, so the GPU branches failed on the real
transition charges against the complex SOC rotation and on Davidson's complex vectors ("expected scalar type
Double but found ComplexDouble"). The torch branches are run here on the CPU: ``is_gpu`` is patched to accept a
torch CPU device, so no GPU is needed.
"""
from unittest import mock

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from qdex.solver import ExcitonSolver  # noqa: E402

_GPU_MODULES = ("qdex.exciton_hamiltonian", "qdex.solver", "qdex.davidson", "qdex.lowdin")


def _as_gpu(device):
    return device is not None and device != "numpy"


@pytest.fixture
def torch_cpu_as_gpu():
    patches = [mock.patch(f"{m}.is_gpu", _as_gpu) for m in _GPU_MODULES]
    for p in patches:
        p.start()
    yield "cpu"
    for p in patches:
        p.stop()


def _system(seed=7, soc=True):
    """4 atoms x 3 AOs, 6 occupied + 6 virtual MOs, a random complex unitary SOC rotation of the 4+4 active MOs."""
    rng = np.random.default_rng(seed)
    syms = ["Cd", "Se", "Cd", "Se"]
    coords = np.array([[0.0, 0.0, 0.0], [2.6, 0.0, 0.0], [0.0, 2.6, 0.0], [2.6, 2.6, 0.4]])
    n_ao = 12
    ranges = [(3 * a, 3 * a + 3) for a in range(4)]
    A = rng.normal(size=(n_ao, n_ao))
    S = np.eye(n_ao) + 0.05 * (A + A.T)
    w, V = np.linalg.eigh(S)
    S_mhalf = V @ np.diag(w ** -0.5) @ V.T
    Q, _ = np.linalg.qr(rng.normal(size=(n_ao, n_ao)))
    C = S_mhalf @ Q  # S-orthonormal MOs
    eps = np.sort(np.concatenate([rng.uniform(-9, -5, 6), rng.uniform(-3, 1, 6)]))
    occ = np.array([2.0] * 6 + [0.0] * 6)
    kw = dict(C=C, eps=eps, occ=occ, overlap=S, atom_symbols=syms, atom_coords=coords, atom_ao_ranges=ranges,
              homo_index=5, n_occ=4, n_virt=4, scissor_ev=0.5, kernel="resta", alpha=1.0, material="CDSE",
              include_exchange=True, spin="singlet")
    if soc:
        k = 8
        M = rng.normal(size=(2 * k, 2 * k)) + 1j * rng.normal(size=(2 * k, 2 * k))
        U, _ = np.linalg.qr(M)
        kw.update(soc_U=U, soc_E=np.sort(np.concatenate([rng.uniform(-9, -5, k), rng.uniform(-3, 1, k)])))
    return kw


@pytest.mark.parametrize("charge_type", ["mulliken", "lowdin"])
def test_soc_dense_torch_matches_numpy(torch_cpu_as_gpu, charge_type):
    ref = ExcitonSolver(**_system(), charge_type=charge_type, device="numpy")
    e_ref, _ = ref.solve(nroots=ref.ham.dim, full_diag=True)
    gpu = ExcitonSolver(**_system(), charge_type=charge_type, device=torch_cpu_as_gpu)
    assert np.iscomplexobj(gpu.ham.q_spinor)
    e, vecs = gpu.solve(nroots=gpu.ham.dim, full_diag=True)
    np.testing.assert_allclose(e, e_ref, atol=1e-10)
    np.testing.assert_allclose(gpu.J_mat, ref.J_mat, atol=1e-10)
    np.testing.assert_allclose(gpu.K_mat, ref.K_mat, atol=1e-10)


@pytest.mark.parametrize("soc", [True, False])
def test_davidson_torch_matches_dense(torch_cpu_as_gpu, soc):
    ref = ExcitonSolver(**_system(soc=soc), charge_type="lowdin", device="numpy")
    e_ref, _ = ref.solve(nroots=ref.ham.dim, full_diag=True)
    gpu = ExcitonSolver(**_system(soc=soc), charge_type="lowdin", device=torch_cpu_as_gpu)
    nroots = 5
    e, _ = gpu.solve(nroots=nroots, full_diag=False, tol=1e-9)
    np.testing.assert_allclose(np.sort(np.real(e))[:nroots], e_ref[:nroots], atol=1e-7)


def test_matvec_torch_matches_numpy_and_keeps_tensors_on_device(torch_cpu_as_gpu):
    ref = ExcitonSolver(**_system(), charge_type="lowdin", device="numpy")
    gpu = ExcitonSolver(**_system(), charge_type="lowdin", device=torch_cpu_as_gpu)
    rng = np.random.default_rng(3)
    x = rng.normal(size=ref.ham.dim) + 1j * rng.normal(size=ref.ham.dim)
    y = gpu.ham.matvec(torch.as_tensor(x))
    np.testing.assert_allclose(y.numpy(), ref.ham.matvec(x), atol=1e-12)
    cached = dict(gpu.ham._device_cache)
    gpu.ham.matvec(torch.as_tensor(x))
    assert all(gpu.ham._device_cache[k][2] is cached[k][2] for k in cached)  # no re-upload on the next call


def test_dense_eigh_out_of_device_memory_falls_back_to_cpu(torch_cpu_as_gpu):
    ref = ExcitonSolver(**_system(), charge_type="lowdin", device="numpy")
    e_ref, _ = ref.solve(nroots=ref.ham.dim, full_diag=True)
    gpu = ExcitonSolver(**_system(), charge_type="lowdin", device=torch_cpu_as_gpu)

    def no_memory(*a, **k):
        raise torch.OutOfMemoryError("HIP out of memory. Tried to allocate 20.19 GiB.")

    with mock.patch("torch.linalg.eigh", no_memory):
        e, _ = gpu.solve(nroots=gpu.ham.dim, full_diag=True)
    np.testing.assert_allclose(e, e_ref, atol=1e-10)
    assert isinstance(gpu.J_mat, np.ndarray) and isinstance(gpu.K_mat, np.ndarray)


@pytest.mark.parametrize("complex_charges", [True, False])
def test_truncated_exchange_blocked_torch_matches_numpy(torch_cpu_as_gpu, complex_charges):
    """The blocked (large-system) branch of the torch K assembly, on a truncated, unordered transition list."""
    from qdex.solver import _assemble_truncated_exchange
    rng = np.random.default_rng(11)
    n_occ, n_virt, n_A = 7, 5, 6
    qh = rng.normal(size=(n_occ, n_occ, n_A)) + (1j * rng.normal(size=(n_occ, n_occ, n_A)) if complex_charges else 0)
    we = rng.normal(size=(n_virt, n_virt, n_A)) + (1j * rng.normal(size=(n_virt, n_virt, n_A)) if complex_charges else 0)
    full = rng.permutation(n_occ * n_virt)[:23]
    vi, va = full // n_virt, full % n_virt
    K_ref = _assemble_truncated_exchange(qh, we, vi, va, device="numpy")
    with mock.patch("qdex.solver._K_BLOCK_BYTES", 16 * n_occ * n_virt * n_virt * 2):  # 2 holes per block
        K = _assemble_truncated_exchange(qh, we, vi, va, max_full_elements=0, device=torch_cpu_as_gpu)
    np.testing.assert_allclose(K.numpy(), K_ref, atol=1e-12)
