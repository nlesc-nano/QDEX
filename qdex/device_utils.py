import platform
import numpy as np

_TORCH_AVAILABLE = False
try:
    import torch
    _TORCH_AVAILABLE = True
except ImportError:
    torch = None


def has_torch():
    return _TORCH_AVAILABLE


def resolve_device(device_str="auto", verbose=False):
    """
    Resolves the compute device.
    Options for device_str: "auto", "cuda", "cpu", "mps", "numpy".

    Returns:
        canonical_name (str): "cuda", "cpu", or "numpy"
        torch_device (torch.device or None): Torch device object if torch is available, else None
    """
    if not _TORCH_AVAILABLE:
        if device_str not in ("auto", "numpy", None):
            if verbose:
                print("  [Device] PyTorch is not available. Falling back to NumPy CPU.")
        return "numpy", None

    device_str = str(device_str).lower() if device_str is not None else "auto"

    if device_str == "numpy":
        return "numpy", None

    if device_str == "auto":
        if torch.cuda.is_available():
            dev = torch.device("cuda")
            if verbose:
                print(f"  [Device] GPU detected: {torch.cuda.get_device_name(0)} (CUDA)")
            return "cuda", dev
        # Note: PyTorch MPS does not support float64/complex128 which is essential for BSE.
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            if verbose:
                print("  [Device] Apple MPS detected, but float64/complex128 is unsupported on MPS.")
                print("           Safely defaulting to CPU NumPy for double precision.")
            return "numpy", None
        return "cpu", torch.device("cpu")

    if device_str == "cuda":
        if not torch.cuda.is_available():
            if verbose:
                print("  [Device] Warning: CUDA requested but not available. Falling back to CPU.")
            return "cpu", torch.device("cpu")
        dev = torch.device("cuda")
        if verbose:
            print(f"  [Device] Using GPU: {torch.cuda.get_device_name(0)} (CUDA)")
        return "cuda", dev

    if device_str == "mps":
        if verbose:
            print("  [Device] Warning: MPS does not support float64/complex128.")
            print("           Falling back to CPU NumPy for numerical precision.")
        return "numpy", None

    if device_str == "cpu":
        return "cpu", torch.device("cpu")

    return "numpy", None


def to_tensor(arr, device, dtype=None):
    """Converts array-like to a torch Tensor on the specified device."""
    if not _TORCH_AVAILABLE or device in ("numpy", None):
        return np.asarray(arr)
    dev = torch.device(device) if isinstance(device, str) else device
    if isinstance(arr, torch.Tensor):
        t = arr.to(device=dev)
        if dtype is not None:
            t = t.to(dtype=dtype)
        return t
    if dtype is None:
        if np.iscomplexobj(arr):
            dtype = torch.complex128
        else:
            dtype = torch.float64
    return torch.as_tensor(arr, dtype=dtype, device=dev)


def to_numpy(t):
    """Converts torch Tensor or NumPy array to contiguous NumPy array."""
    if isinstance(t, np.ndarray):
        return t
    if _TORCH_AVAILABLE and isinstance(t, torch.Tensor):
        return t.detach().cpu().numpy()
    return np.asarray(t)


def is_gpu(device):
    """Returns True if device represents a GPU device."""
    if device is None or device == "numpy":
        return False
    if isinstance(device, str):
        return device.startswith("cuda")
    if _TORCH_AVAILABLE and isinstance(device, torch.device):
        return device.type == "cuda"
    return False
