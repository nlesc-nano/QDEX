"""
miniBSE compatibility shim.
The miniBSE package has been renamed to QDEX (Quantum Dot Excitations & Dynamics).
All functionality is now available under the 'qdex' namespace.
"""

import sys
import warnings
import importlib
import importlib.util

warnings.warn(
    "The 'miniBSE' package has been renamed to 'qdex' (QDEX: Quantum Dot Excitations & Dynamics). "
    "Please update your imports from 'miniBSE' to 'qdex'.",
    DeprecationWarning,
    stacklevel=2,
)

import qdex

# Point sys.modules["miniBSE"] to qdex
sys.modules["miniBSE"] = qdex

# Re-export top-level attributes
for _attr in dir(qdex):
    if not _attr.startswith("__"):
        globals()[_attr] = getattr(qdex, _attr)


class _MiniBSERedirectHook:
    def find_spec(self, fullname, path, target=None):
        if fullname.startswith("miniBSE."):
            qdex_fullname = "qdex." + fullname[len("miniBSE."):]
            try:
                mod = importlib.import_module(qdex_fullname)
                sys.modules[fullname] = mod
                return importlib.util.find_spec(qdex_fullname)
            except Exception:
                pass
        return None


if not any(isinstance(h, _MiniBSERedirectHook) for h in sys.meta_path):
    sys.meta_path.insert(0, _MiniBSERedirectHook())
