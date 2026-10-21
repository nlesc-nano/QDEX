"""Weights of the spin-orbit spinors at a k-point on the spin-free bands around the gap.

    python edge_manifold_weights.py RUN_DIR BASIS_MOLOPT_UZH GTH_SOC_POTENTIALS [kx ky kz]

RUN_DIR is a qdex.bulk_soc run (real-space KS/S, TREXIO, cp2k_job.out); default k = L (0.5, 0.5, 0.5).
Prints the spin-free bands around the gap (relative to the spin-free VBM) and, for every spinor near the
gap, its weight on each of them. When qdex.bulk_soc flags edge_band_order unreliable (the SOC edges mix the
spin-free edge with other bands of its manifold), this shows whether the occupied edge spinors are of valence
or conduction character (PbTe at a_exp: inverted, see summarize_soc_gaps.py).
"""
import sys

import numpy as np

from qdex.bulk_soc import HA_EV, BulkModel, _edge_weights
from qdex.constants import valence_electrons

if __name__ == "__main__":
    run, basis, gth = sys.argv[1:4]
    kf = np.array([float(x) for x in sys.argv[4:7]]) if len(sys.argv) > 6 else np.array([0.5, 0.5, 0.5])
    model = BulkModel(run, basis, gth, "DZVP-MOLOPT-PBE-GTH", "PBE", None, None, 16.0)
    n_occ = sum(valence_electrons[e] for e in model.syms) // 2
    E = model.sf(kf) * HA_EV
    bands = list(range(n_occ - 3, n_occ + 4))
    print("spin-free bands (band - VBM : E - E_VBM):",
          "  ".join(f"{b - n_occ + 1:+d}:{E[b] - E[n_occ - 1]:+.3f}" for b in bands))
    Es, Cs = model.soc(kf, vectors=True)
    Es = Es * HA_EV
    W = np.array([_edge_weights(model, kf, [b], Cs) for b in bands])
    for s in range(2 * n_occ - 6, 2 * n_occ + 8):
        print(f"spinor {s:4d} {'occ' if s < 2 * n_occ else 'emp'} {Es[s] - E[n_occ - 1]:+.3f} eV  weights",
              " ".join(f"{w:.2f}" for w in W[:, s]))
