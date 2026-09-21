import numpy as np
from .constants import HA_TO_EV

def compute_oscillator_strengths(
    eigvals_ev, eigvecs, mu_ia, is_spinor=False, spin_resolved=False,
):
    """
    Computes dimensionless oscillator strengths from excitation energies and transition dipoles.
    
    - is_spinor=False: Assumes spatial MOs (closed-shell singlets). Uses 4/3 prefactor.
    - is_spinor=True or spin_resolved=True: the spin channels are explicit;
      use the 2/3 prefactor.
    """
    eigvals_ha = np.asarray(eigvals_ev) / HA_TO_EV
    prefactor = 2.0 / 3.0 if (is_spinor or spin_resolved) else 4.0 / 3.0

    if hasattr(eigvecs, "order"):
        order = eigvecs.order[:len(eigvals_ha)]
        mu = mu_ia[order]
    else:
        n_roots = len(eigvals_ha)
        if eigvecs.shape[1] >= n_roots:
            mu = eigvecs[:, :n_roots].T @ mu_ia
        else:
            mu = np.zeros((n_roots, mu_ia.shape[1]), dtype=complex if np.iscomplexobj(eigvecs) or np.iscomplexobj(mu_ia) else float)
            mu[:eigvecs.shape[1]] = eigvecs.T @ mu_ia

    f = prefactor * eigvals_ha * np.sum(np.abs(mu)**2, axis=1)
    return np.real(f)
