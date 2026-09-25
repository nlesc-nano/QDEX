"""Dielectric environment of a nanocrystal: one reaction field for QP and BSE.

The nanocrystal is a sphere of electronic permittivity ``eps_in`` embedded in
a medium ``eps_out`` (use the optical value n**2 of the solvent for vertical
excitations).  For two points inside the sphere the electrostatic Green
function is the bulk term plus a reaction (image) term (Böttcher; Brus,
J. Chem. Phys. 80, 4403 (1984))::

    W(r, r') = e^2 / (eps_in |r - r'|)
             + e^2 / eps_in * sum_l c_l (r r')^l / R^(2l+1) P_l(cos theta),
    c_l = (l + 1)(eps_in - eps_out) / (l eps_in + (l + 1) eps_out).

QDEX keeps its atomistic interior kernel (Resta-MNOK) for the first term and
adds the reaction matrix ``W^refl_AB`` built here.  The same matrix enters

* the quasiparticle energies, through the polarization self-energy of each
  orbital, ``Sigma_p = 1/2 q_p^T W^refl q_p`` (occupied levels move down,
  virtual levels up), added to the bulk GW-PBE gap opening; and
* the BSE direct kernel, ``K^d`` with ``W = W^Resta + W^refl``.

Because the charged and neutral excitations see the same operator, the
environmental part of an exciton energy reduces to
``1/2 (q_e - q_h)^T W^refl (q_e - q_h)`` and nearly cancels for overlapping
electron and hole densities (Delerue, Lannoo, Allan, PRL 90, 076803 (2003)).
"""
import numpy as np


COULOMB_EV_ANG = 14.399645  # e^2 / (4 pi eps0) in eV Angstrom


def sphere_cavity(coords_ang, buffer_ang=1.0, center=None):
    """Center and radius of the dielectric cavity enclosing every atom.

    The radius is the largest atomic distance from the center plus
    ``buffer_ang``, so that all charges lie strictly inside the sphere
    (required for convergence of the multipole series).
    """
    coords = np.asarray(coords_ang, dtype=float)
    c = coords.mean(axis=0) if center is None else np.asarray(center, dtype=float)
    r_max = float(np.max(np.linalg.norm(coords - c, axis=1))) if len(coords) else 0.0
    if buffer_ang <= 0.0:
        raise ValueError("cavity buffer must be positive so that all atoms lie inside the sphere")
    return c, r_max + float(buffer_ang)


def reaction_field_matrix(coords_ang, eps_in, eps_out, radius_ang, center=None,
                          tol_ev=1.0e-6, lmax=400, block=2048):
    """Reaction (image) part of the sphere Green function between atoms, in eV.

    Returns a symmetric ``(N, N)`` matrix.  It vanishes identically when
    ``eps_out == eps_in`` and its l = 0 term is the Born charging term
    ``e^2 (1/eps_out - 1/eps_in) / R``.
    """
    eps_in = float(eps_in)
    eps_out = float(eps_out)
    if eps_in <= 0.0 or eps_out <= 0.0:
        raise ValueError("dielectric constants must be positive")
    coords = np.asarray(coords_ang, dtype=float)
    n = len(coords)
    c = coords.mean(axis=0) if center is None else np.asarray(center, dtype=float)
    rel = coords - c
    r = np.linalg.norm(rel, axis=1)
    R = float(radius_ang)
    if n and np.max(r) >= R:
        raise ValueError(f"atom at {np.max(r):.3f} Å lies outside the cavity radius {R:.3f} Å")
    W = np.zeros((n, n))
    if n == 0 or eps_in == eps_out:
        return W

    x = r / R
    unit = np.divide(rel, r[:, None], out=np.zeros_like(rel), where=r[:, None] > 0)
    x_max2 = float(np.max(x)) ** 2
    # Number of multipoles: |P_l| <= 1 and |c_l| <= |eps_in-eps_out|/eps_in bound each term.
    pref0 = COULOMB_EV_ANG / (eps_in * R) * abs(eps_in - eps_out) / min(eps_in, eps_out)
    if x_max2 > 0.0:
        l_needed = int(np.ceil(np.log(max(tol_ev / pref0, 1e-300)) / np.log(x_max2))) + 1
    else:
        l_needed = 1
    L = min(max(l_needed, 1), lmax)

    coef = np.array([(l + 1) * (eps_in - eps_out) / (l * eps_in + (l + 1) * eps_out) for l in range(L + 1)])
    for i0 in range(0, n, block):
        i1 = min(i0 + block, n)
        cos = np.clip(unit[i0:i1] @ unit.T, -1.0, 1.0)
        t = np.outer(x[i0:i1], x)  # (r r'/R^2)
        p_prev = np.ones_like(cos)  # P_0
        p_cur = cos.copy()          # P_1
        tl = np.ones_like(t)
        acc = coef[0] * p_prev
        for l in range(1, L + 1):
            tl *= t
            if l > 1:
                p_prev, p_cur = p_cur, ((2 * l - 1) * cos * p_cur - (l - 1) * p_prev) / l
            acc += coef[l] * tl * p_cur
        W[i0:i1] = acc
    W *= COULOMB_EV_ANG / (eps_in * R)
    return 0.5 * (W + W.T)


def orbital_populations(C, S, atom_ao_ranges, orbital_indices):
    """Mulliken atomic populations q_A(p) for the requested orbitals (columns)."""
    idx = np.asarray(orbital_indices, dtype=int)
    Cp = C[:, idx]
    if hasattr(Cp, "toarray"):
        Cp = Cp.toarray()
    SC = S @ Cp
    if hasattr(SC, "toarray"):
        SC = SC.toarray()
    prod = np.real(np.conj(Cp) * SC)
    q = np.zeros((len(idx), len(atom_ao_ranges)))
    for A, (a0, a1) in enumerate(atom_ao_ranges):
        q[:, A] = prod[a0:a1].sum(axis=0)
    return q


def polarization_self_energies(q, W_refl):
    """Sigma_p = 1/2 q_p^T W^refl q_p for each row of ``q`` (eV)."""
    return 0.5 * np.einsum("pA,AB,pB->p", q, W_refl, q, optimize=True)


def environment_qp_energies(eps, C, S, atom_ao_ranges, homo_index, W_refl, bulk_gap_shift,
                            n_window=100):
    """QP energies = DFT + polarization self-energy + bulk GW opening.

    Occupied level p:  eps_p - Sigma_p ;  virtual level p: eps_p + Delta_bulk + Sigma_p.
    Levels within ``n_window`` of the gap get state-resolved Sigma_p; deeper
    levels take the value of the outermost computed level of their manifold.
    The bulk opening is placed on the virtual manifold (occupied-fixed gauge),
    as for the rigid scissor.
    """
    eps = np.asarray(eps, dtype=float)
    n_mo = len(eps)
    lumo = homo_index + 1
    occ_idx = np.arange(max(0, lumo - n_window), lumo)
    vir_idx = np.arange(lumo, min(n_mo, lumo + n_window))
    sig_occ = polarization_self_energies(orbital_populations(C, S, atom_ao_ranges, occ_idx), W_refl)
    sig_vir = polarization_self_energies(orbital_populations(C, S, atom_ao_ranges, vir_idx), W_refl)

    shift = np.zeros(n_mo)
    shift[:lumo] = -sig_occ[0]
    shift[occ_idx] = -sig_occ
    shift[lumo:] = bulk_gap_shift + sig_vir[-1]
    shift[vir_idx] = bulk_gap_shift + sig_vir
    details = {
        "sigma_pol_homo_ev": float(sig_occ[-1]),
        "sigma_pol_lumo_ev": float(sig_vir[0]),
        "sigma_pol_occ_range_ev": [float(sig_occ.min()), float(sig_occ.max())],
        "sigma_pol_virt_range_ev": [float(sig_vir.min()), float(sig_vir.max())],
        "state_resolved_window": int(n_window),
    }
    return eps + shift, details
