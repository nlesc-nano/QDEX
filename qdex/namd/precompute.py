import os
import glob
import re
import time
import numpy as np
from scipy.optimize import linear_sum_assignment

from qdex.io_utils import (
    read_xyz, parse_basis, build_shell_dicts,
    build_atom_ao_ranges, count_ao_from_shells, read_mos_dense, geometry_source
)
from qdex.integrals import compute_dipole_ao, compute_cross_overlap_ao
import libint_cpp
from qdex.hardness import estimate_gw_qp_gap, estimate_brus_qp_gap, build_resta_mnok, build_gamma
from qdex.constants import BOHR_PER_ANG, HA_TO_EV
import logging

logger = logging.getLogger(__name__)


def natural_sort_key(s):
    """Sort strings with embedded numbers naturally (frame_1, frame_2, ..., frame_10)."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]


def _trivial_crossing_permutation(S_mat, lock_above=0.5):
    """
    Hungarian assignment that relabels only trivial crossings.

    A state whose diagonal overlap is still at least lock_above stays on
    itself. Those are avoided crossings, which surface hopping propagates
    in the adiabatic basis. States with a collapsed diagonal are free to
    match the partner that actually carries their character.
    """
    cost = 1.0 - np.abs(S_mat) ** 2
    diag_mag = np.abs(np.diag(S_mat))
    locked = diag_mag >= lock_above
    if np.any(locked):
        cost = np.array(cost, copy=True)
        cost[locked, :] = 1.0e6
        idx = np.where(locked)[0]
        cost[idx, idx] = 0.0
    _rows, col_ind = linear_sum_assignment(cost)
    return col_ind


def align_phases_and_crossings(S_mat, C_next, track_crossings=True, lock_above=0.5):
    """
    Aligns global signs and tracks trivial orbital crossings.
    S_mat: overlap <phi(t) | phi(t+dt)>, shape (N, N)
    C_next: MO coefficients at t+dt, shape (N_AO, N)
    
    Returns:
      S_aligned: updated overlap matrix
      C_aligned: updated MO coefficients
      perm: permutation vector
    """
    N = S_mat.shape[0]
    perm = np.arange(N)

    if track_crossings:
        perm = _trivial_crossing_permutation(S_mat, lock_above=lock_above)
        S_mat = S_mat[:, perm]
        C_next = C_next[:, perm]

    # Phase correction: ensure diagonal Re(S_ii) >= 0
    diag = np.diagonal(S_mat)
    phases = np.ones(N, dtype=np.float64)
    for i in range(N):
        if np.real(diag[i]) < 0.0:
            phases[i] = -1.0

    C_next = C_next * phases[np.newaxis, :]
    S_mat = S_mat * phases[np.newaxis, :]

    return S_mat, C_next, perm, phases


def align_spinor_phases_and_crossings(S_mat, U_list, track_crossings=True, lock_above=0.5):
    """
    Aligns global U(1) phases and tracks spinor crossings using the Hungarian algorithm.
    S_mat: overlap <psi(t) | psi(t+dt)>, shape (N, N), complex128
    U_list: list of spinor coefficient matrices at t+dt, e.g. [U_alpha, U_beta], each shape (N_mo, N)
    
    Returns:
      S_aligned: updated overlap matrix (N, N)
      U_list_aligned: updated spinor eigenvectors
      perm: permutation vector (N,)
    """
    N = S_mat.shape[0]
    perm = np.arange(N)

    if track_crossings:
        perm = _trivial_crossing_permutation(S_mat, lock_above=lock_above)
        S_mat = S_mat[:, perm]
        for idx in range(len(U_list)):
            U_list[idx] = U_list[idx][:, perm]

    # U(1) phase alignment: ensure a well-defined diagonal is real and positive.
    # A near-zero diagonal has a meaningless argument and is left alone.
    diag = np.diagonal(S_mat)
    phases = np.ones(N, dtype=np.complex128)
    ok = np.abs(diag) > 1e-8
    phases[ok] = np.exp(-1j * np.angle(diag[ok]))

    S_mat = S_mat * phases[np.newaxis, :]
    for idx in range(len(U_list)):
        U_list[idx] = U_list[idx] * phases[np.newaxis, :]

    return S_mat, U_list, perm


DIAGONAL_MODES = ("diagonal_bse", "diagonal_sbse")
INDEPENDENT_MODES = ("independent_qp", "independent_dft")


def frame_kernels(syms, coords, kernel="resta", material=None, alpha=1.0, eps_out=2.0):
    """Direct kernel W and bare exchange interaction gamma for one frame (atom resolution, eV)."""
    import qdex.hardness as _hardness
    from qdex.hardness import build_gamma
    k = str(kernel).lower()
    coords = np.asarray(coords, dtype=float)
    gamma = build_gamma(atom_symbols=syms, coords=coords, alpha=1.0, beta=0.0, exponent=_hardness.MNOK_EXPONENT_K)
    if k == "resta":
        from qdex.hardness import build_resta_mnok
        _, w = build_resta_mnok(atom_symbols=syms, coords=coords, alpha=alpha, material_name=material, eps_out=eps_out)
    elif k in ("dim", "dipole"):
        from qdex.hardness import build_dim_mnok
        w = build_dim_mnok(syms, coords, material_name=material, alpha=alpha)[0]
    elif k in ("bse", "mnok"):
        w = build_gamma(atom_symbols=syms, coords=coords, alpha=alpha, beta=0.0)
    else:
        raise ValueError(f"NAMD precompute: kernel '{kernel}' is not supported (use resta, dim or bse).")
    return w, gamma


def _diag_exchange(C_occ_list, C_virt_list, SC_occ_list, SC_virt_list, atom_ao_ranges, gamma, chunk_elems=4e7):
    """K^x_{ia,ia} = q^{ia} . gamma . q^{ia} with Mulliken transition charges, summed over spin components.

    C_*_list hold one coefficient block per spin component (one for spin-free orbitals, alpha and beta
    for spinors); the charge of each component is added before contraction.
    """
    n_o = C_occ_list[0].shape[1]
    n_v = C_virt_list[0].shape[1]
    n_at = len(atom_ao_ranges)
    kx = np.zeros((n_o, n_v))
    step = max(1, int(chunk_elems // max(1, n_v * n_at)))
    for i0 in range(0, n_o, step):
        i1 = min(n_o, i0 + step)
        q = np.zeros((i1 - i0, n_v, n_at), dtype=np.result_type(*C_occ_list, *C_virt_list))
        for Co, Cv, SCo, SCv in zip(C_occ_list, C_virt_list, SC_occ_list, SC_virt_list):
            for A, (a0, a1) in enumerate(atom_ao_ranges):
                q[:, :, A] += 0.5 * (Co[a0:a1, i0:i1].conj().T @ SCv[a0:a1, :]
                                     + SCo[a0:a1, i0:i1].conj().T @ Cv[a0:a1, :])
        qg = q @ gamma
        kx[i0:i1] = np.real(np.einsum("iam,iam->ia", qg, q.conj()))
    return kx


def compute_frame_diagonal_bse(
    xyz_path,
    mo_path,
    basis_dict,
    n_ao,
    atom_ao_ranges,
    scissor,
    f_homo=0.5,
    f_lumo=0.5,
    w_resta=None,
    nhomos=None,
    nlumos=None,
    excitation_mode="diagonal_sbse",
    include_exchange=True,
    include_direct_eh=True,
    kernel=None,
    material=None,
    alpha=1.0,
    eps_out=2.0,
    energy_window=None,
    fixed_pair_mask=None,
    compute_dipoles=True,
    nthreads=1,
    soc=False,
    gth_file=None,
    device="numpy",
    verbose_soc=False
):
    """
    Computes single-particle QP energies, diagonal BSE exciton energies,
    and transition dipoles / oscillator strengths for a single frame.
    Supports both spin-free (RKS) and 2-component spinor (SOC) representations.

    Diagonal modes: E_ia = e_a^QP - e_i^QP + k_x K^x_ia,ia - K^d_ia,ia with k_x = 2 (spin-free
    singlet) or 1 (spinors). With ``kernel`` given, W and the bare gamma are rebuilt from this
    frame's geometry; otherwise ``w_resta`` (fixed) is used for K^d and K^x is skipped.
    """
    syms, coords = read_xyz(xyz_path)
    shells = build_shell_dicts(syms, coords, basis_dict)
    n_atoms = len(atom_ao_ranges)
    mode = str(excitation_mode).lower()
    use_kernel = mode in DIAGONAL_MODES
    gamma_bare = None
    if use_kernel and kernel is not None:
        w_resta, gamma_bare = frame_kernels(syms, coords, kernel, material, alpha, eps_out)
    do_kd = use_kernel and include_direct_eh and w_resta is not None
    do_kx = use_kernel and include_exchange and gamma_bare is not None

    C_all, eps_all, occ_all = read_mos_dense(mo_path, n_ao)
    eps_all = eps_all * HA_TO_EV
    n_occ_total = int(np.sum(occ_all > 0.5))
    n_mo_total = C_all.shape[1]

    n_occ_act = n_occ_total if nhomos is None else min(nhomos, n_occ_total)
    n_virt_act = (n_mo_total - n_occ_total) if nlumos is None else min(nlumos, n_mo_total - n_occ_total)

    occ_idx = np.arange(n_occ_total - n_occ_act, n_occ_total)
    virt_idx = np.arange(n_occ_total, n_occ_total + n_virt_act)

    if not soc:
        C_occ = C_all[:, occ_idx]
        C_virt = C_all[:, virt_idx]
        eps_occ = eps_all[occ_idx]
        eps_virt = eps_all[virt_idx]

        # Quasiparticle shift (using fixed GW scissor and split fractions)
        dft_gap = float(eps_virt[0] - eps_occ[-1])
        qp_eps_occ = eps_occ - (scissor * f_homo)
        qp_eps_virt = eps_virt + (scissor * f_lumo)
        qp_gap = float(qp_eps_virt[0] - qp_eps_occ[-1])

        # Exciton states: pairs (i, a)
        n_pairs = n_occ_act * n_virt_act
        i_indices = np.repeat(np.arange(n_occ_act), n_virt_act)
        a_indices = np.tile(np.arange(n_virt_act), n_occ_act)
        E_sp = qp_eps_virt[a_indices] - qp_eps_occ[i_indices]

        # Electron-hole interaction for diagonal BSE: direct K^d and exchange K^x
        Kd_mat = None
        Kx_mat = None
        E_diag = E_sp.copy()
        if do_kd or do_kx:
            S_ao_intra = compute_cross_overlap_ao(shells, shells, nthreads=nthreads)
            SC_occ = S_ao_intra @ C_occ
            SC_virt = S_ao_intra @ C_virt
        if do_kx:
            Kx_mat = _diag_exchange([C_occ], [C_virt], [SC_occ], [SC_virt], atom_ao_ranges, gamma_bare)
            E_diag = E_diag + 2.0 * Kx_mat[i_indices, a_indices]
        if do_kd:

            q_occ_diag = np.zeros((n_occ_act, n_atoms), dtype=np.float64)
            q_virt_diag = np.zeros((n_virt_act, n_atoms), dtype=np.float64)
            for A, (a0, a1) in enumerate(atom_ao_ranges):
                q_occ_diag[:, A] = np.sum(C_occ[a0:a1, :] * SC_occ[a0:a1, :], axis=0)
                q_virt_diag[:, A] = np.sum(C_virt[a0:a1, :] * SC_virt[a0:a1, :], axis=0)

            W_virt_diag = q_virt_diag @ w_resta.T
            Kd_mat = q_occ_diag @ W_virt_diag.T
            Kd_pairs = Kd_mat[i_indices, a_indices]
            E_diag = E_diag - Kd_pairs

        if compute_dipoles:
            mu_ao_x, mu_ao_y, mu_ao_z = compute_dipole_ao(shells, nthreads=nthreads)
            mu_mo_x = C_occ.T @ (mu_ao_x @ C_virt)
            mu_mo_y = C_occ.T @ (mu_ao_y @ C_virt)
            mu_mo_z = C_occ.T @ (mu_ao_z @ C_virt)
            mu_sq_grid = mu_mo_x**2 + mu_mo_y**2 + mu_mo_z**2
            mu_sq = mu_sq_grid[i_indices, a_indices]
            # Closed-shell singlet: 4/3, energy in Hartree, dipole in ea0.
            f_osc = (4.0 / 3.0) * (E_diag / HA_TO_EV) * mu_sq
        else:
            mu_sq_grid = None
            f_osc = np.zeros_like(E_diag)

        # Active energy window filtering: lock mask if provided
        if fixed_pair_mask is not None:
            pair_mask = fixed_pair_mask
        elif energy_window is not None:
            e_min, e_max = energy_window
            mask = (E_diag >= e_min) & (E_diag <= e_max)
            pair_mask = np.where(mask)[0]
        else:
            pair_mask = np.arange(n_pairs)

        E_pairs = E_diag[pair_mask]
        f_pairs = f_osc[pair_mask]
        i_pairs = i_indices[pair_mask]
        a_pairs = a_indices[pair_mask]

        return {
            "shells": shells,
            "n_ao": n_ao,
            "C_occ": C_occ,
            "C_virt": C_virt,
            "eps_occ": qp_eps_occ,
            "eps_virt": qp_eps_virt,
            "E_pairs": E_pairs,
            "f_pairs": f_pairs,
            "i_pairs": i_pairs,
            "a_pairs": a_pairs,
            "pair_mask": pair_mask,
            "mu_sq_grid": mu_sq_grid,
            "Kd_mat": Kd_mat,
            "dft_gap": dft_gap,
            "qp_gap": qp_gap,
            "lowest_exc": float(np.min(E_pairs)) if len(E_pairs) > 0 else qp_gap,
            "scissor": scissor,
        }

    else:
        # SOC 2-component spinor calculation
        from qdex.soc_utils import compute_spinor_subspace

        act_idx = np.concatenate([occ_idx, virt_idx])
        C_act = C_all[:, act_idx]
        N_mo = len(act_idx)
        S_ao_intra = compute_cross_overlap_ao(shells, shells, nthreads=nthreads)

        soc_E, soc_U, _ = compute_spinor_subspace(
            atom_symbols=syms,
            coords_ang=coords,
            shells=shells,
            C_AO=C_all,
            eps_Ha=eps_all / HA_TO_EV,
            S_AO=S_ao_intra,
            active_indices=act_idx,
            gth_file=gth_file,
            nthreads=nthreads,
            assume_orthonormal=True,
            device=device,
            verbose=verbose_soc
        )
        soc_E = soc_E * HA_TO_EV

        n_occ_sp = 2 * n_occ_act
        n_virt_sp = 2 * n_virt_act
        eps_occ_sp = soc_E[:n_occ_sp]
        eps_virt_sp = soc_E[n_occ_sp:n_occ_sp + n_virt_sp]

        dft_gap = float(eps_virt_sp[0] - eps_occ_sp[-1])
        qp_eps_occ = eps_occ_sp - (scissor * f_homo)
        qp_eps_virt = eps_virt_sp + (scissor * f_lumo)
        qp_gap = float(qp_eps_virt[0] - qp_eps_occ[-1])

        n_pairs = n_occ_sp * n_virt_sp
        i_indices = np.repeat(np.arange(n_occ_sp), n_virt_sp)
        a_indices = np.tile(np.arange(n_virt_sp), n_occ_sp)
        E_sp = qp_eps_virt[a_indices] - qp_eps_occ[i_indices]

        U_alpha = soc_U[:N_mo, :]
        U_beta  = soc_U[N_mo:, :]
        U_occ_alpha = U_alpha[:, :n_occ_sp]
        U_occ_beta  = U_beta[:, :n_occ_sp]
        U_virt_alpha = U_alpha[:, n_occ_sp:n_occ_sp + n_virt_sp]
        U_virt_beta  = U_beta[:, n_occ_sp:n_occ_sp + n_virt_sp]

        Kd_mat = None
        Kx_mat = None
        E_diag = E_sp.copy()
        if do_kd or do_kx:
            SC_act = S_ao_intra @ C_act
            if do_kx:
                C_sp_a = C_act @ U_alpha
                C_sp_b = C_act @ U_beta
                SC_sp_a = SC_act @ U_alpha
                SC_sp_b = SC_act @ U_beta

                C_occ_sp_a = C_sp_a[:, :n_occ_sp]
                C_virt_sp_a = C_sp_a[:, n_occ_sp:n_occ_sp + n_virt_sp]
                C_occ_sp_b = C_sp_b[:, :n_occ_sp]
                C_virt_sp_b = C_sp_b[:, n_occ_sp:n_occ_sp + n_virt_sp]

                SC_occ_sp_a = SC_sp_a[:, :n_occ_sp]
                SC_virt_sp_a = SC_sp_a[:, n_occ_sp:n_occ_sp + n_virt_sp]
                SC_occ_sp_b = SC_sp_b[:, :n_occ_sp]
                SC_virt_sp_b = SC_sp_b[:, n_occ_sp:n_occ_sp + n_virt_sp]
            else:
                # Optimized low-memory path when exchange is disabled:
                # Transform occupied and virtual subspaces separately to compute diagonal densities.
                C_occ_sp_a = C_act @ U_occ_alpha
                SC_occ_sp_a = SC_act @ U_occ_alpha
                C_occ_sp_b = C_act @ U_occ_beta
                SC_occ_sp_b = SC_act @ U_occ_beta

                dens_occ = np.real(C_occ_sp_a.conj() * SC_occ_sp_a + C_occ_sp_b.conj() * SC_occ_sp_b)
                del C_occ_sp_a, SC_occ_sp_a, C_occ_sp_b, SC_occ_sp_b

                C_virt_sp_a = C_act @ U_virt_alpha
                SC_virt_sp_a = SC_act @ U_virt_alpha
                C_virt_sp_b = C_act @ U_virt_beta
                SC_virt_sp_b = SC_act @ U_virt_beta

                dens_virt = np.real(C_virt_sp_a.conj() * SC_virt_sp_a + C_virt_sp_b.conj() * SC_virt_sp_b)
                del C_virt_sp_a, SC_virt_sp_a, C_virt_sp_b, SC_virt_sp_b

        if do_kx:
            # Spinor transition densities: exchange enters once (K^x - K^d).
            Kx_mat = _diag_exchange([C_occ_sp_a, C_occ_sp_b], [C_virt_sp_a, C_virt_sp_b],
                                    [SC_occ_sp_a, SC_occ_sp_b], [SC_virt_sp_a, SC_virt_sp_b],
                                    atom_ao_ranges, gamma_bare)
            E_diag = E_diag + Kx_mat[i_indices, a_indices]
        if do_kd:
            if do_kx:
                dens_occ = np.real(C_occ_sp_a.conj() * SC_occ_sp_a + C_occ_sp_b.conj() * SC_occ_sp_b)
                dens_virt = np.real(C_virt_sp_a.conj() * SC_virt_sp_a + C_virt_sp_b.conj() * SC_virt_sp_b)

            q_occ_diag = np.empty((n_occ_sp, n_atoms), dtype=np.float64)
            q_virt_diag = np.empty((n_virt_sp, n_atoms), dtype=np.float64)
            for A, (a0, a1) in enumerate(atom_ao_ranges):
                q_occ_diag[:, A] = np.sum(dens_occ[a0:a1, :], axis=0)
                q_virt_diag[:, A] = np.sum(dens_virt[a0:a1, :], axis=0)

            W_virt_diag = q_virt_diag @ w_resta.T
            Kd_mat = q_occ_diag @ W_virt_diag.T
            Kd_pairs = Kd_mat[i_indices, a_indices]
            E_diag = E_diag - Kd_pairs

        if compute_dipoles:
            mu_ao_x, mu_ao_y, mu_ao_z = compute_dipole_ao(shells, nthreads=nthreads)
            M_x = C_act.T @ (mu_ao_x @ C_act)
            M_y = C_act.T @ (mu_ao_y @ C_act)
            M_z = C_act.T @ (mu_ao_z @ C_act)
            mu_sp_x = U_occ_alpha.conj().T @ (M_x @ U_virt_alpha) + U_occ_beta.conj().T @ (M_x @ U_virt_beta)
            mu_sp_y = U_occ_alpha.conj().T @ (M_y @ U_virt_alpha) + U_occ_beta.conj().T @ (M_y @ U_virt_beta)
            mu_sp_z = U_occ_alpha.conj().T @ (M_z @ U_virt_alpha) + U_occ_beta.conj().T @ (M_z @ U_virt_beta)
            mu_sq_grid = np.abs(mu_sp_x)**2 + np.abs(mu_sp_y)**2 + np.abs(mu_sp_z)**2
            mu_sq = mu_sq_grid[i_indices, a_indices]
            # Spinors already carry both components: 2/3, energy in Hartree.
            f_osc = (2.0 / 3.0) * (E_diag / HA_TO_EV) * mu_sq
        else:
            mu_sq_grid = None
            f_osc = np.zeros_like(E_diag)

        if fixed_pair_mask is not None:
            pair_mask = fixed_pair_mask
        elif energy_window is not None:
            e_min, e_max = energy_window
            mask = (E_diag >= e_min) & (E_diag <= e_max)
            pair_mask = np.where(mask)[0]
        else:
            pair_mask = np.arange(n_pairs)

        E_pairs = E_diag[pair_mask]
        f_pairs = f_osc[pair_mask]
        i_pairs = i_indices[pair_mask]
        a_pairs = a_indices[pair_mask]

        return {
            "shells": shells,
            "n_ao": n_ao,
            "C_act": C_act,
            "U_occ_alpha": U_occ_alpha,
            "U_occ_beta": U_occ_beta,
            "U_virt_alpha": U_virt_alpha,
            "U_virt_beta": U_virt_beta,
            "eps_occ": qp_eps_occ,
            "eps_virt": qp_eps_virt,
            "E_pairs": E_pairs,
            "f_pairs": f_pairs,
            "i_pairs": i_pairs,
            "a_pairs": a_pairs,
            "pair_mask": pair_mask,
            "mu_sq_grid": mu_sq_grid,
            "Kd_mat": Kd_mat,
            "dft_gap": dft_gap,
            "qp_gap": qp_gap,
            "lowest_exc": float(np.min(E_pairs)) if len(E_pairs) > 0 else qp_gap,
            "scissor": scissor,
        }


def precompute_namd_data(config):
    """
    Main driver for the NAMD Precomputation stage.
    Iterates through all frames in trajectory.dir, computes AO cross-overlaps,
    contracts to active MOs, tracks phases and crossings, and saves compact NPZ files.
    """
    t0_all = time.time()
    namd_cfg = config.get("namd", {})
    traj_cfg = namd_cfg.get("trajectory", {})
    sys_cfg = config.get("system", {})
    from qdex.config_schema import flatten_config
    phys_cfg = flatten_config(config)
    storage_cfg = namd_cfg.get("storage", {})
    track_cfg = namd_cfg.get("tracking", {})

    traj_dir = traj_cfg.get("dir", ".")
    frame_pattern = traj_cfg.get("frame_pattern", "frame_*")
    mo_name = traj_cfg.get("mo_file", "MOs.mbse")
    # Geometry per frame: xyz_file, or the MO file itself when it is TREXIO HDF5
    xyz_name = geometry_source(traj_cfg.get("xyz_file"), mo_name) or "frame.xyz"
    dt_nuc_fs = float(traj_cfg.get("dt_nuc_fs", 2.0))

    basis_txt = sys_cfg.get("basis_txt", "BASIS_MOLOPT_UZH")
    basis_name = sys_cfg.get("basis_name", "DZVP-MOLOPT-PBE-GTH")
    material = sys_cfg.get("material", "CSPBBR3")
    nthreads = int(sys_cfg.get("nthreads", 1))

    nhomos = phys_cfg.get("nhomos", None)
    nlumos = phys_cfg.get("nlumos", None)
    from qdex.hardness import set_mnok_options, set_bulk_vertex
    set_bulk_vertex(phys_cfg.get("bulk_vertex", "none"), phys_cfg.get("bulk_vertex_factor", 0.8))
    set_mnok_options(phys_cfg.get("mnok_exponent", 2.0), phys_cfg.get("mnok_exponent_exchange"),
                     phys_cfg.get("mnok_onsite", "ip_ea"))

    # Defaults: sBSE on the diagonal (PBE + bulk GW correction, bulk Resta W, K^x and K^d).
    qp_model = str(phys_cfg.get("qp_gap", "bulk"))
    eps_out = float(phys_cfg.get("eps_out", 2.0))
    excitation_mode = str(phys_cfg.get("excitation_mode", "diagonal_sbse")).lower()
    include_exchange = bool(phys_cfg.get("include_exchange", True))                     # K^x
    include_direct_eh = bool(phys_cfg.get("include_direct_eh", phys_cfg.get("exchange", True)))  # K^d
    kernel = phys_cfg.get("kernel", None) or "resta"
    alpha = float(phys_cfg.get("alpha", 1.0))

    qp_key = qp_model.lower()
    if excitation_mode not in DIAGONAL_MODES + INDEPENDENT_MODES:
        raise ValueError(f"NAMD precompute supports diagonal_sbse, diagonal_bse, independent_qp and "
                         f"independent_dft, not '{excitation_mode}'.")
    if qp_key not in ("bulk", "none", "pbe", "dft", "brus", "gw"):
        try:
            float(qp_model)
        except ValueError:
            raise ValueError(f"NAMD precompute: QP model '{qp_model}' is not supported (it needs a QP step "
                             "per frame). Use bulk (default), none, brus, gw or a gap in eV.") from None
    if qp_key == "gw" and excitation_mode in DIAGONAL_MODES and str(kernel).lower() != "resta-sphere":
        logger.warning("  [NAMD Warning] qp_gap 'gw' contains the surface polarization, but the kernel "
              f"'{kernel}' has no matching electron-hole image: the excitons are too high. "
              "Use qp_gap: bulk with diagonal_sbse for a consistent treatment.")
    if qp_key == "gw" and str(kernel).lower() == "resta-sphere":
        raise ValueError("NAMD precompute: the resta-sphere kernel is not available; use qp_gap: bulk.")
    if excitation_mode == "independent_dft":
        qp_model, qp_key = "none", "none"

    precompute_dir = storage_cfg.get("precompute_dir", "namd_precomputed")
    energy_window = storage_cfg.get("active_energy_window_ev", None)
    phase_correction = track_cfg.get("phase_correction", True)
    hungarian_tracking = track_cfg.get("hungarian_tracking", True)
    completeness_thresh = float(track_cfg.get("completeness_threshold", 0.99))

    soc = bool(phys_cfg.get("soc_flag", False) or phys_cfg.get("soc") is True or namd_cfg.get("soc", False))
    gth_file = sys_cfg.get("gth_file", None)
    if soc:
        if not gth_file:
            raise ValueError("system.gth_file must be specified in the configuration when SOC is enabled.")
        if not os.path.isabs(gth_file) and not os.path.exists(gth_file):
            candidate = os.path.join(traj_dir, gth_file)
            if os.path.exists(candidate):
                gth_file = candidate

    os.makedirs(precompute_dir, exist_ok=True)

    # Discover and sort frames
    pattern_path = os.path.join(traj_dir, frame_pattern)
    frame_dirs = sorted(glob.glob(pattern_path), key=natural_sort_key)
    # Filter directories only
    frame_dirs = [d for d in frame_dirs if os.path.isdir(d) and os.path.exists(os.path.join(d, xyz_name))]

    start_idx = traj_cfg.get("start_frame", 1) - 1
    end_idx = traj_cfg.get("end_frame", len(frame_dirs))
    frame_dirs = frame_dirs[start_idx:end_idx]
    n_frames = len(frame_dirs)

    if n_frames < 2:
        raise ValueError(f"Need at least 2 frames for NAMD precomputation, found {n_frames} in {pattern_path}")

    # One-time setup on frame 0: Basis, Geometry, GW Scissor, Resta Matrix
    first_xyz = os.path.join(frame_dirs[0], xyz_name)
    syms0, coords0 = read_xyz(first_xyz)
    try:
        from qdex.cluster_size import cluster_size, format_cluster_size
        size0 = cluster_size(np.array(coords0), syms0, material, sys_cfg.get("inorganic_elements"))
        logger.info(format_cluster_size(size0, None if qp_key in ("bulk", "none", "pbe", "dft") else
                                  f"core hull radius {size0['hull_radius_ang']:.3f} A ('{qp_model}')") + "  [frame 0]")
    except Exception as exc:
        logger.info(f"  [Size] Could not evaluate the cluster size: {exc}")
    basis_dict = parse_basis(basis_txt, basis_name, required_elements=set(syms0))
    shells0 = build_shell_dicts(syms0, coords0, basis_dict)
    n_ao = count_ao_from_shells(shells0)
    atom_ao_ranges = build_atom_ao_ranges(shells0)

    # Quasiparticle shift (computed once on frame 0; constant across trajectory)
    scissor = 0.0
    cluster_radius = None
    if str(qp_model).lower() == "gw":
        from qdex.hardness import estimate_gw_qp_gap
        res = estimate_gw_qp_gap(
            np.array(coords0), syms0, material, eps_out, return_details=True
        )
        if res is not None:
            scissor = float(res[0])
            cluster_radius = res[1].get("cluster_radius_ang", None)
        else:
            raise ValueError(f"GW QP estimation failed for material '{material}'.")
    elif str(qp_model).lower() == "brus":
        from qdex.hardness import estimate_brus_qp_gap
        C0, eps0, occ0 = read_mos_dense(os.path.join(frame_dirs[0], mo_name), n_ao)
        eps0 = eps0 * HA_TO_EV
        n_occ_tot = int(np.sum(occ0 > 0.5))
        dft_gap0 = float(eps0[n_occ_tot] - eps0[n_occ_tot - 1])
        res = estimate_brus_qp_gap(material, np.array(coords0), syms0)
        if res is not None:
            scissor = float(res) - dft_gap0
        else:
            raise ValueError(f"Brus QP estimation failed for material '{material}'.")
    elif str(qp_model).lower() in ("pbe", "none", "dft"):
        scissor = 0.0
    elif str(qp_model).lower() == "bulk":
        from qdex.hardness import bulk_qp_shift
        dft_gap0 = None
        if phys_cfg.get("bulk_vertex", "none") == "scaled":   # the scaled correction needs the DFT gap of frame 0
            C0, eps0, occ0 = read_mos_dense(os.path.join(frame_dirs[0], mo_name), n_ao)
            eps0 = eps0 * HA_TO_EV
            n_occ_tot = int(np.sum(occ0 > 0.5))
            dft_gap0 = float(eps0[n_occ_tot] - eps0[n_occ_tot - 1])
        scissor, _ = bulk_qp_shift(material, dft_gap0)
    else:
        try:
            target_gap = float(qp_model)
            C0, eps0, occ0 = read_mos_dense(os.path.join(frame_dirs[0], mo_name), n_ao)
            eps0 = eps0 * HA_TO_EV
            n_occ_tot = int(np.sum(occ0 > 0.5))
            dft_gap0 = float(eps0[n_occ_tot] - eps0[n_occ_tot - 1])
            scissor = target_gap - dft_gap0
        except ValueError:
            scissor = 0.0

    # Split scissor between valence and conduction bands using anchor asymmetry
    from qdex.hardness import MATERIAL_DB
    entry = MATERIAL_DB.get(str(material).upper(), None) if material else None
    if entry is not None and len(entry) >= 14:
        pbe_h_mono, pbe_l_mono, gw_h_mono, gw_l_mono = entry[10], entry[11], entry[12], entry[13]
        delta_h = gw_h_mono - pbe_h_mono
        delta_l = gw_l_mono - pbe_l_mono
        anchor_gap_opening = delta_l - delta_h
        if anchor_gap_opening > 1.0e-12 and delta_h <= 0.0 and delta_l >= 0.0:
            f_homo = -delta_h / anchor_gap_opening
            f_lumo = delta_l / anchor_gap_opening
        else:
            f_homo = f_lumo = 0.5
    else:
        f_homo, f_lumo = 0.5, 0.5

    # The direct kernel W and the bare exchange gamma are rebuilt for every frame (frame_kernels).
    w_resta = None

    logger.info("=" * 65)
    logger.info(" QDEX - NAMD Precomputation Pipeline")
    logger.info("=" * 65)
    logger.info(f"  Trajectory directory : {traj_dir}")
    logger.info(f"  Frames to process    : {n_frames} (dt = {dt_nuc_fs:.2f} fs)")
    logger.info(f"  Precompute output    : {precompute_dir}")
    logger.info(f"  Active Space         : nhomos={nhomos}, nlumos={nlumos}")
    logger.info(f"  QP Model             : {qp_model}" + ("  (PBE orbitals + bulk GW correction)" if qp_key == "bulk" else ""))
    logger.info(f"  Excitation Framework : {excitation_mode.upper()} (kernel: {kernel}, per frame; "
          f"K^x: {include_exchange}, K^d: {include_direct_eh})")
    if energy_window:
        logger.info(f"  Active Energy Window : [{energy_window[0]:.2f}, {energy_window[1]:.2f}] eV")
    if cluster_radius:
        logger.info(f"  Nanocrystal Radius   : {cluster_radius:.3f} Å (constant across trajectory)")
    logger.info(f"  GW Scissor (Δ_GW)    : {scissor:+.4f} eV (HOMO: {-scissor*f_homo:+.4f} eV, LUMO: {+scissor*f_lumo:+.4f} eV)")
    logger.info(f"  Spin-Orbit Coupling  : soc={soc}" + (f" (GTH: {os.path.basename(gth_file)})" if soc else ""))
    logger.info(f"  Tracking             : phase_correction={phase_correction}, hungarian={hungarian_tracking}")
    logger.info("=" * 65)
    if excitation_mode in DIAGONAL_MODES:
        from qdex.hardness import format_integrals_block
        logger.info(format_integrals_block("mnok", "mulliken", kernel, syms0, include_direct=include_direct_eh,
                                     include_exchange=include_exchange)
              + "\n  (diagonal elements K_ia,ia only; rebuilt from each frame's geometry)")
    logger.info("")

    prev_data = None
    fixed_pair_mask = None
    dft_gaps = []
    qp_gaps = []
    lowest_exc_energies = []

    compress_cache = storage_cfg.get("compress", False)
    save_fn = np.savez_compressed if compress_cache else np.savez

    for k, fdir in enumerate(frame_dirs):
        t0_frame = time.time()
        xyz_path = os.path.join(fdir, xyz_name)
        mo_path = os.path.join(fdir, mo_name)

        print(f"[{k+1}/{n_frames}] Processing {os.path.basename(fdir)} ...", end="", flush=True)

        curr_data = compute_frame_diagonal_bse(
            xyz_path=xyz_path,
            mo_path=mo_path,
            basis_dict=basis_dict,
            n_ao=n_ao,
            atom_ao_ranges=atom_ao_ranges,
            scissor=scissor,
            f_homo=f_homo,
            f_lumo=f_lumo,
            w_resta=w_resta,
            nhomos=nhomos,
            nlumos=nlumos,
            excitation_mode=excitation_mode,
            include_exchange=include_exchange,
            include_direct_eh=include_direct_eh,
            kernel=kernel,
            material=material,
            alpha=alpha,
            eps_out=eps_out,
            energy_window=energy_window if k == 0 else None,
            fixed_pair_mask=fixed_pair_mask,
            compute_dipoles=(k == 0),
            nthreads=nthreads,
            soc=soc,
            gth_file=gth_file,
            verbose_soc=(k == 0)
        )

        dft_gaps.append(curr_data["dft_gap"])
        qp_gaps.append(curr_data["dft_gap"] + curr_data["scissor"])
        lowest_exc_energies.append(float(np.min(curr_data["E_pairs"])))

        if k == 0:
            fixed_pair_mask = curr_data["pair_mask"]

        if prev_data is not None:
            # Compute cross-AO overlap between frame k-1 and frame k
            shells_prev = prev_data["shells"]
            shells_curr = curr_data["shells"]

            S_ao_cross = compute_cross_overlap_ao(shells_prev, shells_curr, nthreads=nthreads)

            if soc:
                M_act = prev_data["C_act"].T @ (S_ao_cross @ curr_data["C_act"])
                S_occ = (
                    prev_data["U_occ_alpha"].conj().T @ (M_act @ curr_data["U_occ_alpha"])
                    + prev_data["U_occ_beta"].conj().T @ (M_act @ curr_data["U_occ_beta"])
                )
                S_virt = (
                    prev_data["U_virt_alpha"].conj().T @ (M_act @ curr_data["U_virt_alpha"])
                    + prev_data["U_virt_beta"].conj().T @ (M_act @ curr_data["U_virt_beta"])
                )

                if phase_correction or hungarian_tracking:
                    U_occ_list = [curr_data["U_occ_alpha"], curr_data["U_occ_beta"]]
                    S_occ, U_occ_list, perm_occ = align_spinor_phases_and_crossings(
                        S_occ, U_occ_list, track_crossings=hungarian_tracking
                    )
                    curr_data["U_occ_alpha"], curr_data["U_occ_beta"] = U_occ_list

                    U_virt_list = [curr_data["U_virt_alpha"], curr_data["U_virt_beta"]]
                    S_virt, U_virt_list, perm_virt = align_spinor_phases_and_crossings(
                        S_virt, U_virt_list, track_crossings=hungarian_tracking
                    )
                    curr_data["U_virt_alpha"], curr_data["U_virt_beta"] = U_virt_list

                    # Permute single-particle QP energies to match tracked spinor basis
                    curr_data["eps_occ"] = curr_data["eps_occ"][perm_occ]
                    curr_data["eps_virt"] = curr_data["eps_virt"][perm_virt]

                    # Update pair energies in the tracked basis
                    if curr_data["Kd_mat"] is not None:
                        curr_data["Kd_mat"] = curr_data["Kd_mat"][np.ix_(perm_occ, perm_virt)]
                        Kd_pairs = curr_data["Kd_mat"][curr_data["i_pairs"], curr_data["a_pairs"]]
                        curr_data["E_pairs"] = (curr_data["eps_virt"][curr_data["a_pairs"]] - curr_data["eps_occ"][curr_data["i_pairs"]]) - Kd_pairs
                    else:
                        curr_data["E_pairs"] = curr_data["eps_virt"][curr_data["a_pairs"]] - curr_data["eps_occ"][curr_data["i_pairs"]]

                    # Update oscillator strengths in the tracked basis if computed
                    if curr_data["mu_sq_grid"] is not None:
                        mu_sq_grid = curr_data["mu_sq_grid"][np.ix_(perm_occ, perm_virt)]
                        curr_data["mu_sq_grid"] = mu_sq_grid
                        curr_data["f_pairs"] = (2.0 / 3.0) * (curr_data["E_pairs"] / HA_TO_EV) * mu_sq_grid[curr_data["i_pairs"], curr_data["a_pairs"]]

            else:
                # Contract to active MOs
                # S_occ = C_occ(t_prev)^T @ S_ao @ C_occ(t_curr)
                # S_virt = C_virt(t_prev)^T @ S_ao @ C_virt(t_curr)
                S_occ = prev_data["C_occ"].T @ (S_ao_cross @ curr_data["C_occ"])
                S_virt = prev_data["C_virt"].T @ (S_ao_cross @ curr_data["C_virt"])

                # Phase and Hungarian tracking
                if phase_correction or hungarian_tracking:
                    S_occ, curr_data["C_occ"], perm_occ, phases_occ = align_phases_and_crossings(
                        S_occ, curr_data["C_occ"], track_crossings=hungarian_tracking
                    )
                    S_virt, curr_data["C_virt"], perm_virt, phases_virt = align_phases_and_crossings(
                        S_virt, curr_data["C_virt"], track_crossings=hungarian_tracking
                    )

                    # Permute single-particle QP energies to match tracked MO basis
                    curr_data["eps_occ"] = curr_data["eps_occ"][perm_occ]
                    curr_data["eps_virt"] = curr_data["eps_virt"][perm_virt]

                    # Update pair energies in the tracked basis
                    if curr_data["Kd_mat"] is not None:
                        curr_data["Kd_mat"] = curr_data["Kd_mat"][np.ix_(perm_occ, perm_virt)]
                        Kd_pairs = curr_data["Kd_mat"][curr_data["i_pairs"], curr_data["a_pairs"]]
                        curr_data["E_pairs"] = (curr_data["eps_virt"][curr_data["a_pairs"]] - curr_data["eps_occ"][curr_data["i_pairs"]]) - Kd_pairs
                    else:
                        curr_data["E_pairs"] = curr_data["eps_virt"][curr_data["a_pairs"]] - curr_data["eps_occ"][curr_data["i_pairs"]]

                    # Update oscillator strengths in the tracked basis if computed
                    if curr_data["mu_sq_grid"] is not None:
                        mu_sq_grid = curr_data["mu_sq_grid"][np.ix_(perm_occ, perm_virt)]
                        curr_data["mu_sq_grid"] = mu_sq_grid
                        curr_data["f_pairs"] = (4.0 / 3.0) * (curr_data["E_pairs"] / HA_TO_EV) * mu_sq_grid[curr_data["i_pairs"], curr_data["a_pairs"]]

            # Completeness and tracking quality check
            min_diag_occ = np.min(np.real(np.diag(S_occ)))
            min_diag_virt = np.min(np.real(np.diag(S_virt)))
            norm_loss_occ = 1.0 - np.sum(np.abs(S_occ)**2, axis=1)
            norm_loss_virt = 1.0 - np.sum(np.abs(S_virt)**2, axis=1)
            max_loss = max(np.max(np.abs(norm_loss_occ)), np.max(np.abs(norm_loss_virt)))

            if min_diag_occ < 0.3 or min_diag_virt < 0.3:
                print(f" [Warn: Low overlap S_diag: occ={min_diag_occ:.3f}, virt={min_diag_virt:.3f}]", end="", flush=True)
            if max_loss > (1.0 - completeness_thresh):
                print(
                    f" [Warn: active-space norm loss {max_loss:.3f} exceeds "
                    f"{1.0 - completeness_thresh:.3f}]",
                    end="", flush=True,
                )

            # Save interval data for step k-1 -> k
            out_file = os.path.join(precompute_dir, f"step_{k-1:05d}_to_{k:05d}.npz")
            save_fn(
                out_file,
                time_prev_fs=(k - 1) * dt_nuc_fs,
                time_curr_fs=k * dt_nuc_fs,
                E_curr=curr_data["E_pairs"],
                eps_occ_prev=prev_data["eps_occ"],
                eps_occ_curr=curr_data["eps_occ"],
                eps_virt_prev=prev_data["eps_virt"],
                eps_virt_curr=curr_data["eps_virt"],
                S_occ=S_occ,
                S_virt=S_virt,
                max_norm_loss=max_loss
            )

        # Save single frame properties (frame 0 contains full pair basis; others are optional)
        save_all_frames = storage_cfg.get("save_all_frames", False)
        if k == 0 or save_all_frames:
            frame_out = os.path.join(precompute_dir, f"frame_{k:05d}.npz")
            save_fn(
                frame_out,
                time_fs=k * dt_nuc_fs,
                E_pairs=curr_data["E_pairs"],
                f_pairs=curr_data["f_pairs"],
                i_pairs=curr_data["i_pairs"],
                a_pairs=curr_data["a_pairs"],
                eps_occ=curr_data["eps_occ"],
                eps_virt=curr_data["eps_virt"],
                dft_gap=curr_data["dft_gap"],
                scissor=curr_data["scissor"],
                f_convention=np.array("au_hartree"),
            )

        dt_f = time.time() - t0_frame
        logger.info(f" done ({dt_f:.2f} s | {len(curr_data['E_pairs'])} active pairs)")
        prev_data = curr_data

    # Save summary metadata with pair indices
    meta_file = os.path.join(precompute_dir, "namd_metadata.npz")
    save_fn(
        meta_file,
        n_frames=n_frames,
        dt_nuc_fs=dt_nuc_fs,
        total_time_fs=(n_frames - 1) * dt_nuc_fs,
        nhomos=nhomos,
        nlumos=nlumos,
        n_occ=curr_data["eps_occ"].shape[0],
        n_virt=curr_data["eps_virt"].shape[0],
        i_pairs=curr_data["i_pairs"],
        a_pairs=curr_data["a_pairs"],
        n_active_pairs=len(fixed_pair_mask) if fixed_pair_mask is not None else len(curr_data["E_pairs"]),
        energy_window=energy_window if energy_window else np.array([]),
        mean_dft_gap=float(np.mean(dft_gaps)),
        mean_qp_gap=float(np.mean(qp_gaps)),
        mean_lowest_exc=float(np.mean(lowest_exc_energies)),
        dft_gaps=np.array(dft_gaps),
        qp_gaps=np.array(qp_gaps),
        lowest_exc_energies=np.array(lowest_exc_energies),
        scissor=scissor,
        soc=soc,
        f_convention=np.array("au_hartree"),
    )

    total_time = time.time() - t0_all
    logger.info(f"\n[NAMD Precompute] Successfully processed {n_frames} frames in {total_time:.2f} s")
    logger.info(f"[NAMD Precompute] Cached data written to: {precompute_dir}\n")


def compact_precomputed_data(precompute_dir, keep_frames=False, verbose=True):
    """
    Compacts an existing precompute directory::

      1. Converts step_*.npz to compressed format with np.savez_compressed.
      2. Removes redundant duplicate arrays (i_pairs, a_pairs, i_pairs_prev, i_pairs_curr,
         a_pairs_prev, a_pairs_curr, E_prev, f_prev, f_curr) from all step files.
      3. Preserves frame_00000.npz and namd_metadata.npz with complete pair mappings.
      4. Deletes redundant frame_00001.npz .. frame_XXXXX.npz (which are not needed
         for dynamics or analysis).

    Reduces disk footprint by 75-95% (e.g. from 22 GB down to ~3-4 GB).
    """
    import glob

    if not os.path.isdir(precompute_dir):
        raise FileNotFoundError(f"Precompute directory '{precompute_dir}' does not exist.")

    def get_dir_size_mb(path):
        total = 0
        for f in glob.glob(os.path.join(path, "*.npz")):
            total += os.path.getsize(f)
        return total / (1024.0 * 1024.0)

    size_before_mb = get_dir_size_mb(precompute_dir)
    if verbose:
        logger.info("=" * 65)
        logger.info(f" QDEX - Compacting Precomputed Data: {precompute_dir}")
        logger.info("=" * 65)
        logger.info(f"  Initial Directory Size : {size_before_mb / 1024.0:.2f} GB ({size_before_mb:.1f} MB)")

    # 1. Compact step files
    step_files = sorted(glob.glob(os.path.join(precompute_dir, "step_*.npz")))
    n_steps = len(step_files)
    if verbose:
        logger.info(f"  Compacting {n_steps} step files (removing duplicate pair indices and compressing)...")

    for k, sf in enumerate(step_files):
        d = np.load(sf)
        # Build clean dict of essential arrays
        clean_dict = {
            "time_prev_fs": d["time_prev_fs"],
            "time_curr_fs": d["time_curr_fs"],
            "E_curr": d["E_curr"],
            "eps_occ_prev": d["eps_occ_prev"],
            "eps_occ_curr": d["eps_occ_curr"],
            "eps_virt_prev": d["eps_virt_prev"],
            "eps_virt_curr": d["eps_virt_curr"],
            "S_occ": d["S_occ"],
            "S_virt": d["S_virt"],
            "max_norm_loss": d["max_norm_loss"]
        }
        # Overwrite in-place with compressed format
        np.savez_compressed(sf, **clean_dict)
        if verbose and (k + 1) % 50 == 0:
            logger.info(f"    [{k+1}/{n_steps}] steps compacted...")

    # 2. Handle frame files
    if not keep_frames:
        frame_files = sorted(glob.glob(os.path.join(precompute_dir, "frame_*.npz")))
        removed_count = 0
        for ff in frame_files:
            basename = os.path.basename(ff)
            if basename != "frame_00000.npz":
                os.remove(ff)
                removed_count += 1
        if verbose:
            logger.info(f"  Removed {removed_count} redundant frame archives (frame_00000.npz retained).")
    else:
        # Compress frame_00000.npz
        f0_path = os.path.join(precompute_dir, "frame_00000.npz")
        if os.path.exists(f0_path):
            f0 = np.load(f0_path)
            f0_dict = {k: f0[k] for k in f0.files}
            np.savez_compressed(f0_path, **f0_dict)

    # 3. Compress namd_metadata.npz
    meta_path = os.path.join(precompute_dir, "namd_metadata.npz")
    if os.path.exists(meta_path):
        m = np.load(meta_path)
        m_dict = {k: m[k] for k in m.files}
        # If i_pairs / a_pairs not in metadata, grab from frame_00000
        if "i_pairs" not in m_dict:
            f0_path = os.path.join(precompute_dir, "frame_00000.npz")
            if os.path.exists(f0_path):
                f0 = np.load(f0_path)
                m_dict["i_pairs"] = f0["i_pairs"]
                m_dict["a_pairs"] = f0["a_pairs"]
        np.savez_compressed(meta_path, **m_dict)

    size_after_mb = get_dir_size_mb(precompute_dir)
    saved_mb = size_before_mb - size_after_mb
    saved_gb = saved_mb / 1024.0
    pct = (saved_mb / max(size_before_mb, 1e-3)) * 100.0

    if verbose:
        logger.info(f"  Final Directory Size   : {size_after_mb / 1024.0:.2f} GB ({size_after_mb:.1f} MB)")
        logger.info(f"  Storage Reclaimed      : {saved_gb:.2f} GB ({pct:.1f}% reduction)")
        logger.info("=" * 65 + "\n")

    return size_before_mb, size_after_mb


def compute_trajectory_decoherence_times(
    precompute_dir,
    out_file="decoherence_times.npz",
    min_tau_fs=1.0,
    max_tau_fs=500.0,
    update_metadata=True,
    verbose=True
):
    """
    Computes state-pair pure-dephasing decoherence times from MD trajectory energy gap fluctuations
    (Prezhdo et al. JCP 111, 8366 (1999); JPCL 5, 4172 (2014); Akimov & Prezhdo, JCP 138, 124102 (2013)).

    The optical dephasing function between states i and j is::

        D_ij(t) = exp( - 1/(2 hbar^2) * var(Delta E_ij) * t^2 )

    yielding the pure-dephasing decoherence time::

        tau_ij = hbar / sigma_ij

    where ``sigma_ij^2 = var(Delta E_ij) = var(E_i) + var(E_j) - 2 * cov(E_i, E_j)``.

    Parameters
    ----------
    precompute_dir : str
        Path to directory containing precomputed NAMD frames/steps.
    out_file : str, optional
        Filename for the saved NPZ archive (default: 'decoherence_times.npz').
    min_tau_fs : float, optional
        Minimum clamp for dephasing time in fs.
    max_tau_fs : float, optional
        Maximum clamp for dephasing time in fs (applied to diagonal and near-degenerate states).
    update_metadata : bool, optional
        Whether to also store tau_occ and tau_virt directly in namd_metadata.npz.
    verbose : bool, optional
        Print summary statistics.

    Returns
    -------
    tau_occ : ndarray
        (n_occ, n_occ) hole dephasing times in fs.
    tau_virt : ndarray
        (n_virt, n_virt) electron dephasing times in fs.
    """
    HBAR_EV_FS = 0.6582119569

    if not os.path.isdir(precompute_dir):
        raise FileNotFoundError(f"Precompute directory '{precompute_dir}' does not exist.")

    step_files = sorted(glob.glob(os.path.join(precompute_dir, "step_*.npz")), key=natural_sort_key)
    f0_path = os.path.join(precompute_dir, "frame_00000.npz")
    if not os.path.exists(f0_path):
        raise FileNotFoundError(f"frame_00000.npz not found in '{precompute_dir}'.")

    f0 = np.load(f0_path)
    eps_occ_0 = f0["eps_occ"]
    eps_virt_0 = f0["eps_virt"]

    occ_history = [eps_occ_0]
    virt_history = [eps_virt_0]

    for sf in step_files:
        d = np.load(sf)
        if "eps_occ_curr" in d and "eps_virt_curr" in d:
            occ_history.append(d["eps_occ_curr"])
            virt_history.append(d["eps_virt_curr"])

    occ_arr = np.array(occ_history, dtype=np.float64)   # shape: (n_frames, n_occ)
    virt_arr = np.array(virt_history, dtype=np.float64) # shape: (n_frames, n_virt)
    n_frames = occ_arr.shape[0]

    if verbose:
        logger.info("=" * 65)
        logger.info(f" QDEX - Computing State-Pair Decoherence Times ({n_frames} frames)")
        logger.info("=" * 65)
        logger.info(f"  Target Directory : {precompute_dir}")
        logger.info(f"  Occupied States  : {occ_arr.shape[1]}")
        logger.info(f"  Virtual States   : {virt_arr.shape[1]}")

    # 1. Hole channel: cov(eps_i, eps_j)
    cov_occ = np.cov(occ_arr, rowvar=False)
    var_occ = np.diag(cov_occ)
    var_diff_occ = np.maximum(var_occ[:, None] + var_occ[None, :] - 2.0 * cov_occ, 0.0)
    sigma_occ = np.sqrt(var_diff_occ)

    tau_occ = np.zeros_like(sigma_occ)
    mask_occ = sigma_occ > 1e-6
    tau_occ[mask_occ] = HBAR_EV_FS / sigma_occ[mask_occ]
    tau_occ[~mask_occ] = max_tau_fs
    tau_occ = np.clip(tau_occ, min_tau_fs, max_tau_fs)
    np.fill_diagonal(tau_occ, max_tau_fs)

    # 2. Electron channel: cov(eps_a, eps_b)
    cov_virt = np.cov(virt_arr, rowvar=False)
    var_virt = np.diag(cov_virt)
    var_diff_virt = np.maximum(var_virt[:, None] + var_virt[None, :] - 2.0 * cov_virt, 0.0)
    sigma_virt = np.sqrt(var_diff_virt)

    tau_virt = np.zeros_like(sigma_virt)
    mask_virt = sigma_virt > 1e-6
    tau_virt[mask_virt] = HBAR_EV_FS / sigma_virt[mask_virt]
    tau_virt[~mask_virt] = max_tau_fs
    tau_virt = np.clip(tau_virt, min_tau_fs, max_tau_fs)
    np.fill_diagonal(tau_virt, max_tau_fs)

    # Save to dedicated npz
    out_path = os.path.join(precompute_dir, out_file)
    np.savez_compressed(
        out_path,
        tau_occ=tau_occ,
        tau_virt=tau_virt,
        sigma_occ=sigma_occ,
        sigma_virt=sigma_virt,
        n_frames=n_frames
    )

    if update_metadata:
        meta_path = os.path.join(precompute_dir, "namd_metadata.npz")
        if os.path.exists(meta_path):
            m = np.load(meta_path)
            m_dict = {k: m[k] for k in m.files}
            m_dict["tau_occ"] = tau_occ
            m_dict["tau_virt"] = tau_virt
            np.savez_compressed(meta_path, **m_dict)

    if verbose:
        offdiag_occ = tau_occ[~np.eye(tau_occ.shape[0], dtype=bool)]
        offdiag_virt = tau_virt[~np.eye(tau_virt.shape[0], dtype=bool)]
        logger.info(f"  Hole Dephasing   (tau_occ)  : min={np.min(offdiag_occ):.2f} fs, median={np.median(offdiag_occ):.2f} fs, max={np.max(offdiag_occ):.2f} fs")
        logger.info(f"  Electron Dephasing (tau_virt): min={np.min(offdiag_virt):.2f} fs, median={np.median(offdiag_virt):.2f} fs, max={np.max(offdiag_virt):.2f} fs")
        logger.info(f"  Cached to: {out_path}")
        logger.info("=" * 65 + "\n")

    return tau_occ, tau_virt

