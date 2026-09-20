import os
import glob
import re
import time
import numpy as np
from scipy.optimize import linear_sum_assignment

from miniBSE.io_utils import (
    read_xyz, parse_basis, build_shell_dicts,
    build_atom_ao_ranges, count_ao_from_shells, read_mos_mbse
)
from miniBSE.integrals import compute_dipole_ao, compute_cross_overlap_ao
import libint_cpp
from miniBSE.hardness import estimate_gw_qp_gap, estimate_brus_qp_gap, build_resta_mnok, build_gamma
from miniBSE.constants import BOHR_PER_ANG, HA_TO_EV


def natural_sort_key(s):
    """Sort strings with embedded numbers naturally (frame_1, frame_2, ..., frame_10)."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]


def align_phases_and_crossings(S_mat, C_next, track_crossings=True):
    """
    Aligns global signs and tracks orbital crossings using the Hungarian algorithm.
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
        # Cost matrix: 1 - |S_ij|^2 (minimizing cost maximizes overlap squared)
        cost = 1.0 - np.abs(S_mat) ** 2
        row_ind, col_ind = linear_sum_assignment(cost)
        perm = col_ind
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


def align_spinor_phases_and_crossings(S_mat, U_list, track_crossings=True):
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
        cost = 1.0 - np.abs(S_mat) ** 2
        row_ind, col_ind = linear_sum_assignment(cost)
        perm = col_ind
        S_mat = S_mat[:, perm]
        for idx in range(len(U_list)):
            U_list[idx] = U_list[idx][:, perm]

    # U(1) phase alignment: ensure diagonal S_kk is real and positive
    diag = np.diagonal(S_mat)
    phases = np.exp(-1j * np.angle(diag))

    S_mat = S_mat * phases[np.newaxis, :]
    for idx in range(len(U_list)):
        U_list[idx] = U_list[idx] * phases[np.newaxis, :]

    return S_mat, U_list, perm


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
    excitation_mode="diagonal_bse",
    include_exchange=True,
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
    """
    syms, coords = read_xyz(xyz_path)
    shells = build_shell_dicts(syms, coords, basis_dict)
    n_atoms = len(atom_ao_ranges)

    C_all, eps_all, occ_all = read_mos_mbse(mo_path, n_ao)
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

        # Electron-hole direct interaction (Kd) for diagonal BSE
        Kd_mat = None
        if excitation_mode == "diagonal_bse" and include_exchange and w_resta is not None:
            S_ao_intra = compute_cross_overlap_ao(shells, shells, nthreads=nthreads)
            SC_occ = S_ao_intra @ C_occ
            SC_virt = S_ao_intra @ C_virt

            q_occ_diag = np.zeros((n_occ_act, n_atoms), dtype=np.float64)
            q_virt_diag = np.zeros((n_virt_act, n_atoms), dtype=np.float64)
            for A, (a0, a1) in enumerate(atom_ao_ranges):
                q_occ_diag[:, A] = np.sum(C_occ[a0:a1, :] * SC_occ[a0:a1, :], axis=0)
                q_virt_diag[:, A] = np.sum(C_virt[a0:a1, :] * SC_virt[a0:a1, :], axis=0)

            W_virt_diag = q_virt_diag @ w_resta.T
            Kd_mat = q_occ_diag @ W_virt_diag.T
            Kd_pairs = Kd_mat[i_indices, a_indices]
            E_diag = E_sp - Kd_pairs
        else:
            E_diag = E_sp

        if compute_dipoles:
            mu_ao_x, mu_ao_y, mu_ao_z = compute_dipole_ao(shells, nthreads=nthreads)
            mu_mo_x = C_occ.T @ (mu_ao_x @ C_virt)
            mu_mo_y = C_occ.T @ (mu_ao_y @ C_virt)
            mu_mo_z = C_occ.T @ (mu_ao_z @ C_virt)
            mu_sq_grid = mu_mo_x**2 + mu_mo_y**2 + mu_mo_z**2
            mu_sq = mu_sq_grid[i_indices, a_indices]
            f_osc = (2.0 / 3.0) * E_diag * mu_sq
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
        from miniBSE.soc_utils import compute_spinor_subspace

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
            assume_orthonormal=False,
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
        if excitation_mode == "diagonal_bse" and include_exchange and w_resta is not None:
            SC_act = S_ao_intra @ C_act
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
            E_diag = E_sp - Kd_pairs
        else:
            E_diag = E_sp

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
            f_osc = (2.0 / 3.0) * E_diag * mu_sq
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
    phys_cfg = config.get("physics", {})
    storage_cfg = namd_cfg.get("storage", {})
    track_cfg = namd_cfg.get("tracking", {})

    traj_dir = traj_cfg.get("dir", ".")
    frame_pattern = traj_cfg.get("frame_pattern", "frame_*")
    xyz_name = traj_cfg.get("xyz_file", "frame.xyz")
    mo_name = traj_cfg.get("mo_file", "MOs.mbse")
    dt_nuc_fs = float(traj_cfg.get("dt_nuc_fs", 2.0))

    basis_txt = sys_cfg.get("basis_txt", "BASIS_MOLOPT_UZH")
    basis_name = sys_cfg.get("basis_name", "DZVP-MOLOPT-PBE-GTH")
    material = sys_cfg.get("material", "CSPBBR3")
    nthreads = int(sys_cfg.get("nthreads", 1))

    nhomos = phys_cfg.get("nhomos", None)
    nlumos = phys_cfg.get("nlumos", None)
    qp_model = phys_cfg.get("qp_gap", "gw")
    eps_out = float(phys_cfg.get("eps_out", 2.4))
    excitation_mode = phys_cfg.get("excitation_mode", "diagonal_bse").lower()
    include_exchange = phys_cfg.get("exchange", True)
    kernel = phys_cfg.get("kernel", "resta")
    alpha = float(phys_cfg.get("alpha", 1.0))

    precompute_dir = storage_cfg.get("precompute_dir", "namd_precomputed")
    energy_window = storage_cfg.get("active_energy_window_ev", None)
    phase_correction = track_cfg.get("phase_correction", True)
    hungarian_tracking = track_cfg.get("hungarian_tracking", True)
    completeness_thresh = float(track_cfg.get("completeness_threshold", 0.99))

    soc = bool(phys_cfg.get("soc", False) or namd_cfg.get("soc", False))
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
    basis_dict = parse_basis(basis_txt, basis_name, required_elements=set(syms0))
    shells0 = build_shell_dicts(syms0, coords0, basis_dict)
    n_ao = count_ao_from_shells(shells0)
    atom_ao_ranges = build_atom_ao_ranges(shells0)

    # Quasiparticle shift (computed once on frame 0; constant across trajectory)
    scissor = 0.0
    cluster_radius = None
    if str(qp_model).lower() == "gw":
        from miniBSE.hardness import estimate_gw_qp_gap
        res = estimate_gw_qp_gap(
            np.array(coords0), syms0, material, eps_out, return_details=True
        )
        if res is not None:
            scissor = float(res[0])
            cluster_radius = res[1].get("cluster_radius_ang", None)
        else:
            raise ValueError(f"GW QP estimation failed for material '{material}'.")
    elif str(qp_model).lower() == "brus":
        from miniBSE.hardness import estimate_brus_qp_gap
        C0, eps0, occ0 = read_mos_mbse(os.path.join(frame_dirs[0], mo_name), n_ao)
        eps0 = eps0 * HA_TO_EV
        n_occ_tot = int(np.sum(occ0 > 0.5))
        dft_gap0 = float(eps0[n_occ_tot] - eps0[n_occ_tot - 1])
        res = estimate_brus_qp_gap(material, np.array(coords0), syms0)
        if res is not None:
            scissor = float(res) - dft_gap0
        else:
            raise ValueError(f"Brus QP estimation failed for material '{material}'.")
    elif str(qp_model).lower() == "pbe":
        scissor = 0.0
    else:
        try:
            target_gap = float(qp_model)
            C0, eps0, occ0 = read_mos_mbse(os.path.join(frame_dirs[0], mo_name), n_ao)
            eps0 = eps0 * HA_TO_EV
            n_occ_tot = int(np.sum(occ0 > 0.5))
            dft_gap0 = float(eps0[n_occ_tot] - eps0[n_occ_tot - 1])
            scissor = target_gap - dft_gap0
        except ValueError:
            scissor = 0.0

    # Split scissor between valence and conduction bands using anchor asymmetry
    from miniBSE.hardness import MATERIAL_DB
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

    # Screened direct kernel W (computed once on frame 0)
    w_resta = None
    if excitation_mode == "diagonal_bse" and include_exchange:
        if str(kernel).lower() == "resta":
            from miniBSE.hardness import build_resta_mnok
            _, w_resta = build_resta_mnok(
                atom_symbols=syms0, coords=coords0, alpha=alpha, material_name=material, eps_out=eps_out
            )
        else:
            from miniBSE.hardness import build_gamma
            w_resta = build_gamma(atom_symbols=syms0, coords=coords0, alpha=alpha, beta=0.0)

    print("=" * 65)
    print(" miniBSE - NAMD Precomputation Pipeline")
    print("=" * 65)
    print(f"  Trajectory directory : {traj_dir}")
    print(f"  Frames to process    : {n_frames} (dt = {dt_nuc_fs:.2f} fs)")
    print(f"  Precompute output    : {precompute_dir}")
    print(f"  Active Space         : nhomos={nhomos}, nlumos={nlumos}")
    print(f"  Excitation Framework : {excitation_mode.upper()} (kernel: {kernel}, exchange: {include_exchange})")
    if energy_window:
        print(f"  Active Energy Window : [{energy_window[0]:.2f}, {energy_window[1]:.2f}] eV")
    if cluster_radius:
        print(f"  Nanocrystal Radius   : {cluster_radius:.3f} Å (constant across trajectory)")
    print(f"  GW Scissor (Δ_GW)    : {scissor:+.4f} eV (HOMO: {-scissor*f_homo:+.4f} eV, LUMO: {+scissor*f_lumo:+.4f} eV)")
    print(f"  Spin-Orbit Coupling  : soc={soc}" + (f" (GTH: {os.path.basename(gth_file)})" if soc else ""))
    print(f"  Tracking             : phase_correction={phase_correction}, hungarian={hungarian_tracking}")
    print("=" * 65 + "\n")

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
                        curr_data["f_pairs"] = (2.0 / 3.0) * curr_data["E_pairs"] * mu_sq_grid[curr_data["i_pairs"], curr_data["a_pairs"]]

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
                        curr_data["f_pairs"] = (2.0 / 3.0) * curr_data["E_pairs"] * mu_sq_grid[curr_data["i_pairs"], curr_data["a_pairs"]]

            # Completeness and tracking quality check
            min_diag_occ = np.min(np.real(np.diag(S_occ)))
            min_diag_virt = np.min(np.real(np.diag(S_virt)))
            norm_loss_occ = 1.0 - np.sum(np.abs(S_occ)**2, axis=1)
            norm_loss_virt = 1.0 - np.sum(np.abs(S_virt)**2, axis=1)
            max_loss = max(np.max(np.abs(norm_loss_occ)), np.max(np.abs(norm_loss_virt)))

            if min_diag_occ < 0.3 or min_diag_virt < 0.3:
                print(f" [Warn: Low overlap S_diag: occ={min_diag_occ:.3f}, virt={min_diag_virt:.3f}]", end="", flush=True)

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
                scissor=curr_data["scissor"]
            )

        dt_f = time.time() - t0_frame
        print(f" done ({dt_f:.2f} s | {len(curr_data['E_pairs'])} active pairs)")
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
        soc=soc
    )

    total_time = time.time() - t0_all
    print(f"\n[NAMD Precompute] Successfully processed {n_frames} frames in {total_time:.2f} s")
    print(f"[NAMD Precompute] Cached data written to: {precompute_dir}\n")


def compact_precomputed_data(precompute_dir, keep_frames=False, verbose=True):
    """
    Compacts an existing precompute directory:
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
        print("=" * 65)
        print(f" miniBSE - Compacting Precomputed Data: {precompute_dir}")
        print("=" * 65)
        print(f"  Initial Directory Size : {size_before_mb / 1024.0:.2f} GB ({size_before_mb:.1f} MB)")

    # 1. Compact step files
    step_files = sorted(glob.glob(os.path.join(precompute_dir, "step_*.npz")))
    n_steps = len(step_files)
    if verbose:
        print(f"  Compacting {n_steps} step files (removing duplicate pair indices and compressing)...")

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
            print(f"    [{k+1}/{n_steps}] steps compacted...")

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
            print(f"  Removed {removed_count} redundant frame archives (frame_00000.npz retained).")
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
        print(f"  Final Directory Size   : {size_after_mb / 1024.0:.2f} GB ({size_after_mb:.1f} MB)")
        print(f"  Storage Reclaimed      : {saved_gb:.2f} GB ({pct:.1f}% reduction)")
        print("=" * 65 + "\n")

    return size_before_mb, size_after_mb

