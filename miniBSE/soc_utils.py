import time

import numpy as np
from scipy.linalg import eigh

import libint_cpp
from miniBSE.constants import BOHR_PER_ANG, HA_TO_EV, valence_electrons
from miniBSE.io_utils import parse_gth_soc_potentials
from miniBSE.device_utils import is_gpu

_L_CACHE = {}


def get_angular_momentum_matrices(l):
    """Returns Lx, Ly, Lz matrices in the complex spherical harmonic basis."""
    if l in _L_CACHE:
        return _L_CACHE[l]

    m = np.arange(-l, l + 1)
    Lz = np.diag(m).astype(complex)
    if l == 0:
        mats = (
            np.zeros((1, 1), dtype=complex),
            np.zeros((1, 1), dtype=complex),
            np.zeros((1, 1), dtype=complex),
        )
    else:
        Lp = np.diag(np.sqrt(l * (l + 1) - m[:-1] * (m[:-1] + 1)), 1).astype(complex)
        Lm = np.diag(np.sqrt(l * (l + 1) - m[1:] * (m[1:] - 1)), -1).astype(complex)
        Lx = 0.5 * (Lp + Lm)
        Ly = -0.5j * (Lp - Lm)
        mats = (Lx, Ly, Lz)

    _L_CACHE[l] = mats
    return mats


def _build_h_soc(k_coeffs, nprj):
    h_soc = np.zeros((nprj, nprj))
    k_idx = 0
    for i in range(nprj):
        for j in range(i, nprj):
            h_soc[i, j] = h_soc[j, i] = k_coeffs[k_idx]
            k_idx += 1
    return h_soc


def _build_soc_projectors(atom_symbols, coords_ang, soc_tbl):
    projectors = []
    proj_groups = {}
    for atom_idx, sym in enumerate(atom_symbols):
        if sym not in soc_tbl or not soc_tbl[sym]['so']:
            continue

        center_bohr = np.array(coords_ang[atom_idx]) * BOHR_PER_ANG
        for block in soc_tbl[sym]['so']:
            if not block.get('k_coeffs'):
                continue

            l = block['l']
            key = (atom_idx, l)
            if key not in proj_groups:
                proj_groups[key] = {
                    'nprj': block['nprj'],
                    'sym': sym,
                    'k_coeffs': block['k_coeffs'],
                    'l': l,
                }

            for i in range(1, block['nprj'] + 1):
                projectors.append({
                    'sym': sym,
                    'atom_idx': atom_idx,
                    'l': l,
                    'i': i,
                    'r_l': block['r'],
                    'center': center_bohr,
                })

    return projectors, proj_groups


def _ortho_active_coeffs(C_act, S_AO, SC_act=None):
    if SC_act is None:
        SC_act = S_AO @ C_act
    S_sub = C_act.conj().T @ SC_act
    S_sub = 0.5 * (S_sub + S_sub.conj().T)
    chol = np.linalg.cholesky(S_sub)
    # C_ortho = C L^{-dagger} for C^dagger S C = L L^dagger.
    return np.linalg.solve(chol.conj(), C_act.T).T


def _project_active_overlaps(
    C_AO, active_indices, B_raw, S_AO, assume_orthonormal=False, SC_AO=None,
):
    C_act = C_AO[:, active_indices]
    if assume_orthonormal:
        return C_act.conj().T @ B_raw
    SC_act = None if SC_AO is None else SC_AO[:, active_indices]
    C_ortho = _ortho_active_coeffs(C_act, S_AO, SC_act=SC_act)
    return C_ortho.conj().T @ B_raw


def _accumulate_soc_component(B_left, B_right, h_soc, L, out):
    """out += 0.5 * B_left @ kron(h_soc, L) @ B_right.conj().T."""
    K = np.kron(h_soc, L)
    BK = B_left @ K
    out += 0.5 * (BK @ B_right.conj().T)



_SPARSE_K_CACHE = {}


def _get_sparse_soc_projector_matrices(proj_groups):
    cache_key = tuple(
        (k, grp['l'], grp['nprj'], tuple(grp['k_coeffs']))
        for k, grp in sorted(proj_groups.items())
    )
    if cache_key in _SPARSE_K_CACHE:
        return _SPARSE_K_CACHE[cache_key]

    from scipy.sparse import block_diag, csr_matrix

    Kx_blocks = []
    Ky_blocks = []
    Kz_blocks = []

    for key in sorted(proj_groups.keys()):
        grp = proj_groups[key]
        l = grp['l']
        nprj = grp['nprj']
        h_soc = _build_h_soc(grp['k_coeffs'], nprj)
        Lx, Ly, Lz = get_angular_momentum_matrices(l)
        Kx_blocks.append(csr_matrix(np.kron(h_soc, Lx).real))
        Ky_blocks.append(csr_matrix(np.kron(h_soc, Ly).imag))  # Ly is purely imaginary (-0.5j*(Lp-Lm))
        Kz_blocks.append(csr_matrix(np.kron(h_soc, Lz).real))

    if not Kx_blocks:
        empty = csr_matrix((0, 0))
        mats = (empty, empty, empty)
        _SPARSE_K_CACHE[cache_key] = mats
        return mats

    Kx_sparse = block_diag(Kx_blocks, format='csr')
    Ky_sparse = block_diag(Ky_blocks, format='csr')
    Kz_sparse = block_diag(Kz_blocks, format='csr')

    mats = (Kx_sparse, Ky_sparse, Kz_sparse)
    _SPARSE_K_CACHE[cache_key] = mats
    return mats


def _assemble_rks_soc_blocks(B_mo, proj_groups, device="numpy"):
    Kx_sparse, Ky_sparse, Kz_sparse = _get_sparse_soc_projector_matrices(proj_groups)

    # Fast sparse projection: B_mo is real (n_mo, n_cols_total)
    is_real = not np.iscomplexobj(B_mo)
    if is_real:
        Mx = Kx_sparse.T.dot(B_mo.T).T
        My_imag = Ky_sparse.T.dot(B_mo.T).T
        Mz = Kz_sparse.T.dot(B_mo.T).T

        Hx = 0.5 * (Mx @ B_mo.T)
        Hy = 0.5j * (My_imag @ B_mo.T)
        Hz = 0.5 * (Mz @ B_mo.T)
    else:
        Mx = Kx_sparse.T.dot(B_mo.T).T
        My_imag = Ky_sparse.T.dot(B_mo.T).T
        Mz = Kz_sparse.T.dot(B_mo.T).T
        B_dag = B_mo.conj().T
        Hx = 0.5 * (Mx @ B_dag)
        Hy = 0.5j * (My_imag @ B_dag)
        Hz = 0.5 * (Mz @ B_dag)

    Hx = 0.5 * (Hx + Hx.conj().T)
    Hy = 0.5 * (Hy + Hy.conj().T)
    Hz = 0.5 * (Hz + Hz.conj().T)

    return Hx, Hy, Hz


def _assemble_uks_soc_chunk(keys, B_alpha, B_beta, proj_groups, col_offsets, device="numpy"):
    Mx_a = np.empty_like(B_alpha, dtype=complex)
    My_a = np.empty_like(B_alpha, dtype=complex)
    Mz_a = np.empty_like(B_alpha, dtype=complex)
    Mx_b = np.empty_like(B_beta, dtype=complex)
    My_b = np.empty_like(B_beta, dtype=complex)
    Mz_b = np.empty_like(B_beta, dtype=complex)

    for key in keys:
        grp = proj_groups[key]
        l = grp['l']
        nprj = grp['nprj']
        num_cols = nprj * (2 * l + 1)
        col_offset = col_offsets[key]
        cols = slice(col_offset, col_offset + num_cols)
        B_a = B_alpha[:, cols]
        B_b = B_beta[:, cols]
        h_soc = _build_h_soc(grp['k_coeffs'], nprj)
        Lx, Ly, Lz = get_angular_momentum_matrices(l)

        Kx = np.kron(h_soc, Lx)
        Ky = np.kron(h_soc, Ly)
        Kz = np.kron(h_soc, Lz)

        Mx_a[:, cols] = B_a @ Kx
        My_a[:, cols] = B_a @ Ky
        Mz_a[:, cols] = B_a @ Kz

        Mx_b[:, cols] = B_b @ Kx
        My_b[:, cols] = B_b @ Ky
        Mz_b[:, cols] = B_b @ Kz

    if is_gpu(device):
        import torch
        dev = torch.device(device)
        Ba_t = torch.from_numpy(B_alpha).to(device=dev, dtype=torch.complex128)
        Bb_t = torch.from_numpy(B_beta).to(device=dev, dtype=torch.complex128)
        Ba_dag_t = Ba_t.conj().T
        Bb_dag_t = Bb_t.conj().T

        Mxa_t = torch.from_numpy(Mx_a).to(device=dev, dtype=torch.complex128)
        Mya_t = torch.from_numpy(My_a).to(device=dev, dtype=torch.complex128)
        Mza_t = torch.from_numpy(Mz_a).to(device=dev, dtype=torch.complex128)
        Mxb_t = torch.from_numpy(Mx_b).to(device=dev, dtype=torch.complex128)
        Myb_t = torch.from_numpy(My_b).to(device=dev, dtype=torch.complex128)
        Mzb_t = torch.from_numpy(Mz_b).to(device=dev, dtype=torch.complex128)

        Hx_aa = (0.5 * torch.matmul(Mxa_t, Ba_dag_t)).cpu().numpy()
        Hy_aa = (0.5 * torch.matmul(Mya_t, Ba_dag_t)).cpu().numpy()
        Hz_aa = (0.5 * torch.matmul(Mza_t, Ba_dag_t)).cpu().numpy()

        Hx_ab = (0.5 * torch.matmul(Mxa_t, Bb_dag_t)).cpu().numpy()
        Hy_ab = (0.5 * torch.matmul(Mya_t, Bb_dag_t)).cpu().numpy()
        Hz_ab = (0.5 * torch.matmul(Mza_t, Bb_dag_t)).cpu().numpy()

        Hx_bb = (0.5 * torch.matmul(Mxb_t, Bb_dag_t)).cpu().numpy()
        Hy_bb = (0.5 * torch.matmul(Myb_t, Bb_dag_t)).cpu().numpy()
        Hz_bb = (0.5 * torch.matmul(Mzb_t, Bb_dag_t)).cpu().numpy()
    else:
        Ba_dag = B_alpha.conj().T
        Bb_dag = B_beta.conj().T

        Hx_aa = 0.5 * (Mx_a @ Ba_dag)
        Hy_aa = 0.5 * (My_a @ Ba_dag)
        Hz_aa = 0.5 * (Mz_a @ Ba_dag)

        Hx_ab = 0.5 * (Mx_a @ Bb_dag)
        Hy_ab = 0.5 * (My_a @ Bb_dag)
        Hz_ab = 0.5 * (Mz_a @ Bb_dag)

        Hx_bb = 0.5 * (Mx_b @ Bb_dag)
        Hy_bb = 0.5 * (My_b @ Bb_dag)
        Hz_bb = 0.5 * (Mz_b @ Bb_dag)

    return Hx_aa, Hy_aa, Hz_aa, Hx_ab, Hy_ab, Hz_ab, Hx_bb, Hy_bb, Hz_bb


def _assemble_uks_soc_blocks(B_alpha, B_beta, proj_groups, device="numpy"):
    keys = sorted(proj_groups.keys())
    col_offsets = {}
    col_offset = 0
    for key in keys:
        grp = proj_groups[key]
        col_offsets[key] = col_offset
        col_offset += grp['nprj'] * (2 * grp['l'] + 1)

    chunk = _assemble_uks_soc_chunk(keys, B_alpha, B_beta, proj_groups, col_offsets, device=device)
    Hx_aa, Hy_aa, Hz_aa, Hx_ab, Hy_ab, Hz_ab, Hx_bb, Hy_bb, Hz_bb = chunk

    Hx_aa = 0.5 * (Hx_aa + Hx_aa.conj().T)
    Hy_aa = 0.5 * (Hy_aa + Hy_aa.conj().T)
    Hz_aa = 0.5 * (Hz_aa + Hz_aa.conj().T)
    Hx_bb = 0.5 * (Hx_bb + Hx_bb.conj().T)
    Hy_bb = 0.5 * (Hy_bb + Hy_bb.conj().T)
    Hz_bb = 0.5 * (Hz_bb + Hz_bb.conj().T)
    Hx_ba = Hx_ab.conj().T
    Hy_ba = Hy_ab.conj().T
    Hz_ba = Hz_ab.conj().T
    return Hx_aa, Hy_aa, Hz_aa, Hx_ab, Hy_ab, Hz_ab, Hx_ba, Hy_ba, Hz_ba, Hx_bb, Hy_bb, Hz_bb


def prepare_soc_overlap_cache(atom_symbols, coords_ang, shells, gth_file, nthreads=1):
    """Compute and cache AO-projector overlaps shared across SOC active windows."""
    elements = {sym: valence_electrons.get(sym) for sym in set(atom_symbols)}
    soc_tbl = parse_gth_soc_potentials(gth_file, elements)
    projectors, proj_groups = _build_soc_projectors(atom_symbols, coords_ang, soc_tbl)
    B_raw = libint_cpp.compute_hgh_overlaps(shells, projectors, nthreads)
    return {
        'soc_tbl': soc_tbl,
        'projectors': projectors,
        'proj_groups': proj_groups,
        'B_raw': B_raw,
    }


def compute_spinor_subspace(
    atom_symbols, coords_ang, shells, C_AO, eps_Ha, S_AO, active_indices, gth_file,
    nthreads=1, soc_cache=None, assume_orthonormal=False, SC_AO=None, device="numpy",
    verbose=True,
):
    if verbose:
        print("\n" + "=" * 60)
        print(" [SOC] Spin-Orbit Coupling Module Initialized")
        print("=" * 60)

    if soc_cache is None:
        if verbose:
            print(f"  -> Reading GTH Potentials from: {gth_file}")
        t0 = time.time()
        soc_cache = prepare_soc_overlap_cache(
            atom_symbols, coords_ang, shells, gth_file, nthreads=nthreads
        )
        if verbose:
            print(f"  -> Parsed potentials and overlaps in {time.time() - t0:.2f}s")
    else:
        if verbose:
            print("  -> Reusing cached AO-projector overlaps")

    proj_groups = soc_cache['proj_groups']
    B_raw = soc_cache['B_raw']
    if verbose:
        print(
            f"  -> Using {len(soc_cache['projectors'])} HGH projectors "
            f"across {len(proj_groups)} angular blocks."
        )
        print(f"  -> Cached overlap matrix shape: {B_raw.shape}")
        print("  -> Projecting overlaps to Active Subspace to accelerate assembly...")
    t0 = time.time()

    if assume_orthonormal and verbose:
        print("  -> Active MOs already S-orthonormal; skipping Cholesky reorthogonalization")
    B_mo_raw = _project_active_overlaps(
        C_AO, active_indices, B_raw, S_AO,
        assume_orthonormal=assume_orthonormal, SC_AO=SC_AO,
    )
    t_proj = time.time()
    Hx_mo, Hy_mo, Hz_mo = _assemble_rks_soc_blocks(B_mo_raw, proj_groups, device=device)
    t_asm = time.time()

    if verbose:
        print(f"  -> MO projection completed in {t_proj - t0:.2f}s")
        print(f"  -> SOC block accumulation completed in {t_asm - t_proj:.2f}s")
        print(f"  -> Hamiltonian assembly completed in {t_asm - t0:.2f}s")

    n_mo = len(active_indices)
    if verbose:
        print(f"  -> Diagonalizing Single-Particle Spinor Hamiltonian (Active Space = {n_mo} MOs)...")
    t0 = time.time()

    eps_act = eps_Ha[active_indices]
    H_total = np.empty((2 * n_mo, 2 * n_mo), dtype=complex)
    H_total[:n_mo, :n_mo] = -0.5 * Hz_mo
    H_total[n_mo:, n_mo:] = 0.5 * Hz_mo
    np.fill_diagonal(H_total[:n_mo, :n_mo], eps_act - 0.5 * np.diag(Hz_mo))
    np.fill_diagonal(H_total[n_mo:, n_mo:], eps_act + 0.5 * np.diag(Hz_mo))

    Hab = -0.5 * (Hx_mo - 1j * Hy_mo)
    H_total[:n_mo, n_mo:] = Hab
    H_total[n_mo:, :n_mo] = Hab.conj().T

    soc_E, soc_U = eigh(H_total)

    if verbose:
        print(f"  -> Spinor diagonalization completed in {time.time() - t0:.2f}s")
    H0_diag = np.sort(np.repeat(eps_act, 2))
    max_shift = np.max(np.abs(soc_E - H0_diag)) * HA_TO_EV
    if verbose:
        print(f"  -> Max SOC-induced energy shift: {max_shift:.3f} eV")
        print("=" * 60 + "\n")

    return soc_E, soc_U, soc_cache


def compute_spinor_subspace_uks(
    atom_symbols, coords_ang, shells, C_alpha_AO, eps_alpha_Ha, active_alpha_indices,
    C_beta_AO, eps_beta_Ha, active_beta_indices, S_AO, gth_file, nthreads=1,
    soc_cache=None, assume_orthonormal=False, SC_alpha_AO=None, SC_beta_AO=None, device="numpy",
):
    print("\n" + "=" * 60)
    print(" [SOC-UKS] Spin-Orbit Coupling Module Initialized")
    print("=" * 60)

    if soc_cache is None:
        print(f"  -> Reading GTH Potentials from: {gth_file}")
        t0 = time.time()
        soc_cache = prepare_soc_overlap_cache(
            atom_symbols, coords_ang, shells, gth_file, nthreads=nthreads
        )
        print(f"  -> Parsed potentials and overlaps in {time.time() - t0:.2f}s")
    else:
        print("  -> Reusing cached AO-projector overlaps")

    proj_groups = soc_cache['proj_groups']
    B_raw = soc_cache['B_raw']
    print(
        f"  -> Using {len(soc_cache['projectors'])} HGH projectors "
        f"across {len(proj_groups)} angular blocks."
    )
    print(f"  -> Cached overlap matrix shape: {B_raw.shape}")

    print("  -> Projecting overlaps to UKS alpha/beta active subspaces...")
    t0 = time.time()

    if assume_orthonormal:
        print("  -> Active MOs already S-orthonormal; skipping Cholesky reorthogonalization")
    B_alpha = _project_active_overlaps(
        C_alpha_AO, active_alpha_indices, B_raw, S_AO,
        assume_orthonormal=assume_orthonormal, SC_AO=SC_alpha_AO,
    )
    B_beta = _project_active_overlaps(
        C_beta_AO, active_beta_indices, B_raw, S_AO,
        assume_orthonormal=assume_orthonormal, SC_AO=SC_beta_AO,
    )
    t_proj = time.time()
    blocks = _assemble_uks_soc_blocks(B_alpha, B_beta, proj_groups, device=device)
    t_asm = time.time()
    Hx_aa, Hy_aa, Hz_aa, Hx_ab, Hy_ab, Hz_ab, Hx_ba, Hy_ba, Hz_ba, Hx_bb, Hy_bb, Hz_bb = blocks

    print(f"  -> MO projection completed in {t_proj - t0:.2f}s")
    print(f"  -> UKS SOC block accumulation completed in {t_asm - t_proj:.2f}s")
    print(f"  -> UKS Hamiltonian assembly completed in {t_asm - t0:.2f}s")

    n_alpha = len(active_alpha_indices)
    n_beta = len(active_beta_indices)
    print(
        f"  -> Diagonalizing UKS Single-Particle Spinor Hamiltonian "
        f"(Alpha={n_alpha}, Beta={n_beta})..."
    )
    t0 = time.time()

    H0 = np.block([
        [np.diag(eps_alpha_Ha[active_alpha_indices]), np.zeros((n_alpha, n_beta))],
        [np.zeros((n_beta, n_alpha)), np.diag(eps_beta_Ha[active_beta_indices])],
    ]).astype(complex)

    H_SO = np.block([
        [0.5 * Hz_aa, 0.5 * (Hx_ab - 1j * Hy_ab)],
        [0.5 * (Hx_ba + 1j * Hy_ba), -0.5 * Hz_bb],
    ])
    H_total = H0 - H_SO
    soc_E, soc_U = eigh(H_total)

    print(f"  -> UKS spinor diagonalization completed in {time.time() - t0:.2f}s")
    H0_diag = np.sort(np.diag(H0).real)
    max_shift = np.max(np.abs(soc_E - H0_diag)) * HA_TO_EV
    print(f"  -> Max SOC-induced energy shift: {max_shift:.3f} eV")
    print("=" * 60 + "\n")

    return soc_E, soc_U, soc_cache
