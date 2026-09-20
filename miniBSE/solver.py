import numpy as np
import time
from miniBSE.davidson import davidson
from miniBSE.exciton_hamiltonian import ExcitonHamiltonian
from miniBSE.device_utils import is_gpu, to_tensor, to_numpy

# --- UPDATED IMPORTS ---
from miniBSE.hardness import build_gamma, build_resta_mnok


def _assemble_truncated_exchange(q_hole, w_elec, vi, va, block_size=128, max_full_elements=12_000_000, device="numpy"):
    """Build K[p,q] = sum_A q_hole[i_p,i_q,A]^* W_elec[a_p,a_q,A]."""
    n_p = len(vi)
    n_occ, _, n_atoms = q_hole.shape
    n_virt = w_elec.shape[0]
    full_dim = n_occ * n_virt

    if is_gpu(device):
        import torch
        dev = torch.device(device) if isinstance(device, str) else device
        qh_t = to_tensor(q_hole, dev)
        we_t = to_tensor(w_elec, dev)
        vi_t = torch.as_tensor(vi, device=dev, dtype=torch.long)
        va_t = torch.as_tensor(va, device=dev, dtype=torch.long)

        if full_dim * full_dim <= max_full_elements:
            K_full = torch.einsum("ijA,abA->iajb", qh_t.conj(), we_t)
            K_2d = K_full.reshape(full_dim, full_dim)
            idx = vi_t * n_virt + va_t
            return K_2d[idx[:, None], idx[None, :]]

        dtype = torch.complex128 if qh_t.is_complex() or we_t.is_complex() else torch.float64
        K = torch.empty((n_p, n_p), dtype=dtype, device=dev)
        for p0 in range(0, n_p, block_size):
            p1 = min(p0 + block_size, n_p)
            qh = qh_t[vi_t[p0:p1, None], vi_t[None, :], :].conj()
            we = we_t[va_t[p0:p1, None], va_t[None, :], :]
            K[p0:p1, :] = torch.einsum("pqA,pqA->pq", qh, we)
        return K

    dtype = np.result_type(q_hole, w_elec)

    # Dense diagonalization already needs an n_p x n_p matrix.  When the full
    # occ-virt grid is modest, assemble the cached exchange tensor with one
    # BLAS-friendly contraction and slice the truncated space from it.  This is
    # much faster than repeatedly gathering screened pairs into temporary
    # blocks, which regressed the 2500x2500 dense path.
    if full_dim * full_dim <= max_full_elements:
        K_full = np.einsum("ijA,abA->iajb", q_hole.conj(), w_elec, optimize=True)
        K_2d = K_full.reshape(full_dim, full_dim)
        idx = vi * n_virt + va
        return np.ascontiguousarray(K_2d[np.ix_(idx, idx)])

    K = np.empty((n_p, n_p), dtype=dtype)

    for p0 in range(0, n_p, block_size):
        p1 = min(p0 + block_size, n_p)
        qh = q_hole[vi[p0:p1, None], vi[None, :], :].conj()
        we = w_elec[va[p0:p1, None], va[None, :], :]
        K[p0:p1, :] = np.einsum("pqA,pqA->pq", qh, we, optimize=True)

    return K


class _TransposedDiagonalBSEVectors:
    def __init__(self, parent):
        self.parent = parent
        self.shape = (parent.n_keep, parent.dim)
        self.ndim = 2
        self.dtype = parent.dtype

    @property
    def T(self):
        return self.parent

    def conj(self):
        return self

    def conjugate(self):
        return self

    def __matmul__(self, other):
        other = np.asarray(other)
        return other[self.parent.order]


class DiagonalBSEVectors:
    """
    Memory-efficient representation of eigenvectors for diagonal BSE / independent transitions.
    In independent-transition mode, each excitation eigenvector is simply a standard basis
    vector e_{order[s]} (a 1.0 at index order[s] and 0 elsewhere).

    Instead of allocating a massive (dim x n_roots) matrix (e.g. 40,000 x 40,000 = 25.6 GB),
    this class stores only the integer permutation `order` and generates 1D slices or
    submatrices on-demand.
    """
    def __init__(self, dim, order, dtype=np.complex128):
        self.dim = int(dim)
        self.order = np.asarray(order, dtype=int)
        self.n_keep = len(self.order)
        self.shape = (self.dim, self.n_keep)
        self.ndim = 2
        self.dtype = np.dtype(dtype)

    def __len__(self):
        return self.dim

    def __getitem__(self, key):
        if isinstance(key, tuple):
            if len(key) != 2:
                raise IndexError(f"Too many indices for 2D array: {key}")
            row_idx, col_idx = key

            # Case 1: vectors[:, n] (single column -> 1D vector of length dim)
            if row_idx == slice(None) and isinstance(col_idx, (int, np.integer)):
                col = int(col_idx)
                if col < -self.n_keep or col >= self.n_keep:
                    raise IndexError(f"Index {col} is out of bounds for axis 1 with size {self.n_keep}")
                col = col if col >= 0 else self.n_keep + col
                vec = np.zeros(self.dim, dtype=self.dtype)
                vec[self.order[col]] = 1.0
                return vec

            # Case 2: vectors[:, :k] or vectors[:, slice] (slice of columns -> 2D array)
            if row_idx == slice(None) and isinstance(col_idx, slice):
                cols = np.arange(self.n_keep)[col_idx]
                k = len(cols)
                res = np.zeros((self.dim, k), dtype=self.dtype)
                res[self.order[cols], np.arange(k)] = 1.0
                return res

            # Case 3: vectors[:, list_or_array] (fancy indexing of columns -> 2D array)
            if row_idx == slice(None) and isinstance(col_idx, (list, np.ndarray)):
                cols = np.asarray(col_idx, dtype=int)
                cols = np.where(cols < 0, self.n_keep + cols, cols)
                k = len(cols)
                res = np.zeros((self.dim, k), dtype=self.dtype)
                res[self.order[cols], np.arange(k)] = 1.0
                return res

            # Case 4: single element vectors[r, c]
            if isinstance(row_idx, (int, np.integer)) and isinstance(col_idx, (int, np.integer)):
                r = int(row_idx) if row_idx >= 0 else self.dim + int(row_idx)
                c = int(col_idx) if col_idx >= 0 else self.n_keep + int(col_idx)
                return self.dtype.type(1.0 if r == self.order[c] else 0.0)

            cols = np.arange(self.n_keep)[col_idx] if isinstance(col_idx, slice) else np.atleast_1d(col_idx)
            full_sub = np.zeros((self.dim, len(cols)), dtype=self.dtype)
            full_sub[self.order[cols], np.arange(len(cols))] = 1.0
            return full_sub[row_idx]

        raise NotImplementedError("1D indexing directly on DiagonalBSEVectors is not supported; use vectors[:, n]")

    @property
    def T(self):
        return _TransposedDiagonalBSEVectors(self)

    def conj(self):
        return self

    def conjugate(self):
        return self

    def __matmul__(self, other):
        other = np.asarray(other)
        res = np.zeros((self.dim, other.shape[1]), dtype=np.result_type(self.dtype, other.dtype))
        res[self.order] = other
        return res

    def __array__(self, dtype=None):
        if self.dim * self.n_keep > 100_000_000:
            raise MemoryError(
                f"Attempting to convert DiagonalBSEVectors of shape {self.shape} to a dense "
                f"numpy array would require {self.dim * self.n_keep * self.dtype.itemsize / (1024**3):.2f} GB! "
                f"Slice columns instead, e.g. vectors[:, :n]."
            )
        arr = np.zeros((self.dim, self.n_keep), dtype=self.dtype if dtype is None else dtype)
        arr[self.order, np.arange(self.n_keep)] = 1.0
        return arr


class ExcitonSolver:
    def __init__(self, C, eps, occ, overlap, atom_symbols, atom_coords, atom_ao_ranges, 
                 homo_index, n_occ, n_virt, scissor_ev, kernel, alpha, beta=0.0, material=None, 
                 include_exchange=False, estimate_qp=False, e_thresh=None, f_thresh=0.0, 
                 mu_ia_x=None, mu_ia_y=None, mu_ia_z=None, eps_out=2.0, 
                 soc_U=None, soc_E=None, device="numpy", precomputed_sigma=None, 
                 vxc_ao_path=None, nthreads=1, spin='singlet', 
                 C_beta=None, eps_beta=None, homo_index_beta=None, charge_type='mulliken',
                  n_occ_beta=None, n_virt_beta=None, include_direct_eh=None,
                  excitation_mode="bse"):

        self.C = C
        self.overlap = overlap
        self.atom_ao_ranges = atom_ao_ranges
        self.n_occ = n_occ
        self.n_virt = n_virt
        self.homo_index = homo_index
        self.soc_flag = (soc_U is not None and soc_E is not None)
        self.spin = spin
        self.device = device
        self.C_beta = C_beta
        self.eps_beta = eps_beta
        self.homo_index_beta = homo_index_beta
        self.excitation_mode = str(excitation_mode).lower()
        # ``exchange`` historically enabled the attractive density-density term.
        # Keep it as an input alias, but use the physical name internally.
        include_direct_eh = include_exchange if include_direct_eh is None else include_direct_eh

        # --- UPDATED KERNEL LOGIC in solver.py ---
        if kernel.lower() == "resta":
            print(f"  [Solver] Using electronic Resta-MNOK direct kernel for material: {material}")
            if alpha != 1.0:
                print("  [Warning] alpha does not tune RESTA and is ignored for this kernel.")
            gamma_qp, w_resta = build_resta_mnok(
                atom_symbols=atom_symbols, coords=atom_coords,
                alpha=alpha, material_name=material, eps_out=eps_out
            )
        else:
            print(f"  [Solver] Using standard Grimme sTDA MNOK kernel.")
            print(f"           -> Screened W/BSE (alpha = {alpha:.3f})")
            g = build_gamma(atom_symbols=atom_symbols, coords=atom_coords, alpha=alpha, beta=0.0)
            gamma_qp, w_resta = g, g

        print(f"  [Solver] Building Bare Kernel V (alpha = 1.000, beta = 0.000)")
        # This MUST be beta=0.0 to preserve your baseline COH polarization!
        gamma_bare = build_gamma(atom_symbols=atom_symbols, coords=atom_coords, alpha=1.0, beta=0.0)

        if beta > 0.0:
            raise ValueError(
                "beta > 0 is disabled until the bare on-site U parameters are defined and validated"
            )
        gamma_penalty = np.zeros_like(gamma_bare)

        self.ham = ExcitonHamiltonian(
            C=C, eps=eps, overlap=overlap, atom_ao_ranges=atom_ao_ranges,
            homo_index=homo_index, n_occ=n_occ, n_virt=n_virt, scissor_ev=scissor_ev,
            gamma_qp=gamma_qp,        
            gamma_bse=w_resta,
            material=material,
            gamma_bare=gamma_bare,
            gamma_penalty=gamma_penalty,
            alpha=alpha,            
            include_exchange=include_direct_eh, estimate_qp=estimate_qp, e_thresh=e_thresh,
            f_thresh=f_thresh, mu_ia_x=mu_ia_x, mu_ia_y=mu_ia_y, mu_ia_z=mu_ia_z, 
            soc_U=soc_U, soc_E=soc_E, device=device, precomputed_sigma=precomputed_sigma,
            vxc_ao_path=vxc_ao_path, nthreads=nthreads, spin=spin,
            C_beta=C_beta, eps_beta=eps_beta, homo_index_beta=homo_index_beta,
            charge_type=charge_type,
            n_occ_beta=n_occ_beta, n_virt_beta=n_virt_beta,
            excitation_mode=excitation_mode
        )

    def solve(self, nroots=10, full_diag=False, tol=1e-5, excitation_mode="bse"):
        if self.ham.dim == 0:
            print("ERROR: Active space dimension is 0! Your energy threshold is filtering out all transitions.")
            import sys; sys.exit(1)

        mode = str(excitation_mode).lower()
        self.excitation_mode = mode
        if mode != "bse":
            energies, kx_diag, kd_diag = self.ham.independent_transition_energies(mode)
            order = np.argsort(energies)
            n_keep = self.ham.dim if full_diag else min(max(1, int(nroots)), self.ham.dim)
            order = order[:n_keep]
            vectors = DiagonalBSEVectors(
                self.ham.dim,
                order,
                dtype=np.result_type(energies, complex if self.soc_flag else float),
            )
            self.diagonal_kx = kx_diag
            self.diagonal_kd = kd_diag
            print(f"  Independent-transition approximation: {mode} (no BSE/TDA diagonalization)")
            return energies[order], vectors

        if full_diag:
            print(f"  Building dense Hamiltonian in truncated space ({self.ham.dim}x{self.ham.dim})...")
            use_gpu = is_gpu(self.device)
            dev = self.device if use_gpu else None
            
            if not self.soc_flag:
                # ==========================================================
                # SPATIAL DENSE BUILDER (Spin-Free Singlets or Triplets)
                # ==========================================================
                is_triplet = getattr(self.ham, 'spin', 'singlet') == 'triplet'
                is_uks_sp  = getattr(self.ham, 'spin', 'singlet') == 'uks_spin_preserving'

                if is_triplet:
                    print("  [Dense] Building triplet BSE/TDA Hamiltonian (Kx_bare absent)...")
                    if use_gpu:
                        import torch
                        H = torch.diag(to_tensor(self.ham.D, dev))
                        J_mat = torch.zeros((self.ham.dim, self.ham.dim), dtype=H.dtype, device=dev)
                        K_mat = torch.zeros_like(J_mat)
                        if self.ham.include_exchange:
                            print("  [Dense] Building screened direct electron-hole matrix (-Kd_screened via GPU)...")
                            t1 = time.time()
                            vi, va = self.ham.valid_i, self.ham.valid_a
                            K_truncated = _assemble_truncated_exchange(
                                self.ham.q_occ, self.ham.W_virt, vi, va, device=dev
                            )
                            H -= K_truncated
                            K_mat = K_truncated
                            print(f"    -> Triplet K built in {time.time()-t1:.2f}s")
                    else:
                        H = np.diag(self.ham.D).copy()
                        J_mat = np.zeros((self.ham.dim, self.ham.dim))
                        K_mat = np.zeros_like(J_mat)
                        if self.ham.include_exchange:
                            print("  [Dense] Building screened direct electron-hole matrix (-Kd_screened)...")
                            t1 = time.time()
                            vi, va = self.ham.valid_i, self.ham.valid_a
                            K_truncated = _assemble_truncated_exchange(
                                self.ham.q_occ, self.ham.W_virt, vi, va
                            )
                            H -= K_truncated
                            K_mat = K_truncated
                            print(f"    -> Triplet K built in {time.time()-t1:.2f}s")

                elif is_uks_sp:
                    # ----------------------------------------------------------
                    # UKS SPIN-PRESERVING (Manifold B) DENSE BUILDER
                    # ----------------------------------------------------------
                    print("  [Dense] Building UKS Spin-Preserving Hamiltonian (Manifold B)...")
                    t0 = time.time()

                    if use_gpu:
                        import torch
                        q_flat_t = to_tensor(self.ham.q_flat, dev)
                        gamma_t = to_tensor(self.ham.gamma, dev)
                        temp = q_flat_t @ gamma_t
                        J_mat = temp @ q_flat_t.T
                        H = torch.diag(to_tensor(self.ham.D, dev)) + 1.0 * J_mat
                        K_mat = torch.zeros_like(J_mat)
                        print(f"    -> Kx_bare built in {time.time()-t0:.2f}s")

                        if self.ham.include_exchange:
                            print("  [Dense] Building block-diagonal screened direct matrix (-Kd_alpha, -Kd_beta via GPU)...")
                            t1 = time.time()
                            dim_a = self.ham.dim_a
                            dim_b = self.ham.dim_b
                            vi_a, va_a = self.ham.vi_a, self.ham.va_a
                            vi_b, va_b = self.ham.vi_b, self.ham.va_b
                            n_a = len(vi_a)
                            n_b = len(vi_b)

                            K_full = torch.zeros((self.ham.dim, self.ham.dim), dtype=H.dtype, device=dev)
                            if n_a > 0:
                                K_alpha = _assemble_truncated_exchange(
                                    self.ham.q_occ_a, self.ham.W_virt_a, vi_a, va_a, device=dev
                                )
                                K_full[:dim_a, :dim_a] = K_alpha

                            if n_b > 0:
                                K_beta = _assemble_truncated_exchange(
                                    self.ham.q_occ_b, self.ham.W_virt_b, vi_b, va_b, device=dev
                                )
                                K_full[dim_a:, dim_a:] = K_beta

                            H -= K_full
                            K_mat = K_full
                            print(f"    -> Exchange built in {time.time()-t1:.2f}s")
                    else:
                        # Bare exchange/local-field Kx couples alpha and beta transitions.
                        temp = self.ham.q_flat @ self.ham.gamma
                        J_mat = temp @ self.ham.q_flat.T
                        H = np.diag(self.ham.D) + 1.0 * J_mat
                        K_mat = np.zeros_like(J_mat)
                        print(f"    -> Kx_bare built in {time.time()-t0:.2f}s")

                        if self.ham.include_exchange:
                            print("  [Dense] Building block-diagonal screened direct matrix (-Kd_alpha, -Kd_beta)...")
                            t1 = time.time()
                            dim_a = self.ham.dim_a
                            dim_b = self.ham.dim_b
                            vi_a, va_a = self.ham.vi_a, self.ham.va_a
                            vi_b, va_b = self.ham.vi_b, self.ham.va_b
                            n_a = len(vi_a)
                            n_b = len(vi_b)

                            K_full = np.zeros((self.ham.dim, self.ham.dim))
                            if n_a > 0:
                                K_alpha = _assemble_truncated_exchange(
                                    self.ham.q_occ_a, self.ham.W_virt_a, vi_a, va_a
                                )
                                K_full[:dim_a, :dim_a] = K_alpha

                            if n_b > 0:
                                K_beta = _assemble_truncated_exchange(
                                    self.ham.q_occ_b, self.ham.W_virt_b, vi_b, va_b
                                )
                                K_full[dim_a:, dim_a:] = K_beta

                            H -= K_full
                            K_mat = K_full
                            print(f"    -> Exchange built in {time.time()-t1:.2f}s")

                else:
                    print("  [Dense] Building bare singlet exchange/local-field term (2Kx_bare)...")
                    t0 = time.time()
                    if use_gpu:
                        import torch
                        q_flat_t = to_tensor(self.ham.q_flat, dev)
                        gamma_t = to_tensor(self.ham.gamma, dev)
                        temp = q_flat_t @ gamma_t
                        J_mat = 2.0 * (temp @ q_flat_t.T)
                        H = torch.diag(to_tensor(self.ham.D, dev)) + J_mat
                        K_mat = torch.zeros_like(J_mat)
                        print(f"    -> Kx_bare built in {time.time()-t0:.2f}s")

                        if self.ham.include_exchange and not is_uks_sp and not is_triplet:
                            print("  [Dense] Building screened direct electron-hole matrix (-Kd_screened via GPU)...")
                            t1 = time.time()
                            c_x = getattr(self.ham, 'c_x', 1.0)
                            vi, va = self.ham.valid_i, self.ham.valid_a
                            K_truncated = _assemble_truncated_exchange(
                                self.ham.q_occ, self.ham.W_virt, vi, va, device=dev
                            )
                            H -= c_x * K_truncated
                            K_mat = c_x * K_truncated
                            print(f"    -> Kd_screened built in {time.time()-t1:.2f}s")
                    else:
                        temp = self.ham.q_flat @ self.ham.gamma
                        J_mat = temp @ self.ham.q_flat.T
                        H = np.diag(self.ham.D) + 2.0 * J_mat
                        K_mat = np.zeros_like(J_mat)
                        print(f"    -> Kx_bare built in {time.time()-t0:.2f}s")

                        if self.ham.include_exchange and not is_uks_sp and not is_triplet:
                            print("  [Dense] Building screened direct electron-hole matrix (-Kd_screened)...")
                            t1 = time.time()
                            c_x = getattr(self.ham, 'c_x', 1.0)
                            vi, va = self.ham.valid_i, self.ham.valid_a
                            K_truncated = _assemble_truncated_exchange(
                                self.ham.q_occ, self.ham.W_virt, vi, va
                            )
                            H -= c_x * K_truncated
                            K_mat = c_x * K_truncated
                            print(f"    -> Kd_screened built in {time.time()-t1:.2f}s")
 
            else:
                # ==========================================================
                # SPINOR DENSE BUILDER (Relativistic Spin-Orbit)
                # ==========================================================
                print("  [Dense-SOC] Building spinor bare exchange/local-field term (Kx_bare)...")
                t0 = time.time()
                if use_gpu:
                    import torch
                    q_sp = to_tensor(self.ham.q_spinor, dev)
                    gamma_t = to_tensor(self.ham.gamma, dev)
                    temp = q_sp.conj() @ gamma_t
                    J_mat = temp @ q_sp.T
                    H = torch.diag(to_tensor(self.ham.D, dev, dtype=torch.complex128)) + J_mat
                    K_mat = torch.zeros_like(J_mat)
                    print(f"    -> Kx_bare built in {time.time()-t0:.2f}s")

                    if self.ham.include_exchange:
                        print("  [Dense-SOC] Building screened direct electron-hole matrix (-Kd_screened via GPU)...")
                        t1 = time.time()
                        n_occ_sp = self.ham.n_occ_spinor
                        n_virt_sp = self.ham.n_virt_spinor
                        if hasattr(self.ham, 'valid_spinor_idx'):
                            v_idx = self.ham.valid_spinor_idx
                            vi_sp = v_idx // n_virt_sp
                            va_sp = v_idx % n_virt_sp
                            W_elec = self.ham.W_elec_spinor
                            K_truncated = _assemble_truncated_exchange(
                                self.ham.q_hole_spinor, W_elec, vi_sp, va_sp, device=dev
                            )
                            H -= K_truncated
                            K_mat = K_truncated
                        else:
                            W = to_tensor(self.ham.W_elec_spinor, dev)
                            qh = to_tensor(self.ham.q_hole_spinor, dev)
                            K_full = torch.tensordot(qh.conj(), W, dims=([2], [2]))
                            K_full_trans = K_full.permute(0, 2, 1, 3)
                            K_2d = K_full_trans.reshape(n_occ_sp * n_virt_sp, n_occ_sp * n_virt_sp)
                            H -= K_2d
                            K_mat = K_2d
                        print(f"    -> Kd_screened built in {time.time()-t1:.2f}s")
                else:
                    # Notice: No factor of 2.0, and requires complex conjugate transpose
                    temp = self.ham.q_spinor.conj() @ self.ham.gamma
                    J_mat = temp @ self.ham.q_spinor.T
                    H = np.diag(self.ham.D).astype(complex) + J_mat
                    K_mat = np.zeros_like(J_mat)
                    print(f"    -> Kx_bare built in {time.time()-t0:.2f}s")

                    if self.ham.include_exchange:
                        print("  [Dense-SOC] Building screened direct electron-hole matrix (-Kd_screened)...")
                        t1 = time.time()
                        n_occ_sp = self.ham.n_occ_spinor
                        n_virt_sp = self.ham.n_virt_spinor

                        if hasattr(self.ham, 'valid_spinor_idx'):
                            v_idx = self.ham.valid_spinor_idx
                            vi_sp = v_idx // n_virt_sp
                            va_sp = v_idx % n_virt_sp
                            W_elec = self.ham.W_elec_spinor
                            K_truncated = _assemble_truncated_exchange(
                                self.ham.q_hole_spinor, W_elec, vi_sp, va_sp
                            )
                            H -= K_truncated
                            K_mat = K_truncated
                        else:
                            W = self.ham.W_elec_spinor
                            K_full = np.tensordot(self.ham.q_hole_spinor.conj(), W, axes=([2], [2]))
                            K_full_trans = np.transpose(K_full, (0, 2, 1, 3))
                            K_2d = K_full_trans.reshape(n_occ_sp * n_virt_sp, n_occ_sp * n_virt_sp)
                            H -= K_2d
                            K_mat = K_2d

                        print(f"    -> Kd_screened built in {time.time()-t1:.2f}s")

            # --- 3. Diagonalization ---
            print(f"  [Dense] Diagonalizing {self.ham.dim}x{self.ham.dim} matrix...")
            t_diag = time.time()
            if use_gpu:
                import torch
                evals_t, evecs_t = torch.linalg.eigh(H)
                evals = to_numpy(evals_t)
                evecs = to_numpy(evecs_t)
                self.J_mat = to_numpy(J_mat)
                self.K_mat = to_numpy(K_mat)
            else:
                evals, evecs = np.linalg.eigh(H)
                self.J_mat = J_mat
                self.K_mat = K_mat
            print(f"    -> Diagonalization complete in {time.time()-t_diag:.2f}s")

            self.Kx_mat = self.J_mat
            self.Kd_mat = self.K_mat
            return evals, evecs

        # Davidson Solver Fallback
        if self.ham.dim <= 2 or nroots >= self.ham.dim:
            return self.solve(nroots=nroots, full_diag=True, tol=tol, excitation_mode="bse")
        nroots = min(nroots, self.ham.dim - 1)
        print(f"  Using Davidson solver on {nroots} roots out of {self.ham.dim} transitions")
        return davidson(self.ham.matvec, self.ham.D, nroots, tol=tol, device=self.device)

    def expectation_components(self, vec):
        """Return <D_QP>, <Kx_bare>, and <-Kd_screened> for one state."""
        mode = getattr(self, "excitation_mode", "bse")
        if mode in {"independent_dft", "dft"}:
            return float(np.real(np.vdot(vec, self.ham.D_dft * vec))), 0.0, 0.0
        d_qp = float(np.real(np.vdot(vec, self.ham.D * vec)))
        if mode in {"independent_qp", "qp"}:
            return d_qp, 0.0, 0.0
        if hasattr(self, "diagonal_kx") and hasattr(self, "diagonal_kd") and self.diagonal_kx is not None and self.diagonal_kd is not None:
            idx = np.argmax(np.abs(vec))
            return d_qp, float(np.real(self.diagonal_kx[idx])), -float(np.real(self.diagonal_kd[idx]))
        if hasattr(self, "Kx_mat"):
            kx = float(np.real(np.vdot(vec, self.Kx_mat @ vec)))
            minus_kd = -float(np.real(np.vdot(vec, self.Kd_mat @ vec)))
            return d_qp, kx, minus_kd
        kx_action, kd_action = self.ham.kernel_actions(vec)
        kx = float(np.real(np.vdot(vec, kx_action)))
        minus_kd = -float(np.real(np.vdot(vec, kd_action)))
        return d_qp, kx, minus_kd

    def main_transition(self, vec):
        """Extracts the dominant hole and electron indices from a state vector."""
        idx = np.argmax(np.abs(vec))
        if getattr(self.ham, 'spin', 'singlet') == 'uks_spin_preserving' and not self.soc_flag:
            # Identify which spin channel the dominant contribution belongs to
            if idx < self.ham.dim_a:
                hole = self.ham.vi_a[idx]
                elec = self.ham.va_a[idx]
            else:
                idx_b = idx - self.ham.dim_a
                hole = self.ham.vi_b[idx_b]
                elec = self.ham.va_b[idx_b]
            return hole, elec, abs(vec[idx])
        if not self.soc_flag:
            hole = self.ham.valid_i[idx]
            elec = self.ham.valid_a[idx]
        else:
            # Map the truncated index back to the full spinor grid index
            full_idx = self.ham.valid_spinor_idx[idx] if hasattr(self.ham, 'valid_spinor_idx') else idx
            hole = full_idx // self.ham.n_virt_spinor
            elec = full_idx % self.ham.n_virt_spinor
        return hole, elec, abs(vec[idx])
