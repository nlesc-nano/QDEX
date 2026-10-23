"""
One-centre extension of the xs representation (integrals.one_centre).

xs keeps only the diagonal AO products |chi_mu|^2 of the Löwdin transition densities (zero differential
overlap). The one-centre model keeps every product of two AOs on the same atom (NDDO-type densities):

    rho_pq(r) ~ sum_A sum_{mu,nu in A} c_{mu p} c_{nu q} lambda_mu(r) lambda_nu(r),     c = S^1/2 C,

with lambda the AOs of each atom orthonormalized among themselves (X_A = S_AA^-1/2), so that the one-centre
charges still sum to delta_pq. Products of AOs on different atoms are dropped, as in xs. The interaction of
two one-centre products is exact, (lambda lambda|lambda lambda), for atom pairs within the cut-off R_c (the
on-site blocks always); farther apart only the diagonal products survive with the xs integrals (mu mu|nu nu).
The model is therefore the xs kernel plus a block-sparse correction,

    (pq|X|rs) = (pq|X|rs)_xs + sum_{A,B: R_AB <= R_c} vec(D^pq_A) . Delta_AB . vec(D^rs_B),
    D^pq_A[mu, nu] = conj(c_{mu p}) c_{nu q}    (mu, nu on A),
    Delta_AB = (lambda_mu lambda_nu | lambda_ka lambda_la) - delta_{mu nu} delta_{ka la} (chi_mu chi_mu | chi_ka chi_ka).

Delta is the bare Coulomb interaction. It is added unchanged to the exchange K^x, where it completes the
bare xs kernel. In the direct term K^d the xs part carries the screened W of the kernel, and every atom-pair
block of the correction is screened with the same ratio xs uses, S_AB = sum_AB W / sum_AB (mu mu|nu nu)
(set_screening): W_1c = W_xs + S_AB Delta_AB. For the Resta-type kernels S_AA = 1, so the on-site blocks stay
bare and the dipole tail between atoms is screened like the xs charges. Only the near blocks are stored,
grouped by element pair, and every BSE product is evaluated without forming pair densities:

    exchange:  P_B = conj(T_B) X V_B^T,   U_A = sum_B Delta_AB : P_B,   (K^x x)_ia += (T^T U conj(V))_ia
    direct:    R_AB = conj(T_A) X V_B^T,  Z_AB = Delta_AB : R_AB,       (K^d x)_ia += sum_AB (T_A^T Z_AB conj(V_B))_ia

with T, V the Löwdin coefficients of the active occupied and virtual orbitals (one pair (T, V) per spin
component for spinors) and X the transition amplitudes. The cost per product is O(n_ao n_occ n_virt) plus the
near blocks; nothing of size n_pairs^2 is stored. Validated against exact four-centre integrals in the
ri4c test (one-centre densities, Löwdin coefficients: S1/T1/dST within ~70 meV for Cd16Se13Cl6 and
Cs20Pb8Br36, against 165-345 meV of missing exchange with xs).
"""
import time
import logging
import numpy as np

from qdex.constants import HA_TO_EV

logger = logging.getLogger(__name__)


class OneCentreKernel:
    """Block-sparse one-centre correction Delta to the xs kernel (energies in eV)."""

    def __init__(self, shells, atom_ao_ranges, coords_ang, overlap, gamma_xs_ev, cutoff_ang=0.0, nthreads=1,
                 tail="dipole"):
        import libint_cpp
        t0 = time.time()
        self.ranges = [tuple(map(int, r)) for r in atom_ao_ranges]
        self.n_atoms = len(self.ranges)
        self.n_ao = self.ranges[-1][1]
        self.cutoff_ang = float(cutoff_ang or 0.0)
        G = np.asarray(gamma_xs_ev)
        if G.shape != (self.n_ao, self.n_ao):
            raise ValueError(f"one-centre correction needs the bare xs matrix (n_ao x n_ao), got {G.shape}")
        self.size = np.array([a1 - a0 for (a0, a1) in self.ranges])
        self._starts = np.array([a0 for (a0, _) in self.ranges])
        self._G_atom = np.add.reduceat(np.add.reduceat(G, self._starts, axis=0), self._starts, axis=1)

        # shells of each atom are contiguous (build_shell_dicts)
        atom_of_shell = np.array([int(sh["atom_idx"]) for sh in shells])
        if np.any(np.diff(atom_of_shell) < 0):
            raise ValueError("one-centre correction expects the shells of each atom to be contiguous")
        start = np.searchsorted(atom_of_shell, np.arange(self.n_atoms + 1))
        element = [str(shells[int(start[A])]["sym"]) for A in range(self.n_atoms)]

        coords = np.asarray(coords_ang, dtype=float)
        pairs = [(A, A) for A in range(self.n_atoms)]
        if self.cutoff_ang > 0.0:
            from scipy.spatial import cKDTree
            near = cKDTree(coords).query_pairs(self.cutoff_ang, output_type="ndarray")
            pairs += [(int(a), int(b)) for a, b in near]

        # intra-atomic orthonormalization X_A = S_AA^-1/2; S_AA depends only on the basis of the element, so
        # kron(X_A, X_A) is built once per element (from the per-atom blocks of S, sparse or dense)
        kron, Xel = {}, {}
        for A, (a0, a1) in enumerate(self.ranges):
            if element[A] in kron:
                continue
            S_AA = overlap[a0:a1, a0:a1]
            S_AA = S_AA.toarray() if hasattr(S_AA, "toarray") else np.asarray(S_AA)
            w, U = np.linalg.eigh(S_AA)
            X = (U / np.sqrt(w)) @ U.T
            Xel[element[A]] = X
            kron[element[A]] = np.kron(X, X)

        # atoms with the same number of AOs ("kinds") are batched: gather index (n_kind, n_A) and position in kind
        kinds = {}
        for A in range(self.n_atoms):
            kinds.setdefault(int(self.size[A]), []).append(A)
        self.kinds = {}
        self.kind_pos = np.zeros(self.n_atoms, dtype=int)
        for nA, atoms in kinds.items():
            self.kinds[nA] = np.stack([np.arange(*self.ranges[A]) for A in atoms])
            self.kind_pos[atoms] = np.arange(len(atoms))

        # corrected blocks, grouped by (n_A, n_B, on-site); computed in chunks of atom pairs so that the raw libint
        # blocks never coexist with all corrected ones. Two layouts: exchange dx [(mu nu), (ka la)] and
        # direct dd [(mu ka), (nu la)].
        items = {}
        chunk = 2048
        for c0 in range(0, len(pairs), chunk):
            sub = pairs[c0:c0 + chunk]
            blocks = libint_cpp.compute_onecentre_coulomb(shells, [int(s) for s in start], sub, int(nthreads))
            for (A, B), raw in zip(sub, blocks):
                nA, nB = int(self.size[A]), int(self.size[B])
                delta = (kron[element[A]].T @ np.asarray(raw) @ kron[element[B]]) * HA_TO_EV
                (a0, a1), (b0, b1) = self.ranges[A], self.ranges[B]
                delta[np.ix_(np.arange(nA) * (nA + 1), np.arange(nB) * (nB + 1))] -= G[a0:a1, b0:b1]
                items.setdefault((nA, nB, A == B), []).append((A, B, delta))
            del blocks
        self.groups = []
        for (nA, nB, onsite), its in items.items():
            dx = np.stack([it[2] for it in its])                                       # (p, nA^2, nB^2)
            its_AB = [(it[0], it[1]) for it in its]
            items[(nA, nB, onsite)] = None
            dd = np.ascontiguousarray(dx.reshape(-1, nA, nA, nB, nB).transpose(0, 1, 3, 2, 4)
                                      ).reshape(-1, nA * nB, nA * nB)                    # (p, nA nB, nA nB)
            A_ = np.array([a for a, _ in its_AB])
            B_ = np.array([b for _, b in its_AB])
            self.groups.append(dict(nA=nA, nB=nB, onsite=onsite, A=A_, B=B_,
                                    posA=self.kind_pos[A_], posB=self.kind_pos[B_], dx=dx, dd=dd,
                                    s=np.ones(len(A_))))
        self.nbytes = sum(g["dx"].nbytes + g["dd"].nbytes for g in self.groups)
        n_near = len(pairs) - self.n_atoms

        # dipole tail: charge-dipole and dipole-dipole interactions of the one-centre densities of atom pairs
        # outside the exact blocks (their charge-charge part is the xs kernel)
        self.tail = str(tail if tail is not None else "dipole").lower()
        if self.tail not in ("dipole", "none"):
            raise ValueError(f"one_centre_tail must be 'dipole' or 'none', not '{tail}'")
        self._tail_cache = {}
        if self.tail == "dipole":
            self._build_dipole_tail(shells, start, element, Xel, coords, pairs, nthreads)
        logger.info(f"  [One-centre] {self.n_atoms} on-site blocks + {n_near} atom pairs within "
                    f"R_c = {self.cutoff_ang:.2f} A; tail: {self.tail}; {self.nbytes / 1e9:.2f} GB; "
                    f"{time.time() - t0:.1f} s")

    def _build_dipole_tail(self, shells, start, element, Xel, coords_ang, pairs, nthreads):
        """Moment matrices m_k[mu, nu] (k = charge, x, y, z) of the lambda products of each atom, and the
        coupling M (4 N_at x 4 N_at, eV) of charges and dipoles of different atoms:
            E = sum_AB  -q_A (t_AB . p_B) + (t_AB . p_A) q_B + p_A . T_AB . p_B,
            t_AB = (R_B - R_A)/R^3,  T_AB = (1 - 3 n n)/R^3   (charge-charge = 0: it is the xs kernel).
        Pairs with exact blocks (A = B, R_AB <= R_c) are excluded."""
        import libint_cpp
        BOHR = 0.529177210903
        mom_el = {}
        for A in range(self.n_atoms):
            el = element[A]
            if el in mom_el:
                continue
            sh = [dict(shells[s], center=np.zeros(3)) for s in range(int(start[A]), int(start[A + 1]))]
            d = np.asarray(libint_cpp.dipole(sh, [0.0, 0.0, 0.0], 1))           # <mu|r - R_A|nu>, bohr
            X = Xel[el]
            nA = X.shape[0]
            mom_el[el] = np.stack([np.eye(nA)] + [X @ d[c] @ X for c in range(3)])   # (4, nA, nA), lambda basis
        self.kmom = {nA: np.stack([mom_el[element[A]] for A in np.where(self.size == nA)[0]])
                     for nA in self.kinds}
        R = np.asarray(coords_ang, dtype=float) / BOHR
        N = self.n_atoms
        Rv = R[None, :, :] - R[:, None, :]                                     # R_B - R_A
        dist = np.linalg.norm(Rv, axis=2)
        mask = np.ones((N, N), dtype=bool)
        np.fill_diagonal(mask, False)
        for (A, B) in pairs:
            mask[A, B] = mask[B, A] = False
        with np.errstate(divide="ignore", invalid="ignore"):
            inv3 = np.where(mask, 1.0 / dist ** 3, 0.0)
            n = np.where(mask[:, :, None], Rv / np.where(dist > 0, dist, 1.0)[:, :, None], 0.0)
        t = Rv * inv3[:, :, None]                                              # (N, N, 3)
        T = (np.eye(3)[None, None] - 3.0 * n[:, :, :, None] * n[:, :, None, :]) * inv3[:, :, None, None]
        M = np.zeros((N, 4, N, 4))
        M[:, 0, :, 1:] = -t
        M[:, 1:, :, 0] = t.transpose(0, 2, 1)
        M[:, 1:, :, 1:] = T.transpose(0, 2, 1, 3)
        self.M = M.reshape(4 * N, 4 * N) * HA_TO_EV
        self.M_direct = self.M
        self.nbytes += self.M.nbytes + sum(m.nbytes for m in self.kmom.values())

    def set_screening(self, W_ev):
        """Screening of the correction in the direct term: S_AB = sum_AB W / sum_AB (mu mu|nu nu), the atom-pair
        ratio of the kernel's screened W to the bare xs integrals (as xs screens its own atom pairs). The exchange
        keeps the bare correction. W_ev = None restores S = 1."""
        if W_ev is None:
            S = np.ones((self.n_atoms, self.n_atoms))
        else:
            W = np.asarray(W_ev)
            if W.shape != (self.n_ao, self.n_ao):
                raise ValueError(f"one-centre screening needs the xs W (n_ao x n_ao), got {W.shape}")
            W_atom = np.add.reduceat(np.add.reduceat(W, self._starts, axis=0), self._starts, axis=1)
            S = W_atom / self._G_atom
            S = 0.5 * (S + S.T)                                      # keep K^d Hermitian
        self.screen_ratio = S
        for g in self.groups:
            g["s"] = S[g["A"], g["B"]]
        if self.tail == "dipole":
            N = self.n_atoms
            self.M_direct = (self.M.reshape(N, 4, N, 4) * S[:, None, :, None]).reshape(4 * N, 4 * N)
        self._tail_cache = {}
        off = ~np.eye(self.n_atoms, dtype=bool)
        logger.info(f"  [One-centre] direct-term screening S_AB = W/v: on-site {np.mean(np.diag(S)):.3f}, "
                    f"between atoms {np.min(S[off]) if off.any() else 1.0:.3f} - {np.max(S[off]) if off.any() else 1.0:.3f}")

    def with_screening(self, W_ev):
        """A view of this kernel (shared blocks and tail coupling) with its own direct-term screening."""
        import copy
        view = copy.copy(self)
        view.groups = [dict(g) for g in self.groups]
        view._tail_cache = {}
        view.set_screening(W_ev)
        return view

    # moments of per-kind densities (n, nA, nA) -> (n, 4) and back
    def _moments(self, P):
        m = np.zeros((self.n_atoms, 4), dtype=np.result_type(*P.values()))
        for nA, idx in self.kinds.items():
            atoms = np.where(self.size == nA)[0]
            m[atoms] = np.einsum("nkab,nab->nk", self.kmom[nA], P[nA])
        return m

    def _tail_direct_data(self, comps):
        """L_k[i, j] = sum_s sum T_mu i conj(T_nu j) m_k[mu nu], E'_k = sum_l M_kl E_l, E_l[a, b] = sum_s sum
        conj(V_ka a) V_la b m_l[ka la] (k = 4A + alpha); cached per set of coefficient arrays."""
        key = tuple((id(T), id(V), T.shape, V.shape) for (T, V) in comps)
        hit = self._tail_cache.get(key)
        # the entry keeps the arrays themselves: same ids are only trusted for the very same objects
        if hit is not None and all(a is b for c0, c1 in zip(hit[0], comps) for a, b in zip(c0, c1)):
            return hit[1]
        n_o, n_v = comps[0][0].shape[1], comps[0][1].shape[1]
        dt = np.result_type(*[a for c in comps for a in c])
        L = np.zeros((self.n_atoms, 4, n_o, n_o), dtype=dt)
        E = np.zeros((self.n_atoms, 4, n_v, n_v), dtype=dt)
        for nA, idx in self.kinds.items():
            atoms = np.where(self.size == nA)[0]
            mk = self.kmom[nA]                                                 # (n, 4, nA, nA)
            for (T, V) in comps:
                TA, VA = T[idx], V[idx]                                        # (n, nA, n_o), (n, nA, n_v)
                L[atoms] += np.einsum("nmi,nkmv,nvj->nkij", TA, mk, np.conj(TA), optimize=True)
                E[atoms] += np.einsum("nma,nkmv,nvb->nkab", np.conj(VA), mk, VA, optimize=True)
        Ep = (self.M_direct @ E.reshape(4 * self.n_atoms, -1)).reshape(4 * self.n_atoms, n_v, n_v)
        data = (L.reshape(4 * self.n_atoms, n_o, n_o), Ep)
        del E
        if len(self._tail_cache) > 1:
            self._tail_cache.clear()
        self._tail_cache[key] = (tuple(tuple(c) for c in comps), data)
        return data

    # ------------------------------------------------------------------ products
    def exchange_action(self, channels):
        """sum_jb (ia|Delta|jb) x_jb for every channel.

        channels: list of (comps, X) with comps a list of (T, V) Löwdin coefficient pairs (n_ao x n_o, n_ao x n_v)
        and X the amplitude matrix (n_o x n_v). The exchange couples all channels (UKS alpha and beta).
        Returns a list of (n_o x n_v) arrays.
        """
        dtype = np.result_type(*[a for comps, X in channels for (T, V) in comps for a in (T, V)],
                               *[X for _, X in channels])
        # on-site densities P_A[mu, nu] = sum_c sum_s conj(T_A) X V_A^T, per kind (n, n_A, n_A)
        P = {nA: np.zeros((len(idx), nA, nA), dtype=dtype) for nA, idx in self.kinds.items()}
        for comps, X in channels:
            for (T, V) in comps:
                Y = (X @ V.T).T                                     # (n_ao, n_o)
                for nA, idx in self.kinds.items():
                    P[nA] += np.conj(T[idx]) @ Y[idx].transpose(0, 2, 1)
        U = {nA: np.zeros_like(P[nA]) for nA in self.kinds}
        for g in self.groups:
            nA, nB = g["nA"], g["nB"]
            PB = P[nB][g["posB"]].reshape(-1, nB * nB, 1)
            np.add.at(U[nA], g["posA"], _rmatmul(g["dx"], PB).reshape(-1, nA, nA))
            if not g["onsite"]:
                PA = P[nA][g["posA"]].reshape(-1, nA * nA, 1)
                np.add.at(U[nB], g["posB"], _rmatmul(g["dx"].transpose(0, 2, 1), PA).reshape(-1, nB, nB))
        if self.tail == "dipole":
            phi = (self.M @ self._moments(P).reshape(-1)).reshape(self.n_atoms, 4)
            for nA, idx in self.kinds.items():
                atoms = np.where(self.size == nA)[0]
                U[nA] += np.einsum("nk,nkab->nab", phi[atoms], self.kmom[nA])
        out = []
        for comps, X in channels:
            n_o, n_v = X.shape
            kx = np.zeros(X.shape, dtype=dtype)
            for (T, V) in comps:
                for nA, idx in self.kinds.items():
                    # sum over atoms and AOs as one contraction: (n_o, n*nA) @ (n*nA, n_v)
                    kx += T[idx].reshape(-1, n_o).T @ (U[nA] @ np.conj(V[idx])).reshape(-1, n_v)
            out.append(kx)
        return out

    def direct_action(self, comps, X):
        """sum_jb (ij|Delta|ab) x_jb for one channel with spin components comps = [(T, V), ...]."""
        dtype = np.result_type(X, *[a for (T, V) in comps for a in (T, V)])
        n_o, n_v = X.shape
        kd = np.zeros(X.shape, dtype=dtype)
        Tk = [{nA: T[idx] for nA, idx in self.kinds.items()} for (T, _) in comps]            # (n, nA, n_o)
        Vk = [{nA: V[idx] for nA, idx in self.kinds.items()} for (_, V) in comps]            # (n, nA, n_v)
        Yk = []
        for (_, V) in comps:
            Y = (X @ V.T).T                                                                  # sum_b x_jb V_lb
            Yk.append({nA: Y[idx] for nA, idx in self.kinds.items()})
        nc = len(comps)
        for g in self.groups:
            nA, nB = g["nA"], g["nB"]
            p = len(g["A"])
            # all spin-component pairs (s, t) in one application of Delta: R_st[nu, la] = conj(T_s,A) X V_t,B^T
            R = np.stack([np.conj(Tk[a][nA][g["posA"]]) @ Yk[b][nB][g["posB"]].transpose(0, 2, 1)
                          for a in range(nc) for b in range(nc)], axis=-1)                   # (p, nA, nB, nc^2)
            Z = _rmatmul(g["dd"], R.reshape(p, nA * nB, -1)).reshape(p, nA, nB, nc, nc)      # Z_st[mu, ka]
            Z = Z * g["s"][:, None, None, None, None]                                        # screened (direct)
            for a in range(nc):
                TA = Tk[a][nA][g["posA"]]
                for b in range(nc):
                    VB = Vk[b][nB][g["posB"]]
                    kd += TA.reshape(-1, n_o).T @ (Z[..., a, b] @ np.conj(VB)).reshape(-1, n_v)
            if not g["onsite"]:
                R2 = np.stack([np.conj(Tk[a][nB][g["posB"]]) @ Yk[b][nA][g["posA"]].transpose(0, 2, 1)
                               for a in range(nc) for b in range(nc)], axis=-1)              # R_BA[la, nu]
                Zt = _rmatmul(g["dd"], R2.transpose(0, 2, 1, 3).reshape(p, nA * nB, -1)).reshape(p, nA, nB, nc, nc)
                Zt = Zt * g["s"][:, None, None, None, None]
                for a in range(nc):
                    TB = Tk[a][nB][g["posB"]]
                    for b in range(nc):
                        VA = Vk[b][nA][g["posA"]]
                        kd += TB.reshape(-1, n_o).T @ (Zt[..., a, b].transpose(0, 2, 1) @ np.conj(VA)).reshape(-1, n_v)
        if self.tail == "dipole":
            L, Ep = self._tail_direct_data(comps)
            kd = kd.astype(np.result_type(kd, L, Ep), copy=False)
            for k0 in range(0, L.shape[0], 512):
                sl = slice(k0, k0 + 512)
                Yk = X @ Ep[sl].transpose(0, 2, 1)                               # (k, n_o, n_v)
                kd += L[sl].transpose(1, 0, 2).reshape(n_o, -1) @ Yk.reshape(-1, n_v)
        return kd

    # ------------------------------------------------------------------ diagonal elements
    def exchange_diag(self, comps, vi, va, chunk=256, batch=128):
        """(ia|Delta|ia) for the transitions (vi, va) of one channel."""
        vi, va = np.asarray(vi), np.asarray(va)
        out = np.zeros(len(vi), dtype=float)
        for t0 in range(0, len(vi), chunk):
            i, a = vi[t0:t0 + chunk], va[t0:t0 + chunk]
            Q = {}
            for nA, idx in self.kinds.items():
                q = 0.0
                for (T, V) in comps:
                    Ti, Va = T[:, i][idx], V[:, a][idx]                           # (n, nA, t): only the chunk
                    # q[n, t, mu, nu] = conj(T_mu i) V_nu a
                    q = q + np.conj(Ti).transpose(0, 2, 1)[:, :, :, None] * Va.transpose(0, 2, 1)[:, :, None, :]
                Q[nA] = q.reshape(len(idx), len(i), nA * nA)
            for g in self.groups:
                nA, nB = g["nA"], g["nB"]
                w = 1.0 if g["onsite"] else 2.0
                for b0 in range(0, len(g["A"]), batch):
                    sl = slice(b0, b0 + batch)
                    QA, QB = Q[nA][g["posA"][sl]], Q[nB][g["posB"][sl]]
                    DQ = _rmatmul(g["dx"][sl], QB.transpose(0, 2, 1)).transpose(0, 2, 1)   # (p, t, nA^2)
                    val = np.sum(np.conj(QA) * DQ, axis=(0, 2)).real
                    out[t0:t0 + chunk] += w * val
            if self.tail == "dipole":
                m = np.zeros((self.n_atoms, 4, len(i)), dtype=np.result_type(*[a for c in comps for a in c]))
                for nA, idx in self.kinds.items():
                    atoms = np.where(self.size == nA)[0]
                    m[atoms] += np.einsum("ntab,nkab->nkt", Q[nA].reshape(len(idx), len(i), nA, nA),
                                          self.kmom[nA], optimize=True)
                m = m.reshape(4 * self.n_atoms, len(i))
                out[t0:t0 + chunk] += np.sum(np.conj(m) * (self.M @ m), axis=0).real
        return out

    def direct_diag(self, comps, ochunk=256, vchunk=256, batch=128):
        """(ii|Delta|aa) as an (n_o x n_v) matrix for one channel, in chunks of occupied and virtual orbitals."""
        n_o, n_v = comps[0][0].shape[1], comps[0][1].shape[1]
        out = np.zeros((n_o, n_v), dtype=float)

        def moments(cols, which, nA, idx):
            m = 0.0
            for (T, V) in comps:
                C = (T if which == "o" else V)[:, cols][idx].transpose(0, 2, 1)   # (n, k, nA)
                if which == "o":
                    m = m + C[:, :, :, None] * np.conj(C)[:, :, None, :]           # T_mu i conj(T_nu i)
                else:
                    m = m + np.conj(C)[:, :, :, None] * C[:, :, None, :]           # conj(V_ka a) V_la a
            return m.reshape(len(idx), len(cols), nA * nA)

        for o0 in range(0, n_o, ochunk):
            oc = np.arange(o0, min(n_o, o0 + ochunk))
            H = {nA: moments(oc, "o", nA, idx) for nA, idx in self.kinds.items()}
            for v0 in range(0, n_v, vchunk):
                vc = np.arange(v0, min(n_v, v0 + vchunk))
                E = {nA: moments(vc, "v", nA, idx) for nA, idx in self.kinds.items()}
                blk = np.zeros((len(oc), len(vc)), dtype=float)
                for g in self.groups:
                    nA, nB = g["nA"], g["nB"]
                    for b0 in range(0, len(g["A"]), batch):
                        sl = slice(b0, b0 + batch)
                        dx = g["dx"][sl]
                        HA, EB = H[nA][g["posA"][sl]], E[nB][g["posB"][sl]]
                        L = _rmatmul(dx.transpose(0, 2, 1), HA.transpose(0, 2, 1)).transpose(0, 2, 1)  # (p, k_o, nB^2)
                        L = L * g["s"][sl][:, None, None]
                        blk += (L.transpose(1, 0, 2).reshape(len(oc), -1)
                                @ EB.transpose(1, 0, 2).reshape(len(vc), -1).T).real
                        if not g["onsite"]:
                            HB, EA = H[nB][g["posB"][sl]], E[nA][g["posA"][sl]]
                            L2 = _rmatmul(dx, HB.transpose(0, 2, 1)).transpose(0, 2, 1)                 # (p, k_o, nA^2)
                            L2 = L2 * g["s"][sl][:, None, None]
                            blk += (L2.transpose(1, 0, 2).reshape(len(oc), -1)
                                    @ EA.transpose(1, 0, 2).reshape(len(vc), -1).T).real
                out[o0:o0 + len(oc), v0:v0 + len(vc)] = blk
        if self.tail == "dipole":
            # diagonal moments only: L_k[i, i] and E_l[a, a] (no n_o^2 or n_v^2 arrays)
            Ld = np.zeros((self.n_atoms, 4, n_o), dtype=complex)
            Ed = np.zeros((self.n_atoms, 4, n_v), dtype=complex)
            for nA, idx in self.kinds.items():
                atoms = np.where(self.size == nA)[0]
                mk = self.kmom[nA]
                for (T, V) in comps:
                    TA, VA = T[idx], V[idx]
                    Ld[atoms] += np.einsum("nmi,nkmv,nvi->nki", TA, mk, np.conj(TA), optimize=True)
                    Ed[atoms] += np.einsum("nma,nkmv,nva->nka", np.conj(VA), mk, VA, optimize=True)
            Ld = Ld.reshape(4 * self.n_atoms, n_o)
            Ed = (self.M_direct @ Ed.reshape(4 * self.n_atoms, n_v))
            out += (Ld.T @ Ed).real
        return out


def _rmatmul(D, Y):
    """D @ Y for real D and real or complex Y, without converting D to complex."""
    if not np.iscomplexobj(Y):
        return D @ Y
    k = Y.shape[-1]
    r = D @ np.concatenate([Y.real, Y.imag], axis=-1)
    return r[..., :k] + 1j * r[..., k:]


_CACHE = {}


def get_onecentre_kernel(shells, atom_ao_ranges, coords_ang, overlap, gamma_xs_ev, cutoff_ang=0.0, nthreads=1,
                         tail="dipole"):
    """OneCentreKernel, reused by the spin-free and SOC solvers of one run (same geometry, basis and xs matrix)."""
    import hashlib
    G = np.asarray(gamma_xs_ev)
    h = hashlib.sha1()
    h.update(np.ascontiguousarray(np.asarray(coords_ang, dtype=float)).tobytes())
    # the diagonal and one strided sample of G identify the xs matrix without hashing all n_ao^2 numbers
    h.update(np.ascontiguousarray(np.diag(G), dtype=float).tobytes())
    h.update(np.ascontiguousarray(G.ravel()[:: max(1, G.size // 4096)], dtype=float).tobytes())
    key = (G.shape[0], round(float(cutoff_ang or 0.0), 6), str(tail), h.hexdigest())
    if key not in _CACHE:
        _CACHE.clear()
        _CACHE[key] = OneCentreKernel(shells, atom_ao_ranges, coords_ang, overlap, G, cutoff_ang or 0.0, nthreads,
                                      tail=tail)
    return _CACHE[key]
