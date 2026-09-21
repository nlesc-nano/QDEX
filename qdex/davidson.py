import numpy as np
from qdex.device_utils import is_gpu, to_tensor, to_numpy

def davidson(matvec, diag, nroots, max_iter=500, tol=1e-6, max_subspace=None, device="numpy"):
    
    n = len(diag)
    # Dynamically set max subspace if not provided (5x roots is standard)
    if max_subspace is None:
        max_subspace = max(50, 5 * nroots)

    if is_gpu(device):
        import torch
        dev = torch.device(device) if isinstance(device, str) else device
        diag_t = to_tensor(diag, dev)

        # Complex arithmetic is required for SOC/spinor Hamiltonians.
        rng = np.random.default_rng(0)
        V_np = rng.standard_normal((n, nroots)) + 1j * rng.standard_normal((n, nroots))
        V = to_tensor(V_np, dev, dtype=torch.complex128)
        V, _ = torch.linalg.qr(V)

        def eval_mv(vec):
            out = matvec(vec)
            if isinstance(out, torch.Tensor):
                return out.to(device=dev, dtype=torch.complex128)
            return to_tensor(out, dev, dtype=torch.complex128)

        AV = torch.column_stack([eval_mv(V[:, i]) for i in range(V.shape[1])])

        for it in range(max_iter):
            Hsub = V.conj().T @ AV
            Hsub = 0.5 * (Hsub + Hsub.conj().T)

            evals, evecs = torch.linalg.eigh(Hsub)
            Ritz = V @ evecs[:, :nroots]
            Ritz_AV = AV @ evecs[:, :nroots]

            residuals = Ritz_AV - evals[:nroots] * Ritz
            norms = torch.linalg.norm(residuals, dim=0)
            print(f"[DAV] Iter {it:3d} residuals:", np.round(to_numpy(norms), 6))

            if bool(torch.all(norms < tol)):
                print(f"[DAV] Converged in {it} iterations.")
                return to_numpy(evals[:nroots]), to_numpy(Ritz)

            if V.shape[1] >= max_subspace:
                V = Ritz
                V, _ = torch.linalg.qr(V)
                AV = torch.column_stack([eval_mv(V[:, i]) for i in range(V.shape[1])])
                continue

            new_vecs = []
            for i in range(nroots):
                if norms[i] > tol:
                    diff = diag_t.to(dtype=evals.dtype).clone() - evals[i]
                    diff[torch.abs(diff) < 1e-4] = 1e-4
                    delta = residuals[:, i] / diff
                    new_vecs.append(delta)

            for delta in new_vecs:
                for j in range(V.shape[1]):
                    overlap = torch.vdot(V[:, j], delta)
                    delta = delta - overlap * V[:, j]

                norm = torch.linalg.norm(delta)
                if norm > 1e-5:
                    delta = delta / norm
                    V = torch.column_stack((V, delta))
                    AV = torch.column_stack((AV, eval_mv(delta)))

        raise RuntimeError(
            f"Davidson did not converge after {max_iter} iterations. Max residual: {torch.max(norms).item():.2e}"
        )

    # Complex arithmetic is required for SOC/spinor Hamiltonians.  It is also
    # harmless for a real Hermitian problem and keeps one reference path.
    rng = np.random.default_rng(0)
    V = rng.standard_normal((n, nroots)) + 1j * rng.standard_normal((n, nroots))
    V, _ = np.linalg.qr(V)

    def eval_mv_np(vec):
        out = matvec(vec)
        return to_numpy(out) if hasattr(out, 'detach') else out

    # Pre-compute AV to save matvec calls in the loop
    AV = np.column_stack([eval_mv_np(V[:, i]) for i in range(V.shape[1])])

    for it in range(max_iter):

        # Subspace Hamiltonian
        Hsub = V.conj().T @ AV
        Hsub = 0.5 * (Hsub + Hsub.conj().T)

        # Diagonalize the small subspace Hamiltonian
        evals, evecs = np.linalg.eigh(Hsub)
        
        # Get Ritz vectors and their matvecs for the lowest nroots
        Ritz = V @ evecs[:, :nroots]
        Ritz_AV = AV @ evecs[:, :nroots]

        # Compute residuals
        residuals = np.zeros((n, nroots), dtype=np.result_type(V, AV))
        for i in range(nroots):
            residuals[:, i] = Ritz_AV[:, i] - evals[i] * Ritz[:, i]

        norms = np.linalg.norm(residuals, axis=0)
        print(f"[DAV] Iter {it:3d} residuals:", np.round(norms, 6))

        # Check convergence
        if np.all(norms < tol):
            print(f"[DAV] Converged in {it} iterations.")
            return evals[:nroots], Ritz

        # Check if subspace needs to collapse
        if V.shape[1] >= max_subspace:
            # Collapse back down to just the best current Ritz vectors
            V = Ritz
            V, _ = np.linalg.qr(V)
            # Must recompute AV after a collapse
            AV = np.column_stack([eval_mv_np(V[:, i]) for i in range(V.shape[1])])
            continue

        # Generate new correction vectors
        new_vecs = []
        for i in range(nroots):
            if norms[i] > tol:
                # Preconditioner
                diff = np.asarray(diag, dtype=np.result_type(diag, evals)).copy() - evals[i]
                
                # Safeguard against denominator singularity!
                diff[np.abs(diff) < 1e-4] = 1e-4
                
                delta = residuals[:, i] / diff
                new_vecs.append(delta)

        # Modified Gram-Schmidt orthogonalization
        for delta in new_vecs:
            # Orthogonalize against all existing vectors in V
            for j in range(V.shape[1]):
                overlap = np.vdot(V[:, j], delta)
                delta -= overlap * V[:, j]
            
            norm = np.linalg.norm(delta)
            # Only add to subspace if it is sufficiently linearly independent
            if norm > 1e-5:
                delta /= norm
                V = np.column_stack((V, delta))
                # Compute matvec ONLY for the new vector and append to AV
                AV = np.column_stack((AV, eval_mv_np(delta)))

    raise RuntimeError(f"Davidson did not converge after {max_iter} iterations. Max residual: {np.max(norms):.2e}")
