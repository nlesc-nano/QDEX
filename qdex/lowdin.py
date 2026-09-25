import numpy as np
from qdex.device_utils import is_gpu, to_tensor, to_numpy

def lowdin_sqrt(S, device="numpy"):
    if is_gpu(device):
        import torch
        dev = torch.device(device) if isinstance(device, str) else device
        S_t = to_tensor(S, dev, dtype=torch.float64)
        eigvals, eigvecs = torch.linalg.eigh(S_t)
        eigvals = torch.clamp(eigvals, min=1e-15)
        S_half = eigvecs @ torch.diag(torch.sqrt(eigvals)) @ eigvecs.T
        return to_numpy(S_half)

    eigvals, eigvecs = np.linalg.eigh(S)
    eigvals = np.clip(eigvals, a_min=1e-15, a_max=None)
    return eigvecs @ np.diag(np.sqrt(eigvals)) @ eigvecs.T

def transform_mos(C, S, device="numpy"):
    S_half = lowdin_sqrt(S, device=device)
    return S_half @ C

def build_lowdin_transition_charges_flat(C_occ_act, C_virt_act, S, atom_ao_ranges, valid_i, valid_a, device="numpy"):
    """
    Computes transition charges using Löwdin symmetric orthogonalization:
    C^L = S^{1/2} * C.
    q^L_{ia, A} = sum_{mu in A} C^L_{mu i} * C^L_{mu a}.
    """
    if is_gpu(device):
        import torch
        dev = torch.device(device) if isinstance(device, str) else device
        S_t = to_tensor(S, dev, dtype=torch.float64)
        eigvals, eigvecs = torch.linalg.eigh(S_t)
        eigvals = torch.clamp(eigvals, min=1e-15)
        S_half = eigvecs @ torch.diag(torch.sqrt(eigvals)) @ eigvecs.T

        C_occ_t = to_tensor(C_occ_act, dev, dtype=torch.float64)
        C_virt_t = to_tensor(C_virt_act, dev, dtype=torch.float64)
        C_occ_lowdin = S_half @ C_occ_t
        C_virt_lowdin = S_half @ C_virt_t

        dim = len(valid_i)
        n_atoms = len(atom_ao_ranges)
        q_flat_t = torch.zeros((dim, n_atoms), dtype=torch.float64, device=dev)
        vi_t = torch.as_tensor(valid_i, device=dev, dtype=torch.long)
        va_t = torch.as_tensor(valid_a, device=dev, dtype=torch.long)

        for A, (a0, a1) in enumerate(atom_ao_ranges):
            Ci_A = C_occ_lowdin[a0:a1, :][:, vi_t]
            Ca_A = C_virt_lowdin[a0:a1, :][:, va_t]
            q_flat_t[:, A] = torch.sum(Ci_A * Ca_A, dim=0)

        return to_numpy(q_flat_t)

    S_half = lowdin_sqrt(S, device=device)
    C_occ_lowdin = S_half @ C_occ_act
    C_virt_lowdin = S_half @ C_virt_act
    
    dim = len(valid_i)
    n_atoms = len(atom_ao_ranges)
    q_flat = np.zeros((dim, n_atoms))
    
    for A, (a0, a1) in enumerate(atom_ao_ranges):
        Ci_A = C_occ_lowdin[a0:a1, :][:, valid_i]
        Ca_A = C_virt_lowdin[a0:a1, :][:, valid_a]
        q_flat[:, A] = np.sum(Ci_A * Ca_A, axis=0)
        
    return q_flat


def build_xs_transition_densities_flat(C_occ_act, C_virt_act, S, valid_i, valid_a, device="numpy"):
    """
    Computes AO-resolved transition densities for XsTD-DFT:
    C^L = S^{1/2} * C
    Q_{ia}(\mu) = C^L_{\mu i} * C^L_{\mu a}
    Returns shape (dim, n_ao).
    """
    if is_gpu(device):
        import torch
        dev = torch.device(device) if isinstance(device, str) else device
        S_t = to_tensor(S, dev, dtype=torch.float64)
        eigvals, eigvecs = torch.linalg.eigh(S_t)
        eigvals = torch.clamp(eigvals, min=1e-15)
        S_half = eigvecs @ torch.diag(torch.sqrt(eigvals)) @ eigvecs.T

        C_occ_t = to_tensor(C_occ_act, dev, dtype=torch.float64)
        C_virt_t = to_tensor(C_virt_act, dev, dtype=torch.float64)
        C_occ_lowdin = S_half @ C_occ_t
        C_virt_lowdin = S_half @ C_virt_t

        vi_t = torch.as_tensor(valid_i, device=dev, dtype=torch.long)
        va_t = torch.as_tensor(valid_a, device=dev, dtype=torch.long)
        Q_t = (C_occ_lowdin[:, vi_t] * C_virt_lowdin[:, va_t]).T
        return to_numpy(Q_t)

    S_half = lowdin_sqrt(S, device=device)
    C_occ_lowdin = S_half @ C_occ_act
    C_virt_lowdin = S_half @ C_virt_act
    Q = (C_occ_lowdin[:, valid_i] * C_virt_lowdin[:, valid_a]).T
    return Q


def build_xs_state_densities(C_occ_act, C_virt_act, S, device="numpy"):
    """
    Computes AO-resolved pair densities for occupied and virtual MOs for XsTD-DFT:
    Q_occ[i, j, mu] = C^L_{mu i} * C^L_{mu j}
    Q_virt[a, b, mu] = C^L_{mu a} * C^L_{mu b}
    Q_ov[i, a, mu] = C^L_{mu i} * C^L_{mu a}
    """
    S_half = lowdin_sqrt(S, device=device)
    C_occ_lowdin = S_half @ C_occ_act
    C_virt_lowdin = S_half @ C_virt_act

    Q_occ = np.einsum("mi,mj->ijm", C_occ_lowdin, C_occ_lowdin, optimize=True)
    Q_virt = np.einsum("ma,mb->abm", C_virt_lowdin, C_virt_lowdin, optimize=True)
    Q_ov = np.einsum("mi,ma->iam", C_occ_lowdin, C_virt_lowdin, optimize=True)
    return Q_occ, Q_virt, Q_ov
