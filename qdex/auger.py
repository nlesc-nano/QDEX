"""
QDEX Auger Recombination Module.

Implements non-radiative multi-carrier Auger recombination rates in colloidal
quantum dots and nanoclusters using atom-centered transition charges and the
microscopic Resta-screened dielectric kernel.

Channels:
  - Negative Trion / Biexciton eeh: e2 + h -> recombine, ejecting e1 -> e' (conduction continuum)
  - Positive Trion / Biexciton hhe: h2 + e -> recombine, ejecting h1 -> h' (valence continuum)
  - Neutral Biexciton XX: Gamma_XX = 2 * Gamma_eeh + 2 * Gamma_hhe, tau_XX = 1 / Gamma_XX
"""

import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

# Physical constants
HBAR_EV_FS = 0.6582119569  # hbar in eV * fs
KB_EV = 8.617333262e-5      # Boltzmann constant in eV / K
TWO_PI_OVER_HBAR = 2.0 * np.pi / HBAR_EV_FS  # ~9.546387 fs^-1 eV^-1


@dataclass
class AugerChannelDetail:
    """Detailed record of a single Auger transition pathway."""
    channel_type: str        # 'eeh' or 'hhe'
    carrier_e1: int          # 0-indexed orbital index of spectator carrier
    carrier_recomb1: int     # 0-indexed orbital index of recombining carrier 1
    carrier_recomb2: int     # 0-indexed orbital index of recombining carrier 2
    continuum_state: int     # 0-indexed orbital index of ejected carrier
    energy_recomb_ev: float  # Recombination energy released (eV)
    energy_eject_ev: float   # Excitation energy of ejected carrier (eV)
    energy_mismatch_ev: float # E_i - E_f (eV)
    v_direct_ev: float       # Direct Coulomb matrix element (eV)
    v_exchange_ev: float     # Exchange Coulomb matrix element (eV)
    m_eff_ev: float          # Net effective matrix element |M_if| (eV)
    rate_contrib_fs: float   # Rate contribution in fs^-1


@dataclass
class AugerResult:
    """Comprehensive container for Auger recombination calculations."""
    # Rates in fs^-1
    rate_eeh_fs: float = 0.0
    rate_hhe_fs: float = 0.0
    rate_biexciton_fs: float = 0.0

    # Rates in ps^-1
    rate_eeh_ps: float = 0.0
    rate_hhe_ps: float = 0.0
    rate_biexciton_ps: float = 0.0

    # Rates in s^-1
    rate_eeh_s: float = 0.0
    rate_hhe_s: float = 0.0
    rate_biexciton_s: float = 0.0

    # Lifetimes in ps
    tau_eeh_ps: float = np.inf
    tau_hhe_ps: float = np.inf
    tau_biexciton_ps: float = np.inf

    # Diagnostics & metadata
    fundamental_gap_ev: float = 0.0
    broadening_mode: str = "gaussian"
    sigma_ev: float = 0.05
    n_eeh_pathways: int = 0
    n_hhe_pathways: int = 0
    dominant_eeh: List[AugerChannelDetail] = field(default_factory=list)
    dominant_hhe: List[AugerChannelDetail] = field(default_factory=list)
    computation_time_s: float = 0.0

    def summary_table(self) -> str:
        """Returns a formatted ASCII summary of the Auger recombination results."""
        lines = [
            "================================================================================",
            "                     QDEX AUGER RECOMBINATION REPORT                           ",
            "================================================================================",
            f"  Fundamental Bandgap (E_g)       : {self.fundamental_gap_ev:8.4f} eV",
            f"  Energy Conservation Line Shape  : {self.broadening_mode.upper()} (sigma = {self.sigma_ev*1e3:.1f} meV)",
            f"  Evaluated Active Pathways       : {self.n_eeh_pathways} (eeh) | {self.n_hhe_pathways} (hhe)",
            "--------------------------------------------------------------------------------",
            "  Channel                 Rate (fs^-1)        Rate (ps^-1)         Lifetime (ps)",
            "--------------------------------------------------------------------------------",
            f"  Negative Trion (eeh)   {self.rate_eeh_fs:14.4e}    {self.rate_eeh_ps:14.4e}    {self.tau_eeh_ps:14.3f} ps",
            f"  Positive Trion (hhe)   {self.rate_hhe_fs:14.4e}    {self.rate_hhe_ps:14.4e}    {self.tau_hhe_ps:14.3f} ps",
            f"  Biexciton (XX)         {self.rate_biexciton_fs:14.4e}    {self.rate_biexciton_ps:14.4e}    {self.tau_biexciton_ps:14.3f} ps",
            "--------------------------------------------------------------------------------",
        ]

        if self.dominant_eeh:
            lines.append("  Dominant eeh Pathways (Spectator e1 + [e2 -> h] -> e'):")
            lines.append("    e1   e2    h ->   e' |  Delta E_recomb |  E_mismatch |   |M_eff|   |  tau_path (ps)")
            for d in self.dominant_eeh[:5]:
                tau_path = 1.0 / (d.rate_contrib_fs * 1e3) if d.rate_contrib_fs > 0 else np.inf
                lines.append(
                    f"   {d.carrier_e1:3d}  {d.carrier_recomb1:3d}  {d.carrier_recomb2:3d} -> {d.continuum_state:4d} | "
                    f"   {d.energy_recomb_ev:7.3f} eV  |  {d.energy_mismatch_ev:8.4f} eV | {d.m_eff_ev*1e3:7.2f} meV | {tau_path:10.2f} ps"
                )

        if self.dominant_hhe:
            lines.append("  Dominant hhe Pathways (Spectator h1 + [h2 -> e] -> h'):")
            lines.append("    h1   h2    e ->   h' |  Delta E_recomb |  E_mismatch |   |M_eff|   |  tau_path (ps)")
            for d in self.dominant_hhe[:5]:
                tau_path = 1.0 / (d.rate_contrib_fs * 1e3) if d.rate_contrib_fs > 0 else np.inf
                lines.append(
                    f"   {d.carrier_e1:3d}  {d.carrier_recomb1:3d}  {d.carrier_recomb2:3d} -> {d.continuum_state:4d} | "
                    f"   {d.energy_recomb_ev:7.3f} eV  |  {d.energy_mismatch_ev:8.4f} eV | {d.m_eff_ev*1e3:7.2f} meV | {tau_path:10.2f} ps"
                )

        lines.append("================================================================================")
        return "\n".join(lines)


def compute_transition_density_charges(
    C: np.ndarray,
    S: np.ndarray,
    atom_ao_ranges: List[Tuple[int, int]],
    i_indices: Union[List[int], np.ndarray],
    j_indices: Union[List[int], np.ndarray],
    spinor: bool = False,
    U_spinor_alpha: Optional[np.ndarray] = None,
    U_spinor_beta: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Computes atom-centered transition charges using Mulliken population analysis.

    For spatial orbitals:
      q_A^{ij} = sum_{mu in A, nu} C_{mu i} S_{mu nu} C_{nu j}

    For two-component spinors:
      q_A^{ij} = sum_{mu in A, nu} [ (U_alpha)_{mu i}^* S_{mu nu} (U_alpha)_{nu j} +
                                      (U_beta)_{mu i}^* S_{mu nu} (U_beta)_{nu j} ]

    Parameters
    ----------
    C : ndarray of shape (N_ao, N_mo)
        Molecular orbital coefficients (real or complex).
    S : ndarray of shape (N_ao, N_ao)
        Atomic orbital overlap matrix.
    atom_ao_ranges : list of (start, end) tuples
        AO slice ranges corresponding to each atomic center.
    i_indices : list or array of ints
        Left state indices.
    j_indices : list or array of ints
        Right state indices.
    spinor : bool
        Whether to compute spinor transition charges using U_spinor_alpha/beta.
    U_spinor_alpha, U_spinor_beta : ndarray, optional
        Spinor expansion coefficients in AO basis for alpha/beta components.

    Returns
    -------
    q : ndarray of shape (len(i_indices), len(j_indices), N_atoms)
        Atom-partitioned transition charges.
    """
    i_idx = np.asarray(i_indices, dtype=int)
    j_idx = np.asarray(j_indices, dtype=int)
    n_i, n_j = len(i_idx), len(j_idx)
    n_atoms = len(atom_ao_ranges)

    if spinor and U_spinor_alpha is not None and U_spinor_beta is not None:
        dtype = np.complex128
        Ui_a = U_spinor_alpha[:, i_idx]
        Uj_a = U_spinor_alpha[:, j_idx]
        Ui_b = U_spinor_beta[:, i_idx]
        Uj_b = U_spinor_beta[:, j_idx]

        SUj_a = S @ Uj_a
        SUj_b = S @ Uj_b

        q = np.zeros((n_i, n_j, n_atoms), dtype=dtype)
        for A, (a0, a1) in enumerate(atom_ao_ranges):
            block_a = Ui_a[a0:a1, :].conj().T @ SUj_a[a0:a1, :]
            block_b = Ui_b[a0:a1, :].conj().T @ SUj_b[a0:a1, :]
            q[:, :, A] = block_a + block_b
        return q

    # Spatial orbitals (vectorized per-atom DGEMM)
    dtype = C.dtype
    Ci = C[:, i_idx]
    Cj = C[:, j_idx]
    SCj = S @ Cj

    q = np.zeros((n_i, n_j, n_atoms), dtype=dtype)
    for A, (a0, a1) in enumerate(atom_ao_ranges):
        # (N_i, N_ao_atom) @ (N_ao_atom, N_j) -> (N_i, N_j)
        q[:, :, A] = Ci[a0:a1, :].conj().T @ SCj[a0:a1, :]

    return q


def compute_auger_matrix_element(
    q_recomb: np.ndarray,
    q_eject_1: np.ndarray,
    q_eject_2: np.ndarray,
    W_resta: np.ndarray,
) -> Tuple[float, float, float]:
    """
    Computes direct, exchange, and effective net Auger matrix elements:
      V_dir  = (q_recomb)^T @ W_resta @ q_eject_1
      V_exch = (q_recomb_alt)^T @ W_resta @ q_eject_2
      M_eff  = V_dir - V_exch

    Parameters
    ----------
    q_recomb : ndarray of shape (N_atoms,)
        Transition charge for carrier annihilation (e.g. e2 -> h).
    q_eject_1 : ndarray of shape (N_atoms,)
        Transition charge for spectator carrier excitation (e.g. e1 -> e').
    q_eject_2 : ndarray of shape (N_atoms,)
        Alternative transition charge for exchange channel (e.g. e2 -> e').
    W_resta : ndarray of shape (N_atoms, N_atoms)
        Resta screened Coulomb kernel in eV.

    Returns
    -------
    v_dir, v_exch, m_eff : float (in eV)
    """
    W_q_recomb = W_resta @ q_recomb
    v_dir = float(np.vdot(q_eject_1, W_q_recomb).real)

    W_q_alt = W_resta @ q_eject_2
    v_exch = float(np.vdot(q_recomb, W_q_alt).real)

    m_eff = v_dir - v_exch
    return v_dir, v_exch, m_eff


def calculate_auger_rates(
    C: np.ndarray,
    eps: np.ndarray,
    S: np.ndarray,
    atom_ao_ranges: List[Tuple[int, int]],
    coords: np.ndarray,
    atom_symbols: List[str],
    homo_idx: int,
    W_resta: Optional[np.ndarray] = None,
    material_name: Optional[str] = "DEFAULT",
    eps_out: float = 2.0,
    sigma_ev: float = 0.05,
    broadening_mode: str = "gaussian",
    channel: str = "all",
    n_initial_elec: int = 1,
    n_initial_hole: int = 1,
    e_search_sigma_factor: float = 4.0,
    spinor: bool = False,
    U_spinor_alpha: Optional[np.ndarray] = None,
    U_spinor_beta: Optional[np.ndarray] = None,
    lambda_reorg_ev: Optional[float] = None,
    temperature_k: float = 300.0,
    verbose: bool = True,
) -> AugerResult:
    """
    Computes non-radiative Auger recombination rates and lifetimes using Fermi's Golden Rule
    and the atom-centered Resta-screened Coulomb interaction.

    Parameters
    ----------
    C : ndarray of shape (N_ao, N_mo)
        Molecular orbital coefficients (spatial or spinor basis).
    eps : ndarray of shape (N_mo,)
        Quasiparticle or DFT orbital eigenvalues in eV.
    S : ndarray of shape (N_ao, N_ao)
        AO overlap matrix.
    atom_ao_ranges : list of (start, end) tuples
        AO slice ranges for each atomic center.
    coords : ndarray of shape (N_atoms, 3)
        Atomic coordinates in Angstroms.
    atom_symbols : list of str
        Element symbol of each atom.
    homo_idx : int
        0-indexed position of the highest occupied orbital (HOMO).
    W_resta : ndarray of shape (N_atoms, N_atoms), optional
        Precomputed Resta kernel matrix in eV. If None, built on-the-fly.
    material_name : str, optional
        Material identifier for dielectric parameters (default: 'DEFAULT').
    eps_out : float, optional
        Solvent/external dielectric constant.
    sigma_ev : float, optional
        Energy conservation Gaussian broadening parameter in eV (default: 0.05 eV).
    broadening_mode : str, optional
        'gaussian' (default) or 'fcwd' (Marcus multi-phonon Franck-Condon line shape).
    channel : str, optional
        'all' (both eeh and hhe), 'eeh', or 'hhe'.
    n_initial_elec : int, optional
        Number of conduction band-edge states to consider as initial electron carriers (default: 1, LUMO).
    n_initial_hole : int, optional
        Number of valence band-edge states to consider as initial hole carriers (default: 1, HOMO).
    e_search_sigma_factor : float, optional
        Search window multiplier for energy conservation: abs(E_i - E_f) <= factor * sigma_ev.
    spinor : bool, optional
        Whether wavefunctions are two-component spinors with SOC enabled.
    U_spinor_alpha, U_spinor_beta : ndarray, optional
        Alpha and beta spinor AO expansion matrices if spinor=True.
    lambda_reorg_ev : float, optional
        Nuclear reorganization energy in eV for FCWD line shape.
    temperature_k : float, optional
        Temperature in Kelvin (default: 300 K).
    verbose : bool, optional
        Print diagnostic progress and summary table.

    Returns
    -------
    AugerResult
        Dataclass containing total rates, lifetimes, and channel breakdown.
    """
    t0 = time.time()
    n_mo = len(eps)
    n_atoms = len(atom_ao_ranges)

    # 1. Build Resta Kernel if not provided
    if W_resta is None:
        from qdex.hardness import build_resta_mnok
        W_resta, _ = build_resta_mnok(
            atom_symbols=atom_symbols,
            coords=coords,
            alpha=0.0,
            material_name=material_name or "DEFAULT",
            eps_out=eps_out,
        )

    # 2. Identify Active Frontier States
    lumo_idx = homo_idx + 1
    e_homo = eps[homo_idx]
    e_lumo = eps[lumo_idx]
    e_gap = e_lumo - e_homo

    initial_elecs = [lumo_idx + k for k in range(n_initial_elec) if lumo_idx + k < n_mo]
    initial_holes = [homo_idx - k for k in range(n_initial_hole) if homo_idx - k >= 0]

    # Pre-filter virtual candidates e' in conduction band (for eeh)
    # Expected energy: eps[e'] ~ eps[e1] + e_gap
    # Pre-filter deep valence candidates h' (for hhe)
    # Expected energy: eps[h'] ~ eps[h1] - e_gap
    all_virt_indices = np.arange(lumo_idx, n_mo)
    all_occ_indices = np.arange(0, homo_idx + 1)

    eeh_details: List[AugerChannelDetail] = []
    hhe_details: List[AugerChannelDetail] = []

    rate_eeh_fs = 0.0
    rate_hhe_fs = 0.0

    # Line-shape evaluator
    kb_t = KB_EV * temperature_k
    two_sigma_sq = 2.0 * (sigma_ev ** 2)
    norm_gauss = 1.0 / (np.sqrt(2.0 * np.pi) * sigma_ev)

    def eval_lineshape(dE: float) -> float:
        if broadening_mode.lower() == "fcwd" and lambda_reorg_ev is not None and lambda_reorg_ev > 0:
            # Marcus multi-phonon FCWD line shape:
            denom = np.sqrt(4.0 * np.pi * lambda_reorg_ev * kb_t)
            exponent = -((dE - lambda_reorg_ev) ** 2) / (4.0 * lambda_reorg_ev * kb_t)
            return float((1.0 / denom) * np.exp(exponent))
        # Standard Gaussian line shape
        return float(norm_gauss * np.exp(-(dE ** 2) / two_sigma_sq))

    # =========================================================================
    # CHANNEL 1: NEGATIVE TRION / BIEXCITON (eeh)
    # Spectator: e1 (conduction), Recombining pair: (e2, h), Ejected: e' (virt)
    # =========================================================================
    if channel.lower() in ["all", "eeh"]:
        for e1 in initial_elecs:
            for e2 in initial_elecs:
                for h in initial_holes:
                    # Energy released by (e2, h) recombination
                    dE_recomb = eps[e2] - eps[h]
                    # Required final state energy: eps[e'] = eps[e1] + dE_recomb
                    target_e_prime = eps[e1] + dE_recomb

                    # Candidate virtual states within energy window
                    window = e_search_sigma_factor * sigma_ev
                    cand_e_prime = [
                        ep for ep in all_virt_indices
                        if abs(eps[ep] - target_e_prime) <= window and ep not in (e1, e2)
                    ]

                    if not cand_e_prime:
                        continue

                    # Precompute transition charges for this set of states
                    # 1. q^{e2, h}: recombination charge
                    q_recomb_arr = compute_transition_density_charges(
                        C=C, S=S, atom_ao_ranges=atom_ao_ranges,
                        i_indices=[e2], j_indices=[h],
                        spinor=spinor, U_spinor_alpha=U_spinor_alpha, U_spinor_beta=U_spinor_beta
                    )
                    q_recomb = q_recomb_arr[0, 0, :].real  # (N_atoms,)

                    # 2. q^{e1, h}: alternative recombination charge for exchange
                    q_recomb_alt_arr = compute_transition_density_charges(
                        C=C, S=S, atom_ao_ranges=atom_ao_ranges,
                        i_indices=[e1], j_indices=[h],
                        spinor=spinor, U_spinor_alpha=U_spinor_alpha, U_spinor_beta=U_spinor_beta
                    )
                    q_recomb_alt = q_recomb_alt_arr[0, 0, :].real  # (N_atoms,)

                    # Pre-multiply kernel with recombination vectors:
                    W_q_recomb = W_resta @ q_recomb        # (N_atoms,)
                    W_q_recomb_alt = W_resta @ q_recomb_alt # (N_atoms,)

                    # 3. Transition charges from e1 and e2 to all candidate e'
                    q_eject_1_arr = compute_transition_density_charges(
                        C=C, S=S, atom_ao_ranges=atom_ao_ranges,
                        i_indices=cand_e_prime, j_indices=[e1],
                        spinor=spinor, U_spinor_alpha=U_spinor_alpha, U_spinor_beta=U_spinor_beta
                    )[:, 0, :].real  # (len(cand_e_prime), N_atoms)

                    q_eject_2_arr = compute_transition_density_charges(
                        C=C, S=S, atom_ao_ranges=atom_ao_ranges,
                        i_indices=cand_e_prime, j_indices=[e2],
                        spinor=spinor, U_spinor_alpha=U_spinor_alpha, U_spinor_beta=U_spinor_beta
                    )[:, 0, :].real  # (len(cand_e_prime), N_atoms)

                    # Vectorized matrix elements
                    v_dir_vec = q_eject_1_arr @ W_q_recomb       # (len(cand),)
                    v_exch_vec = q_eject_2_arr @ W_q_recomb_alt   # (len(cand),)

                    # Effective net amplitude (spin-averaged singlet biexciton factor)
                    # For singlet electrons: |V_dir - V_exch|^2 + |V_dir|^2 = 2|V_dir|^2 + |V_exch|^2 - 2 Re(V_dir* V_exch)
                    if spinor:
                        # Full spinor anti-symmetrization directly includes spin factors
                        m_eff_vec = v_dir_vec - v_exch_vec
                        m_sq_vec = np.abs(m_eff_vec) ** 2
                    else:
                        # Spatial closed-shell spin statistics for 2e in conduction edge:
                        # 1 channel with parallel spins (V_dir - V_exch) and 1 with anti-parallel (V_dir)
                        m_sq_vec = 0.5 * ((v_dir_vec - v_exch_vec) ** 2 + v_dir_vec ** 2)

                    for idx, ep in enumerate(cand_e_prime):
                        dE_mismatch = (eps[ep] - eps[e1]) - dE_recomb
                        rho = eval_lineshape(dE_mismatch)
                        d_rate_fs = TWO_PI_OVER_HBAR * m_sq_vec[idx] * rho
                        rate_eeh_fs += d_rate_fs

                        eeh_details.append(
                            AugerChannelDetail(
                                channel_type="eeh",
                                carrier_e1=e1,
                                carrier_recomb1=e2,
                                carrier_recomb2=h,
                                continuum_state=ep,
                                energy_recomb_ev=float(dE_recomb),
                                energy_eject_ev=float(eps[ep] - eps[e1]),
                                energy_mismatch_ev=float(dE_mismatch),
                                v_direct_ev=float(v_dir_vec[idx]),
                                v_exchange_ev=float(v_exch_vec[idx]),
                                m_eff_ev=float(np.sqrt(m_sq_vec[idx])),
                                rate_contrib_fs=float(d_rate_fs),
                            )
                        )

    # =========================================================================
    # CHANNEL 2: POSITIVE TRION / BIEXCITON (hhe)
    # Spectator: h1 (valence), Recombining pair: (h2, e), Ejected: h' (deep valence)
    # =========================================================================
    if channel.lower() in ["all", "hhe"]:
        for h1 in initial_holes:
            for h2 in initial_holes:
                for e in initial_elecs:
                    dE_recomb = eps[e] - eps[h2]
                    # Target deep valence energy: eps[h'] = eps[h1] - dE_recomb
                    target_h_prime = eps[h1] - dE_recomb

                    window = e_search_sigma_factor * sigma_ev
                    cand_h_prime = [
                        hp for hp in all_occ_indices
                        if abs(eps[hp] - target_h_prime) <= window and hp not in (h1, h2)
                    ]

                    if not cand_h_prime:
                        continue

                    # 1. q^{e, h2}: recombination charge
                    q_recomb_arr = compute_transition_density_charges(
                        C=C, S=S, atom_ao_ranges=atom_ao_ranges,
                        i_indices=[e], j_indices=[h2],
                        spinor=spinor, U_spinor_alpha=U_spinor_alpha, U_spinor_beta=U_spinor_beta
                    )
                    q_recomb = q_recomb_arr[0, 0, :].real

                    # 2. q^{e, h1}: alternative recombination charge for exchange
                    q_recomb_alt_arr = compute_transition_density_charges(
                        C=C, S=S, atom_ao_ranges=atom_ao_ranges,
                        i_indices=[e], j_indices=[h1],
                        spinor=spinor, U_spinor_alpha=U_spinor_alpha, U_spinor_beta=U_spinor_beta
                    )
                    q_recomb_alt = q_recomb_alt_arr[0, 0, :].real

                    W_q_recomb = W_resta @ q_recomb
                    W_q_recomb_alt = W_resta @ q_recomb_alt

                    # 3. Transition charges from h1 and h2 to all candidate deep valence h'
                    q_eject_1_arr = compute_transition_density_charges(
                        C=C, S=S, atom_ao_ranges=atom_ao_ranges,
                        i_indices=[h1], j_indices=cand_h_prime,
                        spinor=spinor, U_spinor_alpha=U_spinor_alpha, U_spinor_beta=U_spinor_beta
                    )[0, :, :].real  # (len(cand_h_prime), N_atoms)

                    q_eject_2_arr = compute_transition_density_charges(
                        C=C, S=S, atom_ao_ranges=atom_ao_ranges,
                        i_indices=[h2], j_indices=cand_h_prime,
                        spinor=spinor, U_spinor_alpha=U_spinor_alpha, U_spinor_beta=U_spinor_beta
                    )[0, :, :].real

                    v_dir_vec = q_eject_1_arr @ W_q_recomb
                    v_exch_vec = q_eject_2_arr @ W_q_recomb_alt

                    if spinor:
                        m_eff_vec = v_dir_vec - v_exch_vec
                        m_sq_vec = np.abs(m_eff_vec) ** 2
                    else:
                        m_sq_vec = 0.5 * ((v_dir_vec - v_exch_vec) ** 2 + v_dir_vec ** 2)

                    for idx, hp in enumerate(cand_h_prime):
                        dE_mismatch = (eps[h1] - eps[hp]) - dE_recomb
                        rho = eval_lineshape(dE_mismatch)
                        d_rate_fs = TWO_PI_OVER_HBAR * m_sq_vec[idx] * rho
                        rate_hhe_fs += d_rate_fs

                        hhe_details.append(
                            AugerChannelDetail(
                                channel_type="hhe",
                                carrier_e1=h1,
                                carrier_recomb1=h2,
                                carrier_recomb2=e,
                                continuum_state=hp,
                                energy_recomb_ev=float(dE_recomb),
                                energy_eject_ev=float(eps[h1] - eps[hp]),
                                energy_mismatch_ev=float(dE_mismatch),
                                v_direct_ev=float(v_dir_vec[idx]),
                                v_exchange_ev=float(v_exch_vec[idx]),
                                m_eff_ev=float(np.sqrt(m_sq_vec[idx])),
                                rate_contrib_fs=float(d_rate_fs),
                            )
                        )

    # 3. Sort dominant channels by rate contribution
    eeh_details.sort(key=lambda d: d.rate_contrib_fs, reverse=True)
    hhe_details.sort(key=lambda d: d.rate_contrib_fs, reverse=True)

    # 4. Total Biexciton Auger Rate & Lifetimes
    # In a neutral biexciton (2e + 2h): Gamma_XX = 2 * Gamma_eeh + 2 * Gamma_hhe
    rate_xx_fs = 2.0 * rate_eeh_fs + 2.0 * rate_hhe_fs

    rate_eeh_ps = rate_eeh_fs * 1.0e3
    rate_hhe_ps = rate_hhe_fs * 1.0e3
    rate_xx_ps = rate_xx_fs * 1.0e3

    rate_eeh_s = rate_eeh_fs * 1.0e15
    rate_hhe_s = rate_hhe_fs * 1.0e15
    rate_xx_s = rate_xx_fs * 1.0e15

    tau_eeh_ps = 1.0 / rate_eeh_ps if rate_eeh_ps > 0.0 else np.inf
    tau_hhe_ps = 1.0 / rate_hhe_ps if rate_hhe_ps > 0.0 else np.inf
    tau_xx_ps = 1.0 / rate_xx_ps if rate_xx_ps > 0.0 else np.inf

    res = AugerResult(
        rate_eeh_fs=rate_eeh_fs,
        rate_hhe_fs=rate_hhe_fs,
        rate_biexciton_fs=rate_xx_fs,
        rate_eeh_ps=rate_eeh_ps,
        rate_hhe_ps=rate_hhe_ps,
        rate_biexciton_ps=rate_xx_ps,
        rate_eeh_s=rate_eeh_s,
        rate_hhe_s=rate_hhe_s,
        rate_biexciton_s=rate_xx_s,
        tau_eeh_ps=tau_eeh_ps,
        tau_hhe_ps=tau_hhe_ps,
        tau_biexciton_ps=tau_xx_ps,
        fundamental_gap_ev=float(e_gap),
        broadening_mode=broadening_mode,
        sigma_ev=float(sigma_ev),
        n_eeh_pathways=len(eeh_details),
        n_hhe_pathways=len(hhe_details),
        dominant_eeh=eeh_details[:10],
        dominant_hhe=hhe_details[:10],
        computation_time_s=float(time.time() - t0),
    )

    if verbose:
        print(res.summary_table())

    return res


def compute_trajectory_auger_rates(
    step_files: List[str],
    atom_ao_ranges: List[Tuple[int, int]],
    coords: np.ndarray,
    atom_symbols: List[str],
    homo_idx: int,
    material_name: Optional[str] = "DEFAULT",
    sigma_ev: float = 0.05,
    broadening_mode: str = "gaussian",
    max_frames: Optional[int] = None,
    verbose: bool = True,
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Evaluates dynamic Auger recombination rates across an MD trajectory of NAMD step files.

    Parameters
    ----------
    step_files : list of str
        List of paths to precomputed step_*.npz files containing MO eigenvalues and coefficients.
    atom_ao_ranges : list of (start, end) tuples
        AO slice ranges for each atomic center.
    coords : ndarray of shape (N_atoms, 3)
        Atomic coordinates.
    atom_symbols : list of str
        Element symbols.
    homo_idx : int
        0-indexed HOMO position.
    material_name : str, optional
        Material identifier for Resta kernel.
    sigma_ev : float, optional
        Broadening width in eV.
    broadening_mode : str, optional
        'gaussian' or 'fcwd'.
    max_frames : int, optional
        Maximum number of frames to sample along trajectory.
    verbose : bool, optional
        Print progress.

    Returns
    -------
    dict
        Mean and standard deviation of Auger rates and lifetimes.
    """
    from qdex.hardness import build_resta_mnok

    W_resta, _ = build_resta_mnok(
        atom_symbols=atom_symbols,
        coords=coords,
        alpha=0.0,
        material_name=material_name or "DEFAULT",
    )

    files_to_process = step_files[:max_frames] if max_frames else step_files
    n_frames = len(files_to_process)
    if n_frames == 0:
        return {}

    rates_xx_ps = np.zeros(n_frames, dtype=np.float64)
    rates_eeh_ps = np.zeros(n_frames, dtype=np.float64)
    rates_hhe_ps = np.zeros(n_frames, dtype=np.float64)

    if verbose:
        print(f"\nEvaluating Trajectory-Averaged Auger Recombination across {n_frames} frames...")

    for k, sfile in enumerate(files_to_process):
        data = np.load(sfile, allow_pickle=True)
        # Extract MO eigenvalues and coefficients if cached
        if "C_mo" in data and "eps_mo" in data and "overlap" in data:
            C = data["C_mo"]
            eps = data["eps_mo"]
            S = data["overlap"]
        elif "eigenvalues" in data:
            eps = data["eigenvalues"]
            C = data.get("C", None)
            S = data.get("S", None)
        else:
            continue

        if C is None or S is None:
            continue

        res = calculate_auger_rates(
            C=C, eps=eps, S=S,
            atom_ao_ranges=atom_ao_ranges,
            coords=coords, atom_symbols=atom_symbols,
            homo_idx=homo_idx,
            W_resta=W_resta,
            sigma_ev=sigma_ev,
            broadening_mode=broadening_mode,
            verbose=False,
        )
        rates_xx_ps[k] = res.rate_biexciton_ps
        rates_eeh_ps[k] = res.rate_eeh_ps
        rates_hhe_ps[k] = res.rate_hhe_ps

    mean_xx_ps = float(np.mean(rates_xx_ps))
    std_xx_ps = float(np.std(rates_xx_ps))
    mean_tau_xx_ps = float(1.0 / mean_xx_ps) if mean_xx_ps > 0 else np.inf

    summary = {
        "mean_rate_xx_ps": mean_xx_ps,
        "std_rate_xx_ps": std_xx_ps,
        "mean_tau_xx_ps": mean_tau_xx_ps,
        "mean_rate_eeh_ps": float(np.mean(rates_eeh_ps)),
        "mean_rate_hhe_ps": float(np.mean(rates_hhe_ps)),
        "rates_xx_ps": rates_xx_ps,
    }

    if verbose:
        print(f"  Trajectory-Averaged Biexciton Auger Rate : {mean_xx_ps:.3e} +/- {std_xx_ps:.3e} ps^-1")
        print(f"  Trajectory-Averaged Biexciton Lifetime  : {mean_tau_xx_ps:.2f} ps")

    return summary
