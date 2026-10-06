import numpy as np
import logging

logger = logging.getLogger(__name__)

def compute_spin_character(vec, soc_U, n_occ_sp, n_virt_sp, valid_mask=None):
    """
    Projects the Spinor BSE eigenvector back onto the spatial basis 
    to calculate the total Singlet and Triplet weights.
    """
    X_IA = np.zeros((n_occ_sp, n_virt_sp), dtype=vec.dtype)
    if valid_mask is not None:
        if valid_mask.ndim == 1:
            X_IA.reshape(-1)[valid_mask] = vec
        else:
            X_IA[valid_mask] = vec
    else:
        X_IA = vec.reshape((n_occ_sp, n_virt_sp))
    
    n_mo = soc_U.shape[0] // 2
    U_occ_a = soc_U[:n_mo, :n_occ_sp]
    U_occ_b = soc_U[n_mo:, :n_occ_sp]
    
    U_virt_a = soc_U[:n_mo, n_occ_sp:]
    U_virt_b = soc_U[n_mo:, n_occ_sp:]
    
    rho_aa = U_occ_a.conj() @ X_IA @ U_virt_a.T
    rho_bb = U_occ_b.conj() @ X_IA @ U_virt_b.T
    
    S_mat = (rho_aa + rho_bb) / np.sqrt(2.0)
    singlet_weight = np.sum(np.abs(S_mat)**2)
    
    singlet_weight = min(1.0, max(0.0, singlet_weight))
    triplet_weight = 1.0 - singlet_weight
    
    return singlet_weight * 100, triplet_weight * 100


def spin_label_from_S(S):
    labels = {
        0: "S",
        1: "T",
        2: "Q",
        3: "7",
        4: "9",
    }
    S_int = int(round(S))
    return labels.get(S_int, f"{2 * S_int + 1}")


def spin_multiplicity_name(S):
    multiplicity = int(round(2.0 * float(S) + 1.0))
    names = {
        1: "singlet",
        2: "doublet",
        3: "triplet",
        4: "quartet",
        5: "quintet",
        6: "sextet",
        7: "septet",
    }
    return names.get(multiplicity, f"{multiplicity}-plet")


def infer_reference_spin(n_alpha, n_beta):
    return abs(float(n_alpha) - float(n_beta)) / 2.0


def compute_uks_soc_spin_character(vec, ham, soc_U, n_alpha_ref, n_beta_ref):
    """
    Approximate UKS-SOC spin-sector character.

    For a closed-shell reference this reduces to the usual singlet/triplet
    alpha/beta transition decomposition. For high-spin UKS references, the
    scalar excitation component is assigned to S_ref and the spin-vector
    component is distributed over S_ref-1, S_ref, S_ref+1 by product-space
    dimensions. This is a labeling diagnostic, not an exact <S^2> projection.
    """
    if not hasattr(ham, "n_occ_act_b"):
        s_pct, t_pct = compute_spin_character(vec, soc_U, ham.n_occ_spinor, ham.n_virt_spinor)
        return {"S": s_pct, "T": t_pct}, 0.0

    X_spinor = np.zeros(ham.dim_spinor_full, dtype=complex)
    if hasattr(ham, "valid_spinor_idx"):
        X_spinor[ham.valid_spinor_idx] = vec
    else:
        X_spinor[:len(vec)] = vec
    X_spinor = X_spinor.reshape(ham.n_occ_spinor, ham.n_virt_spinor)

    n_alpha = ham.n_occ_act + ham.n_virt_act
    occ_cols = slice(0, ham.n_occ_spinor)
    virt_cols = slice(ham.n_occ_spinor, ham.n_occ_spinor + ham.n_virt_spinor)
    U_occ_a = soc_U[0:ham.n_occ_act, occ_cols]
    U_virt_a = soc_U[ham.n_occ_act:n_alpha, virt_cols]
    U_occ_b = soc_U[n_alpha:n_alpha + ham.n_occ_act_b, occ_cols]
    U_virt_b = soc_U[n_alpha + ham.n_occ_act_b:n_alpha + ham.n_occ_act_b + ham.n_virt_act_b, virt_cols]

    X_a = U_occ_a.conj() @ X_spinor @ U_virt_a.T
    X_b = U_occ_b.conj() @ X_spinor @ U_virt_b.T

    n_i = min(X_a.shape[0], X_b.shape[0])
    n_a = min(X_a.shape[1], X_b.shape[1])
    paired_a = X_a[:n_i, :n_a]
    paired_b = X_b[:n_i, :n_a]

    scalar = np.sum(np.abs((paired_a + paired_b) / np.sqrt(2.0)) ** 2)
    vector = np.sum(np.abs((paired_a - paired_b) / np.sqrt(2.0)) ** 2)
    vector += np.sum(np.abs(X_a[n_i:, :]) ** 2) + np.sum(np.abs(X_b[n_i:, :]) ** 2)
    vector += np.sum(np.abs(X_a[:n_i, n_a:]) ** 2) + np.sum(np.abs(X_b[:n_i, n_a:]) ** 2)

    S_ref = infer_reference_spin(n_alpha_ref, n_beta_ref)
    sectors = {}
    if S_ref < 1e-8:
        sectors["S"] = scalar
        sectors["T"] = vector
    else:
        sectors[spin_label_from_S(S_ref)] = sectors.get(spin_label_from_S(S_ref), 0.0) + scalar
        allowed = [S for S in (S_ref - 1.0, S_ref, S_ref + 1.0) if S >= 0.0]
        dim_sum = sum(2.0 * S + 1.0 for S in allowed)
        for S in allowed:
            label = spin_label_from_S(S)
            sectors[label] = sectors.get(label, 0.0) + vector * ((2.0 * S + 1.0) / dim_sum)

    total = sum(sectors.values())
    if total <= 1e-14:
        return {k: 0.0 for k in sectors}, S_ref
    return {k: 100.0 * v / total for k, v in sectors.items()}, S_ref


def format_spin_character(sectors):
    order = ["S", "T", "Q", "7", "9"]
    keys = [k for k in order if k in sectors] + [k for k in sectors if k not in order]
    return " / ".join(f"{sectors[k]:4.1f}% {k}" for k in keys)


def compute_uks_soc_spin_free_channels(vec, ham, soc_U, n_alpha_ref, n_beta_ref):
    """
    Project a UKS-SOC exciton from spinor-transition space back to spin-free
    UKS transition channels.

    Channels are labeled by the spin of the removed occupied electron and the
    spin of the created virtual electron:
      aa: alpha -> alpha, bb: beta -> beta,
      ab: alpha -> beta, ba: beta -> alpha.
    """
    X_spinor = np.zeros(ham.dim_spinor_full, dtype=complex)
    if hasattr(ham, "valid_spinor_idx"):
        X_spinor[ham.valid_spinor_idx] = vec
    else:
        X_spinor[:len(vec)] = vec
    X_spinor = X_spinor.reshape(ham.n_occ_spinor, ham.n_virt_spinor)

    n_alpha = ham.n_occ_act + ham.n_virt_act
    occ_cols = slice(0, ham.n_occ_spinor)
    virt_cols = slice(ham.n_occ_spinor, ham.n_occ_spinor + ham.n_virt_spinor)

    U_occ_a = soc_U[0:ham.n_occ_act, occ_cols]
    U_virt_a = soc_U[ham.n_occ_act:n_alpha, virt_cols]
    U_occ_b = soc_U[n_alpha:n_alpha + ham.n_occ_act_b, occ_cols]
    U_virt_b = soc_U[n_alpha + ham.n_occ_act_b:n_alpha + ham.n_occ_act_b + ham.n_virt_act_b, virt_cols]

    X_aa = U_occ_a.conj() @ X_spinor @ U_virt_a.T
    X_ab = U_occ_a.conj() @ X_spinor @ U_virt_b.T
    X_ba = U_occ_b.conj() @ X_spinor @ U_virt_a.T
    X_bb = U_occ_b.conj() @ X_spinor @ U_virt_b.T

    weights = {
        "aa": float(np.sum(np.abs(X_aa) ** 2)),
        "bb": float(np.sum(np.abs(X_bb) ** 2)),
        "ab": float(np.sum(np.abs(X_ab) ** 2)),
        "ba": float(np.sum(np.abs(X_ba) ** 2)),
    }
    total = sum(weights.values())
    if total <= 1e-14:
        total = 1.0
    channels = {k: 100.0 * v / total for k, v in weights.items()}
    delta_ms = {
        "-1": channels["ab"],
        "0": channels["aa"] + channels["bb"],
        "+1": channels["ba"],
    }

    S_ref = infer_reference_spin(n_alpha_ref, n_beta_ref)
    closed_shell = None
    if S_ref < 1e-8:
        n_i = min(X_aa.shape[0], X_bb.shape[0])
        n_a = min(X_aa.shape[1], X_bb.shape[1])
        paired_a = X_aa[:n_i, :n_a]
        paired_b = X_bb[:n_i, :n_a]
        singlet = float(np.sum(np.abs((paired_a + paired_b) / np.sqrt(2.0)) ** 2))
        triplet = float(np.sum(np.abs((paired_a - paired_b) / np.sqrt(2.0)) ** 2))
        triplet += weights["ab"] + weights["ba"]
        st_total = singlet + triplet
        if st_total <= 1e-14:
            st_total = 1.0
        closed_shell = {"S": 100.0 * singlet / st_total, "T": 100.0 * triplet / st_total}

    return {
        "S_ref": S_ref,
        "ref_label": spin_label_from_S(S_ref),
        "ref_name": spin_multiplicity_name(S_ref),
        "channels": channels,
        "delta_ms": delta_ms,
        "closed_shell": closed_shell,
    }


def compute_uks_spin_free_channels(vec, ham, n_alpha_ref, n_beta_ref):
    """
    Decompose a non-SOC UKS BSE eigenvector into alpha-alpha and beta-beta
    spin-free transition channels.

    This is a channel diagnostic. For open-shell UKS references it should not
    be read as an exact many-electron S^2 decomposition.
    """
    X_a = np.zeros((ham.n_occ_act, ham.n_virt_act), dtype=complex)
    X_b = np.zeros((ham.n_occ_act_b, ham.n_virt_act_b), dtype=complex)

    if ham.dim_a:
        X_a[ham.vi_a, ham.va_a] = vec[:ham.dim_a]
    if ham.dim_b:
        X_b[ham.vi_b, ham.va_b] = vec[ham.dim_a:ham.dim_a + ham.dim_b]

    weights = {
        "aa": float(np.sum(np.abs(X_a) ** 2)),
        "bb": float(np.sum(np.abs(X_b) ** 2)),
    }
    total = sum(weights.values())
    if total <= 1e-14:
        total = 1.0
    channels = {k: 100.0 * v / total for k, v in weights.items()}

    S_ref = infer_reference_spin(n_alpha_ref, n_beta_ref)
    closed_shell = None
    if S_ref < 1e-8:
        n_i = min(X_a.shape[0], X_b.shape[0])
        n_a = min(X_a.shape[1], X_b.shape[1])
        paired_a = X_a[:n_i, :n_a]
        paired_b = X_b[:n_i, :n_a]
        singlet = float(np.sum(np.abs((paired_a + paired_b) / np.sqrt(2.0)) ** 2))
        triplet = float(np.sum(np.abs((paired_a - paired_b) / np.sqrt(2.0)) ** 2))
        st_total = singlet + triplet
        if st_total <= 1e-14:
            st_total = 1.0
        closed_shell = {"S": 100.0 * singlet / st_total, "T": 100.0 * triplet / st_total}

    return {
        "S_ref": S_ref,
        "ref_label": spin_label_from_S(S_ref),
        "ref_name": spin_multiplicity_name(S_ref),
        "channels": channels,
        "delta_ms": {"0": 100.0},
        "closed_shell": closed_shell,
    }


def format_uks_soc_spin_free_character(character):
    ch = character["channels"]
    dm = character["delta_ms"]
    closed_shell = character.get("closed_shell")
    ref = f"ref {character.get('ref_label', '?')}(S={character['S_ref']:.1f})"
    if closed_shell is not None:
        return (
            f"{ref} | S {closed_shell['S']:4.1f}% / T {closed_shell['T']:4.1f}% | "
            f"aa {ch['aa']:4.1f} bb {ch['bb']:4.1f} ab {ch['ab']:4.1f} ba {ch['ba']:4.1f}"
        )
    return (
        f"{ref} | dMs0 {dm['0']:4.1f}% -1 {dm['-1']:4.1f}% +1 {dm['+1']:4.1f}% | "
        f"aa {ch['aa']:4.1f} bb {ch['bb']:4.1f} ab {ch['ab']:4.1f} ba {ch['ba']:4.1f}"
    )


def format_uks_spin_free_character(character):
    ch = character["channels"]
    closed_shell = character.get("closed_shell")
    ref = f"ref {character.get('ref_label', '?')}(S={character['S_ref']:.1f})"
    if closed_shell is not None:
        return (
            f"{ref} | S {closed_shell['S']:4.1f}% / T {closed_shell['T']:4.1f}% | "
            f"aa {ch['aa']:4.1f} bb {ch['bb']:4.1f}"
        )
    return f"{ref} | dMs0 100.0% | aa {ch['aa']:4.1f} bb {ch['bb']:4.1f}"


def _normalize_population_tags(population_bars, n_atoms):
    if not population_bars:
        return []

    atom_index_base = int(population_bars.get("atom_index_base", 1))
    tagged_atoms = population_bars.get("tagged_atoms", [])
    normalized = []
    used_atoms = set()

    for entry in tagged_atoms:
        if not isinstance(entry, dict):
            continue

        raw_indices = entry.get("indices", None)
        if raw_indices is None:
            raw_indices = [entry.get("index", None)]
        elif np.isscalar(raw_indices):
            raw_indices = [raw_indices]

        atom_indices = []
        for raw in raw_indices:
            if raw is None:
                continue
            atom_idx = int(raw) - atom_index_base
            if atom_idx < 0 or atom_idx >= n_atoms:
                continue
            if atom_idx in used_atoms:
                continue
            used_atoms.add(atom_idx)
            atom_indices.append(atom_idx)

        if not atom_indices:
            continue

        normalized.append({
            "label": entry.get("label", f"atom[{atom_indices[0] + atom_index_base}]"),
            "atom_indices": atom_indices,
        })

    return normalized


def compute_orbital_centroids(
    pops, shells, coords_ang, syms,
    core_elements=None,
    xi_core=0.35,
    xi_trap=0.55,
    f_core_min=0.60,
    f_core_trap=0.40,
):
    """
    Computes spatial centroid <r> = sum_A p_A R_A, distance from QD core center d_COM,
    normalized centroid displacement xi = d_COM / R_core, radius of gyration sigma,
    and localization classification ([Core], [Mixed], [Surf/Trap]) for each orbital.

    Parameters:
        pops: np.ndarray, shape (n_ao, n_states)
            Mulliken populations on AOs.
        shells: list of dict
            Shell metadata with 'l' and 'atom_idx'.
        coords_ang: np.ndarray or list of shape (n_atoms, 3)
            Cartesian coordinates of atoms in Angstrom.
        syms: list of str
            Element symbols of all atoms.
        core_elements: list or set of str, optional
            Core elements (e.g. ['Cd', 'Se']). If None, auto-detected by excluding
            typical ligand/capping elements (Cl, Br, I, F, H, C, N, O, P).
        xi_core: float
            Maximum xi for delocalized core envelope states (default 0.35).
        xi_trap: float
            Minimum xi for surface-localized / trap states (default 0.55).
        f_core_min: float
            Minimum core population fraction for [Core] classification (default 0.60).
        f_core_trap: float
            Maximum core population fraction below which state is [Surf/Trap] (default 0.40).

    Returns:
        dict containing:
            'r_cent': (n_states, 3), centroid coordinates
            'd_com': (n_states,), distance from core center in Angstrom
            'xi': (n_states,), normalized centroid displacement d_com / R_core
            'sigma': (n_states,), RMS orbital spread in Angstrom
            'sigma_tilde': (n_states,), normalized spread sigma / R_core
            'f_core': (n_states,), core population fraction
            'labels': list of str, '[Core]', '[Mixed]', or '[Surf/Trap]'
            'com_core': (3,), core center of mass coordinates
            'r_core': float, core radius in Angstrom
            'core_mask': (n_atoms,) boolean array
    """
    coords_ang = np.asarray(coords_ang, dtype=np.float64)
    n_atoms = len(syms)
    n_states = pops.shape[1]

    # Map AOs to atoms
    ao_atom_indices = []
    for sh in shells:
        l_int = int(sh['l'])
        atom_idx = int(sh.get("atom_idx", 0))
        nbf = 2 * l_int + 1
        ao_atom_indices.extend([atom_idx] * nbf)

    atom_pops = np.zeros((n_atoms, n_states), dtype=np.float64)
    np.add.at(atom_pops, ao_atom_indices, pops)
    tot_pop = np.sum(atom_pops, axis=0)
    tot_pop_safe = np.where(tot_pop == 0.0, 1.0, tot_pop)

    # Core elements detection
    if core_elements is not None:
        core_set = {str(e).strip().capitalize() for e in core_elements}
    else:
        # Common organic ligands & passivating halides
        passivants = {"H", "C", "N", "O", "F", "Cl", "Br", "I", "P"}
        candidate_core = {str(s).strip().capitalize() for s in syms if str(s).strip().capitalize() not in passivants}
        if candidate_core:
            core_set = candidate_core
        else:
            core_set = {str(s).strip().capitalize() for s in syms}

    syms_cap = [str(s).strip().capitalize() for s in syms]
    core_mask = np.array([s in core_set for s in syms_cap], dtype=bool)
    if not np.any(core_mask):
        core_mask = np.ones(n_atoms, dtype=bool)

    coords_core = coords_ang[core_mask]
    com_core = np.mean(coords_core, axis=0)
    dists_core = np.linalg.norm(coords_core - com_core, axis=1)
    r_core = float(np.max(dists_core)) if len(dists_core) > 0 else 1.0
    if r_core < 1e-3:
        r_core = 1.0

    # Electronic centroids
    r_cent = (atom_pops.T @ coords_ang) / tot_pop_safe[:, None]
    d_com = np.linalg.norm(r_cent - com_core, axis=1)
    xi = d_com / r_core

    # Orbital spread (radius of gyration)
    d2 = np.sum((coords_ang[None, :, :] - r_cent[:, None, :]) ** 2, axis=2)
    var_r = np.sum((atom_pops.T / tot_pop_safe[:, None]) * d2, axis=1)
    sigma = np.sqrt(np.maximum(0.0, var_r))
    sigma_tilde = sigma / r_core

    # Core fraction
    f_core = np.sum(atom_pops[core_mask, :], axis=0) / tot_pop_safe

    # Classification
    labels = []
    for k in range(n_states):
        if xi[k] <= xi_core and f_core[k] >= f_core_min:
            labels.append("[Core]")
        elif xi[k] >= xi_trap or f_core[k] < f_core_trap:
            labels.append("[Surf/Trap]")
        else:
            labels.append("[Mixed]")

    return {
        "r_cent": r_cent,
        "d_com": d_com,
        "xi": xi,
        "sigma": sigma,
        "sigma_tilde": sigma_tilde,
        "f_core": f_core,
        "labels": labels,
        "com_core": com_core,
        "r_core": r_core,
        "core_mask": core_mask,
    }


def print_orbital_summary(
    energies_eV, occ, homo_idx, pops, syms, shells, is_soc=False, offset=0,
    print_range=15, population_bars=None, qp_breakdown=None, coords_ang=None,
    core_elements=None, xi_core=0.35, xi_trap=0.55, f_core_min=0.60, f_core_trap=0.40
):
    """
    Fast Mulliken population analysis broken down by Element and Angular Momentum (s, p, d).
    Expects precomputed 'pops' matrix to avoid duplicating S @ C multiplications.
    If qp_breakdown is provided, also displays the microscopic breakdown of each orbital's QP shift.
    If coords_ang is provided, computes electronic centroids <r>, displacement xi = d_COM / R_core,
    orbital spread sigma/R, and core/trap classification.
    """
    has_qp = (qp_breakdown is not None) and (not is_soc)
    has_coords = coords_ang is not None

    centroid_info = None
    if has_coords:
        centroid_info = compute_orbital_centroids(
            pops, shells, coords_ang, syms, core_elements=core_elements,
            xi_core=xi_core, xi_trap=xi_trap, f_core_min=f_core_min, f_core_trap=f_core_trap
        )

    if has_coords:
        table_width = 145 if has_qp else 136
    else:
        table_width = 115 if has_qp else 115

    logger.info("\n" + "=" * table_width)
    if has_qp:
        if has_coords:
            logger.info(f"{'Orbital':>14} | {'Index':>6} | {'DFT (eV)':>10} | {'QP (eV)':>10} | {'Occ':>5} | {'d_COM (Å)':>9} | {'ξ':>6} | {'σ/R':>5} | {'Loc':>11} | {'Main Contributions':>40}")
        else:
            logger.info(f"{'Orbital':>14} | {'Index':>6} | {'DFT (eV)':>10} | {'QP (eV)':>10} | {'Occ':>5} | {'Main Contributions':>45}")
    else:
        if has_coords:
            logger.info(f"{'Orbital':>14} | {'Index':>6} | {'Energy (eV)':>12} | {'Occ':>5} | {'d_COM (Å)':>9} | {'ξ':>6} | {'σ/R':>5} | {'Loc':>11} | {'Main Contributions':>40}")
        else:
            logger.info(f"{'Orbital':>14} | {'Index':>6} | {'Energy (eV)':>12} | {'Occ':>5} | {'Main Contributions':>45}")
    logger.info("-" * table_width)
    
    n_states = pops.shape[1]
    l_char = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g'}
    ao_labels = []
    ao_atom_indices = []
    
    for sh in shells:
        sym = sh['sym']
        atom_idx = int(sh.get("atom_idx", 0))
        l_int = int(sh['l'])
        l_str = l_char.get(l_int, str(l_int))
        label = f"{sym}({l_str})"
        nbf = 2 * l_int + 1
        ao_labels.extend([label] * nbf)
        ao_atom_indices.extend([atom_idx] * nbf)
        
    tags = _normalize_population_tags(population_bars, len(syms))
    tagged_atoms = {}
    for tag in tags:
        for atom_idx in tag["atom_indices"]:
            tagged_atoms[atom_idx] = tag["label"]

    unique_labels = list(dict.fromkeys(
        [f"{tag['label']}" for tag in tags] + [
            f"{label}" for label in dict.fromkeys(ao_labels)
        ]
    ))
    label_pops = np.zeros((len(unique_labels), n_states))

    label_to_idx = {label: i for i, label in enumerate(unique_labels)}
    for i_ao, base_label in enumerate(ao_labels):
        atom_idx = ao_atom_indices[i_ao]
        tag_label = tagged_atoms.get(atom_idx)
        label = tag_label if tag_label is not None else base_label
        label_pops[label_to_idx[label], :] += pops[i_ao, :]
        
    col_sums = np.sum(label_pops, axis=0)
    col_sums[col_sums == 0] = 1.0
    label_pops = label_pops / col_sums[np.newaxis, :]
    
    start_idx = max(0, homo_idx - print_range + 1)
    end_idx = min(n_states, homo_idx + 1 + print_range)
    
    for idx in range(end_idx - 1, start_idx - 1, -1):
        if not is_soc:
            rel = idx - homo_idx
            label = ("LUMO" if rel == 1 else f"LUMO+{rel-1}") if rel > 0 else ("HOMO" if rel == 0 else f"HOMO{rel}")
        else:
            rel = idx - homo_idx
            label = ("spL" if rel == 1 else f"spL+{rel-1}") if rel > 0 else ("spH" if rel == 0 else f"spH{rel}")
                
        state_pops = label_pops[:, idx]
        top_indices = np.argsort(-state_pops)[:5]
        contrib_str = ", ".join([f"{unique_labels[i]} ({state_pops[i]*100:.0f}%)" for i in top_indices if state_pops[i] > 0.05])
        
        if has_qp:
            e_dft = float(qp_breakdown["eps_dft"][idx])
            e_qp = float(qp_breakdown["eps_qp"][idx])
            if has_coords:
                d_c = centroid_info["d_com"][idx]
                x_c = centroid_info["xi"][idx]
                s_c = centroid_info["sigma_tilde"][idx]
                loc_c = centroid_info["labels"][idx]
                logger.info(f"{label:>14} | {idx + offset:6d} | {e_dft:10.4f} | {e_qp:10.4f} | {occ[idx]:5.1f} | {d_c:9.2f} | {x_c:6.2f} | {s_c:5.2f} | {loc_c:>11} | {contrib_str}")
            else:
                logger.info(f"{label:>14} | {idx + offset:6d} | {e_dft:10.4f} | {e_qp:10.4f} | {occ[idx]:5.1f} | {contrib_str}")
        else:
            if has_coords:
                d_c = centroid_info["d_com"][idx]
                x_c = centroid_info["xi"][idx]
                s_c = centroid_info["sigma_tilde"][idx]
                loc_c = centroid_info["labels"][idx]
                logger.info(f"{label:>14} | {idx + offset:6d} | {energies_eV[idx]:12.4f} | {occ[idx]:5.1f} | {d_c:9.2f} | {x_c:6.2f} | {s_c:5.2f} | {loc_c:>11} | {contrib_str}")
            else:
                logger.info(f"{label:>14} | {idx + offset:6d} | {energies_eV[idx]:12.4f} | {occ[idx]:5.1f} | {contrib_str}")
        
        if idx == homo_idx + 1:
            if has_qp:
                if has_coords:
                    logger.info(f"   {'-- FERMI --':>11} | {'------':>6} | {'----------':>10} | {'----------':>10} | {'-----':>5} | {'---------':>9} | {'------':>6} | {'-----':>5} | {'-----------':>11} | {'-'*40}")
                else:
                    logger.info(f"   {'-- FERMI --':>11} | {'------':>6} | {'----------':>10} | {'----------':>10} | {'-----':>5} | {'-'*45}")
            else:
                if has_coords:
                    logger.info(f"   {'-- FERMI --':>11} | {'------':>6} | {'------------':>12} | {'-----':>5} | {'---------':>9} | {'------':>6} | {'-----':>5} | {'-----------':>11} | {'-'*40}")
                else:
                    logger.info(f"   {'-- FERMI --':>11} | {'------':>6} | {'------------':>12} | {'-----':>5} | {'-'*45}")
    logger.info("=" * table_width + "\n")
    return centroid_info

