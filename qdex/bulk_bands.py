import os
import re
import logging
import numpy as np

logger = logging.getLogger(__name__)

DEFAULT_BULK_DIR = os.path.join(os.path.dirname(__file__), "data", "bulk_bands")


def parse_cp2k_bs(filepath):
    """
    Parses a CP2K band structure (.bs) output file.
    Stitches multiple k-point branches seamlessly into a single array (n_k, n_bands).
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"CP2K band structure file not found: {filepath}")

    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()

    sets = []
    current_set = None
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        m_set = re.match(r"^# Set\s+(\d+):\s+(\d+)\s+special points,\s+(\d+)\s+k-points,\s+(\d+)\s+bands", line)
        if m_set:
            current_set = {"points": [], "n_bands": int(m_set.group(4))}
            sets.append(current_set)
        elif line.startswith("#  Point"):
            parts = line.split()
            pt_idx = int(parts[2])
            kx, ky, kz = float(parts[5]), float(parts[6]), float(parts[7])
            i += 1  # skip header
            bands_e = []
            bands_occ = []
            for _ in range(current_set["n_bands"]):
                i += 1
                bparts = lines[i].split()
                bands_e.append(float(bparts[1]))
                bands_occ.append(float(bparts[2]))
            current_set["points"].append({
                "point": pt_idx,
                "k": (kx, ky, kz),
                "energies": np.array(bands_e, dtype=np.float64),
                "occ": np.array(bands_occ, dtype=np.float64),
            })
        i += 1

    if not sets:
        raise ValueError(f"No band data found in {filepath}")

    all_bands = []
    all_k = []
    all_occ = []
    for s_idx, s in enumerate(sets):
        pts = s["points"]
        # Skip the first point of subsequent branches because CP2K duplicates the boundary k-point
        pts_to_add = pts if s_idx == 0 else pts[1:]
        for p in pts_to_add:
            all_bands.append(p["energies"])
            all_k.append(p["k"])
            all_occ.append(p["occ"])

    bands = np.array(all_bands, dtype=np.float64)  # (n_k, n_bands)
    occ = np.array(all_occ, dtype=np.float64)
    n_k, n_bands = bands.shape

    # Determine VBM, CBM, and band gap
    # Find occupied vs unoccupied bands from occupations
    mean_occ = np.mean(occ, axis=0)
    occ_band_indices = np.where(mean_occ > 0.5)[0]
    virt_band_indices = np.where(mean_occ <= 0.5)[0]

    if len(occ_band_indices) > 0 and len(virt_band_indices) > 0:
        vbm_band = occ_band_indices[-1]
        cbm_band = virt_band_indices[0]
        vbm = float(np.max(bands[:, vbm_band]))
        cbm = float(np.min(bands[:, cbm_band]))
        direct_gap_gamma = float(bands[0, cbm_band] - bands[0, vbm_band])
        gap = float(cbm - vbm)
    else:
        vbm_band = n_bands // 2 - 1
        cbm_band = n_bands // 2
        vbm = float(np.max(bands[:, vbm_band]))
        cbm = float(np.min(bands[:, cbm_band]))
        direct_gap_gamma = float(bands[0, cbm_band] - bands[0, vbm_band])
        gap = float(cbm - vbm)

    midgap = 0.5 * (vbm + cbm)

    return {
        "bands": bands,
        "occ": occ,
        "k_points": all_k,
        "n_k": n_k,
        "n_bands": n_bands,
        "vbm_band": int(vbm_band),
        "cbm_band": int(cbm_band),
        "vbm": vbm,
        "cbm": cbm,
        "gap": gap,
        "direct_gap_gamma": direct_gap_gamma,
        "midgap": midgap,
    }


def find_bulk_bs(material="CdSe", custom_path=None):
    """
    Finds CP2K .bs file for a given material in data/bulk_bands or custom path.
    """
    if custom_path and os.path.exists(custom_path):
        return custom_path

    if not material:
        return None

    mat_clean = str(material).strip().upper()
    candidates = [
        f"{mat_clean}_bulk.bs",
        f"{mat_clean}.bs",
        f"{mat_clean.lower()}_bulk.bs",
        f"{mat_clean.lower()}.bs",
        f"CdSe_bulk.bs" if "CDSE" in mat_clean else None,
    ]

    for cand in candidates:
        if not cand:
            continue
        p = os.path.join(DEFAULT_BULK_DIR, cand)
        if os.path.exists(p):
            return p

    return None


def get_aligned_bulk_bands(
    material="CdSe",
    bs_path=None,
    alignment_mode="core_level",
    ewin=(-5.0, 5.0),
    qd_homo_rel=None,
    qd_lumo_rel=None,
    core_vbm_rel=None,
    core_cbm_rel=None,
    qd_mo_abs_ev=None,
    qd_mo_rel_ev=None,
    qd_elem_fraction=None,
    elements=None,
    qd_h5_path=None,
):
    """
    Loads bulk bands and aligns them to the QD energy axis for overlay.

    Parameters:
        material: str
            Material name (e.g. 'CdSe').
        bs_path: str, optional
            Path to CP2K .bs file.
        alignment_mode: str
            'core_level' (default): Align bulk band edges to the nanocrystal using
                the chemically inert, deep semicore manifold (e.g. Cd 4d for CdSe).
                This core-level anchor is immune to surface traps, quantum confinement,
                and asymmetric band shifts.
            'core_vbm': Align bulk VBM with the QD's True Core VBM.
            'vbm': Align bulk VBM with nominal QD HOMO.
            'midgap': Align bulk midgap with QD midgap (E = 0).
        ewin: tuple of (float, float)
            Energy window in eV to retain bands.
        qd_homo_rel: float, optional
            Nominal QD HOMO energy relative to midgap.
        qd_lumo_rel: float, optional
            Nominal QD LUMO energy relative to midgap.
        core_vbm_rel: float, optional
            True Core VBM energy relative to midgap.
        core_cbm_rel: float, optional
            True Core CBM energy relative to midgap.
        qd_mo_abs_ev: array-like, optional
            Absolute DFT MO eigenvalues (eV) of the nanocrystal.
        qd_mo_rel_ev: array-like, optional
            Relative DFT MO eigenvalues (eV) of the nanocrystal (midgap at 0).
        qd_elem_fraction: ndarray, optional
            (n_mo, n_elements) element fraction per MO.
        elements: list of str, optional
            Element symbols corresponding to qd_elem_fraction columns.
        qd_h5_path: str, optional
            Path to qdex_electronic.h5 to load MO eigenvalues and element fractions
            if not directly provided.

    Returns:
        dict with:
            'bands_aligned': (n_k, n_selected_bands)
            'band_indices': list of original band indices
            'vbm_aligned': float
            'cbm_aligned': float
            'gap': float
            'n_k': int
            'shift': float, energy shift applied to bulk eigenvalues
            'source_file': str
            'alignment_mode': str
            'anchor_detail': dict or str
    """
    path = find_bulk_bs(material=material, custom_path=bs_path)
    if not path:
        return None

    data = parse_cp2k_bs(path)
    bands_raw = data["bands"]
    vbm_raw = data["vbm"]
    cbm_raw = data["cbm"]
    mid_raw = data["midgap"]

    anchor_detail = {}
    shift = None

    # Attempt core_level (semicore) alignment
    if alignment_mode in ("core_level", "semicore", "auto"):
        # Load from HDF5 if arrays not explicitly provided
        if qd_mo_abs_ev is None and qd_h5_path and os.path.exists(qd_h5_path):
            try:
                import h5py
                with h5py.File(qd_h5_path, "r") as f:
                    if "sf/mo/energy_dft_abs_ev" in f:
                        qd_mo_abs_ev = f["sf/mo/energy_dft_abs_ev"][:]
                    if "sf/mo/energy_ev" in f:
                        qd_mo_rel_ev = f["sf/mo/energy_ev"][:]
                    if "sf/mo/element_fraction" in f:
                        qd_elem_fraction = f["sf/mo/element_fraction"][:]
                    if "sf/mo/elements" in f:
                        elements = [e.decode() if isinstance(e, bytes) else str(e) for e in f["sf/mo/elements"][:]]
            except Exception as e_h5:
                logger.debug(f"Failed to read MO data for core_level alignment from HDF5: {e_h5}")

        # Core-level anchoring for CdSe: Cd 4d semicore manifold (bulk bands 2-6)
        mat_upper = str(material).strip().upper()
        if "CDSE" in mat_upper and qd_mo_abs_ev is not None and qd_mo_rel_ev is not None:
            # 1. Bulk Cd 4d reference (bands 2 through 6 in CdSe CP2K calculation, 0-indexed: 1..5)
            # Bands 1..5 are the narrow 4d manifold
            bulk_d_bands = bands_raw[:, 1:6]
            bulk_d_mean = float(np.mean(bulk_d_bands))
            bulk_vbm_to_d = vbm_raw - bulk_d_mean  # invariant separation: ~7.587 eV

            # 2. Nanocrystal Cd 4d semicore identification
            # Cd 4d states lie between -16.0 eV and -13.0 eV in absolute energy scale
            if qd_elem_fraction is not None and elements and "Cd" in elements:
                cd_idx = elements.index("Cd")
                mask_d = (
                    (qd_mo_abs_ev >= -16.5)
                    & (qd_mo_abs_ev <= -13.0)
                    & (qd_elem_fraction[:, cd_idx] > 0.85)
                )
            else:
                mask_d = (qd_mo_abs_ev >= -16.0) & (qd_mo_abs_ev <= -13.5)

            if np.sum(mask_d) >= 10:
                qd_d_mean = float(np.mean(qd_mo_abs_ev[mask_d]))
                # Predicted bulk VBM from Cd 4d semicore core level
                pred_bulk_vbm_abs = qd_d_mean + bulk_vbm_to_d
                midgap_offset = float(qd_mo_abs_ev[0] - qd_mo_rel_ev[0])
                pred_bulk_vbm_rel_semicore = pred_bulk_vbm_abs - midgap_offset

                # 3. Dense DOS Cross-Correlation (Deep Valence Manifold)
                # Compute continuous DOS of bulk valence bands (Bands 2-9) vs QD core-occupied MOs
                dos_corr_delta = 0.0
                dos_corr_r = None
                try:
                    all_bulk_vals = bands_raw[:, 1:9].flatten() - bulk_d_mean
                    sigma_dos = 0.15
                    dE_dos = 0.02
                    grid_dos = np.arange(-1.5, 6.5, dE_dos)
                    bulk_dos = np.sum(np.exp(-0.5 * ((grid_dos[:, None] - all_bulk_vals[None, :]) / sigma_dos) ** 2), axis=1)
                    b_sub = bulk_dos - np.mean(bulk_dos)

                    # QD core weights
                    if qd_elem_fraction is not None and elements:
                        cd_idx = elements.index("Cd") if "Cd" in elements else None
                        se_idx = elements.index("Se") if "Se" in elements else None
                        core_w = np.zeros(len(qd_mo_abs_ev))
                        if cd_idx is not None:
                            core_w += qd_elem_fraction[:, cd_idx]
                        if se_idx is not None:
                            core_w += qd_elem_fraction[:, se_idx]
                    else:
                        core_w = np.ones(len(qd_mo_abs_ev))

                    mask_deep = (qd_mo_abs_ev >= -16.5) & (qd_mo_abs_ev <= -6.0)
                    qd_vals_shifted = qd_mo_abs_ev[mask_deep] - qd_d_mean
                    qd_w_deep = core_w[mask_deep]

                    offsets = np.arange(-0.35, 0.35, 0.005)
                    corrs = []
                    norm_b = np.linalg.norm(b_sub)
                    for off in offsets:
                        q_dos = np.sum(qd_w_deep[None, :] * np.exp(-0.5 * ((grid_dos[:, None] - (qd_vals_shifted[None, :] + off)) / sigma_dos) ** 2), axis=1)
                        q_dos_sub = q_dos - np.mean(q_dos)
                        norm_q = np.linalg.norm(q_dos_sub)
                        c = np.dot(b_sub, q_dos_sub) / (norm_b * norm_q) if norm_q > 1e-12 else 0.0
                        corrs.append(c)

                    best_idx = int(np.argmax(corrs))
                    dos_corr_delta = float(offsets[best_idx])
                    dos_corr_r = float(corrs[best_idx])
                except Exception as e_dos:
                    logger.debug(f"DOS cross-correlation refinement skipped: {e_dos}")

                # Both anchors computed:
                # Semicore anchor:
                #   pred_bulk_vbm_rel_semicore
                # DOS-refined anchor:
                #   pred_bulk_vbm_rel_dos = pred_bulk_vbm_rel_semicore + dos_corr_delta
                # Hybrid consensus: weight 50% semicore + 50% DOS correlation (or select via alignment_mode)
                if alignment_mode == "semicore":
                    final_vbm_rel = pred_bulk_vbm_rel_semicore
                elif alignment_mode == "dos_anchor":
                    final_vbm_rel = pred_bulk_vbm_rel_semicore + dos_corr_delta
                else:
                    # Consensus dual anchor: average of semicore and deep DOS cross-correlation
                    final_vbm_rel = pred_bulk_vbm_rel_semicore + 0.5 * dos_corr_delta

                shift = final_vbm_rel - vbm_raw

                anchor_detail = {
                    "anchor_type": "Dual (Cd 4d semicore + deep valence DOS correlation)",
                    "n_semicore_mos": int(np.sum(mask_d)),
                    "qd_cd_4d_mean_ev": qd_d_mean,
                    "bulk_cd_4d_mean_ev": bulk_d_mean,
                    "vbm_semicore_rel_ev": pred_bulk_vbm_rel_semicore,
                    "vbm_dos_rel_ev": pred_bulk_vbm_rel_semicore + dos_corr_delta,
                    "dos_corr_delta_ev": dos_corr_delta,
                    "dos_corr_pearson_r": dos_corr_r,
                    "final_anchored_vbm_rel_ev": final_vbm_rel,
                }
                logger.info(
                    f"  [Bulk Bands] Dual Anchoring Consensus: "
                    f"Semicore VBM = {pred_bulk_vbm_rel_semicore:.3f} eV, "
                    f"DOS-Correlation VBM = {pred_bulk_vbm_rel_semicore + dos_corr_delta:.3f} eV (r={dos_corr_r:.3f}) "
                    f"-> Final Anchored Bulk VBM = {final_vbm_rel:.3f} eV relative to midgap."
                )

    # Fallback options
    if shift is None:
        if alignment_mode == "core_vbm" and core_vbm_rel is not None:
            shift = core_vbm_rel - vbm_raw
        elif alignment_mode == "vbm" and qd_homo_rel is not None:
            shift = qd_homo_rel - vbm_raw
        else:
            # Default fallback: midgap alignment (bulk midgap = 0.0)
            shift = -mid_raw

    bands_shifted = bands_raw + shift
    vbm_aligned = vbm_raw + shift
    cbm_aligned = cbm_raw + shift

    # Filter bands that cross or lie within ewin
    e_lo, e_hi = ewin[0] - 0.5, ewin[1] + 0.5
    selected_bands = []
    selected_indices = []

    for b in range(data["n_bands"]):
        b_min = np.min(bands_shifted[:, b])
        b_max = np.max(bands_shifted[:, b])
        if b_max >= e_lo and b_min <= e_hi:
            selected_bands.append(bands_shifted[:, b])
            selected_indices.append(b)

    if not selected_bands:
        return None

    bands_aligned = np.column_stack(selected_bands)

    return {
        "bands_aligned": bands_aligned,
        "band_indices": selected_indices,
        "vbm_aligned": vbm_aligned,
        "cbm_aligned": cbm_aligned,
        "gap": data["gap"],
        "n_k": data["n_k"],
        "shift": shift,
        "source_file": path,
        "alignment_mode": alignment_mode,
        "anchor_detail": anchor_detail,
    }
