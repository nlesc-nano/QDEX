import os
import re
import logging
import numpy as np

logger = logging.getLogger(__name__)

DEFAULT_BULK_DIR = os.path.join(os.path.dirname(__file__), "data", "bulk_bands")

# Each material in data/bulk_bands has <name>.bs (spin-free), <name>_soc.bs (spin-orbit) and
# <name>.json (band edges, gaps, Gamma band characters and the semicore manifold), written by
# qdex.bulk_soc; <name> is the stem of the CIF the fuzzy bands are computed with (e.g. CdSe_zb).
STRUCTURE_PREFERENCE = ("zb", "rs", "cubic", "wz")


def parse_cp2k_bs(filepath):
    """
    Parses a CP2K band structure (.bs) output file.

    Returns every k-point set as a segment (fractional k-points and bands), and the sets
    stitched into one array (n_k, n_bands). A set's first point is dropped from the stitched
    array only when it repeats the previous set's last point (not across a path break).
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"CP2K band structure file not found: {filepath}")

    if str(filepath).endswith(".gz"):
        import gzip
        with gzip.open(filepath, "rt", encoding="utf-8") as f:
            lines = f.readlines()
    else:
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

    segments = []
    all_bands, all_k, all_occ = [], [], []
    for s in sets:
        pts = s["points"]
        kfrac = np.array([p["k"] for p in pts], dtype=np.float64)
        segments.append({
            "kfrac": kfrac,
            "bands": np.array([p["energies"] for p in pts], dtype=np.float64),
        })
        skip_first = bool(all_k) and np.allclose(kfrac[0], all_k[-1], atol=1e-6)
        for p in (pts[1:] if skip_first else pts):
            all_bands.append(p["energies"])
            all_k.append(p["k"])
            all_occ.append(p["occ"])

    bands = np.array(all_bands, dtype=np.float64)  # (n_k, n_bands)
    occ = np.array(all_occ, dtype=np.float64)
    n_k, n_bands = bands.shape

    # Occupied vs unoccupied bands from occupations
    mean_occ = np.mean(occ, axis=0)
    occ_band_indices = np.where(mean_occ > 0.5)[0]
    virt_band_indices = np.where(mean_occ <= 0.5)[0]

    if len(occ_band_indices) > 0 and len(virt_band_indices) > 0:
        vbm_band = occ_band_indices[-1]
        cbm_band = virt_band_indices[0]
    else:
        vbm_band = n_bands // 2 - 1
        cbm_band = n_bands // 2
    vbm = float(np.max(bands[:, vbm_band]))
    cbm = float(np.min(bands[:, cbm_band]))
    direct_gap_gamma = float(bands[0, cbm_band] - bands[0, vbm_band])
    gap = float(cbm - vbm)

    return {
        "bands": bands,
        "occ": occ,
        "k_points": all_k,
        "segments": segments,
        "n_k": n_k,
        "n_bands": n_bands,
        "vbm_band": int(vbm_band),
        "cbm_band": int(cbm_band),
        "vbm": vbm,
        "cbm": cbm,
        "gap": gap,
        "direct_gap_gamma": direct_gap_gamma,
        "midgap": 0.5 * (vbm + cbm),
    }


def _stem(path):
    return os.path.splitext(os.path.basename(str(path)))[0]


def find_bulk_bs(material="CdSe", custom_path=None, soc=False, cif=None):
    """
    Bulk band file for a run: the bands of its CIF (data/bulk_bands/<cif stem>[_soc].bs) or, failing
    that, of the material (<Material>_<structure>[_soc].bs, preferring zb, rs, cubic, wz).
    custom_path: an explicit spin-free .bs file.
    """
    tag = "_soc" if soc else ""
    if custom_path and not soc and os.path.exists(custom_path):
        return custom_path
    if cif:
        for ext in (".bs", ".bs.gz"):
            p = os.path.join(DEFAULT_BULK_DIR, f"{_stem(cif)}{tag}{ext}")
            if os.path.exists(p):
                return p
    if not material:
        return None
    mat = str(material).strip().upper()
    found = {}
    for f in os.listdir(DEFAULT_BULK_DIR):
        g = f[:-3] if f.endswith(".gz") else f
        if not g.endswith(f"{tag}.bs") or (not soc and g.endswith("_soc.bs")):
            continue
        name = g[: -len(f"{tag}.bs")]
        parts = name.split("_")
        if parts[0].upper() == mat:
            found[parts[1].lower() if len(parts) > 1 else ""] = os.path.join(DEFAULT_BULK_DIR, f)
    for key in STRUCTURE_PREFERENCE:
        if key in found:
            return found[key]
    return next(iter(sorted(found.values())), None)


def load_bulk_meta(bs_path):
    """The <name>.json written by qdex.bulk_soc next to a band file, or None."""
    import json
    if not bs_path:
        return None
    base = str(bs_path)
    base = base[:-3] if base.endswith(".gz") else base
    base = base[:-3]
    base = base[:-4] if base.endswith("_soc") else base
    path = base + ".json"
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def _semicore_candidates(meta):
    meta = meta or {}
    cands = meta.get("semicore_candidates")
    if cands is None:
        cands = [meta["semicore"]] if meta.get("semicore") else []
    return cands


def bulk_semicore_level(meta, bs_data, spinor=False, label=None):
    """Mean energy of the bulk semicore bands along the k-path (bulk eV scale), or None.
    spinor: the bands are spin-orbit bands (two per spin-free band); label: the manifold
    (e.g. 'Cd-d'), by default the preferred one."""
    cands = _semicore_candidates(meta)
    semi = next((c for c in cands if label is None or c["label"] == label), None)
    if not semi:
        return None
    lo, hi = semi["bands_soc" if spinor else "bands_sf"]
    return float(np.mean(bs_data["bands"][:, lo:hi]))


def qd_semicore_level(material, C, S, shells, syms, coords_ang, energies, homo_index,
                      coordination, bond_ang, bs_path=None, window_ev=2.5, cif=None, anchor="auto", bond_pairs=None):
    """
    Semicore (e.g. Cd 4d) level of the dot's bulk-like interior atoms, on the axis of `energies`.

    For each atom of the semicore element, the level is the Mulliken-weighted mean energy of its
    semicore-l population over the orbitals near the expected semicore energy (bulk VBM - semicore
    separation below the HOMO, +- window_ev). Atoms bonded to ligands carry a chemical shift
    (Cd 4d of Cd-Cl lies ~0.4 eV deeper than in Cd-Se), so only bulk-like atoms are used: full
    crystal coordination by the crystal's own elements, and among those the inner half by radius.

    The manifolds of the material's metadata are tried in order of preference (cation d, Cs 5p,
    anion s); the first one whose energy range is covered by the orbitals is used (MO files often
    hold only the orbitals near the gap). With anchor='auto', manifolds of elements outside the
    bonded framework (coordination 0 in the crystal bond star: the perovskite A-site Cs) are tried
    last: their level follows the size and shape of the cage, which the octahedral tilts change
    (Cs 5p lies 0.6-0.9 eV closer to the VBM in orthorhombic than in cubic CsPbX3, and drifts by
    0.6 eV from the core to the surface of a dot), while the halide s level is flat to 0.1 eV.
    Any other anchor value is the label of the manifold to use (e.g. 'Br-s', 'Cs-p').

    bond_pairs: element pairs bonded in the crystal; the bulk-like test then counts only bonded
    neighbours (a Cs next to a perovskite halide does not disqualify it).

    Returns a dict (level_ev, label, n_atoms, spread_ev, all_atoms_level_ev, n_bulk_like) or None.
    """
    path = find_bulk_bs(material=material, custom_path=bs_path, cif=cif)
    meta = load_bulk_meta(path)
    cands = _semicore_candidates(meta)
    if path is None or not cands:
        return None
    if anchor and anchor != "auto":
        cands = [c for c in cands if c["label"].lower() == str(anchor).lower()]
        if not cands:
            logger.warning(f"  [Bulk Bands] No '{anchor}' manifold in the bulk data of {material}: no semicore anchor.")
            return None
    else:   # framework manifolds first, then the A-site cation (coordination 0)
        cands = sorted(cands, key=lambda c: coordination.get(c["element"], 1) == 0)
    bs = parse_cp2k_bs(path)
    energies = np.asarray(energies, dtype=float)
    cfg = None
    for c in cands:
        separation = bs["vbm"] - bulk_semicore_level(meta, bs, label=c["label"])
        target = energies[homo_index] - separation
        if energies.min() <= target - 1.0:
            cfg = c
            break
        logger.info(f"  [Bulk Bands] {c['label']} manifold ({separation:.1f} eV below the VBM) is below the "
                    f"lowest orbital ({energies[homo_index] - energies.min():.1f} eV below the HOMO).")
    if cfg is None:
        return None

    syms = list(syms)
    coords = np.asarray(coords_ang, dtype=float)
    sel = np.where(np.abs(energies - target) <= window_ev)[0]
    elem = cfg["element"]
    atoms = np.array([i for i, s in enumerate(syms) if s == elem])
    if sel.size == 0 or atoms.size == 0:
        return None

    # semicore-l AO rows of every atom of the element
    rows = {int(a): [] for a in atoms}
    ao = 0
    for sh in shells:
        n = 2 * int(sh["l"]) + 1
        if int(sh["l"]) == int(cfg["l"]) and int(sh["atom_idx"]) in rows:
            rows[int(sh["atom_idx"])].extend(range(ao, ao + n))
        ao += n

    C_sel = np.asarray(C[:, sel])
    P = C_sel * (S @ C_sel)                      # Mulliken AO populations (n_ao, n_sel)
    levels = np.full(len(atoms), np.nan)
    for i, a in enumerate(atoms):
        w = P[rows[int(a)]].sum(axis=0)
        if w.sum() > 1e-8:
            levels[i] = float(w @ energies[sel] / w.sum())

    # bulk-like: fully coordinated, only by the crystal's own elements
    from scipy.spatial import cKDTree
    crystal = set(coordination)
    tree = cKDTree(coords)
    n_full = coordination.get(elem, 4)
    bulk_like = []
    for i, a in enumerate(atoms):
        nbrs = [j for j in tree.query_ball_point(coords[a], 1.25 * bond_ang) if j != a]
        bonded = nbrs if bond_pairs is None else [j for j in nbrs if (elem, syms[j]) in bond_pairs]
        if len(bonded) == n_full and all(syms[j] in crystal for j in nbrs) and np.isfinite(levels[i]):
            bulk_like.append(i)
    bulk_like = np.array(bulk_like, dtype=int)

    finite = np.isfinite(levels)
    result = dict(all_atoms_level_ev=float(np.mean(levels[finite])), n_bulk_like=int(bulk_like.size))
    if bulk_like.size == 0:
        logger.warning(f"  [Bulk Bands] No bulk-like {elem} atoms: semicore anchor uses every {elem} atom.")
        use = np.where(finite)[0]
    else:
        r = np.linalg.norm(coords[atoms[bulk_like]] - coords[atoms].mean(axis=0), axis=1)
        inner = bulk_like[r <= np.median(r)]
        use = inner if inner.size >= 4 else bulk_like
    bulk_vbm_rel = float(np.mean(levels[use])) + separation
    if bulk_vbm_rel < energies[homo_index] - 0.2:
        logger.warning(f"  [Bulk Bands] The anchored bulk VBM ({bulk_vbm_rel:.2f} eV) lies "
                       f"{energies[homo_index] - bulk_vbm_rel:.2f} eV below the dot's HOMO; a confined hole lies "
                       f"below the bulk VBM. Was the dot computed at another level of theory (functional, basis, "
                       f"pseudopotential) than the bulk bands (PBE, DZVP-MOLOPT-PBE-GTH)?")
    result["bulk_vbm_rel_ev"] = bulk_vbm_rel
    result.update(level_ev=float(np.mean(levels[use])), spread_ev=float(np.std(levels[use])),
                  n_atoms=int(use.size), element=elem, label=cfg.get("label", elem), bulk_file=os.path.basename(path))
    logger.info(f"  [Bulk Bands] {cfg.get('label', elem)} semicore level: {result['level_ev']:.3f} eV from {use.size} interior bulk-like "
                f"atoms (spread {result['spread_ev']:.3f} eV); all {finite.sum()} {elem} atoms: "
                f"{result['all_atoms_level_ev']:.3f} eV.")
    return result


def map_bulk_to_path(segments, path_frac, atol=1e-4):
    """
    Positions of the bulk k-points on the fuzzy path axis (k-point index).

    Each bulk segment is matched by its end points to the straight stretch of the fuzzy path
    that joins them; points in between are placed by their fractional position along it.
    Returns a list of (x, bands) per matched segment.
    """
    path_frac = np.asarray(path_frac, dtype=float)
    mapped = []
    for seg in segments:
        ks, ke = seg["kfrac"][0], seg["kfrac"][-1]
        length = np.linalg.norm(ke - ks)
        if length < atol:
            continue
        starts = np.where(np.linalg.norm(path_frac - ks, axis=1) < atol)[0]
        ends = np.where(np.linalg.norm(path_frac - ke, axis=1) < atol)[0]
        best = None
        for i in starts:
            for j in ends:
                lo, hi = min(i, j), max(i, j)
                # a stretch needs interior points: two adjacent path points are a path break (X|M)
                if hi - lo < 2 or (best is not None and hi - lo >= abs(best[1] - best[0])):
                    continue
                inner = path_frac[lo:hi + 1]
                t = (inner - ks) @ (ke - ks) / length ** 2
                off = np.linalg.norm(inner - (ks + np.outer(t, ke - ks)), axis=1)
                if np.all(off < 1e-3) and np.all(np.abs(np.diff(t)) > 0):
                    best = (i, j)
        if best is None:
            continue
        i, j = best
        t = (seg["kfrac"] - ks) @ (ke - ks) / length ** 2
        mapped.append((i + t * (j - i), seg["bands"]))
    return mapped


def find_unfolded(material=None, cif=None):
    """Unfolded supercell bands for a run: data/bulk_bands/<cif stem>_unfolded.npz, else
    <Material>_<structure>_unfolded.npz (qdex.bulk_unfold), or None."""
    if cif:
        p = os.path.join(DEFAULT_BULK_DIR, f"{_stem(cif)}_unfolded.npz")
        if os.path.exists(p):
            return p
    if material:
        mat = str(material).strip().upper()
        for f in sorted(os.listdir(DEFAULT_BULK_DIR)):
            if f.endswith("_unfolded.npz") and f.split("_")[0].upper() == mat:
                return os.path.join(DEFAULT_BULK_DIR, f)
    return None


def get_aligned_unfolded_bands(npz_path, ewin=(-5.0, 5.0), qd_semicore_rel=None, path_frac=None, soc=False,
                               qd_semicore_label=None, min_weight=0.03, alignment_mode="core_level", qd_homo_rel=None):
    """Unfolded bulk bands (qdex.bulk_unfold) on the dot's energy axis, for the overlay.

    Same anchor as get_aligned_bulk_bands (the semicore level of the dot's interior atoms); every
    (k, band) point carries its unfolding weight. Returns the dict of get_aligned_bulk_bands with
    'segments' [(x, bands, weights)] and 'unfolded': True, or None."""
    import json
    d = np.load(npz_path)
    with open(os.path.splitext(npz_path)[0] + ".json") as f:
        meta = json.load(f)
    tag = "soc" if soc else "sf"
    E, P = d[f"E_{tag}"], d[f"P_{tag}"]
    cands = meta.get("semicore_candidates") or []
    semi = next((c for c in cands if qd_semicore_label is None or c["label"] == qd_semicore_label), None)
    if alignment_mode == "vbm" and qd_homo_rel is not None:
        shift, mode = float(qd_homo_rel) - meta[tag]["vbm"], "vbm"
    elif alignment_mode in ("core_level", "semicore", "auto") and semi is not None and qd_semicore_rel is not None:
        shift, mode = float(qd_semicore_rel) - float(semi[f"level_{tag}"]), "core_level"
        logger.info(f"  [Bulk Bands] Unfolded {os.path.basename(npz_path)}: semicore anchor ({semi['label']}), "
                    f"bulk VBM at {meta[tag]['vbm'] + shift:.3f} eV, CBM at {meta[tag]['cbm'] + shift:.3f} eV.")
    else:
        shift, mode = -0.5 * (meta[tag]["vbm"] + meta[tag]["cbm"]), "midgap"
    Es = E + shift
    keep = np.where((Es.max(axis=0) >= ewin[0] - 0.5) & (Es.min(axis=0) <= ewin[1] + 0.5))[0]
    if keep.size == 0:
        return None
    bounds = np.concatenate([[0], np.cumsum(d["seg_len"])])
    segs = [dict(kfrac=d["kfrac"][a:b], bands=np.hstack([Es[a:b][:, keep], P[a:b][:, keep]]))
            for a, b in zip(bounds[:-1], bounds[1:])]
    segments = None
    if path_frac is not None:
        mapped = map_bulk_to_path(segs, path_frac)
        if mapped:
            n = len(keep)
            segments = [(x, b[:, :n], b[:, n:]) for x, b in mapped]
    if segments is None:
        logger.warning("  [Bulk Bands] Unfolded bulk k-path does not match the fuzzy path: overlay skipped.")
        return None
    return {"bands_aligned": Es[:, keep], "segments": segments, "band_indices": keep.tolist(),
            "vbm_aligned": meta[tag]["vbm"] + shift, "cbm_aligned": meta[tag]["cbm"] + shift,
            "gap": meta[tag]["gap"], "n_k": len(E), "shift": shift, "source_file": npz_path,
            "alignment_mode": mode, "anchor_detail": {}, "soc": soc, "unfolded": True,
            "min_weight": min_weight}


def get_aligned_bulk_bands(
    material="CdSe",
    bs_path=None,
    alignment_mode="core_level",
    ewin=(-5.0, 5.0),
    qd_homo_rel=None,
    qd_lumo_rel=None,
    core_vbm_rel=None,
    qd_semicore_rel=None,
    path_frac=None,
    soc=False,
    cif=None,
    qd_semicore_label=None,
    unfolded="auto",
):
    """
    Loads bulk bands and aligns them to the QD energy axis for overlay.

    Parameters:
        material: str
            Material name (e.g. 'CdSe').
        bs_path: str, optional
            Path to CP2K .bs file.
        alignment_mode: str
            'core_level' (default): the bulk semicore bands (e.g. Cd 4d) are placed on the
                semicore level of the dot's interior bulk-like atoms (qd_semicore_rel, from
                qd_semicore_level). Falls back to 'midgap' when that level is not available.
            'core_vbm': Align bulk VBM with the QD's True Core VBM (core_vbm_rel).
            'vbm': Align bulk VBM with nominal QD HOMO.
            'midgap': Align bulk midgap with QD midgap (E = 0).
        ewin: tuple of (float, float)
            Energy window in eV to retain bands.
        qd_homo_rel, qd_lumo_rel: float, optional
            Nominal QD HOMO / LUMO energies on the plot axis.
        core_vbm_rel: float, optional
            True Core VBM energy on the plot axis.
        qd_semicore_rel: float, optional
            Semicore level of the dot's interior atoms on the plot axis.
        path_frac: (n_k, 3) array, optional
            Fractional k-points of the fuzzy path; bulk segments are then placed by coordinates.
        soc: bool
            Use the spin-orbit bulk bands (<name>_soc.bs) when available; 'soc' in the returned
            dict says whether they were found (otherwise the spin-free bands are used).
        cif: str, optional
            CIF of the run; its stem selects the band files (e.g. CdSe_wz.cif -> CdSe_wz.bs).
        unfolded: 'auto' (default) or 'off'
            With 'auto', unfolded bands of the measured supercell structure (<cif stem>_unfolded.npz,
            qdex.bulk_unfold; the tilted perovskites) replace the primitive-cell bands when present.
        qd_semicore_label: str, optional
            Manifold qd_semicore_rel was measured on (e.g. 'Cd-d'; default: the preferred one).

    Returns:
        dict with 'bands_aligned' (n_k, n_selected), 'segments' [(x, bands_selected)] or None,
        'band_indices', 'vbm_aligned', 'cbm_aligned', 'gap', 'n_k', 'shift', 'source_file',
        'alignment_mode' (the mode actually used), 'anchor_detail'.
    """
    if unfolded and unfolded != "off":
        # tilted perovskites: the measured (orthorhombic) structure unfolded onto the cubic path
        u = find_unfolded(material=material, cif=cif)
        if u:
            res = get_aligned_unfolded_bands(u, ewin=ewin, qd_semicore_rel=qd_semicore_rel, path_frac=path_frac,
                                             soc=soc, qd_semicore_label=qd_semicore_label,
                                             alignment_mode=alignment_mode, qd_homo_rel=qd_homo_rel)
            if res is not None:
                return res
    path = find_bulk_bs(material=material, soc=True, cif=cif) if soc else None
    soc_bands = path is not None
    if not soc_bands:
        path = find_bulk_bs(material=material, custom_path=bs_path, cif=cif)
    if not path:
        return None

    data = parse_cp2k_bs(path)
    bands_raw = data["bands"]
    vbm_raw = data["vbm"]
    cbm_raw = data["cbm"]

    anchor_detail = {}
    shift = None
    mode_used = alignment_mode

    if alignment_mode in ("core_level", "semicore", "auto"):
        bulk_level = bulk_semicore_level(load_bulk_meta(path), data, spinor=soc_bands, label=qd_semicore_label)
        if bulk_level is not None and qd_semicore_rel is not None:
            shift = float(qd_semicore_rel) - bulk_level
            mode_used = "core_level"
            anchor_detail = {
                "anchor_type": "semicore level of interior bulk-like atoms",
                "bulk_semicore_ev": bulk_level,
                "qd_semicore_rel_ev": float(qd_semicore_rel),
                "vbm_aligned_rel_ev": vbm_raw + shift,
            }
            logger.info(f"  [Bulk Bands] Semicore anchor: bulk VBM at {vbm_raw + shift:.3f} eV, "
                        f"CBM at {cbm_raw + shift:.3f} eV on the plot axis.")
        else:
            logger.info("  [Bulk Bands] No semicore level for this material/run: bulk midgap aligned to E = 0.")

    if shift is None:
        if alignment_mode == "core_vbm" and core_vbm_rel is not None:
            shift = core_vbm_rel - vbm_raw
        elif alignment_mode == "vbm" and qd_homo_rel is not None:
            shift = qd_homo_rel - vbm_raw
        else:
            shift = -data["midgap"]
            mode_used = "midgap"

    bands_shifted = bands_raw + shift

    # Bands that cross or lie within ewin
    e_lo, e_hi = ewin[0] - 0.5, ewin[1] + 0.5
    selected_indices = [b for b in range(data["n_bands"])
                        if np.max(bands_shifted[:, b]) >= e_lo and np.min(bands_shifted[:, b]) <= e_hi]
    if not selected_indices:
        return None

    segments = None
    if path_frac is not None:
        mapped = map_bulk_to_path(data["segments"], path_frac)
        if mapped:
            segments = [(x, b[:, selected_indices] + shift) for x, b in mapped]
            if len(mapped) < len(data["segments"]):
                logger.info(f"  [Bulk Bands] {len(mapped)}/{len(data['segments'])} bulk k-segments lie on the fuzzy path.")
        else:
            logger.warning("  [Bulk Bands] Bulk k-path does not match the fuzzy path: bulk bands stretched by index.")

    return {
        "bands_aligned": bands_shifted[:, selected_indices],
        "segments": segments,
        "band_indices": selected_indices,
        "vbm_aligned": vbm_raw + shift,
        "cbm_aligned": cbm_raw + shift,
        "gap": data["gap"],
        "n_k": data["n_k"],
        "shift": shift,
        "source_file": path,
        "alignment_mode": mode_used,
        "anchor_detail": anchor_detail,
        "soc": soc_bands,
    }
