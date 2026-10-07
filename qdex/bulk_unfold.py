"""Bulk bands of a supercell structure unfolded onto the k-path of its primitive (pseudo-cubic) cell.

Tilted lead halide perovskites (Pnma, four formula units) are supercells of the cubic cell: their
bands fold, and the cubic R and M points both land on the orthorhombic Gamma. The fuzzy bands of a
tilted nanocube are drawn on the cubic path, so the bulk reference drawn over them is the
orthorhombic band structure unfolded onto that path (spectral weight, Popescu and Zunger, Phys. Rev.
B 85, 085201 (2012)):

    P_Kn(k) = sum_{g in primitive lattice} |<k+g|psi_Kn>|^2 / sum_{G in supercell lattice} |<k+G|psi_Kn>|^2,

for the supercell eigenstates at K = k. The plane-wave amplitudes <q|psi> = sum_mu c_mu F_mu(q) use the
analytic Fourier transforms of the Gaussian basis (libint_cpp.ao_ft_complex, the fuzzy-band
machinery), with the Bloch coefficients of qdex.bulk_soc.BulkModel (spin-free or two-component). The
primitive lattice is the average pseudo-cubic lattice of the supercell: the nearest B-B (Pb-Pb)
vectors give the integer supercell matrix M, and A_p = M^-1 A_s exactly.

    python -m qdex.bulk_unfold RUN_DIR --basis BASIS_MOLOPT_UZH --gth GTH_SOC_POTENTIALS \
        --path-bs qdex/data/bulk_bands/CsPbBr3_cubic.bs.gz --b-site Pb -o qdex/data/bulk_bands/CsPbBr3_cubic_unfolded.npz
"""
import argparse
import json
import logging
import os

import numpy as np

logger = logging.getLogger(__name__)
HA_EV = 27.211386245988
BOHR_PER_ANG = 1.8897259886


def pseudocubic_cell(cell, syms, coords, b_site="Pb"):
    """Primitive (pseudo-cubic) cell A_p (rows, A) of a supercell and the integer matrix M, A_s = M A_p."""
    cell = np.asarray(cell, float)
    X = np.asarray(coords, float)[[i for i, s in enumerate(syms) if s == b_site]]
    vecs = []
    rng = range(-1, 2)
    for i in range(len(X)):
        for j in range(len(X)):
            for n in ((a, b, c) for a in rng for b in rng for c in rng):
                v = X[j] + np.dot(n, cell) - X[i]
                d = np.linalg.norm(v)
                if d > 1e-3:
                    vecs.append(v)
    vecs = np.array(vecs)
    d = np.linalg.norm(vecs, axis=1)
    nn = vecs[d < 1.15 * d.min()]                      # the six nearest B-B directions
    basis = [nn[0]]
    for v in nn[1:]:                                   # three non-collinear ones
        if all(abs(np.dot(v, b)) / np.linalg.norm(v) / np.linalg.norm(b) < 0.5 for b in basis):
            basis.append(v)
        if len(basis) == 3:
            break
    P = np.array(basis)
    if np.linalg.det(P) < 0:
        P[2] *= -1
    M = np.rint(cell @ np.linalg.inv(P)).astype(int)
    if abs(round(np.linalg.det(M))) < 1:
        raise ValueError("could not find the primitive cell of the supercell")
    A_p = np.linalg.inv(M) @ cell
    return A_p, M


def _g_vectors(B, q_max):
    """Reciprocal lattice vectors (rows of B combined with integers) shorter than q_max (1/A)."""
    n = int(np.ceil(q_max / min(np.linalg.norm(B, axis=1)))) + 1
    r = np.arange(-n, n + 1)
    m = np.array(np.meshgrid(r, r, r, indexing="ij")).reshape(3, -1).T
    G = m @ B
    return G[np.linalg.norm(G, axis=1) < q_max]


def unfold(model, kfrac_prim, A_p, soc=False, q_max=5.0, nthreads=1):
    """Energies (eV) and unfolding weights of the supercell bands at the primitive k-points.

    kfrac_prim: (n_k, 3) fractional coordinates in the reciprocal basis of A_p. Returns (E, P),
    each (n_k, n_bands)."""
    import libint_cpp
    A_s = np.asarray(model.cell, float)
    B_s = 2 * np.pi * np.linalg.inv(A_s).T
    B_p = 2 * np.pi * np.linalg.inv(A_p).T
    G = _g_vectors(B_s, q_max)
    prim = np.all(np.abs((G @ A_p.T) / (2 * np.pi) - np.rint((G @ A_p.T) / (2 * np.pi))) < 1e-6, axis=1)
    n = model.n_ao
    E_out, P_out = [], []
    for kf in np.asarray(kfrac_prim, float):
        k = kf @ B_p
        K_frac_s = k @ A_s.T / (2 * np.pi)
        E, C = (model.soc if soc else model.sf)(K_frac_s, vectors=True)
        F = libint_cpp.ao_ft_complex(model.shells, (k + G) / BOHR_PER_ANG, nthreads)      # (n_ao, n_G)
        if soc:
            w = np.abs(C[:n].T @ F) ** 2 + np.abs(C[n:].T @ F) ** 2                         # (n_b, n_G)
        else:
            w = np.abs(C.T @ F) ** 2
        tot = w.sum(axis=1)
        P = np.divide(w[:, prim].sum(axis=1), tot, out=np.zeros_like(tot), where=tot > 0)
        E_out.append(E * HA_EV)
        P_out.append(P)
    return np.array(E_out), np.array(P_out)


def unfold_to_file(run_dir, basis_file, gth_file, path_bs, out_npz, b_site="Pb", q_max=5.0, e_window=(-6.0, 6.0),
                   basis_name="DZVP-MOLOPT-PBE-GTH", functional="PBE", nthreads=1):
    """Unfold a supercell run onto the k-path of a primitive-cell band file (its segments) and write
    out_npz (+ .json metadata: semicore anchor levels, edges, primitive cell)."""
    from qdex.bulk_bands import parse_cp2k_bs
    from qdex.bulk_soc import BulkModel, choose_semicore
    from qdex.constants import valence_electrons
    model = BulkModel(run_dir, basis_file, gth_file, basis_name, functional)
    A_p, M = pseudocubic_cell(model.cell, model.syms, model.coords, b_site)
    logger.info(f"  [Unfold] supercell = M x primitive, M = {M.tolist()} (det {round(np.linalg.det(M))}); "
                f"primitive |a| = {np.round(np.linalg.norm(A_p, axis=1), 4).tolist()} A")
    segs = parse_cp2k_bs(path_bs)["segments"]
    kf = np.concatenate([np.asarray(s["kfrac"]) for s in segs])
    seg_len = np.array([len(s["kfrac"]) for s in segs])
    n_occ = sum(valence_electrons[s] for s in model.syms) // 2
    out = dict(kfrac=kf, seg_len=seg_len, primitive_cell=A_p, supercell_matrix=M)
    meta = dict(source=os.path.abspath(run_dir), path_bs=os.path.basename(path_bs), b_site=b_site,
                n_occ_sf=int(n_occ), q_max_inv_ang=q_max, supercell_matrix=M.tolist(),
                primitive_a_ang=np.linalg.norm(A_p, axis=1).round(4).tolist())
    for tag, soc in (("sf", False), ("soc", True)):
        E, P = unfold(model, kf, A_p, soc=soc, q_max=q_max, nthreads=nthreads)
        nv = n_occ * (2 if soc else 1)
        vbm, cbm = float(E[:, nv - 1].max()), float(E[:, nv].min())
        keep = np.where((E.max(axis=0) >= vbm + e_window[0] - 8.0) & (E.min(axis=0) <= cbm + e_window[1]))[0]
        out[f"E_{tag}"], out[f"P_{tag}"], out[f"bands_{tag}"] = E[:, keep], P[:, keep], keep
        meta[tag] = dict(vbm=vbm, cbm=cbm, gap=cbm - vbm)
        logger.info(f"  [Unfold] {tag}: gap {cbm - vbm:.4f} eV; mean weight of the edge bands "
                    f"{P[:, nv - 1].mean():.2f} / {P[:, nv].mean():.2f}")
    # semicore anchors: the manifolds of the primitive-cell file (Cs 5p, anion s), each the n bands with the
    # largest weight on it below the valence band (not necessarily contiguous: in the tilted cells other
    # bands can fall inside the manifold), n = atoms x (2l + 1); the level is their mean energy along the path
    from qdex.bulk_bands import load_bulk_meta
    g = np.zeros(3)
    Eg, labels, frac, _ = model.characters(g)
    counts = {e: model.syms.count(e) for e in set(model.syms)}
    prim_meta = load_bulk_meta(path_bs) or {}
    wanted = [c["label"] for c in (prim_meta.get("semicore_candidates") or [])] or \
        [c["label"] for c in choose_semicore(Eg, labels, frac, counts, meta["sf"]["vbm"])]
    cands = []
    for lab in wanted:
        el, l = lab.split("-")[0], "spdf".index(lab.split("-")[1])
        if lab not in labels or el not in counts:
            continue
        nb = counts[el] * (2 * l + 1)
        w = frac[:, labels.index(lab)]
        below = np.where(Eg < meta["sf"]["vbm"] - 2.0)[0]
        best = np.sort(below[np.argsort(-w[below])[:nb]])
        if len(best) < nb or w[best].min() < 0.5:
            continue
        rows_sf = [list(out["bands_sf"]).index(b) for b in best if b in out["bands_sf"]]
        rows_soc = [list(out["bands_soc"]).index(b) for b in np.concatenate([2 * best, 2 * best + 1])
                    if b in out["bands_soc"]]
        if len(rows_sf) < nb or len(rows_soc) < 2 * nb:
            continue
        cands.append(dict(element=el, l=l, label=lab, bands_sf=best.tolist(), min_fraction=float(w[best].min()),
                          level_sf=float(out["E_sf"][:, rows_sf].mean()), level_soc=float(out["E_soc"][:, rows_soc].mean())))
    meta["semicore_candidates"] = cands
    meta["semicore"] = cands[0] if cands else None
    logger.info(f"  [Unfold] semicore anchors: {[(c['label'], round(c['level_sf'], 3)) for c in cands]}")
    np.savez_compressed(out_npz, **out)
    with open(os.path.splitext(out_npz)[0] + ".json", "w") as f:
        json.dump(meta, f, indent=1)
    return out, meta


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--basis", required=True)
    ap.add_argument("--gth", required=True)
    ap.add_argument("--path-bs", required=True, help="band file of the primitive cell whose k-path is used")
    ap.add_argument("--b-site", default="Pb")
    ap.add_argument("--q-max", type=float, default=5.0, help="plane-wave cutoff |k+G| (1/A)")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--nthreads", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    unfold_to_file(a.run_dir, a.basis, a.gth, a.path_bs, a.out, a.b_site, a.q_max, nthreads=a.nthreads)


if __name__ == "__main__":
    main()
