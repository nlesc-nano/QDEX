"""
Machine-readable results for databases: two HDF5 files per run.

  qdex_electronic.h5    structure, QP correction, spin-free MOs and SOC spinors with
                        their element / shell / surface projections, IPR, COOP,
                        PDOS on a grid and the raw fuzzy-band weights
  qdex_excitations.h5   every exciton up to `output.excitations_emax` (exciton energy,
                        after the BSE): energy, oscillator strength, transition
                        dipole, hole/electron orbitals, energy decomposition and
                        Dreuw-Plasser descriptors, spin-free and SOC

Switched on with `output: {h5: true}`. Arrays are collected while QDEX runs
(`ResultStore.put`) and written at the end (`ResultStore.write`), so a failed
run leaves no half-written file. Energies are in eV, lengths in Å, transition
dipoles in atomic units (e bohr). Eigenvectors are not stored: in the diagonal
(independent-transition) modes each state is one hole-electron pair, and any
state can be recomputed from the orbitals.
"""
from __future__ import annotations

import json
import logging
import subprocess
import time
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

FILES = {"electronic": "qdex_electronic.h5", "excitations": "qdex_excitations.h5"}
L_LABELS = "spdfghi"


class ResultStore:
    def __init__(self):
        self.data = {name: {} for name in FILES}      # file -> {path: array}
        self.attrs = {name: {} for name in FILES}     # file -> {group path: {key: value}}
        self.cache = {}                               # in-run helpers, never written

    def put(self, file: str, path: str, value) -> None:
        self.data[file][path] = value

    def attr(self, file: str, group: str, **values) -> None:
        self.attrs[file].setdefault(group, {}).update(values)

    def write(self, args, out_dir: str = ".") -> list:
        import h5py
        written = []
        header = run_header(args)
        for name, fname in FILES.items():
            if not self.data[name] and not self.attrs[name]:
                continue
            path = Path(out_dir) / fname
            tmp = path.with_suffix(".h5.tmp")
            with h5py.File(tmp, "w") as h5:
                for k, v in header.items():
                    h5.attrs[k] = v
                for p, value in sorted(self.data[name].items()):
                    _dataset(h5, p, value)
                for group, values in self.attrs[name].items():
                    g = h5.require_group(group) if group not in ("", "/") else h5
                    for k, v in values.items():
                        g.attrs[k] = _attr_value(v)
            tmp.replace(path)
            written.append(str(path))
            logger.info(f"  [Store] Wrote {path} ({path.stat().st_size / 1e6:.2f} MB)")
        return written


def _dataset(h5, path, value):
    if value is None:
        return
    if isinstance(value, (list, tuple)) and value and isinstance(value[0], str):
        h5.create_dataset(path, data=np.array(value, dtype=object), dtype=h5py_str())
        return
    arr = np.asarray(value)
    if arr.dtype.kind in "OU":
        h5.create_dataset(path, data=arr.astype(object), dtype=h5py_str())
        return
    if arr.ndim == 0:
        h5.create_dataset(path, data=arr)
        return
    if arr.dtype == np.float64 and arr.size > 64:
        arr = arr.astype(np.float32) if path.endswith("intensity") else arr
    kw = {"compression": "gzip", "compression_opts": 4, "shuffle": True} if arr.size > 64 else {}
    h5.create_dataset(path, data=arr, **kw)


def h5py_str():
    import h5py
    return h5py.string_dtype(encoding="utf-8")


def _attr_value(v):
    if isinstance(v, (dict, list)) and not (isinstance(v, list) and v and isinstance(v[0], (int, float))):
        return json.dumps(v, default=_json_default)
    if v is None:
        return "none"
    if isinstance(v, np.generic):
        return v.item()
    return v


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.generic):
        return o.item()
    return str(o)


def run_header(args) -> dict:
    """Version, commit and the run settings, for every file."""
    here = Path(__file__).resolve().parent
    try:
        commit = subprocess.run(["git", "-C", str(here), "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, timeout=5).stdout.strip() or "unknown"
        dirty = bool(subprocess.run(["git", "-C", str(here), "status", "--porcelain", "--untracked-files=no"],
                                    capture_output=True, text=True, timeout=5).stdout.strip())
    except Exception:
        commit, dirty = "unknown", False
    rev = here.parent / "REVISION"              # source exports without .git carry their commit here
    if commit == "unknown" and rev.is_file():
        commit = rev.read_text().strip() or "unknown"
    settings = {k: v for k, v in vars(args).items()
                if not k.startswith("_") and isinstance(v, (str, int, float, bool, list, type(None)))}
    return {"program": "QDEX", "commit": commit, "dirty": dirty,
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "settings": json.dumps(settings, default=_json_default),
            "units": "energies eV, lengths Angstrom, transition dipoles e*bohr"}


# ---------------------------------------------------------------------------
# Orbital projections
# ---------------------------------------------------------------------------

def ao_labels(shells):
    """Element and element(l) label of every AO, in AO order."""
    el, lab = [], []
    for sh in shells:
        n = 2 * int(sh["l"]) + 1
        sym = sh.get("sym", "X")
        el += [sym] * n
        lab += [f"{sym}({L_LABELS[int(sh['l'])]})"] * n
    return np.array(el), np.array(lab)


def projections(P, shells, surface_mask):
    """
    Per-orbital fractions from Mulliken AO populations P (n_ao, n_states):
    element, element(l) and surface; normalised per orbital.
    """
    el, lab = ao_labels(shells)
    tot = P.sum(axis=0)
    tot = np.where(np.abs(tot) > 1e-12, tot, 1.0)
    elements = list(dict.fromkeys(el))
    shells_l = list(dict.fromkeys(lab))
    f_el = np.stack([P[el == e].sum(axis=0) / tot for e in elements], axis=1)
    f_lab = np.stack([P[lab == x].sum(axis=0) / tot for x in shells_l], axis=1)
    f_surf = P[np.asarray(surface_mask, bool)].sum(axis=0) / tot
    return elements, f_el, shells_l, f_lab, f_surf


def pdos_grid(energies, fractions, ewin, sigma, spin_factor, n=1000):
    """Gaussian-broadened PDOS per column of `fractions` (not stacked)."""
    grid = np.linspace(ewin[0], ewin[1], n)
    g = np.exp(-0.5 * ((grid[:, None] - np.asarray(energies)[None, :]) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    return grid, spin_factor * g @ np.asarray(fractions)


def put_orbitals(store, group, energies, occupation, analysis, shells, coop_pairs, ewin, pdos_sigma,
                 spin_factor, indices=None, extra=None):
    """One orbital set (spin-free MOs or SOC spinors) with projections, COOP and PDOS."""
    s = store
    s.put("electronic", f"{group}/energy_ev", np.asarray(energies, float))
    s.put("electronic", f"{group}/occupation", np.asarray(occupation, float))
    if indices is not None:
        s.put("electronic", f"{group}/index", np.asarray(indices, int))
    for k, v in (extra or {}).items():
        s.put("electronic", f"{group}/{k}", v)
    if analysis is None:
        return
    elements, f_el, labels, f_lab, f_surf = projections(analysis["P_weights"], shells, analysis["surface_ao_mask"])
    s.put("electronic", f"{group}/elements", elements)
    s.put("electronic", f"{group}/element_fraction", f_el)
    s.put("electronic", f"{group}/shells", labels)
    s.put("electronic", f"{group}/shell_fraction", f_lab)
    s.put("electronic", f"{group}/surface_fraction", f_surf)
    s.put("electronic", f"{group}/ipr", np.asarray(analysis["IPR"], float))
    pairs = [p for p in coop_pairs if p in analysis["coop_results"]]
    if pairs:
        s.put("electronic", f"{group}/coop_pairs", pairs)
        s.put("electronic", f"{group}/coop", np.stack([analysis["coop_results"][p] for p in pairs], axis=1))
    grid, pd = pdos_grid(energies, f_el, ewin, pdos_sigma, spin_factor)
    s.put("electronic", f"{group}/pdos/energy_ev", grid)
    s.put("electronic", f"{group}/pdos/elements", elements)
    s.put("electronic", f"{group}/pdos/pdos", pd)
    s.attr("electronic", f"{group}/pdos", sigma_ev=float(pdos_sigma), spin_factor=float(spin_factor),
           note="Gaussian-broadened, per element, not stacked; energies as in ../energy_ev")
    s.attr("electronic", group, surface_definition="AOs at >= 0.75 R_max from the centroid",
           ipr="sum over AOs of the squared Mulliken populations",
           coop="2 sum_{A-B} C_A S_AB C_B per orbital (non-zero within ewin +- 1 eV)")


def put_fuzzy(store, group, kpts, labels, energies, intensity, sigma, ewin, indices=None):
    """Raw fuzzy-band weights |<phi_n|k>|^2 along the k-path (smear with `sigma` to plot)."""
    ticks = [(i, str(l).replace("\\Gamma", "Γ").replace("GAMMA", "Γ").replace("$", ""))
             for i, l in enumerate(labels) if l]
    s = store
    s.put("electronic", f"{group}/fuzzy/kpoints_inv_ang", np.asarray(kpts, float))
    s.put("electronic", f"{group}/fuzzy/tick_index", np.array([t[0] for t in ticks], int))
    s.put("electronic", f"{group}/fuzzy/tick_label", [t[1] for t in ticks] or [""])
    s.put("electronic", f"{group}/fuzzy/energy_ev", np.asarray(energies, float))
    s.put("electronic", f"{group}/fuzzy/intensity", np.asarray(intensity, np.float32))
    if indices is not None:
        s.put("electronic", f"{group}/fuzzy/index", np.asarray(indices, int))
    s.attr("electronic", f"{group}/fuzzy", sigma_ev=float(sigma), ewin_ev=list(map(float, ewin)),
           note="intensity[state, k] = |<phi_state|k>|^2; plot sum_state intensity * Gauss(E - energy, sigma)")


# ---------------------------------------------------------------------------
# Excitations
# ---------------------------------------------------------------------------

def _orbital_moments(q, coords):
    """Centroid (n, 3) and variance (n,) of normalised atomic populations q (n_atoms, n)."""
    q = q / (q.sum(axis=0, keepdims=True) + 1e-12)
    r = q.T @ coords
    d2 = ((coords[None, :, :] - r[:, None, :]) ** 2).sum(axis=2)       # (n, n_atoms)
    return r, (q.T * d2).sum(axis=1)


def _atom_pops(P, atom_ao_ranges):
    return np.stack([P[a:b].sum(axis=0) for a, b in atom_ao_ranges], axis=0)


def put_excitations(store, solver, energies_ev, vectors, f_strengths, mu_ia, analyzer, analysis_results,
                    args, suffix, soc_U=None):
    """Every state up to excitations_emax (at least excitations_min_states) of one solver run."""
    group = "soc" if solver.soc_flag else "sf"
    emax = float(getattr(args, "excitations_emax", 4.5))
    n_min = int(getattr(args, "excitations_min_states", 20))
    n_all = len(energies_ev)
    n = max(int(np.sum(np.asarray(energies_ev) <= emax)), min(n_min, n_all))
    ham = solver.ham
    coords = np.asarray(analyzer.coords, float)
    E = np.asarray(energies_ev[:n], float)
    out = {"energy_ev": E, "f_osc": np.asarray(f_strengths[:n], float)}
    diagonal = hasattr(vectors, "order")
    if n_all and E[-1] < emax and not diagonal and n == n_all:
        logger.warning(f"  [Store] {group}: the {n_all} computed states end at {E[-1]:.3f} eV < "
                       f"excitations_emax {emax} eV; raise nroots to cover it")
    if diagonal:
        p = vectors.order[:n]
        mu = np.asarray(mu_ia)[p]
        out["d_qp_ev"] = np.real(np.asarray(ham.D)[p])
        if getattr(solver, "diagonal_kx", None) is not None:
            out["kx_ev"] = np.real(np.asarray(solver.diagonal_kx)[p])
            out["minus_kd_ev"] = -np.real(np.asarray(solver.diagonal_kd)[p])
        out["pr"] = np.ones(n)
        if not solver.soc_flag:
            i, a = np.asarray(ham.valid_i)[p], np.asarray(ham.valid_a)[p]
            out["hole_mo"] = solver.homo_index - ham.n_occ_act + 1 + i
            out["electron_mo"] = solver.homo_index + 1 + a
            r_o, v_o = _orbital_moments(analyzer.q_occ, coords)
            r_v, v_v = _orbital_moments(analyzer.q_virt, coords)
            rh, re, vh, ve = r_o[i], r_v[a], v_o[i], v_v[a]
        else:
            full = np.asarray(ham.valid_spinor_idx)[p] if hasattr(ham, "valid_spinor_idx") else p
            h, e = full // ham.n_virt_spinor, full % ham.n_virt_spinor
            out["hole_spinor"] = h
            out["electron_spinor"] = ham.n_occ_spinor + e
            q_sp = store.cache.get("soc_spinor_atom_pops")
            if q_sp is not None:
                r_s, v_s = _orbital_moments(q_sp, coords)
                rh, re = r_s[h], r_s[ham.n_occ_spinor + e]
                vh, ve = v_s[h], v_s[ham.n_occ_spinor + e]
            else:
                rh = re = vh = ve = None
            n_mo = soc_U.shape[0] // 2
            n_os = ham.n_occ_spinor
            Uoa, Uob, Uva, Uvb = soc_U[:n_mo, :n_os], soc_U[n_mo:, :n_os], soc_U[:n_mo, n_os:], soc_U[n_mo:, n_os:]
            x = np.sum(Uoa.conj() * Uob, axis=0)         # per occupied spinor
            y = np.sum(Uva * Uvb.conj(), axis=0)         # per virtual spinor
            na, nb = np.sum(np.abs(Uoa) ** 2, 0), np.sum(np.abs(Uob) ** 2, 0)
            ma, mb = np.sum(np.abs(Uva) ** 2, 0), np.sum(np.abs(Uvb) ** 2, 0)
            sw = 0.5 * (na[h] * ma[e] + nb[h] * mb[e] + 2.0 * np.real(x[h] * y[e]))
            out["singlet_fraction"] = np.clip(sw, 0.0, 1.0)
        if rh is not None:
            d_ct = np.linalg.norm(re - rh, axis=1)
            d_eh = np.sqrt(vh + ve + d_ct ** 2)            # covariance 0 for a single pair
            out.update(d_ct_ang=d_ct, d_eh_ang=d_eh, sigma_h_ang=np.sqrt(vh), sigma_e_ang=np.sqrt(ve),
                       ct_character=d_ct / (d_eh + 1e-6))
    else:
        V = vectors[:, :n] if not isinstance(vectors, np.ndarray) else vectors[:, :n]
        V = np.asarray(V)
        mu = (V.T[:, :, None] * np.asarray(mu_ia)[None, :, :]).sum(axis=1)
        w = np.abs(V) ** 2
        out["pr"] = 1.0 / np.sum((w / w.sum(axis=0, keepdims=True)) ** 2, axis=0)
        m = min(n, len(analysis_results))
        for key, src in (("d_ct_ang", "d_CT"), ("d_eh_ang", "d_eh"), ("sigma_h_ang", "sigma_h"),
                         ("sigma_e_ang", "sigma_e"), ("ct_character", "CT_Character")):
            col = np.full(n, np.nan)
            col[:m] = [analysis_results[k][src] for k in range(m)]
            out[key] = col
    if np.iscomplexobj(mu):
        out["mu_re_au"], out["mu_im_au"] = np.real(mu), np.imag(mu)
    else:
        out["mu_au"] = np.asarray(mu, float)
    if "d_eh_ang" in out:
        ct, deh = out["ct_character"], out["d_eh_ang"]
        out["type"] = np.where(ct > 0.6, "CT", np.where((deh < 3.0) & (ct < 0.2), "Frenkel", "Wannier"))
    for k, v in out.items():
        store.put("excitations", f"{group}/{k}", v)
    # spectrum of the stored states, on a fixed grid up to emax + 0.5 eV
    from qdex.spectrum import generate_spectrum
    x, yv = generate_spectrum(E, out["f_osc"], e_min=0.0, e_max=max(emax, E.max() if n else emax) + 0.5,
                              sigma=args.sigma, profile=args.broadening if args.broadening != "none" else "gaussian")
    store.put("excitations", f"{group}/spectrum/energy_ev", x)
    store.put("excitations", f"{group}/spectrum/intensity", yv)
    store.attr("excitations", f"{group}/spectrum", sigma_ev=float(args.sigma), profile=str(args.broadening))
    store.attr("excitations", group, n_states_stored=int(n), n_states_computed=int(n_all),
               n_transitions=int(ham.dim), excitations_emax_ev=emax, mode=str(args.excitation_mode),
               kernel=str(args.kernel),
               indices=("hole_mo / electron_mo: MO index in the orbital file" if not solver.soc_flag else
                        "hole_spinor / electron_spinor: index into soc/bse_spinor/energy_ev of qdex_electronic.h5"),
               descriptors=("diagonal mode: one hole-electron pair per state, so PR = 1 and the e-h covariance "
                            "is 0; SOC descriptors from the spinors' own Mulliken populations"
                            if diagonal else "Dreuw-Plasser analysis of the first 100 states (NaN beyond)"))
