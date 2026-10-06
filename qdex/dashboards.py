"""
Rebuild QDEX's HTML dashboards from a finished run's HDF5 files (output.h5: true),
without recomputing anything:

    qdex dashboards <run dir> [-o <out dir>] [--keep-data]

writes fuzzy_dashboard_sf.html and fuzzy_dashboard_soc.html (fuzzy bands, PDOS,
IPR, surface/core, COOP, and the MO cubes found in the run directory),
exciton_analysis_sf.html and exciton_analysis_soc.html, and the spectrum plots.
The intermediate CSV / NPZ files the plotting code reads are recreated from
qdex_electronic.h5 and qdex_excitations.h5 in a scratch folder and removed
afterwards (kept with --keep-data). Database runs store only the HDF5 files; the
pages for a webapp are made from them whenever they are needed.
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)
N_ANALYSIS = 100          # states on the exciton analysis page (as in a QDEX run)


def _str(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def _attrs(group):
    out = {}
    for k, v in group.attrs.items():
        out[k] = v.item() if isinstance(v, np.generic) else v
    return out


def _orbital_files(g, prefix, ewin, sigma, spin_factor):
    """ipr / surf_core / pdos / coop CSVs of one orbital group, in QDEX's formats."""
    E = g["energy_ev"][:]
    mask = (E >= ewin[0]) & (E <= ewin[1])
    with open(f"ipr_data_{prefix}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["MO_Energy_eV", "IPR"])
        w.writerows(zip(E[mask], g["ipr"][:][mask]))
    surf = g["surface_fraction"][:]
    with open(f"surf_core_data_{prefix}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["MO_Energy_eV", "Surface", "Core"])
        w.writerows(zip(E[mask], surf[mask], 1.0 - surf[mask]))
    elements = [_str(e) for e in g["elements"][:]]
    frac = g["element_fraction"][:]
    grid = np.linspace(ewin[0], ewin[1], 1000)
    G = np.exp(-0.5 * ((grid[:, None] - E[None, :]) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    Ycum = np.cumsum(spin_factor * G @ frac, axis=1)
    with open(f"pdos_data_{prefix}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Energy_eV"] + elements)
        for e, row in zip(grid, Ycum):
            w.writerow([e] + list(row))
    if "coop" in g:
        pairs = [_str(p) for p in g["coop_pairs"][:]]
        coop = g["coop"][:]
        with open(f"coop_data_{prefix}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["MO_Energy_eV"] + pairs)
            for e, row in zip(E[mask], coop[mask]):
                w.writerow([e] + list(row))
    return E


def _fuzzy_file(fz, prefix):
    from qdex.fuzzy_bands import smear_and_export_fuzzy
    n_k = fz["kpoints_inv_ang"].shape[0]
    labels = [""] * n_k
    for i, lab in zip(fz["tick_index"][:], fz["tick_label"][:]):
        labels[int(i)] = _str(lab)
    a = _attrs(fz)
    ewin = [float(x) for x in (a["ewin_ev"] if not isinstance(a["ewin_ev"], str) else
                               a["ewin_ev"].strip("[]").split())]
    smear_and_export_fuzzy(fz["intensity"][:].astype(float), fz["energy_ev"][:], labels, ewin,
                           float(a["sigma_ev"]), prefix=prefix)
    return ewin


def _edges(E, occupation=None):
    if occupation is not None:
        occ = np.flatnonzero(np.asarray(occupation) > 0.5)
        h = int(occ.max())
        return float(E[h]), float(E[h + 1])
    return float(E[E <= 0].max()), float(E[E > 0].min())


def fuzzy_dashboards(el, run_dir: Path, out: Path, cubes: bool = True):
    from qdex.plot_fuzzy import generate_interactive_plot
    material = _str(el["structure"].attrs.get("material", "DEFAULT"))
    if cubes:                     # the 3D panel embeds every cube found (several MB per page)
        for cube in list(run_dir.glob("spatial_*.cube")) + list(run_dir.glob("spinor_*.cube")):
            shutil.copy(cube, cube.name)
    made = []
    for prefix, orb, spin in (("sf", "sf/mo", 2.0), ("soc", "soc/spinor", 1.0)):
        if f"{prefix}/fuzzy" not in el or orb not in el:
            continue
        ewin = _fuzzy_file(el[f"{prefix}/fuzzy"], prefix)
        sigma = float(el[f"{orb}/pdos"].attrs["sigma_ev"])
        E = _orbital_files(el[orb], prefix, ewin, sigma, spin)
        if prefix == "sf":
            h = int(el[orb].attrs["homo_index"])
            e_homo, e_lumo = float(E[h]), float(E[h + 1])
        else:
            e_homo, e_lumo = _edges(E)
        html = out / f"fuzzy_dashboard_{prefix}.html"
        generate_interactive_plot(prefix=prefix, material=material, ef=0.0, e_homo=e_homo, e_lumo=e_lumo,
                                  normalize_coop=False, energy_label="DFT MO energy (eV)", output_html=str(html))
        made.append(html)
    return made


def exciton_pages(el, ex, out: Path):
    from qdex.exciton_analysis import plot_analysis_summary
    from qdex.spectrum import plot_spectrum
    qp = _attrs(el["qp"])
    made = []
    for g in ("sf", "soc"):
        if g not in ex:
            continue
        grp = ex[g]
        n = min(N_ANALYSIS, grp["energy_ev"].shape[0])
        E, f = grp["energy_ev"][:n], grp["f_osc"][:n]
        col = lambda k, default=0.0: grp[k][:n] if k in grp else np.full(n, default)
        rows = [{"energy": float(E[i]), "f_osc": float(f[i]), "PR": float(col("pr", 1.0)[i]),
                 "d_eh": float(col("d_eh_ang")[i]), "d_CT": float(col("d_ct_ang")[i]),
                 "sigma_h": float(col("sigma_h_ang")[i]), "sigma_e": float(col("sigma_e_ang")[i]),
                 "corr_eh": 0.0, "CT_Character": float(col("ct_character")[i])} for i in range(n)]
        is_soc = g == "soc"
        ref_gap = float(qp["soc_qp_gap_ev"]) if is_soc and "soc_qp_gap_ev" in qp else float(qp["qp_gap_ev"])
        metrics = {"dft_gap": float(qp["dft_gap_ev"]), "qp_correction": float(qp["scissor_ev"]),
                   "confinement_energy": float(qp.get("confinement_energy_ev", 0.0)),
                   "binding_energy": ref_gap - float(E[0]) if n else 0.0,
                   "first_exc_energy": float(E[0]) if n else 0.0, "is_soc": is_soc}
        if is_soc and "soc_dft_gap_ev" in qp:
            metrics["soc_gap"] = float(qp["soc_dft_gap_ev"])
        sp = grp["spectrum"]
        sigmas = sorted(k for k in sp if k.startswith("sigma_"))
        sigma = float(sigmas[0].split("_")[1]) if sigmas else 0.03
        html = out / f"exciton_analysis_{g}.html"
        plot_analysis_summary(rows, physics_metrics=metrics, filename=str(html), broadening="gaussian", sigma=sigma)
        made.append(html)
        if sigmas:
            x, y = sp["energy_ev"][:], sp[sigmas[0]][:]
            all_E, all_f = grp["energy_ev"][:], grp["f_osc"][:]
            png = out / f"spectrum_{g}.png"
            plot_spectrum(x, y, all_E, all_f, filename=str(png))
            made += [png, out / f"spectrum_{g}_nm.png"]
    return made


def rebuild(run_dir: str, out_dir: str | None = None, keep_data: bool = False, cubes: bool = True) -> list:
    import h5py
    run = Path(run_dir).resolve()
    out = Path(out_dir or run_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="qdex_dash_", dir=out))
    here = os.getcwd()
    made = []
    try:
        os.chdir(work)
        with h5py.File(run / "qdex_electronic.h5", "r") as el:
            made += fuzzy_dashboards(el, run, out, cubes)
            ex_path = run / "qdex_excitations.h5"
            if ex_path.is_file():
                with h5py.File(ex_path, "r") as ex:
                    made += exciton_pages(el, ex, out)
    finally:
        os.chdir(here)
        if keep_data:
            for p in work.iterdir():
                if p.suffix in (".csv", ".npz"):
                    shutil.move(str(p), out / p.name)
        shutil.rmtree(work, ignore_errors=True)
    return [str(p) for p in made if Path(p).exists()]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="qdex dashboards", description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir", help="QDEX run directory with qdex_electronic.h5 (and qdex_excitations.h5)")
    ap.add_argument("-o", "--out", help="where the pages go (default: the run directory)")
    ap.add_argument("--keep-data", action="store_true", help="also keep the CSV / NPZ files the pages are drawn from")
    ap.add_argument("--no-cubes", action="store_true",
                    help="leave the 3D MO panel out of the fuzzy pages (much smaller; show the .cube files elsewhere)")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    for p in rebuild(a.run_dir, a.out, a.keep_data, cubes=not a.no_cubes):
        print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
