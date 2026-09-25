#!/usr/bin/env python
"""Full sweep of consistent QP x BSE combinations for one quantum dot.

Copy this file into a system folder that contains ``config.yaml`` and the
files it references (MO file or MO .gz, xyz, basis, GTH SOC file), then run::

    python qp_bse_sweep.py --list              # show the cases
    python qp_bse_sweep.py --dry-run           # only write sweep/<case>/config.yaml
    python qp_bse_sweep.py --jobs 2 --nthreads 4
    python qp_bse_sweep.py --groups A,E        # a subset
    python qp_bse_sweep.py --collect           # rebuild sweep/summary.* from finished runs

Every case gets its own folder ``sweep/<case>/`` with a complete
``config.yaml`` (so each combination can be inspected and rerun by hand with
``qdex --config config.yaml``), its log ``run.out`` and the usual QDEX outputs.
Finished cases are skipped unless ``--force`` is given.

Consistency rule: one W for QP and BSE
--------------------------------------
* ``two_electron_integrals`` (mnok | xs) is the representation of W. The QP
  correction and the BSE kernel use the same representation.
* Delta-W QP models (sgw-*, evgw-*, qsgw-*; Resta or DIM) define W; the BSE uses
  that W (``kernel: qp``). Resta and DIM are never mixed.
* ``gw`` (two-anchor) with the sphere polarization uses ``resta-sphere``:
  bulk Resta W plus the reaction field of the same dielectric sphere.
* ``gw`` with the legacy curve and ``brus`` define no W. They use the bulk
  Resta kernel; this is the legacy (not consistent) reference.

Groups
------
A  core matrix: integrals {mnok, xs} x QP model x eps_out {1, solvent}
   QP models: gw (sphere), gw (legacy), brus, sgw/evgw/qsgw x {resta, dim}
B  quasiparticle weight: Z = 1 and Z = 0.8 against the derived default
C  QP levels: rigid scissor against orbital-resolved levels
D  two-anchor sensitivity: residual power p = 1.5, 3; gw-sphere energies
   with the bulk-only kernel (shows what the reaction field does)
E  spin-orbit coupling for the main models, in the solvent
F  active-space convergence (50 x 50, 100 x 100, Davidson)
G  transition charges: Loewdin instead of Mulliken
H  singlet-triplet (exchange) splitting, spin-free
"""
import argparse
import copy
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml

HERE = Path.cwd()

QP_MODELS = {
    # name: (physics overrides)
    "gw-sphere": {"qp_gap": "gw", "qp_polarization": "sphere", "kernel": "resta-sphere"},
    "gw-legacy": {"qp_gap": "gw", "qp_polarization": "legacy", "kernel": "resta"},
    "brus": {"qp_gap": "brus", "kernel": "resta"},
    "sgw-resta": {"qp_gap": "sgw-resta", "kernel": "qp"},
    "evgw-resta": {"qp_gap": "evgw-resta", "kernel": "qp"},
    "qsgw-resta": {"qp_gap": "qsgw-resta", "kernel": "qp"},
    "sgw-dim": {"qp_gap": "sgw-dim", "kernel": "qp"},
    "evgw-dim": {"qp_gap": "evgw-dim", "kernel": "qp"},
    "qsgw-dim": {"qp_gap": "qsgw-dim", "kernel": "qp"},
}


def build_cases(eps_solvent):
    es = f"{eps_solvent:g}"
    cases = []

    def add(group, name, phys, soc=False, solver=None):
        cases.append({"group": group, "name": name, "physics": phys, "soc": soc, "solver": solver or {}})

    # A: core matrix
    for ints in ("mnok", "xs"):
        for m, ph in QP_MODELS.items():
            for eo in (1.0, eps_solvent):
                tag = "vac" if eo == 1.0 else f"eps{es}"
                add("A", f"{ints}_{m}_{tag}", {**ph, "two_electron_integrals": ints, "eps_out": eo})
    # B: Z
    for m in ("sgw-resta", "sgw-dim", "qsgw-dim"):
        for z in ("1.0", "0.8"):
            add("B", f"mnok_{m}_Z{z}_vac", {**QP_MODELS[m], "two_electron_integrals": "mnok", "eps_out": 1.0,
                                            "qp_z": z})
    add("B", f"mnok_sgw-resta_Z1.0_eps{es}", {**QP_MODELS["sgw-resta"], "two_electron_integrals": "mnok",
                                              "eps_out": eps_solvent, "qp_z": "1.0"})
    # C: rigid scissor vs orbital-resolved levels
    for m in ("gw-sphere", "sgw-resta", "sgw-dim"):
        for ints in ("mnok", "xs"):
            add("C", f"{ints}_{m}_rigid_vac", {**QP_MODELS[m], "two_electron_integrals": ints, "eps_out": 1.0,
                                               "qp_levels": "rigid"})
    # D: two-anchor sensitivity
    for p in (1.5, 3.0):
        add("D", f"mnok_gw-sphere_p{p:g}_vac", {**QP_MODELS["gw-sphere"], "two_electron_integrals": "mnok",
                                               "eps_out": 1.0, "qp_residual_power": p})
    for eo in (1.0, eps_solvent):
        tag = "vac" if eo == 1.0 else f"eps{es}"
        add("D", f"mnok_gw-sphere+bulkkernel_{tag}", {**QP_MODELS["gw-sphere"], "kernel": "resta",
                                                      "two_electron_integrals": "mnok", "eps_out": eo})
    # E: SOC in the solvent
    for ints in ("mnok", "xs"):
        for m in ("gw-sphere", "sgw-resta", "sgw-dim", "qsgw-dim"):
            add("E", f"{ints}_{m}_soc_eps{es}", {**QP_MODELS[m], "two_electron_integrals": ints,
                                                 "eps_out": eps_solvent}, soc=True)
    # F: active space
    for n in (50, 100):
        for m in ("gw-sphere", "sgw-resta"):
            add("F", f"mnok_{m}_as{n}_vac", {**QP_MODELS[m], "two_electron_integrals": "mnok", "eps_out": 1.0,
                                             "nhomos": n, "nlumos": n, "nroots": 10},
                solver={"full_diag": False})
    # G: charges
    for m in ("gw-sphere", "sgw-resta"):
        add("G", f"mnok_{m}_lowdin_vac", {**QP_MODELS[m], "two_electron_integrals": "mnok", "eps_out": 1.0,
                                          "charge_type": "lowdin"})
    # H: triplets
    for m in ("gw-sphere", "sgw-resta"):
        add("H", f"mnok_{m}_triplet_vac", {**QP_MODELS[m], "two_electron_integrals": "mnok", "eps_out": 1.0,
                                           "triplet": True})
    return cases


# ---------------------------------------------------------------------------
def input_files(cfg):
    files = []
    for key in ("mo_file", "xyz", "basis_txt", "mo_file_beta", "vxc_ao"):
        val = cfg.get("system", {}).get(key) or cfg.get("physics", {}).get(key)
        if val:
            files.append(val)
    if cfg.get("soc", {}).get("gth_file"):
        files.append(cfg["soc"]["gth_file"])
    for f in files:
        p = HERE / f
        if not p.exists() and (HERE / (f + ".gz")).exists():
            print(f"decompressing {f}.gz ...")
            with gzip.open(HERE / (f + ".gz"), "rb") as src, open(p, "wb") as dst:
                shutil.copyfileobj(src, dst)
        if not p.exists():
            raise FileNotFoundError(f"{p} (referenced by config.yaml) is missing")
    return files


def case_config(base, case, nthreads):
    c = copy.deepcopy(base)
    c.setdefault("system", {})["nthreads"] = nthreads
    phys = c.setdefault("physics", {})
    if "exchange" in phys:  # deprecated alias of include_direct_eh
        phys.setdefault("include_direct_eh", phys.pop("exchange"))
    for k in ("qp_z", "qp_levels", "qp_residual_power", "charge_type", "triplet", "nroots", "qp_polarization",
              "dynamic_z"):
        phys.pop(k, None)
    phys.update(case["physics"])
    c.setdefault("soc", {})["soc_flag"] = bool(case["soc"])
    c.setdefault("solver", {}).update({"full_diag": True, **case["solver"]})
    out = c.setdefault("output", {})
    out["plot"] = False
    out["cube"] = False
    c.setdefault("fuzzy", {})["run_fuzzy"] = False
    return c


def is_done(d):
    log = d / "run.out"
    return log.exists() and "All calculations finished successfully." in log.read_text(errors="replace")


def run_case(case, base, files, nthreads, force, dry):
    d = HERE / "sweep" / case["name"]
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "config.yaml", "w") as fh:
        yaml.safe_dump(case_config(base, case, nthreads), fh, sort_keys=False)
    for f in files:
        link = d / Path(f).name
        if not link.exists():
            os.symlink(HERE / f, link)
    if dry or (is_done(d) and not force):
        return case["name"], "skipped" if not dry else "written", 0.0
    t0 = time.time()
    with open(d / "run.out", "w") as log:
        rc = subprocess.run(["qdex", "--config", "config.yaml"], cwd=d, stdout=log, stderr=subprocess.STDOUT).returncode
    wall = time.time() - t0
    (d / "wall_s.txt").write_text(f"{wall:.1f}\n")
    return case["name"], "ok" if rc == 0 else f"FAILED (rc={rc})", wall


# ---------------------------------------------------------------------------
def _csv(path):
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh))


def main_peak(rows, frac=0.3, sigma=0.05):
    """First absorption peak: first local maximum of the Gaussian-broadened
    stick spectrum (sigma in eV) that reaches frac x its global maximum.

    With SOC the band-edge oscillator strength spreads over many weak lines, so
    a threshold on single lines (f > 0.05) can pick a weak side line; the
    broadened spectrum picks the first band that actually dominates the edge.
    """
    import numpy as np
    E = np.array([float(x["Energy_eV"]) for x in rows])
    f = np.array([float(x["f_osc"]) for x in rows])
    if len(E) == 0 or f.max() <= 0.0:
        return None
    grid = np.arange(E.min() - 0.2, E.max() + 0.2, 0.002)
    spec = (f[None, :] * np.exp(-0.5 * ((grid[:, None] - E[None, :]) / sigma) ** 2)).sum(axis=1)
    top = spec.max()
    for i in range(1, len(grid) - 1):
        if spec[i] >= frac * top and spec[i] >= spec[i - 1] and spec[i] >= spec[i + 1]:
            return float(grid[i])
    return float(grid[int(spec.argmax())])


def parse_case(case):
    d = HERE / "sweep" / case["name"]
    r = {"group": case["group"], "case": case["name"], "ok": is_done(d)}
    log = (d / "run.out").read_text(errors="replace") if (d / "run.out").exists() else ""
    m = re.findall(r"\[QP\] Final QP gap: ([-\d.]+)", log)
    r["qp_gap"] = float(m[-1]) if m else None
    m = re.search(r"\[DFT\] Initial Gap\s*:\s*([-\d.]+)", log)
    r["dft_gap"] = float(m.group(1)) if m else None
    prov = {}
    if (d / "qp_provenance.json").exists():
        prov = json.loads((d / "qp_provenance.json").read_text())
        prov = prov.get("details", prov)
    for key, pk in [("qp_homo", "modeled_qp_homo_ev"), ("qp_lumo", "modeled_qp_lumo_ev"),
                    ("z_homo", "z_homo_orbital"), ("z_lumo", "z_lumo_orbital"),
                    ("z_min", "z_min_window"), ("z_max", "z_max_window"),
                    ("eps_eff", "eps_eff_qd"), ("qp_levels", "qp_levels"),
                    ("spread_occ", "qp_shift_spread_occ_ev"), ("spread_virt", "qp_shift_spread_virt_ev"),
                    ("f_homo", "f_homo"), ("radius", "cluster_radius_ang")]:
        r[key] = prov.get(pk)
    if r["z_homo"] is None:
        r["z_homo"], r["z_lumo"] = prov.get("z_homo"), prov.get("z_lumo")
    sf = _csv(d / "exciton_results.csv") or _csv(d / "exciton_results_sf.csv")
    soc = _csv(d / "exciton_results_soc.csv")
    if sf:
        r["s1"] = float(sf[0]["Energy_eV"])
        r["f1"] = float(sf[0]["f_osc"])
        b = [x for x in sf if float(x["f_osc"]) > 0.05]
        r["bright"] = float(b[0]["Energy_eV"]) if b else None
        r["peak"] = main_peak(sf)
    if soc:
        r["s1_soc"] = float(soc[0]["Energy_eV"])
        b = [x for x in soc if float(x["f_osc"]) > 0.05]
        r["bright_soc"] = float(b[0]["Energy_eV"]) if b else None
        r["peak_soc"] = main_peak(soc)
    w = d / "wall_s.txt"
    r["wall_s"] = float(w.read_text()) if w.exists() else None
    err = re.findall(r"^(\w*Error: .*)$", log, flags=re.M)
    r["error"] = err[-1][:200] if err else None
    return r


def collect(cases):
    rows = [parse_case(c) for c in cases]
    cols = ["group", "case", "ok", "dft_gap", "qp_gap", "qp_homo", "qp_lumo", "f_homo", "z_homo", "z_lumo",
            "z_min", "z_max", "eps_eff", "qp_levels", "spread_occ", "spread_virt", "s1", "f1", "bright", "peak",
            "s1_soc", "bright_soc", "peak_soc", "radius", "wall_s", "error"]
    out = HERE / "sweep"
    with open(out / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    def f(v, n=3):
        return "—" if v is None else (f"{v:.{n}f}" if isinstance(v, float) else str(v))

    lines = [f"# QP x BSE sweep: {HERE}", "",
             "| group | case | QP gap | HOMO | LUMO | Z_H / Z_L | S1 | bright | peak | S1 SOC | peak SOC | status |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        z = "—" if r.get("z_homo") is None else f"{r['z_homo']:.3f} / {r['z_lumo']:.3f}"
        st = "ok" if r["ok"] else (r.get("error") or "missing")
        lines.append(f"| {r['group']} | `{r['case']}` | {f(r.get('qp_gap'))} | {f(r.get('qp_homo'))} | "
                     f"{f(r.get('qp_lumo'))} | {z} | {f(r.get('s1'))} | {f(r.get('bright'))} | {f(r.get('peak'))} | "
                     f"{f(r.get('s1_soc'))} | {f(r.get('peak_soc'))} | {st} |")
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    n_ok = sum(r["ok"] for r in rows)
    print(f"{n_ok}/{len(rows)} cases finished. Summary: {out / 'summary.md'} and summary.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--groups", default="ABCDEFGH", help="groups to run, e.g. A or A,E (default: all)")
    ap.add_argument("--only", default=None, help="regular expression on case names")
    ap.add_argument("--eps-solvent", type=float, default=2.24, help="solvent eps_out (optical, n^2); toluene 2.24")
    ap.add_argument("--nthreads", type=int, default=4, help="threads per qdex run")
    ap.add_argument("--jobs", type=int, default=1, help="qdex runs in parallel")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="write the case folders and YAMLs only")
    ap.add_argument("--collect", action="store_true", help="only rebuild the summary")
    ap.add_argument("--force", action="store_true", help="rerun finished cases")
    a = ap.parse_args()

    base = yaml.safe_load((HERE / "config.yaml").read_text())
    groups = set(a.groups.replace(",", ""))
    cases = [c for c in build_cases(a.eps_solvent) if c["group"] in groups]
    if a.only:
        cases = [c for c in cases if re.search(a.only, c["name"])]
    if a.list:
        for c in cases:
            extra = " soc" if c["soc"] else ""
            print(f"{c['group']}  {c['name']:40s} {json.dumps(c['physics'])}{extra}")
        print(f"{len(cases)} cases")
        return
    if a.collect:
        collect(cases)
        return
    files = input_files(base)
    (HERE / "sweep").mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as pool:
        futs = [pool.submit(run_case, c, base, files, a.nthreads, a.force, a.dry_run) for c in cases]
        for k, fu in enumerate(futs, 1):
            name, status, wall = fu.result()
            print(f"[{k:3d}/{len(cases)}] {name:40s} {status}  ({wall:.0f} s)", flush=True)
    if not a.dry_run:
        collect(cases)


if __name__ == "__main__":
    main()
