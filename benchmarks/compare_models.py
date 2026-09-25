#!/usr/bin/env python
"""Run and cross-check QDEX quasiparticle and BSE model combinations.

The runs are organised around five questions, so the many possible
combinations are compared one factor at a time:

  A  QP model       Which quasiparticle gap does each ``qp_gap`` option give
                    (same Resta BSE kernel)?
  B  Environment    Does the optical gap stay (nearly) constant when the
                    solvent changes, while the QP gap moves?
  C  BSE kernel     How much exciton binding does each screened W give
                    (same QP gap)?
  D  Framework      Do the excitation modes, spin and charge partition behave
                    as they must (exact inequalities)?
  E  Convergence    Is the lowest exciton converged in the active space?

Every case runs the normal ``qdex`` CLI in its own directory, so the results
are exactly what a user run would give.  Automatic checks are reported as
PASS / WARN / FAIL with the reason.  A WARN is not necessarily a bug: some mark
known limitations of a model (for example the environment check for the
scissor-type QP models).

Usage::

    python benchmarks/compare_models.py --system tests/CdSe --profile standard
    python benchmarks/compare_models.py --system tests/CdSe --profile quick --exp-gap 2.70 2.95
    python benchmarks/compare_models.py --system my/dir --only A,B --eps-solvent 2.02

``--system`` must contain the ``config.yaml`` and the files it references.
Finished cases are reused unless ``--force`` is given, so an interrupted run
can simply be restarted.  Results: ``<out>/summary.md``, ``summary.csv`` and
``summary.png``.
"""
import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

# ---------------------------------------------------------------------------
# Case definitions: (group, name, extra CLI args, profiles)
# ---------------------------------------------------------------------------
Q, S, F = "quick", "standard", "full"
ALL = (Q, S, F)
STD = (S, F)


def build_cases(eps_solvent, eps_inf):
    es = f"{eps_solvent:g}"
    ei = f"{eps_inf:g}"
    cases = [
        # A: QP models, vacuum, Resta BSE
        ("A", "pbe", ["--qp_gap", "pbe"], ALL),
        ("A", "brus", ["--qp_gap", "brus"], STD),
        ("A", "gw", ["--qp_gap", "gw"], ALL),
        ("A", "sgw-resta", ["--qp_gap", "sgw-resta"], STD),
        ("A", "sgw-resta-pure", ["--qp_gap", "sgw-resta-pure"], F),
        ("A", "sgw-dim", ["--qp_gap", "sgw-dim"], STD),
        ("A", "evgw-resta", ["--qp_gap", "evgw-resta"], F),
        ("A", "evgw-dim", ["--qp_gap", "evgw-dim"], F),
        ("A", "qsgw-resta", ["--qp_gap", "qsgw-resta"], F),
        ("A", "qsgw-dim", ["--qp_gap", "qsgw-dim"], STD),
        ("A", "sgw+sbse", ["--qp_gap", "sgw", "--kernel", "sbse"], F),
        ("A", "env", ["--qp_gap", "env"], ALL),
        ("A", "env+residual", ["--qp_gap", "env", "--env-anchor-residual"], STD),
        # B: environment consistency (vacuum run is the group-A case)
        ("B", f"gw@{es}", ["--qp_gap", "gw", "--eps-out", es], ALL),
        ("B", f"sgw-resta@{es}", ["--qp_gap", "sgw-resta", "--eps-out", es], STD),
        ("B", f"sgw-dim@{es}", ["--qp_gap", "sgw-dim", "--eps-out", es], F),
        ("B", f"env@{es}", ["--qp_gap", "env", "--eps-out", es], ALL),
        ("B", f"env@{ei}", ["--qp_gap", "env", "--eps-out", ei], STD),
        ("B", f"env+residual@{es}", ["--qp_gap", "env", "--env-anchor-residual", "--eps-out", es], STD),
        # C: BSE kernels at the gw QP gap (resta is the group-A 'gw' case)
        ("C", "kernel:dim", ["--qp_gap", "gw", "--kernel", "dim"], STD),
        ("C", "kernel:xs-resta", ["--qp_gap", "gw", "--kernel", "xs-resta"], F),
        ("C", "kernel:sbse", ["--qp_gap", "gw", "--kernel", "sbse"], STD),
        ("C", "kernel:mnok-bare", ["--qp_gap", "gw", "--kernel", "bse"], STD),
        ("C", "kernel:xs-resta+env", ["--qp_gap", "env", "--kernel", "xs-resta"], F),
        # D: frameworks, spin, charges (all at gw + Resta)
        ("D", "mode:independent_dft", ["--qp_gap", "gw", "--excitation-mode", "independent_dft"], ALL),
        ("D", "mode:independent_qp", ["--qp_gap", "gw", "--excitation-mode", "independent_qp"], ALL),
        ("D", "mode:diagonal_bse", ["--qp_gap", "gw", "--excitation-mode", "diagonal_bse"], ALL),
        ("D", "spin:triplet", ["--qp_gap", "gw", "--triplet"], STD),
        ("D", "charges:lowdin", ["--qp_gap", "gw", "--charge_type", "lowdin"], F),
        ("D", "no-exchange", ["--qp_gap", "gw", "--no-exchange"], F),
        # E: active-space convergence (Davidson)
        ("E", "as:50x50", ["--qp_gap", "gw", "--nhomos", "50", "--nlumos", "50"], STD),
        ("E", "as:100x100", ["--qp_gap", "gw", "--nhomos", "100", "--nlumos", "100", "--nroots", "10", "__davidson__"], F),
    ]
    return cases


# ---------------------------------------------------------------------------
# Running and parsing
# ---------------------------------------------------------------------------
def prepare_inputs(system, cfg):
    """Return the list of input files to link, decompressing a .gz MO file if needed."""
    files = []
    for key in ("mo_file", "xyz", "basis_txt", "mo_file_beta", "vxc_ao"):
        val = cfg.get("system", {}).get(key) or cfg.get("physics", {}).get(key)
        if val:
            files.append(val)
    soc = cfg.get("soc", {})
    if soc.get("gth_file"):
        files.append(soc["gth_file"])
    for f in files:
        p = system / f
        if not p.exists() and (system / (f + ".gz")).exists():
            print(f"  decompressing {f}.gz ...")
            with gzip.open(system / (f + ".gz"), "rb") as src, open(p, "wb") as dst:
                shutil.copyfileobj(src, dst)
        if not p.exists():
            raise FileNotFoundError(f"input file {p} referenced by config.yaml is missing")
    return files


def write_config(cfg, dest, nthreads, davidson, soc):
    c = json.loads(json.dumps(cfg))  # deep copy
    c.setdefault("system", {})["nthreads"] = nthreads
    c.setdefault("soc", {})["soc_flag"] = bool(soc)
    c.setdefault("output", {})["plot"] = False
    c["output"]["cube"] = False
    phys = c.setdefault("physics", {})
    if "exchange" in phys:  # deprecated alias of include_direct_eh
        phys.setdefault("include_direct_eh", phys.pop("exchange"))
    phys.setdefault("kernel", "resta")
    c.setdefault("solver", {})["full_diag"] = not davidson
    with open(dest, "w") as fh:
        yaml.safe_dump(c, fh, sort_keys=False)


def parse_run(workdir):
    log = (workdir / "run.out").read_text(errors="replace")
    res = {"ok": "All calculations finished successfully." in log}
    m = re.findall(r"\[QP\] Final QP gap: ([-\d.]+)", log)
    res["qp_gap"] = float(m[-1]) if m else None
    m = re.search(r"\[DFT\] Initial Gap\s*:\s*([-\d.]+)", log)
    res["dft_gap"] = float(m.group(1)) if m else None
    csvs = sorted(workdir.glob("exciton_results*.csv"))
    rows = []
    if csvs:
        with open(csvs[0]) as fh:
            rows = list(csv.DictReader(fh))
    if rows:
        res["s1"] = float(rows[0]["Energy_eV"])
        res["f1"] = float(rows[0]["f_osc"])
        bright = [r for r in rows if float(r["f_osc"]) > 0.05]
        res["bright"] = float(bright[0]["Energy_eV"]) if bright else None
    res["qsgw_converged"] = None
    if "qsGW" in log:
        res["qsgw_converged"] = bool(re.search(r"qsGW-\w+ converged", log))
    err = re.findall(r"^(\w*Error: .*)$", log, flags=re.M)
    res["error"] = err[-1][:200] if err else None
    return res


def run_case(system, files, cfg, out, name, extra, nthreads, force, soc):
    d = out / re.sub(r"[^A-Za-z0-9_.@+-]", "_", name)
    done = d / "result.json"
    if done.exists() and not force:
        return json.loads(done.read_text())
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    for f in files:
        os.symlink((system / f).resolve(), d / Path(f).name)
    davidson = "__davidson__" in extra
    args = [a for a in extra if a != "__davidson__"]
    write_config(cfg, d / "config.yaml", nthreads, davidson, soc)
    cmd = [sys.executable, "-c", "from qdex.cli import main; main()", "--config", "config.yaml", *args]
    t0 = time.time()
    with open(d / "run.out", "w") as fh:
        rc = subprocess.run(cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT).returncode
    res = parse_run(d)
    res.update({"name": name, "returncode": rc, "wall_s": round(time.time() - t0, 1), "args": " ".join(args)})
    done.write_text(json.dumps(res, indent=1))
    return res


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
def checks(R, groups, eps_solvent, eps_inf, bulk_shift, exp_gap):
    out = []

    def add(status, what, detail):
        out.append((status, what, detail))

    def get(name, key):
        r = R.get(name)
        return None if r is None else r.get(key)

    es, ei = f"{eps_solvent:g}", f"{eps_inf:g}"
    for name, r in R.items():
        if not r.get("ok"):
            add("FAIL", f"{name} ran", r.get("error") or f"return code {r.get('returncode')}; see {name}/run.out")
    for name, r in R.items():
        if r.get("qsgw_converged") is False:
            add("WARN", f"{name} converged", "qsGW loop reached max iterations")

    # D: exact relations
    if get("mode:independent_dft", "s1") is not None and get("gw", "dft_gap") is not None:
        d = get("mode:independent_dft", "s1") - get("gw", "dft_gap")
        add("PASS" if abs(d) < 2e-3 else "FAIL", "independent_dft lowest = DFT gap", f"difference {d*1000:+.1f} meV")
    if get("mode:independent_qp", "s1") is not None and get("gw", "qp_gap") is not None:
        d = get("mode:independent_qp", "s1") - get("gw", "qp_gap")
        add("PASS" if abs(d) < 2e-3 else "FAIL", "independent_qp lowest = QP gap", f"difference {d*1000:+.1f} meV")
    if get("mode:diagonal_bse", "s1") is not None and get("gw", "s1") is not None:
        d = get("mode:diagonal_bse", "s1") - get("gw", "s1")
        add("PASS" if d >= -1e-4 else "FAIL", "full BSE S1 <= lowest diagonal element",
            f"diagonal - full = {d*1000:+.1f} meV (must be >= 0 by the variational principle)")
    if get("spin:triplet", "s1") is not None and get("gw", "s1") is not None:
        d = get("gw", "s1") - get("spin:triplet", "s1")
        add("PASS" if d >= -1e-4 else "FAIL", "T1 <= S1 (exchange is repulsive)", f"S1 - T1 = {d*1000:+.1f} meV")
    if get("no-exchange", "s1") is not None and get("spin:triplet", "s1") is not None:
        d = get("no-exchange", "s1") - get("spin:triplet", "s1")
        add("PASS" if abs(d) < 2e-3 else "FAIL", "singlet without exchange = triplet", f"difference {d*1000:+.1f} meV")
    for name, r in R.items():
        if r.get("s1") is not None and r.get("qp_gap") is not None and "independent" not in name:
            eb = r["qp_gap"] - r["s1"]
            if eb < -1e-3:
                add("FAIL", f"{name}: binding > 0", f"QP - S1 = {eb:.3f} eV")

    # B: environment consistency
    for model in ["gw", "sgw-resta", "sgw-dim", "env", "env+residual"]:
        vac, sol = R.get(model), R.get(f"{model}@{es}")
        if not (vac and sol and vac.get("s1") and sol.get("s1")):
            continue
        dqp = sol["qp_gap"] - vac["qp_gap"]
        ds1 = sol["s1"] - vac["s1"]
        ratio = abs(ds1) / max(abs(dqp), 1e-6)
        status = "PASS" if ratio < 0.2 else "WARN"
        why = "" if status == "PASS" else "  (optical gap follows the QP gap: BSE kernel lacks the solvent term)"
        add(status, f"{model}: optical gap insensitive to solvent",
            f"eps_out 1 -> {es}: dQP = {dqp:+.3f} eV, dS1 = {ds1:+.3f} eV, ratio {ratio:.2f}{why}")
    r = R.get(f"env@{ei}")
    if r and r.get("qp_gap") is not None and r.get("dft_gap") is not None:
        d = r["qp_gap"] - (r["dft_gap"] + bulk_shift)
        add("PASS" if abs(d) < 2e-3 else "FAIL", "env, matched medium: QP gap = DFT + bulk shift",
            f"difference {d*1000:+.1f} meV (reaction field must vanish)")

    # C: kernel spread
    kern = {"resta": get("gw", "s1"), "dim": get("kernel:dim", "s1"), "xs-resta": get("kernel:xs-resta", "s1"),
            "sbse": get("kernel:sbse", "s1"), "mnok-bare": get("kernel:mnok-bare", "s1")}
    qp = get("gw", "qp_gap")
    kb = {k: qp - v for k, v in kern.items() if v is not None and qp is not None}
    if len(kb) > 1:
        spread = max(kb.values()) - min(kb.values())
        add("PASS" if spread < 0.2 else "WARN", "binding energy robust to the W model",
            "; ".join(f"{k} {v:.2f} eV" for k, v in kb.items()) + f"  (spread {spread:.2f} eV)")

    # E: convergence
    seq = [("25", get("gw", "s1")), ("50", get("as:50x50", "s1")), ("100", get("as:100x100", "s1"))]
    seq = [(n, v) for n, v in seq if v is not None]
    if len(seq) >= 2:
        d = seq[-1][1] - seq[-2][1]
        add("PASS" if abs(d) < 0.01 else "WARN", "S1 converged in the active space",
            " -> ".join(f"{n}x{n}: {v:.3f}" for n, v in seq) + f" eV (last change {d*1000:+.1f} meV)")

    # experiment
    if exp_gap:
        lo, hi = exp_gap
        for name in [f"gw@{es}", f"sgw-resta@{es}", f"env@{es}", f"env+residual@{es}", "env", "env+residual"]:
            v = get(name, "s1")
            if v is None:
                continue
            dev = 0.0 if lo <= v <= hi else (v - hi if v > hi else v - lo)
            add("PASS" if dev == 0 else "WARN", f"{name}: S1 within experiment [{lo}, {hi}] eV",
                f"S1 = {v:.3f} eV" + ("" if dev == 0 else f" ({dev:+.2f} eV outside)"))
    return out


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
TITLES = {
    "A": "A. Quasiparticle model (Resta BSE kernel, vacuum)",
    "B": "B. Environment: solvent run vs. vacuum run of the same model",
    "C": "C. BSE kernel at fixed QP gap (gw)",
    "D": "D. Excitation framework, spin and charge partition (gw + Resta)",
    "E": "E. Active-space convergence (gw + Resta)",
}


def report(out, cases, R, chk, meta):
    lines = [f"# QDEX model comparison: {meta['system']}", "",
             f"profile `{meta['profile']}`, eps_solvent = {meta['eps_solvent']}, eps_inf = {meta['eps_inf']}, "
             f"bulk GW-PBE opening = {meta['bulk_shift']:.3f} eV, spin-free", ""]
    rows_csv = []
    for g in "ABCDE":
        names = [c[1] for c in cases if c[0] == g]
        if g == "B":
            names = [n for n in ["gw", "sgw-resta", "sgw-dim", "env", "env+residual"] if n in R] + names
        if g in "CDE":
            names = ["gw"] + names
        names = [n for n in dict.fromkeys(names) if n in R]
        if not names:
            continue
        lines += [f"## {TITLES[g]}", "", "| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | wall (s) | status |",
                  "|---|---|---|---|---|---|---|"]
        for n in names:
            r = R[n]
            qp, s1 = r.get("qp_gap"), r.get("s1")
            eb = f"{qp - s1:.3f}" if (qp is not None and s1 is not None) else "—"
            st = "ok" if r.get("ok") else "FAILED"
            qp_s = f"{qp:.3f}" if qp is not None else "—"
            s1_s = f"{s1:.3f}" if s1 is not None else "—"
            f1_s = f"{r['f1']:.3f}" if r.get("f1") is not None else "—"
            lines.append(f"| `{n}` | {qp_s} | {s1_s} | {f1_s} | {eb} | {r.get('wall_s')} | {st} |")
            rows_csv.append({"group": g, "case": n, **{k: r.get(k) for k in ["qp_gap", "s1", "f1", "bright", "dft_gap", "wall_s", "ok", "args"]}})
        lines.append("")
    lines += ["## Checks", "", "| status | check | detail |", "|---|---|---|"]
    order = {"FAIL": 0, "WARN": 1, "PASS": 2}
    for st, what, det in sorted(chk, key=lambda x: order[x[0]]):
        lines.append(f"| **{st}** | {what} | {det} |")
    lines += ["", "How to read this report: `docs/validation/model_comparison.rst`.", ""]
    (out / "summary.md").write_text("\n".join(lines))
    with open(out / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_csv[0].keys()) if rows_csv else ["case"])
        w.writeheader()
        w.writerows(rows_csv)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        names = [n for n in R if R[n].get("s1") is not None and R[n].get("qp_gap") is not None]
        fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(names)), 4))
        x = range(len(names))
        ax.bar(x, [R[n]["qp_gap"] for n in names], color="#c4b5fd", label="QP gap")
        ax.bar(x, [R[n]["s1"] for n in names], width=0.45, color="#059669", label="S1")
        if meta.get("exp_gap"):
            ax.axhspan(*meta["exp_gap"], color="#fde68a", alpha=0.6, lw=0, label="experiment")
        ax.set_xticks(list(x))
        ax.set_xticklabels(names, rotation=70, ha="right", fontsize=7.5)
        ax.set_ylabel("eV")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(out / "summary.png", dpi=130)
    except Exception as exc:  # plotting is optional
        print("  (plot skipped:", exc, ")")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", required=True, help="directory with config.yaml and its input files")
    ap.add_argument("--profile", choices=ALL, default=S)
    ap.add_argument("--out", default=None, help="output directory (default: <system>/compare_<profile>)")
    ap.add_argument("--only", default=None, help="comma-separated groups, e.g. A,B")
    ap.add_argument("--eps-solvent", type=float, default=2.24, help="optical eps_out of the solvent (n^2); toluene 2.24, hexane 1.89")
    ap.add_argument("--exp-gap", type=float, nargs=2, default=None, metavar=("LO", "HI"),
                    help="experimental first-exciton window in eV for the plausibility check")
    ap.add_argument("--nthreads", type=int, default=4)
    ap.add_argument("--soc", action="store_true", help="keep SOC on (slower); default spin-free")
    ap.add_argument("--force", action="store_true", help="rerun cases that already have results")
    ap.add_argument("--list", action="store_true", help="only list the cases of the profile")
    a = ap.parse_args()

    system = Path(a.system).resolve()
    cfg = yaml.safe_load((system / "config.yaml").read_text())
    material = str(cfg.get("system", {}).get("material", "DEFAULT")).upper()
    from qdex.hardness import MATERIAL_DB
    entry = MATERIAL_DB.get(material, MATERIAL_DB["DEFAULT"])
    eps_inf = float(entry[0])
    bulk_shift = float(entry[8] - entry[7]) if len(entry) >= 9 else 0.0

    groups = set(a.only.split(",")) if a.only else set("ABCDE")
    cases = [c for c in build_cases(a.eps_solvent, eps_inf) if a.profile in c[3] and c[0] in groups]
    # cases other groups need as references
    needed = {"gw"} | ({"gw", "sgw-resta", "sgw-dim", "env", "env+residual"} if "B" in groups else set())
    ref = [c for c in build_cases(a.eps_solvent, eps_inf) if c[1] in needed and c not in cases and a.profile in c[3]]
    cases = ref + cases
    if a.list:
        for g, n, args, _ in cases:
            print(f"{g}  {n:28s} {' '.join(x for x in args if x != '__davidson__')}")
        return

    out = Path(a.out) if a.out else system / f"compare_{a.profile}"
    out.mkdir(parents=True, exist_ok=True)
    files = prepare_inputs(system, cfg)
    R = {}
    for k, (g, name, extra, _) in enumerate(cases, 1):
        print(f"[{k:2d}/{len(cases)}] {g} {name:28s}", end=" ", flush=True)
        r = run_case(system, files, cfg, out, name, extra, a.nthreads, a.force, a.soc)
        R[name] = r
        s1 = r.get("s1")
        print(f"QP={r.get('qp_gap')}  S1={s1 if s1 is None else round(s1, 3)}  "
              f"{'ok' if r.get('ok') else 'FAILED'}  ({r.get('wall_s')} s)")
    meta = {"system": str(system), "profile": a.profile, "eps_solvent": a.eps_solvent,
            "eps_inf": eps_inf, "bulk_shift": bulk_shift, "exp_gap": a.exp_gap}
    chk = checks(R, groups, a.eps_solvent, eps_inf, bulk_shift, a.exp_gap)
    report(out, cases, R, chk, meta)
    n_fail = sum(1 for c in chk if c[0] == "FAIL")
    n_warn = sum(1 for c in chk if c[0] == "WARN")
    print(f"\n{len(chk)} checks: {n_fail} FAIL, {n_warn} WARN.  Report: {out / 'summary.md'}")


if __name__ == "__main__":
    main()
