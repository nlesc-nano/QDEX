"""Numerical regression of the CdSe validation numbers (docs/validation/cdse_experiment.rst).

Guards the default settings (MNOK with gamma_AA = IP - EA and exponent 2, Mulliken charges,
SAXS radius, bulk QSGW correction) against silent drifts. Cd16Se13Cl6 (1.2 nm), spin-free,
25 x 25, dense diagonalization; about 7 s per case::

    QDEX_RUN_REGRESSION=1 python -m pytest tests/test_cdse_regression.py

The 2.0 nm cases take about a minute each and also need QDEX_RUN_CDSE=1.
References: benchmarks/results/cdse_validation_2026-09-27.csv. If a change is meant to move
these numbers, rerun the validation, update the docs and the CSV, and then the values here.
"""
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent / "CdSe"
REPO = Path(__file__).resolve().parent.parent
RUN = os.environ.get("QDEX_RUN_REGRESSION") == "1"
RUN_LARGE = RUN and os.environ.get("QDEX_RUN_CDSE") == "1"
TOL = 0.005  # eV

# (name, size, flags, eps_out, QP gap, S1, first bright state)
CASES = [
    ("sgw_resta_vacuum", "1.2nm", ["--qp_gap", "sgw-resta"], 1.0, 6.767, 3.454, 3.819),
    ("sgw_resta_toluene", "1.2nm", ["--qp_gap", "sgw-resta"], 2.24, 5.407, 3.406, 3.746),
    ("sgw_dim_vacuum", "1.2nm", ["--qp_gap", "sgw-dim"], 1.0, 6.635, 3.476, 3.845),
    ("qsgw_dim_vacuum", "1.2nm", ["--qp_gap", "qsgw-dim"], 1.0, 6.626, 3.437, 3.809),
    ("sbse_resta", "1.2nm", ["--qp_gap", "bulk", "--excitation-mode", "sbse", "--kernel", "resta"],
     1.0, 4.209, 3.323, 3.617),
    ("sbse_dim", "1.2nm", ["--qp_gap", "bulk", "--excitation-mode", "sbse", "--kernel", "dim"],
     1.0, 4.209, 3.137, 3.442),
    ("stda_dielectric", "1.2nm", ["--qp_gap", "bulk", "--excitation-mode", "stda", "--stda-ax", "dielectric"],
     1.0, 4.209, 3.808, None),
    ("stda_pbe", "1.2nm", ["--qp_gap", "none", "--excitation-mode", "stda", "--stda-functional", "pbe"],
     1.0, 2.639, 2.748, None),
    ("sgw_resta_toluene_2nm", "2.0nm", ["--qp_gap", "sgw-resta"], 2.24, 3.707, 2.964, 2.964),
    ("sbse_resta_2nm", "2.0nm", ["--qp_gap", "bulk", "--excitation-mode", "sbse", "--kernel", "resta"],
     1.0, 3.027, 2.779, 2.779),
]
MO_FILES = {"1.2nm": "MOs_cleaned_12ang.txt", "2.0nm": "MOs_cleaned_20ang.txt"}


def _states(workdir):
    with open(Path(workdir) / "exciton_results.csv") as fh:
        rows = list(csv.DictReader(fh))
    energies = [float(r["Energy_eV"]) for r in rows]
    f_osc = [float(r["f_osc"]) for r in rows]
    fmax = max(f_osc)
    bright = next(e for e, f in zip(energies, f_osc) if f >= 0.1 * fmax)
    return energies[0], bright


def _qp_gap(log):
    lines = [l for l in log.splitlines() if l.startswith("[QP] Final QP gap")]
    if lines:
        return float(lines[-1].split()[-2])
    lines = [l for l in log.splitlines() if "[QP]  Target Gap" in l]
    return float(lines[-1].split()[-2])


@unittest.skipUnless(RUN, "set QDEX_RUN_REGRESSION=1 to run the CdSe regression test")
class CdSeRegressionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="qdex_regression_")
        cls.dirs = {}
        for size, mo_name in MO_FILES.items():
            if size == "2.0nm" and not RUN_LARGE:
                continue
            src, dst = ROOT / size, Path(cls.tmp) / size
            dst.mkdir()
            mo = src / mo_name
            if mo.exists():
                os.symlink(mo, dst / mo_name)
            else:
                with gzip.open(src / (mo_name + ".gz"), "rb") as fin, open(dst / mo_name, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
            for name in ("geom.xyz", "BASIS_MOLOPT_UZH", "GTH_SOC_POTENTIALS.txt"):
                os.symlink(src / name, dst / name)
            text = (src / "config.yaml").read_text()
            text = (text.replace("enabled: true", "enabled: false").replace("nthreads: 12", "nthreads: 2")
                    .replace("plot: true", "plot: false"))
            (dst / "config.yaml").write_text(text)
            cls.dirs[size] = dst

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _check(self, name, size, flags, eps, qp_ref, s1_ref, bright_ref):
        if size not in self.dirs:
            self.skipTest("2.0 nm cases need QDEX_RUN_CDSE=1")
        base = self.dirs[size]
        work = Path(tempfile.mkdtemp(dir=base))
        for f in base.iterdir():
            if f.is_file() or f.is_symlink():
                os.symlink(f.resolve(), work / f.name)
        cmd = [sys.executable, "-c", "from qdex.cli import main; main()", "--config", "config.yaml",
               "--eps-out", str(eps), *flags]
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(REPO), env.get("PYTHONPATH")]))  # qdex, libint_cpp
        res = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=1800, env=env)
        self.assertEqual(res.returncode, 0, res.stdout[-3000:] + res.stderr[-3000:])
        s1, bright = _states(work)
        self.assertAlmostEqual(_qp_gap(res.stdout), qp_ref, delta=TOL, msg=f"{name}: QP gap")
        self.assertAlmostEqual(s1, s1_ref, delta=TOL, msg=f"{name}: S1")
        if bright_ref is not None:
            self.assertAlmostEqual(bright, bright_ref, delta=TOL, msg=f"{name}: first bright state")


def _make_test(case):
    def test(self):
        self._check(*case)
    return test


for _case in CASES:
    setattr(CdSeRegressionTest, f"test_{_case[0]}", _make_test(_case))
