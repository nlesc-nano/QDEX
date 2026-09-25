"""End-to-end regression on the CdSe 2.0 nm example in ``tests/CdSe/2.0nm``.

The MO file is large, so this test is opt-in::

    QDEX_RUN_CDSE=1 python -m pytest tests/test_cdse_integration.py

It decompresses ``MOs_cleaned_20ang.txt.gz`` into a temporary directory when the
plain-text file is absent.  Reference values were produced on the spin-free
path with the configuration in ``tests/CdSe/2.0nm/config.yaml`` (25 x 25 active
space, Resta-MNOK direct kernel, dense diagonalization).
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

HERE = Path(__file__).resolve().parent / "CdSe" / "2.0nm"
RUN = os.environ.get("QDEX_RUN_CDSE") == "1"


def _first_state(workdir):
    with open(Path(workdir) / "exciton_results.csv") as fh:
        row = next(csv.DictReader(fh))
    return float(row["Energy_eV"]), float(row["f_osc"])


def _qp_gap(log):
    lines = [l for l in log.splitlines() if l.startswith("[QP] Final QP gap")]
    return float(lines[-1].split()[-2])


@unittest.skipUnless(RUN, "set QDEX_RUN_CDSE=1 to run the CdSe integration test")
class CdSeIntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="qdex_cdse_")
        mo = HERE / "MOs_cleaned_20ang.txt"
        if mo.exists():
            os.symlink(mo, Path(cls.tmp) / mo.name)
        else:
            with gzip.open(HERE / "MOs_cleaned_20ang.txt.gz", "rb") as src, \
                    open(Path(cls.tmp) / mo.name, "wb") as dst:
                shutil.copyfileobj(src, dst)
        for name in ("geom.xyz", "BASIS_MOLOPT_UZH", "GTH_SOC_POTENTIALS.txt"):
            os.symlink(HERE / name, Path(cls.tmp) / name)
        text = (HERE / "config.yaml").read_text()
        text = text.replace("soc_flag: true", "soc_flag: false").replace("nthreads: 12", "nthreads: 4")
        (Path(cls.tmp) / "config.yaml").write_text(text)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _run(self, *extra):
        work = tempfile.mkdtemp(dir=self.tmp)
        for f in Path(self.tmp).iterdir():
            if f.is_file() or f.is_symlink():
                os.symlink(f, Path(work) / f.name)
        cmd = [sys.executable, "-c", "from qdex.cli import main; main()", "--config", "config.yaml", *extra]
        res = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=900)
        self.assertEqual(res.returncode, 0, res.stdout[-3000:] + res.stderr[-3000:])
        return work, res.stdout

    def test_gw_vacuum_reference(self):
        # Sphere polarization + resta-sphere kernel (config default), residual scaled by E_conf.
        work, log = self._run("--eps-out", "1.0")
        self.assertAlmostEqual(_qp_gap(log), 4.020, places=2)
        e1, f1 = _first_state(work)
        self.assertAlmostEqual(e1, 2.443, places=2)
        self.assertGreater(f1, 0.1)

    def test_gw_legacy_reference(self):
        # kappa/(R + ell) curve with the bulk Resta kernel.  3.869 eV with the
        # corrected anchor radius R0 = 5.3133 A (3.861 eV with the old 5.258 A).
        work, log = self._run("--eps-out", "1.0", "--qp-polarization", "legacy", "--kernel", "resta",
                              "--qp-residual-scaling", "power")
        self.assertAlmostEqual(_qp_gap(log), 3.869, places=2)
        e1, _ = _first_state(work)
        self.assertAlmostEqual(e1, 3.642, places=2)

    def test_brus_model_runs(self):
        # Regression: the CLI passed (coords, symbols, material) to a
        # function whose signature is (material, coords, symbols).
        work, log = self._run("--qp_gap", "brus", "--kernel", "resta")
        self.assertIn("[Brus Model]", log)
        e1, _ = _first_state(work)
        self.assertTrue(1.5 < e1 < 5.0)


if __name__ == "__main__":
    unittest.main()
