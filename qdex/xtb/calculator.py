"""ASE calculator for g-xTB through the xtb command line (the g-xTB "bleed" binary)."""
import os
import shutil
import subprocess

import numpy as np
from ase.calculators.calculator import Calculator, CalculationFailed, all_changes
from ase.units import Bohr, Hartree

METHOD_ARGS = {
    "gxtb": ["--gxtb"],
    "gfn2": ["--gfn", "2"],
    "gfn1": ["--gfn", "1"],
}


def read_engrad(path):
    """Energy (Eh) and gradient (natoms x 3, Eh/bohr) from an xtb ``.engrad`` file."""
    values = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
    n = int(values[0])
    energy = float(values[1])
    grad = np.array([float(v) for v in values[2:2 + 3 * n]]).reshape(n, 3)
    return energy, grad


class GXTB(Calculator):
    """
    Energy and forces from ``xtb <geom> --gxtb --grad``.

    All calculations run in ``directory``, so xtb reuses its ``xtbrestart`` file as the
    SCF guess of the next step. Set ``molden_dest`` to a path before a calculation to keep
    that calculation's ``molden.input`` there.

    The binary is taken from ``binary`` or ``$GXTB``. ``libdir`` (or ``$XTB_LIBDIR``) is put on
    the dynamic library path: the macOS bleed binary needs the conda xtb runtime libraries.
    """

    implemented_properties = ["energy", "forces"]
    default_parameters = dict(method="gxtb", charge=0, uhf=0, acc=1.0, etemp=300.0, nthreads=4)

    def __init__(self, binary=None, libdir=None, xtbpath=None, directory="gxtb_calc", **kwargs):
        super().__init__(directory=directory, **kwargs)
        self.binary = binary or os.environ.get("GXTB") or shutil.which("xtb")
        if not self.binary:
            raise FileNotFoundError("No xtb binary: pass binary= or set $GXTB.")
        self.libdir = libdir or os.environ.get("XTB_LIBDIR")
        self.xtbpath = xtbpath or os.environ.get("XTBPATH")
        self.molden_dest = None
        self.last_output = None

    def _env(self):
        env = dict(os.environ)
        env["OMP_NUM_THREADS"] = str(self.parameters.nthreads)
        env.setdefault("OMP_STACKSIZE", "1G")  # 2G x 8 threads crashed xtb at startup
        if self.xtbpath:
            env["XTBPATH"] = self.xtbpath
        if self.libdir:
            for key in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
                env[key] = self.libdir + (os.pathsep + env[key] if env.get(key) else "")
        return env

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        p = self.parameters
        os.makedirs(self.directory, exist_ok=True)
        with open(os.path.join(self.directory, "geom.xyz"), "w") as f:
            f.write(f"{len(self.atoms)}\n\n")
            for s, (x, y, z) in zip(self.atoms.get_chemical_symbols(), self.atoms.positions):
                f.write(f"{s:<3s} {x:20.12f} {y:20.12f} {z:20.12f}\n")
        for stale in ("geom.engrad", "molden.input"):
            path = os.path.join(self.directory, stale)
            if os.path.exists(path):
                os.remove(path)

        args = [self.binary, "geom.xyz", *METHOD_ARGS[p.method], "--grad",
                "--chrg", str(p.charge), "--acc", str(p.acc), "--etemp", str(p.etemp)]
        if p.uhf:
            # an explicit --uhf (even 0) makes the tblite path spin-unrestricted
            args += ["--uhf", str(p.uhf)]
        if self.molden_dest:
            args.append("--molden")
        res = subprocess.run(args, cwd=self.directory, env=self._env(), capture_output=True, text=True)
        self.last_output = res.stdout
        with open(os.path.join(self.directory, "xtb.out"), "w") as f:
            f.write(res.stdout)
        engrad = os.path.join(self.directory, "geom.engrad")
        if res.returncode != 0 or "normal termination" not in res.stderr + res.stdout or not os.path.exists(engrad):
            tail = "\n".join((res.stdout + res.stderr).splitlines()[-20:])
            raise CalculationFailed(f"xtb failed in {self.directory} (exit {res.returncode}):\n{tail}")

        energy, grad = read_engrad(engrad)
        self.results["energy"] = energy * Hartree
        self.results["forces"] = -grad * (Hartree / Bohr)
        if self.molden_dest:
            os.makedirs(os.path.dirname(os.path.abspath(self.molden_dest)), exist_ok=True)
            shutil.move(os.path.join(self.directory, "molden.input"), self.molden_dest)
            self.molden_dest = None
