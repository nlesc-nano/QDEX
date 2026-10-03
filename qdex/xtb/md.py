"""
g-xTB molecular dynamics with ASE, writing one QDEX NAMD frame per production step.

Equilibration: Langevin thermostat. Production: velocity Verlet (NVE, classical path
approximation). In production every force call also writes the molden file, which is
converted at once into ``frame_XXXXXX/{frame.xyz, BASIS_GXTB, MOs.mbse}``, so forces and
orbitals of a frame come from the same calculation.

Usage:
    python -m qdex.xtb.md start.xyz --outdir run --temp 300 --dt 2 --n-equil 500 --n-prod 50
"""
import argparse
import csv
import os
import time

import numpy as np
from ase import units
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from ase.md.verlet import VelocityVerlet

from qdex.xtb.calculator import GXTB
from qdex.xtb.molden import convert_molden
import logging

logger = logging.getLogger(__name__)


class _Logger:
    """Per-step CSV log: time, Epot, Ekin, Etot (eV), T (K), wall time of the step (s)."""

    def __init__(self, atoms, path, dt_fs, t0_fs=0.0):
        self.atoms, self.dt, self.t0 = atoms, dt_fs, t0_fs
        self.step = 0
        self.wall = time.time()
        new = not os.path.exists(path)
        self.f = open(path, "a", newline="")
        self.w = csv.writer(self.f)
        if new:
            self.w.writerow(["time_fs", "epot_ev", "ekin_ev", "etot_ev", "temp_k", "wall_s"])

    def __call__(self):
        a = self.atoms
        epot, ekin = a.get_potential_energy(), a.get_kinetic_energy()
        now = time.time()
        self.w.writerow([f"{self.t0 + self.step * self.dt:.3f}", f"{epot:.8f}", f"{ekin:.8f}",
                         f"{epot + ekin:.8f}", f"{a.get_temperature():.3f}", f"{now - self.wall:.2f}"])
        self.f.flush()
        self.wall = now
        self.step += 1


def make_calculator(workdir, nthreads=4, binary=None, libdir=None, xtbpath=None, acc=1.0, etemp=300.0, charge=0):
    return GXTB(binary=binary, libdir=libdir, xtbpath=xtbpath, directory=workdir, nthreads=nthreads,
                acc=acc, etemp=etemp, charge=charge)


def equilibrate(atoms, outdir, n_steps, temp_k=300.0, dt_fs=2.0, friction_per_fs=0.01, seed=1):
    """Langevin equilibration from Maxwell-Boltzmann velocities. Writes equil.traj and equil_log.csv."""
    rng = np.random.default_rng(seed)
    if not atoms.get_momenta().any():
        MaxwellBoltzmannDistribution(atoms, temperature_K=temp_k, rng=rng)
        Stationary(atoms)
        ZeroRotation(atoms)
    dyn = Langevin(atoms, timestep=dt_fs * units.fs, temperature_K=temp_k,
                   friction=friction_per_fs / units.fs, fixcm=True, rng=rng)
    traj = Trajectory(os.path.join(outdir, "equil.traj"), "w", atoms)
    dyn.attach(_Logger(atoms, os.path.join(outdir, "equil_log.csv"), dt_fs), interval=1)
    dyn.attach(traj.write, interval=10)
    dyn.run(n_steps)
    traj.close()
    write(os.path.join(outdir, "equil_final.traj"), atoms)
    return atoms


def production(atoms, outdir, n_frames, dt_fs=2.0, first_frame=1, nthreads=1, keep_molden=False, tol=1e-4):
    """
    NVE velocity Verlet; frame k (1-based) is written to ``outdir/frames/frame_{k:06d}``.

    The first frame is the starting geometry. Each force call writes the molden file of the
    new positions, which is converted and (unless ``keep_molden``) deleted.
    """
    calc = atoms.calc
    frames = os.path.join(outdir, "frames")
    os.makedirs(frames, exist_ok=True)
    log = _Logger(atoms, os.path.join(outdir, "prod_log.csv"), dt_fs)
    traj = Trajectory(os.path.join(outdir, "prod.traj"), "a", atoms)

    def frame_dir(k):
        return os.path.join(frames, f"frame_{k:06d}")

    def store(k, t_fs):
        d = frame_dir(k)
        molden = os.path.join(d, "molden.input")
        info = convert_molden(molden, d, nthreads=nthreads, tol=tol,
                              comment=f"t = {t_fs:.3f} fs  E = {atoms.get_potential_energy():.8f} eV")
        if not keep_molden:
            os.remove(molden)
        np.save(os.path.join(d, "velocities.npy"), atoms.get_velocities())
        traj.write(atoms)
        log()
        logger.info(f"  frame {k:6d}  t = {t_fs:8.2f} fs  T = {atoms.get_temperature():7.1f} K  "
                    f"gap = {info['gap_ev']:.4f} eV  |C^T S C - 1| = {info['orthonormality_error']:.1e}")

    # Frame 1: single point with molden at the starting geometry
    calc.molden_dest = os.path.join(frame_dir(first_frame), "molden.input")
    calc.reset()
    atoms.get_forces()
    store(first_frame, 0.0)

    dyn = VelocityVerlet(atoms, timestep=dt_fs * units.fs)
    for i in range(1, n_frames):
        calc.molden_dest = os.path.join(frame_dir(first_frame + i), "molden.input")
        dyn.run(1)
        store(first_frame + i, i * dt_fs)
    traj.close()
    write(os.path.join(outdir, "prod_final.traj"), atoms)
    return atoms


def main(argv=None):
    ap = argparse.ArgumentParser(description="g-xTB MD with ASE that writes QDEX NAMD frames.")
    ap.add_argument("start", help="Starting geometry (xyz) or an ASE .traj with velocities (restart)")
    ap.add_argument("--outdir", default="gxtb_md")
    ap.add_argument("--temp", type=float, default=300.0, help="Temperature (K)")
    ap.add_argument("--dt", type=float, default=2.0, help="Time step (fs)")
    ap.add_argument("--n-equil", type=int, default=500, help="Langevin equilibration steps (0 to skip)")
    ap.add_argument("--n-prod", type=int, default=50, help="NVE production frames")
    ap.add_argument("--friction", type=float, default=0.01, help="Langevin friction (1/fs)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--acc", type=float, default=1.0, help="xtb SCF accuracy (--acc)")
    ap.add_argument("--etemp", type=float, default=300.0, help="Electronic temperature (K)")
    ap.add_argument("--nthreads", type=int, default=8)
    ap.add_argument("--binary", default=None, help="xtb binary with g-xTB (default $GXTB)")
    ap.add_argument("--libdir", default=None, help="Runtime library dir for the binary (default $XTB_LIBDIR)")
    ap.add_argument("--keep-molden", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    os.makedirs(args.outdir, exist_ok=True)
    atoms = read(args.start)
    atoms.calc = make_calculator(os.path.join(args.outdir, "calc"), args.nthreads, args.binary, args.libdir,
                                 acc=args.acc, etemp=args.etemp, charge=args.charge)
    if args.n_equil > 0:
        logger.info(f"Equilibration: {args.n_equil} Langevin steps of {args.dt} fs at {args.temp} K")
        equilibrate(atoms, args.outdir, args.n_equil, args.temp, args.dt, args.friction, args.seed)
    if args.n_prod > 0:
        logger.info(f"Production: {args.n_prod} NVE frames of {args.dt} fs")
        production(atoms, args.outdir, args.n_prod, args.dt, nthreads=args.nthreads, keep_molden=args.keep_molden)


if __name__ == "__main__":
    main()
