"""Bulk zinc-blende CdSe with PBE0 as a Gamma-point supercell (no RI): 4x4x4 conventional cell, 512 atoms, ADMM2 with
4-centre truncated-Coulomb exchange (the dot setup), GPW 400 Ry, DZVP-MOLOPT-PBE-GTH. The supercell Gamma point folds
in Gamma, X and L of the primitive cell. pbe.in (PBE OT, atomic guess) -> pbe0.in (PBE0 OT from it) -> diag.in (one
diagonalization, 40 virtuals, eigenvalues).     python mk_bulk_sc.py NAME A_ANGSTROM RC_ANGSTROM [N=4]"""
import os, sys
import numpy as np

name, a, rc = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
N = int(sys.argv[4]) if len(sys.argv) > 4 else 4
W = os.environ["WORK"]
B = f"{W}/CP2K/CdSe_dots"
fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
atoms = [(s, (f + off + np.array(c)) * a) for c in np.ndindex(N, N, N)
         for s, off in (("Cd", np.zeros(3)), ("Se", np.full(3, .25))) for f in fcc]
L = N * a
assert rc < L / 2, "CUTOFF_RADIUS must stay below half the supercell"
coord = "".join(f"      {s} {x:.8f} {y:.8f} {z:.8f}\n" for s, (x, y, z) in atoms)


def inp(project, pbe0, scf, prnt, restart):
    kinds = "".join(f"""    &KIND {e}
      BASIS_SET ORB DZVP-MOLOPT-PBE-GTH-{q}
{f'      BASIS_SET AUX_FIT admm-dzp-{q}' + chr(10) if pbe0 else ''}      POTENTIAL GTH-PBE-{q}
    &END KIND
""" for e, q in (("Cd", "q12"), ("Se", "q6")))
    xc = f"""    &XC
      &XC_FUNCTIONAL
        &PBE
          SCALE_X 0.75
          SCALE_C 1.0
        &END PBE
      &END XC_FUNCTIONAL
      &HF
        FRACTION 0.25
        &SCREENING
          EPS_SCHWARZ 1.0E-6
          SCREEN_ON_INITIAL_P TRUE
        &END SCREENING
        &MEMORY
          MAX_MEMORY 20000
          EPS_STORAGE_SCALING 0.1
        &END MEMORY
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE TRUNCATED
          CUTOFF_RADIUS {rc}
          T_C_G_DATA {W}/sw/cp2k-2026.2/data/t_c_g.dat
        &END INTERACTION_POTENTIAL
      &END HF
    &END XC
    &AUXILIARY_DENSITY_MATRIX_METHOD
      ADMM_TYPE ADMM2
    &END AUXILIARY_DENSITY_MATRIX_METHOD
""" if pbe0 else "    &XC\n      &XC_FUNCTIONAL PBE\n      &END XC_FUNCTIONAL\n    &END XC\n"
    rst = f"    WFN_RESTART_FILE_NAME {restart}\n" if restart else ""
    return f"""&GLOBAL
  PROJECT {project}
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME {B}/basis_2024.3/BASIS_MOLOPT_UZH
    BASIS_SET_FILE_NAME {W}/cp2k_basis/BASIS_ADMM_UZH
    POTENTIAL_FILE_NAME {B}/basis_2024.3/POTENTIAL_UZH
{rst}    &MGRID
      CUTOFF 400
      NGRIDS 4
    &END MGRID
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-12
    &END QS
    &POISSON
      PERIODIC XYZ
      POISSON_SOLVER PERIODIC
    &END POISSON
{scf}{xc}{prnt}  &END DFT
  &SUBSYS
    &CELL
      ABC {L:.8f} {L:.8f} {L:.8f}
      PERIODIC XYZ
    &END CELL
    &COORD
{coord}    &END COORD
{kinds}  &END SUBSYS
&END FORCE_EVAL
"""


def ot(guess, precond):
    return f"""    &SCF
      SCF_GUESS {guess}
      EPS_SCF 1.0E-6
      MAX_SCF 200
      IGNORE_CONVERGENCE_FAILURE
      &OT
        MINIMIZER DIIS
        N_DIIS 7
        PRECONDITIONER {precond}
      &END OT
    &END SCF
"""


DIAG = """    &SCF
      SCF_GUESS RESTART
      ADDED_MOS 40
      EPS_SCF 1.0E-6
      MAX_SCF 1
      IGNORE_CONVERGENCE_FAILURE
    &END SCF
"""
PRINT = """    &PRINT
      &MO
        EIGENVALUES T
        OCCUPATION_NUMBERS F
        COEFFICIENTS F
        NDIGITS 8
      &END MO
    &END PRINT
"""
os.makedirs(name, exist_ok=True)
open(f"{name}/pbe.in", "w").write(inp("pbe", False, ot("ATOMIC", "FULL_KINETIC"), "", None))
open(f"{name}/pbe0.in", "w").write(inp("pbe0", True, ot("RESTART", "FULL_SINGLE_INVERSE"), "", "pbe-RESTART.wfn"))
open(f"{name}/diag.in", "w").write(inp("pbe0_diag", True, DIAG, PRINT, "pbe0-RESTART.wfn"))
print(f"{name}: {len(atoms)} atoms, L = {L:.3f} A, R_c = {rc} A")
