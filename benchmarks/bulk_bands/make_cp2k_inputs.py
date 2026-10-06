"""CP2K bulk band-structure inputs for the bulk CIFs of the QDSpaceWebApp.

Same primitive cell and k-path as qdex.fuzzy_bands (pymatgen standard primitive cell, HighSymmKpath,
line density 50); prints the spin-free bands, the real-space KS/S matrices and a TREXIO file for
qdex.bulk_soc. Usage:

    python make_cp2k_inputs.py CIF_ROOT OUT_DIR [BASIS_DIR]

CIF_ROOT holds II-VI/, III-V/, IV-VI/, ABX3/ with bulk_cifs/*.cif; BASIS_DIR is where CP2K finds
BASIS_MOLOPT_UZH and POTENTIAL_UZH on the cluster.
"""
import json, os, sys
import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

CIF_ROOT = sys.argv[1]
OUT = sys.argv[2]
BASIS_DIR = sys.argv[3] if len(sys.argv) > 3 else "/scratch/iinfante/cp2k_basis"
CIF_DIRS = [os.path.join(CIF_ROOT, f, "bulk_cifs") for f in ("II-VI", "III-V", "IV-VI", "ABX3")]
from qdex.constants import valence_electrons as Q
KSPACING = 0.22          # 1/A (2pi included); gives 8x8x8 for zinc-blende CdSe
# PBE semimetals whose SCF needs Fermi smearing (300 K) and slower Broyden mixing
SMEARED = {"InSb_zb"}
SMEAR = """      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.1
        NBROYDEN 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR"""

TEMPLATE = """&GLOBAL
  PROJECT {name}
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME {basis_dir}/BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME {basis_dir}/POTENTIAL_UZH
    &MGRID
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
    &KPOINTS
      SCHEME MONKHORST-PACK {kgrid}
      SYMMETRY ON
      EPS_GEO 1.0E-5
      FULL_GRID .TRUE.
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-7
      MAX_SCF 150
      ADDED_MOS {added}
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.3
      &END MIXING
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &PRINT
      &BAND_STRUCTURE
        FILE_NAME {name}.bs
{ksets}
      &END BAND_STRUCTURE
      &KS_CSR_WRITE
        REAL_SPACE .TRUE.
        THRESHOLD 1.0E-14
        UPPER_TRIANGULAR .FALSE.
        BINARY .FALSE.
      &END KS_CSR_WRITE
      &S_CSR_WRITE
        REAL_SPACE .TRUE.
        THRESHOLD 1.0E-14
        UPPER_TRIANGULAR .FALSE.
        BINARY .FALSE.
      &END S_CSR_WRITE
      &TREXIO
        FILENAME {name}_trexio
      &END TREXIO
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      A {A}
      B {B}
      C {C}
      PERIODIC XYZ
    &END CELL
    &COORD
{coords}
    &END COORD
{kinds}
  &END SUBSYS
&END FORCE_EVAL
"""

SLURM = """#!/bin/bash
#SBATCH -J bs_{name}
#SBATCH -t 00:09:00
#SBATCH -n 48
#SBATCH --mem-per-cpu=2g
#SBATCH --qos=regular

module purge
module load CP2K/2026.2-foss-2025a
export OMP_NUM_THREADS=1
mpirun -np ${{SLURM_NTASKS}} cp2k.popt -i cp2k_job.in -o cp2k_job.out
"""


def kpoint_sets(prim):
    kpts, labels = HighSymmKpath(prim).get_kpoints(line_density=50, coords_are_cartesian=False)
    kpts = np.asarray(kpts)
    ticks = [i for i, l in enumerate(labels) if l]
    sets = []
    for a, b in zip(ticks[:-1], ticks[1:]):
        if b - a < 2:                       # duplicated end point or a path break
            continue
        sets.append((kpts[a], kpts[b], b - a, labels[a], labels[b]))
    out = []
    for k0, k1, n, l0, l1 in sets:
        out.append(f"""        &KPOINT_SET
          UNITS B_VECTOR
          NPOINTS {n}
          SPECIAL_POINT {k0[0]:.8f} {k0[1]:.8f} {k0[2]:.8f}  # {l0}
          SPECIAL_POINT {k1[0]:.8f} {k1[1]:.8f} {k1[2]:.8f}  # {l1}
        &END KPOINT_SET""")
    return "\n".join(out), len(kpts)


meta = {}
import glob
for path in sorted(p for d in CIF_DIRS for p in glob.glob(os.path.join(d, "*.cif"))):
    name = os.path.basename(path)[:-4]
    s = Structure.from_file(path)
    sga = SpacegroupAnalyzer(s)
    prim = sga.get_primitive_standard_structure()
    lat = prim.lattice.matrix
    rec = prim.lattice.reciprocal_lattice.matrix
    kgrid = [max(2, int(round(np.linalg.norm(b) / KSPACING))) for b in rec]
    syms = [site.specie.symbol for site in prim]
    n_el = sum(Q[e] for e in syms)
    ksets, nk = kpoint_sets(prim)
    coords = "\n".join(f"      {e} {x:.8f} {y:.8f} {z:.8f}" for e, (x, y, z) in zip(syms, prim.cart_coords))
    kinds = "\n".join(f"    &KIND {e}\n      BASIS_SET DZVP-MOLOPT-PBE-GTH-q{Q[e]}\n      POTENTIAL GTH-PBE-q{Q[e]}\n    &END KIND"
                      for e in sorted(set(syms)))
    d = os.path.join(OUT, name)
    os.makedirs(d, exist_ok=True)
    fmt = lambda v: " ".join(f"{x:.8f}" for x in v)
    inp = TEMPLATE.format(
        name=name, basis_dir=BASIS_DIR, kgrid=" ".join(map(str, kgrid)), added=max(20, n_el // 2), ksets=ksets,
        A=fmt(lat[0]), B=fmt(lat[1]), C=fmt(lat[2]), coords=coords, kinds=kinds)
    if name in SMEARED:
        i = inp.index("      &MIXING")
        j = inp.index("      &END MIXING") + len("      &END MIXING")
        inp = inp[:i] + SMEAR + inp[j:]
    open(os.path.join(d, "cp2k_job.in"), "w").write(inp)
    open(os.path.join(d, "cp2k.slurm"), "w").write(SLURM.format(name=name))
    meta[name] = dict(cif=path, a_conv=round(s.lattice.a, 5), spacegroup=sga.get_space_group_symbol(), n_atoms=len(prim), formula=prim.composition.reduced_formula,
                      kgrid=kgrid, n_electrons=n_el, n_path=nk)
    print(name, meta[name])
json.dump(meta, open(os.path.join(OUT, "materials.json"), "w"), indent=2)
