"""CP2K inputs for the bulk PBE gaps of MATERIAL_DB at the experimental and at the PBE lattice constant.

For every material of MATERIAL_DB with a bulk CIF in inputs/<name>/, the input there (CIF lattice) is
scaled isotropically to
  - the experimental lattice constant (MATERIAL_DB index 2), tag "exp", and
  - a_CIF x (0.97, 0.985, 1.00, 1.015, 1.03), tags "s0.970" ... , an energy-volume scan whose
    Birch-Murnaghan minimum is the PBE lattice constant at the level of theory of the dots.
Materials whose PBE band order at Gamma is inverted or close to it also print the real-space KS/S
matrices and a TREXIO file (qdex.bulk_soc gives the s-p band order); the others only print the bands.

    python make_lattice_scan.py OUT_DIR
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from qdex.hardness import MATERIAL_DB  # noqa: E402

NAMES = {"CSPBCL3": "CsPbCl3_cubic", "CSPBBR3": "CsPbBr3_cubic", "CSPBI3": "CsPbI3_cubic",
         "ZNS": "ZnS_zb", "ZNSE": "ZnSe_zb", "ZNTE": "ZnTe_zb", "CDS": "CdS_zb", "CDSE": "CdSe_zb",
         "CDTE": "CdTe_zb", "HGS": "HgS_zb", "HGSE": "HgSe_zb", "HGTE": "HgTe_zb", "ALP": "AlP_zb",
         "ALAS": "AlAs_zb", "ALSB": "AlSb_zb", "GAP": "GaP_zb", "GAAS": "GaAs_zb", "GASB": "GaSb_zb",
         "INP": "InP_zb", "INAS": "InAs_zb", "INSB": "InSb_zb", "PBS": "PbS_rs", "PBSE": "PbSe_rs"}
# in data/bulk_bands but not in MATERIAL_DB: energy-volume scan only (isotropic; wurtzite keeps c/a and u)
EXTRA = {"PBTE": "PbTe_rs", "CDSE_WZ": "CdSe_wz"}
BAND_ORDER = {"HgS_zb", "HgSe_zb", "HgTe_zb", "InAs_zb", "InSb_zb", "GaAs_zb", "GaSb_zb"}
SCAN = (0.970, 0.985, 1.000, 1.015, 1.030)

SLURM = """#!/bin/bash
#SBATCH -J {job}
#SBATCH -t 00:09:00
#SBATCH -n 48
#SBATCH --mem-per-cpu=2g
#SBATCH --qos=regular

module purge
module load CP2K/2026.2-foss-2025a
export OMP_NUM_THREADS=1
mpirun -np ${{SLURM_NTASKS}} cp2k.popt -i cp2k_job.in -o cp2k_job.out
"""


def scaled_input(src, name, tag, s, band_order):
    out = []
    block = None
    for line in src.splitlines():
        t = line.strip()
        if t in ("&CELL", "&COORD"):
            block = t
        elif t in ("&END CELL", "&END COORD"):
            block = None
        elif block == "&CELL" and re.match(r"^[ABC] ", t):
            k, *v = t.split()
            line = "      " + k + " " + " ".join(f"{float(x) * s:.8f}" for x in v)
        elif block == "&COORD" and re.match(r"^[A-Z][a-z]? ", t):
            e, *v = t.split()
            line = "      " + e + " " + " ".join(f"{float(x) * s:.8f}" for x in v)
        out.append(line)
    inp = "\n".join(out) + "\n"
    inp = inp.replace(f"PROJECT {name}", f"PROJECT {name}_{tag}").replace(f"FILE_NAME {name}.bs", f"FILE_NAME {name}_{tag}.bs")
    if not band_order:                                          # bands only
        inp = re.sub(r"      &KS_CSR_WRITE.*?&END TREXIO\n", "", inp, flags=re.S)
    return inp


if __name__ == "__main__":
    dest = sys.argv[1]
    only = sys.argv[2:]
    for key, name in {**NAMES, **EXTRA}.items():
        if only and name not in only:
            continue
        src = open(os.path.join(HERE, "inputs", name, "cp2k_job.in")).read()
        a = re.search(r"\n\s+A ([-\d. ]+)\n", src).group(1).split()
        a_cif = (sum(float(x) ** 2 for x in a) ** 0.5)
        # conventional a: fcc primitive vector length a/sqrt2, simple cubic a
        a_conv = a_cif * 2 ** 0.5 if not name.endswith("_cubic") else a_cif
        runs = {"exp": MATERIAL_DB[key][2] / a_conv} if key in MATERIAL_DB else {}
        runs.update({f"s{s:.3f}": s for s in SCAN})
        for tag, s in runs.items():
            d = os.path.join(dest, name, tag)
            os.makedirs(d, exist_ok=True)
            open(os.path.join(d, "cp2k_job.in"), "w").write(scaled_input(src, name, tag, s, name in BAND_ORDER))
            open(os.path.join(d, "cp2k.slurm"), "w").write(SLURM.format(job=f"lat_{name[:6]}_{tag}"))
        print(f"{name:14s} a_CIF {a_conv:.4f}  a_exp {MATERIAL_DB[key][2] if key in MATERIAL_DB else float('nan'):.3f}"
              f"{'  + KS/S matrices' if name in BAND_ORDER else ''}")
