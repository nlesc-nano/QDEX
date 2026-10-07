"""CP2K band-structure inputs of the fuzzy-band overlay at the PBE lattice constant.

The inputs/<name>/cp2k_job.in (CIF lattice) are scaled isotropically by scale_pbe of lattice_scan.json
(energy-volume minimum at the level of theory of the dots; wurtzite keeps c/a and u of the CIF).

    python make_pbe_lattice_inputs.py OUT_DIR
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_lattice_scan import SLURM, scaled_input  # noqa: E402
SMEAR = """      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.1
        NBROYDEN 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR"""

# gapless or inverted in PBE at their PBE lattice: Fermi smearing (InSb has it in inputs/ already)
SMEARED_PBE = {"GaAs_zb", "GaSb_zb"}

if __name__ == "__main__":
    dest = sys.argv[1]
    scan = json.load(open(os.path.join(HERE, "lattice_scan.json")))
    for key, t in scan.items():
        name = t["name"]
        src = open(os.path.join(HERE, "inputs", name, "cp2k_job.in")).read()
        inp = scaled_input(src, name, "pbe", t["scale_pbe"], band_order=True)
        inp = inp.replace(f"PROJECT {name}_pbe", f"PROJECT {name}").replace(f"FILE_NAME {name}_pbe.bs", f"FILE_NAME {name}.bs")
        if name in SMEARED_PBE:
            i = inp.index("      &MIXING")
            j = inp.index("      &END MIXING") + len("      &END MIXING")
            inp = inp[:i] + SMEAR + inp[j:]
        d = os.path.join(dest, name)
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, "cp2k_job.in"), "w").write(inp)
        open(os.path.join(d, "cp2k.slurm"), "w").write(SLURM.format(job=f"bsp_{name[:10]}"))
        print(f"{name:14s} scale {t['scale_pbe']:.5f}  a_PBE {t['a_pbe']:.4f}")
