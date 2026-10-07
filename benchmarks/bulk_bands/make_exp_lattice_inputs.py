"""CP2K inputs at the experimental lattice constant with the real-space KS/S matrices (spin-orbit gaps).

The inputs/<name>/cp2k_job.in (CIF lattice) of the MATERIAL_DB materials are scaled isotropically to
a_exp (MATERIAL_DB index 2, as in lattice_scan.json) and print the real-space matrices and the TREXIO
file, so that qdex.bulk_soc gives the spin-free and the PBE+SOC gaps at a_exp: the lattice of the
literature QSGW gaps (MATERIAL_DB_BULK_PBE gap_soc_exp_lattice).

    python make_exp_lattice_inputs.py OUT_DIR
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_lattice_scan import SLURM, scaled_input  # noqa: E402

if __name__ == "__main__":
    dest = sys.argv[1]
    scan = json.load(open(os.path.join(HERE, "lattice_scan.json")))
    for key, t in scan.items():
        if t.get("a_exp") is None:
            continue
        name = t["name"]
        src = open(os.path.join(HERE, "inputs", name, "cp2k_job.in")).read()
        s = t["a_exp"] / t["a_cif"]
        inp = scaled_input(src, name, "exp", s, band_order=True)
        inp = inp.replace(f"PROJECT {name}_exp", f"PROJECT {name}").replace(f"FILE_NAME {name}_exp.bs", f"FILE_NAME {name}.bs")
        d = os.path.join(dest, name)
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, "cp2k_job.in"), "w").write(inp)
        open(os.path.join(d, "cp2k.slurm"), "w").write(SLURM.format(job=f"bse_{name[:10]}"))
        print(f"{name:14s} scale {s:.5f}  a_exp {t['a_exp']:.4f}")
