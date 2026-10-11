"""HLE17 inputs for zinc-blende InAs from the PBE bands input (../../inputs/InAs_zb/cp2k_job.in).

scan/<s>: CIF lattice x s (s = 0.970 ... 1.030) and scan/exp (a_exp = 6.060 A); InAs_zb_hle17/ is written by
`python make_inputs.py bands <a>` once the scan has been fitted. Run on NHR from $WORK/bulk_hle17/InAs_zb.
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PBE = os.path.join(HERE, "cp2k_job_pbe.in")
A_CIF = 2 * 3.05356166          # cubic lattice of the PBE input cell (fcc primitive vectors a/2 (0,1,1) ...)
A_EXP = 6.060           # MATERIAL_DB a_exp, the lattice of the QSGW reference


def hle17_input(scale, project):
    s = open(PBE).read()
    s = s.replace("/scratch/iinfante/cp2k_basis", "/scratch/usr/becivinf/cp2k_basis")
    s = s.replace("PROJECT InAs_zb", f"PROJECT {project}").replace("FILE_NAME InAs_zb.bs", f"FILE_NAME {project}.bs")
    s = s.replace("&XC_FUNCTIONAL PBE", "&XC_FUNCTIONAL\n        &MGGA_XC_HLE17\n        &END MGGA_XC_HLE17")
    head, sub = s.split("&SUBSYS", 1)
    num = lambda m: f"{float(m.group(0)) * scale:.8f}"
    cellcoord = re.sub(r"(?s)(&CELL.*?&END COORD)", lambda m: re.sub(r"-?\d+\.\d+", num, m.group(1)), sub)
    return head + "&SUBSYS" + cellcoord


def write(d, scale, project):
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "cp2k_job.in"), "w").write(hle17_input(scale, project))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "bands":
        a = float(sys.argv[2])
        write("InAs_zb_hle17", a / A_CIF, "InAs_zb_hle17")
        print(f"InAs_zb_hle17: a = {a:.4f} A (CIF x {a / A_CIF:.5f})")
    else:
        for s in (0.970, 0.985, 1.000, 1.015, 1.030):
            write(f"scan/s{s:.3f}", s, f"InAs_zb_s{s:.3f}")
        write("scan/exp", A_EXP / A_CIF, "InAs_zb_exp")
