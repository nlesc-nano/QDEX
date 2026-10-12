"""Bulk zinc-blende CdSe with PBE0 (k-point RI-HFX + ADMM2, CP2K 2026.2; settings of the CP2K manual page "HFX-RI with
k-Points": EPS_PGF_ORB 1e-5 everywhere (a tight value multiplies the periodic images), default RI metric, R_c ~ 5 A,
KP_NGROUPS 16, PBE k-point restart, progress bar and memory estimate), from the PBE bulk-band input (GPW 400 Ry,
EPS_DEFAULT 1e-12, DZVP-MOLOPT-PBE-GTH, 8x8x8, same band path, real-space KS/S for qdex.bulk_soc). Writes pbe.in
(PBE k-point SCF, the restart CP2K asks for) and cp2k_job.in (PBE0 from it).
    python mk_bulk.py NAME A_ANGSTROM RC_ANGSTROM     (in $WORK/bulk_pbe0/CdSe_zb)"""
import os, re, sys

name, a, rc = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
W = os.environ["WORK"]
s = open("bulk_pbe_template.in").read()
s = s.replace("/scratch/iinfante/cp2k_basis", f"{W}/cp2k_basis")
h = a / 2
s = re.sub(r"      A .*\n      B .*\n      C .*\n",
           f"      A -0.00000000 {h:.8f} {h:.8f}\n      B {h:.8f} 0.00000000 {h:.8f}\n      C {h:.8f} {h:.8f} 0.00000000\n", s)
s = re.sub(r"      Se .*\n", f"      Se {1.5 * h:.8f} {1.5 * h:.8f} {1.5 * h:.8f}\n", s)
s = s.replace("      EPS_DEFAULT 1.0E-12", "      EPS_DEFAULT 1.0E-12\n      EPS_PGF_ORB 1.0E-5")
s = s.replace("  PRINT_LEVEL LOW", "  PRINT_LEVEL MEDIUM")
os.makedirs(name, exist_ok=True)

# PBE k-point SCF without the print section (restart file only)
pbe = re.sub(r"    &PRINT\n.*?    &END PRINT\n", "", s, flags=re.S).replace("PROJECT CdSe_zb", "PROJECT CdSe_zb_pbe")
open(f"{name}/pbe.in", "w").write(pbe)

p0 = s.replace("PROJECT CdSe_zb", "PROJECT CdSe_zb_pbe0").replace("FILE_NAME CdSe_zb.bs", "FILE_NAME CdSe_zb_pbe0.bs")
p0 = p0.replace(f"    BASIS_SET_FILE_NAME {W}/cp2k_basis/BASIS_MOLOPT_UZH",
                f"    BASIS_SET_FILE_NAME {W}/cp2k_basis/BASIS_MOLOPT_UZH\n    BASIS_SET_FILE_NAME {W}/cp2k_basis/BASIS_ADMM_UZH\n"
                "    AUTO_BASIS RI_HFX MEDIUM\n    SORT_BASIS EXP\n    WFN_RESTART_FILE_NAME CdSe_zb_pbe-RESTART.kp")
p0 = p0.replace("      SCF_GUESS ATOMIC", "      SCF_GUESS RESTART")
xc = f"""    &XC
      &XC_FUNCTIONAL
        &PBE
          SCALE_X 0.75
          SCALE_C 1.0
        &END PBE
      &END XC_FUNCTIONAL
      &HF
        FRACTION 0.25
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE TRUNCATED
          CUTOFF_RADIUS {rc}
        &END INTERACTION_POTENTIAL
        &RI
          KP_NGROUPS 16
          EPS_PGF_ORB 1.0E-5
          &PRINT
            KP_RI_MEMORY_ESTIMATE
            KP_RI_PROGRESS_BAR
          &END PRINT
        &END RI
      &END HF
    &END XC
    &AUXILIARY_DENSITY_MATRIX_METHOD
      ADMM_TYPE ADMM2
    &END AUXILIARY_DENSITY_MATRIX_METHOD
"""
p0 = re.sub(r"    &XC\n.*?    &END XC\n", xc, p0, flags=re.S)
for e, q in (("Cd", "q12"), ("Se", "q6")):
    p0 = p0.replace(f"      BASIS_SET DZVP-MOLOPT-PBE-GTH-{q}\n",
                    f"      BASIS_SET ORB DZVP-MOLOPT-PBE-GTH-{q}\n      BASIS_SET AUX_FIT admm-dzp-{q}\n")
open(f"{name}/cp2k_job.in", "w").write(p0)
print(name, a, rc)
