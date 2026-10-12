# PBE0 bulk reference of zinc-blende CdSe

Bands and gaps behind `qdex/data/bulk_bands/CdSe_zb_pbe0*` and `MATERIAL_DB_BULK_DFT["pbe0"]` in
`qdex/hardness.py`. Used by QDEX runs with `system.functional: pbe0` (orbitals from PBE0 single points at the PBE
geometries of the dots). Same level of theory as the PBE files (CP2K 2026.2 built in NHR `$WORK/sw`, GPW 400 Ry,
EPS_DEFAULT 1e-12, DZVP-MOLOPT-PBE-GTH, GTH-PBE, 8x8x8), with PBE0 through ADMM2 (`admm-dzp`) and k-point RI-HFX.
Runs: NHR `$WORK/bulk_pbe0/CdSe_zb`.

1. `mk_bulk.py NAME A RC` writes `pbe.in` (PBE k-point SCF) and `cp2k_job.in` (PBE0 restarted from it) from
   `bulk_pbe_template.in`; `run.slurm` runs both on one genoa node (96 ranks x 2 threads) and then `qdex.bulk_soc`.
2. `CdSe_zb_pbe0/`: bands at the PBE lattice (6.217 A), truncated Coulomb R_c 8 A. `CdSe_zb_pbe0_aexp/`: at a_exp.

Settings that keep k-point RI-HFX from stalling (CP2K manual, "HFX-RI with k-Points"): `EPS_PGF_ORB 1e-5` in `&QS`
and `&RI` (a tight value multiplies the periodic images, and the MOLOPT basis is diffuse), the default RI metric,
`KP_NGROUPS 16`, ADMM, a restart from the converged PBE k-point run, and `KP_RI_PROGRESS_BAR` /
`KP_RI_MEMORY_ESTIMATE` to see progress. A first attempt without ADMM, from the atomic guess and with
`EPS_PGF_ORB 1e-6` and the overlap metric, made no progress in 35 min. One PBE0 SCF at a_PBE takes ~1.5 h (20
steps, ~240 s each), at a_exp ~2.3 h (~400 s per step, more images).

| | spin-free gap | +SOC gap | Cd 4d below VBM |
|---|---|---|---|
| PBE0 at a_PBE, R_c 5 A | 1.586 | 1.466 | |
| PBE0 at a_PBE, R_c 8 A | 1.617 | 1.498 | 8.15 |
| PBE0 at a_exp, R_c 8 A | 1.823 | 1.704 | 8.34 |
| PBE at a_PBE | 0.474 | 0.352 | 7.51 |

Delta_Sigma(PBE0) = QSGW+SOC 2.164 - PBE0+SOC 1.704 = 0.460 eV (PBE: 1.642 eV).

Cross-check without RI (`mk_bulk_sc.py`): 4x4x4 conventional supercell (512 atoms) at Gamma, ADMM2 with 4-centre
truncated-Coulomb exchange: spin-free Gamma gap 1.623 (R_c 6 A) / 1.652 eV (R_c 9 A) at a_PBE, 1.863 eV at a_exp,
within 35 meV of the k-point runs.
