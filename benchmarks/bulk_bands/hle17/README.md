# HLE17 bulk reference of zinc-blende CdSe

Bands and gaps behind `qdex/data/bulk_bands/CdSe_zb_hle17*` and `MATERIAL_DB_BULK_DFT["hle17"]` in
`qdex/hardness.py`. Used by QDEX runs with `system.functional: hle17`. Same level of theory as the PBE files (CP2K
2026.2 built in NHR `$WORK/sw`, GPW 400 Ry, EPS_DEFAULT 1e-12, DZVP-MOLOPT-PBE-GTH, GTH-PBE, 8x8x8), with
`&MGGA_XC_HLE17` in place of PBE. Runs: NHR `$WORK/bulk_hle17/CdSe_zb`.

1. `scan/`: energy-volume scan, CIF lattice × 0.970 … 1.030 plus a_exp = 6.05 Å (all six in one job, 32 ranks
   each). The Birch-Murnaghan fit gives a_HLE17 = 6.108 Å, B0 = 51.4 GPa (rms 0.03 meV), against PBE 6.217 Å.
2. `CdSe_zb_hle17/`: bands at a_HLE17 (CIF × 0.99472), real-space KS/S matrices and TREXIO. `run.slurm` runs it,
   then `qdex.bulk_soc` on it and on `scan/exp`.

   Use `srun --exact -n 32 -c 1 --cpu-bind=cores`. A plain `srun -n 64` made the first SCF step take 1856 s
   instead of 3 s.

| | spin-free gap | +SOC gap | Cd 4d below VBM |
|---|---|---|---|
| HLE17 at a_HLE17 | 1.452 | 1.322 | 8.63 |
| HLE17 at a_exp | 1.535 | 1.405 | |
| PBE at a_PBE | 0.474 | 0.352 | 7.51 |
| PBE at a_exp | 0.644 | 0.522 | |

Delta_Sigma(HLE17) = QSGW+SOC 2.164 − HLE17+SOC 1.405 = 0.759 eV (PBE: 1.642 eV).
