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

## Zinc-blende InAs (`InAs_zb/`, 10 Oct 2026)

Same protocol for the HLE17 InAs dots of Hyun (In 4d 15.7 eV below the dot HOMO, PBE bulk 14.31 eV below the VBM).
`make_inputs.py` writes the scan (CIF lattice 6.1071 Å × 0.970 … 1.030, plus a_exp = 6.060 Å, the lattice of the
QSGW reference) from the PBE input; `run.slurm` runs the scan, fits it (`fit_scan.py`), runs the bands at a_HLE17 and
`qdex.bulk_soc` on both. One test-queue job, 7 min. Runs: NHR `$WORK/bulk_hle17/InAs_zb`.

a_HLE17 = 6.053 Å (B0 57.4 GPa, rms 0.11 meV; PBE 6.204, exp 6.06).

| | Γ6 − Γ8 spin-free | +SOC | In 4d below VBM |
|---|---|---|---|
| HLE17 at a_HLE17 | +0.044 | −0.073 | 16.08 |
| HLE17 at a_exp | +0.023 | −0.094 | 16.07 |
| PBE at a_PBE | −0.603 | −0.713 | 14.31 |
| PBE at a_exp | −0.247 | −0.360 | |

Delta_Sigma(HLE17) = QSGW+SOC 0.675 − HLE17+SOC (−0.094) = 0.769 eV (PBE: 1.035 eV).
