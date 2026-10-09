#!/bin/bash
#SBATCH -J hle17_bulk_bands
#SBATCH -p cpu-genoa:test
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH -o bands-%j.out
# HLE17 bands of zinc-blende CdSe at the HLE17 lattice (6.108 A), then qdex.bulk_soc for it and for a_exp (scan/exp)
source $WORK/sw/cp2k-2026.2/install/cp2k_env
export OMP_NUM_THREADS=1
srun --exact -n 32 -c 1 --cpu-bind=cores cp2k.psmp -i cp2k_job.in -o cp2k_job.out
source $HOME/miniforge3/etc/profile.d/conda.sh && conda activate qdex
B=$WORK/cp2k_basis; C=$WORK/bulk_hle17/CdSe_zb
cd $WORK/QDEX
python -m qdex.bulk_soc $C/CdSe_zb_hle17 --basis $B/BASIS_MOLOPT_UZH --gth $B/GTH_SOC_POTENTIALS --name CdSe_zb_hle17 --cif $C/CdSe_zb.cif -d $C/out --gzip
python -m qdex.bulk_soc $C/scan/exp --basis $B/BASIS_MOLOPT_UZH --gth $B/GTH_SOC_POTENTIALS --name CdSe_zb_hle17_exp --cif $C/CdSe_zb.cif -d $C/out_exp
