#!/bin/bash
#SBATCH -J hle17_InAs
#SBATCH -p cpu-genoa:test
#SBATCH -N 1
#SBATCH --ntasks-per-node=192
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH -o run-%j.out
# HLE17 bulk reference of zinc-blende InAs: energy-volume scan (6 lattices in parallel, 32 ranks each), BM fit,
# bands at a_HLE17, then qdex.bulk_soc at a_HLE17 and at a_exp (scan/exp).
source $WORK/sw/cp2k-2026.2/install/cp2k_env
export OMP_NUM_THREADS=1
C=$WORK/bulk_hle17/InAs_zb; cd $C
for d in exp s0.970 s0.985 s1.000 s1.015 s1.030; do
  (cd scan/$d && srun --exact -n 32 -c 1 --cpu-bind=cores cp2k.psmp -i cp2k_job.in -o cp2k_job.out; echo "$(date +%T) $d done") &
done
wait
source $HOME/miniforge3/etc/profile.d/conda.sh && conda activate qdex
python fit_scan.py | tee fit.txt
a=$(tail -1 fit.txt)
python make_inputs.py bands $a
(cd InAs_zb_hle17 && srun --exact -n 32 -c 1 --cpu-bind=cores cp2k.psmp -i cp2k_job.in -o cp2k_job.out; echo "$(date +%T) bands done")
B=$WORK/cp2k_basis
cd $WORK/QDEX
python -m qdex.bulk_soc $C/InAs_zb_hle17 --basis $B/BASIS_MOLOPT_UZH --gth $B/GTH_SOC_POTENTIALS --name InAs_zb_hle17 --cif $C/InAs_zb.cif -d $C/out --gzip
python -m qdex.bulk_soc $C/scan/exp --basis $B/BASIS_MOLOPT_UZH --gth $B/GTH_SOC_POTENTIALS --name InAs_zb_hle17_exp --cif $C/InAs_zb.cif -d $C/out_exp
