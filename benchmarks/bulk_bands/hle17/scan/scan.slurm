#!/bin/bash
#SBATCH -J hle17_bulk_scan
#SBATCH -p cpu-genoa:test
#SBATCH -N 1
#SBATCH --ntasks-per-node=192
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH -o scan-%j.out
# HLE17 energy-volume scan of zinc-blende CdSe (6 lattices in parallel, 32 ranks each)
source $WORK/sw/cp2k-2026.2/install/cp2k_env
export OMP_NUM_THREADS=1
for d in exp s0.970 s0.985 s1.000 s1.015 s1.030; do
  (cd $d && srun --exact -n 32 -c 1 --cpu-bind=cores cp2k.psmp -i cp2k_job.in -o cp2k_job.out; echo "$d done") &
done
wait
