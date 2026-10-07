#!/bin/bash
#SBATCH -J gw_CdSe
#SBATCH -t 04:00:00
#SBATCH -p preemption
#SBATCH -N 2
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G

module purge
module load CP2K/2026.2-foss-2025a
export OMP_NUM_THREADS=8
mpirun -np ${SLURM_NTASKS} --map-by slot:PE=8 cp2k.psmp -i cp2k_job.in -o cp2k_job.out
