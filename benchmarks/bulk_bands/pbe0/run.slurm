#!/bin/bash
#SBATCH -J pbe0_bulk
#SBATCH -p cpu-genoa
#SBATCH -N 1
#SBATCH --ntasks-per-node=96
#SBATCH -c 2
#SBATCH -t 12:00:00
#SBATCH -o bulk-%j.out
# PBE0 bands of zinc-blende CdSe (k-point RI-HFX + ADMM2) on the full node: PBE k-point SCF -> PBE0 restart ->
# qdex.bulk_soc.   sbatch run_bulk.slurm DIR [DIR ...]
source $WORK/sw/cp2k-2026.2/install/cp2k_env
export OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close OMP_STACKSIZE=256M
ulimit -s unlimited
C=$WORK/bulk_pbe0/CdSe_zb; B=$WORK/cp2k_basis
for d in "$@"; do
  cd $C/$d
  [ -f CdSe_zb_pbe-RESTART.kp ] || { srun --cpu-bind=cores cp2k.psmp -i pbe.in -o pbe.out; echo "$(date +%T) $d PBE done"; }
  t0=$(date +%s); srun --cpu-bind=cores cp2k.psmp -i cp2k_job.in -o cp2k_job.out
  echo "$(date +%T) $d PBE0: $(( $(date +%s) - t0 )) s, converged $(grep -c 'SCF run converged' cp2k_job.out)"
done
source $HOME/miniforge3/etc/profile.d/conda.sh && conda activate qdex
cd $WORK/QDEX
for d in "$@"; do
  python -m qdex.bulk_soc $C/$d --basis $B/BASIS_MOLOPT_UZH --gth $B/GTH_SOC_POTENTIALS --name CdSe_zb_pbe0_$d \
         --cif $WORK/bulk_hle17/CdSe_zb/CdSe_zb.cif -d $C/$d/out > $C/$d/bulk_soc.log 2>&1; echo "$d bulk_soc exit $?"
done
