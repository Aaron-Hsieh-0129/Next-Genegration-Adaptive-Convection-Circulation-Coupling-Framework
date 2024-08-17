#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --nodelist=node01
#SBATCH -o log/SP_new/prof_RKM_dt1800_2_csswm_2_vvm_2E5diff_cross.o
#SBATCH -e log/SP_new/prof_RKM_dt1800_2_csswm_2_vvm_2E5diff_cross.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b08209006@ntu.edu.tw

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi
export OMP_NUM_THREADS=$omp_threads
echo $OMP_NUM_THREADS

rm -rf build
mkdir build
cd build/ && cmake ../ && make -j 4 && ./TMIF
