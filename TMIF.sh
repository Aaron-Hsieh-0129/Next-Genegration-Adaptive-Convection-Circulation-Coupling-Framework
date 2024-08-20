#!/bin/bash
#SBATCH -N 1
#SBATCH -n 22
#SBATCH --nodelist=node01
#SBATCH -o log/SP_old_CSSWM/AB2_prof_RKM_dt900_3_csswm_3_vvm_2E5diff_lineless_3B.o
#SBATCH -e log/SP_old_CSSWM/error_log/AB2_prof_RKM_dt900_3_csswm_3_vvm_2E5diff_lineless_3B.e
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
