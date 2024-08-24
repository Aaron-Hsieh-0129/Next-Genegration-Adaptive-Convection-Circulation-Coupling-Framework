#!/bin/bash
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --nodelist=mogamd
#SBATCH -o log/nudge/dt600_1_csswm_1_vvm_2E5diff_5vvm_4B_1non_10kmcouple_Q1inter_30K.o
#SBATCH -e log/nudge/error_log/dt600_1_csswm_1_vvm_2E5diff_5vvm_4B_1non_10kmcouple_Q1inter_30K.e
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
