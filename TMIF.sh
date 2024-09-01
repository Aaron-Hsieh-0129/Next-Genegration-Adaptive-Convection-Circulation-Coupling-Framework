#!/bin/bash
#SBATCH -N 1
#SBATCH -n 7
#SBATCH --nodelist=mogamd
#SBATCH -o log/0901_stable_version/dt600_1_csswm_1_vvm_2E5diff_7vvm_3B_4non_10kmcouple.o
#SBATCH -e log/0901_stable_version/error_log/dt600_1_csswm_1_vvm_2E5diff_7vvm_3B_4non_10kmcouple.e
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
