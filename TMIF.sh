#!/bin/bash
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --nodelist=node01
#SBATCH -o log/0901_stable_version/600_600_12kmcouple_7vvm_3B_4non.o
#SBATCH -e log/0901_stable_version/error_log/600_600_12kmcouple_7vvm_3B_4non.e
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
