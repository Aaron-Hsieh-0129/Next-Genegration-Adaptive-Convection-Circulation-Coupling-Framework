#!/bin/bash
#SBATCH -N 1
#SBATCH -n 7
#SBATCH --nodelist=node01
#SBATCH -o log/0904_couple_time/200_600_7vvm_3B_4non.o
#SBATCH -e log/0904_couple_time/error_log/200_600_7vvm_3B_4non.e
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
