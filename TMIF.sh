#!/bin/bash
#SBATCH -N 1
#SBATCH -n 7
#SBATCH --nodelist=mogamd
#SBATCH -o log/0912_ensemble_smaller_perturb/200_3600_7vvm_3B_random1s_seed60_4non.o
#SBATCH -e log/0912_ensemble_smaller_perturb/error_log/200_3600_7vvm_3B_random1s_seed60_4non.e
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
