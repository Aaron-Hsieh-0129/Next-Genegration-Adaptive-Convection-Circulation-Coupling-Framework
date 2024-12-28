#!/bin/bash
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --nodelist=mogamd
#SBATCH -o log/EquatorialWave.o
#SBATCH -e log/EquatorialWave.e
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
