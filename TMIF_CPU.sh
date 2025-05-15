#!/bin/bash
#SBATCH --account=MST113255
#SBATCH --partition=n3atm_280
#SBATCH --job-name=NextACC
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=log/0514/job-%j.out
#SBATCH --error=log/0514/job-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zz85721@gmail.com

#if [ -n "$SLURM_CPUS_PER_TASK" ]; then
#  omp_threads=$SLURM_CPUS_PER_TASK
#else
#  omp_threads=1
#fi
# export OMP_NUM_THREADS=$omp_threads
# echo $OMP_NUM_THREADS
export OMP_NUM_THREADS=1
export OMP_MAX_ACTIVE_LEVELS=0

# rm -rf build
# mkdir build

# wiht GPU
# cd build/ && cmake .. && make -j 8 && mpirun -np 2 -mca btl_base_warn_component_unused 0 -np 1 -x CUDA_VISIBLE_DEVICES=0 ./TMIF : -np 1 -x CUDA_VISIBLE_DEVICES=0 ./TMIF

# wiht CPU
cd build/ && cmake .. && make -j 8 && mpirun -np 2 ./TMIF