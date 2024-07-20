#!/bin/bash
#SBATCH --nodelist=mogamd,node01
#SBATCH --ntasks-per-node=40
#SBATCH -o log/newTur2_dt_4.7_cloud_2_csswm_1E6diff_60p1.o
#SBATCH -e log/newTur2_dt_4.7_cloud_2_csswm_1E6diff_60p1.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b08209006@ntu.edu.tw

# Get the node names and save them to hostfile_mpi
srun hostname -s > hostfile_mpi

# Clean build directory, create a new one, configure and build the project
rm -rf build
mkdir build
cd build/
cmake ../
make -j 4

# Run the application using the custom hostfile
mpirun --hostfile ../hostfile_mpi -np 81 ./TMIF