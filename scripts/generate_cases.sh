#!/usr/bin/bash

project_root="/home/Aaron/TMIF_VVM_CSSWM"

# Define model parameters
timesteps=(600 3600 5400 7200 9000 10800 12600 14400 16200 18000)   # Different time steps
# timesteps=(9000)   # Different time steps

# Loop over seeds from 0 to 190, incrementing by 10
for seed in $(seq 90 10 90); do
    for timestep in "${timesteps[@]}"; do
        # Define unique case name
        case_name="200_${timestep}_7vvm_3B_random1s_seed${seed}_4non"
        case_dir="${project_root}/RUN/1024_2speed/$case_name"
        mkdir -p $case_dir

        cat <<EOL > $case_dir/config.txt
OUTPUTPATH=/data/Aaron/TMIF/1024_2speed/$case_name/
SEED=$seed
COUPLETIME=$timestep
Bubble_p_i_j=[(1,46,47), (1,47,47), (1,48,47)]
NotBubble_p_i_j=[(1,45,47), (1,49,47)]
EOL

        # Write the job script (TMIF.sh) for the case
        cat <<EOL > $case_dir/TMIF.sh
#!/usr/bin/bash
#SBATCH -J $case_name
#SBATCH -N 1
#SBATCH -c 5
#SBATCH --nodelist=mogamd
#SBATCH -o $case_name_%j.o
#SBATCH -e $case_name_%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aaronfinland0129@gmail.com

# Handle OpenMP threads
if [ -n "\$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=\$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi
export OMP_NUM_THREADS=\$omp_threads
echo \$OMP_NUM_THREADS

# Clean and rebuild the model
rm -rf build
mkdir build
cd build/ && cmake $project_root && make -j 4

# Run the atmospheric model with dynamic parameters
./TMIF

EOL

        # Make the script executable
        chmod +x $case_dir/TMIF.sh

    done
done
