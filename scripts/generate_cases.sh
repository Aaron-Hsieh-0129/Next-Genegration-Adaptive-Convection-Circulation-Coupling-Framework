#!/usr/bin/bash

project_root="/home/Aaron/TMIF_VVM_CSSWM"
expname="1202_2speed_eva"
csswm_gravity=0.4
time_values=(5400)

# Loop over seeds from 0 to 190, incrementing by 10
for seed in $(seq 0 10 30); do
    for timestep in "${time_values[@]}"; do
        # Define unique case name
        case_name="200_${timestep}_7vvm_3B_random1s_seed${seed}_4non"
        case_dir="${project_root}/RUN/${expname}/${case_name}"
        rm -rf $case_dir
        mkdir -p $case_dir

        rsync -aq --exclude='build' --exclude='docs' --exclude='log' --exclude='*.md' --exclude='*.rst' $project_root/2DVVM $case_dir
        rsync -aq --exclude='build' --exclude='docs' --exclude='log' ${project_root}/CSSWM ${case_dir}
        cp -r ${project_root}/src ${case_dir}/src
        cp ${project_root}/CMakeLists.txt ${case_dir}/CMakeLists.txt

        cat <<EOL > ${case_dir}/config.txt
OUTPUTPATH=/data/Aaron/TMIF/${expname}/${case_name}/
SEED=${seed}
COUPLETIME=${timestep}
Bubble_p_i_j=[(1,46,47),(1,47,47),(1,48,47)]
NotBubble_p_i_j=[(1,45,47),(1,49,47)]
BubbleCase=1
CSSWM_GRAVITY=${csswm_gravity}
EOL

        # Write the job script (TMIF.sh) for the case
        cat <<EOL > ${case_dir}/TMIF.sh
#!/usr/bin/bash
#SBATCH -J ${case_name}
#SBATCH -N 1
#SBATCH -c 5
#SBATCH --nodelist=mogamd
#SBATCH -o ${case_name}_%j.o
#SBATCH -e ${case_name}_%j.e
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
cd build/ && cmake $case_dir && make -j 4

# Run the atmospheric model with dynamic parameters
./TMIF

EOL

        # Make the script executable
        chmod +x ${case_dir}/TMIF.sh

    done
done
