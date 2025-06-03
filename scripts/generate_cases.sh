#!/usr/bin/bash

project_root="/home/Aaron/NextACC_twnia3"
expname="0516_p3_newprof_grav20"
csswm_gravity=0.0383
time_values=(600 1800 3600 5400 7200 10800 14400 18000 21600 25200 28800 32400 36000)
# time_values=(600)

# Initialize core counter starting from 5
core_counter=100

# Loop over seeds from 0 to 190, incrementing by 10
# for seed in $(seq 200 10 280); do
for seed in 290; do
    for timestep in "${time_values[@]}"; do
        # Define unique case name
        case_name="${timestep}_seed${seed}"
        case_dir="${project_root}/RUN/${expname}/${case_name}"
        rm -rf $case_dir
        mkdir -p $case_dir

        rsync -aq --exclude='build' --exclude='docs' --exclude='log' --exclude='DATA' --exclude='*.git' --exclude='*.md' --exclude='*.rst' $project_root/2DVVM $case_dir
        rsync -aq --exclude='build' --exclude='docs' --exclude='log' --exclude='DATA' --exclude='*.git' ${project_root}/CSSWM ${case_dir}
        cp -r ${project_root}/src ${case_dir}/src
        cp ${project_root}/CMakeLists.txt ${case_dir}/CMakeLists.txt

        cat <<EOL > ${case_dir}/config.txt
OUTPUTPATH=/data/Aaron/NextACC_twnia3/${expname}/${case_name}/
SEED=${seed}                                          # random seed for initial perturbation in the CRM bottom
COUPLETIME=${timestep}                                # Coupling time for NextGCC [s]
Bubble_p_i_j=[(1,47,47)]          # CRMs with bubble inside
NotBubble_p_i_j=[]                 # CRMs with nothing inside
BubbleCase=1                                          # Case0: Nothing, Case1: Bubble, Case2: Bubble+wind shear
CSSWM_GRAVITY=${csswm_gravity}                        # gravity wave speed for CSSWM

CSSWM_DT=200
CSSWM_TIMEEND=86400                                   # Integration Time [s]
CSSWM_OUTPUTSTEP=1                                    # Output frequency
CSSWM_DIFFUSION_KX=200000
CSSWM_DIFFUSION_KY=200000
CSSWM_DIFFUSION_TS=0.06
CSSWM_ADDFORCING_TIME=0                               # If the user specifies adding forcing, the adding time can be specified here
CSSWM_H_NUDGE_TIME=0                                  # CSSWM h nudging time scale, if it is 0, the nudge will be closed.

VVM_XRANGE=100000                                     # Domain for x [m]
VVM_ZRANGE=20000                                      # Domain for z [m]
VVM_DT=2
VVM_DX=200
VVM_DZ=200                                            # Should be same as dx
VVM_TIMEEND=86400                                     # Integration Time [s]                                
VVM_OUTPUTSTEP=50                                     # Output frequency                 
VVM_MOISTURE_NUDGE_TIME=0                             # VVM moisture nudging time scale, if it is 0, the nudge will be closed.
VVM_YEAR=2025                                       # Simulation start year
VVM_MONTH=3                                         # Simulation start year
VVM_DAY=21                                          # Simulation start year
VVM_HOUR=12                                         # Simulation start year
VVM_MINUTE=0                                        # Simulation start year
VVM_SECOND=0                                        # Simulation start year
VVM_LON=0                                           # VVM Longitude
VVM_LAT=0                                           # VVM Latitude
EOL

        # Write the job script (TMIF.sh) for the case
        cat <<EOL > ${case_dir}/TMIF.sh
#!/usr/bin/bash
#SBATCH -J ${case_name}
#SBATCH -N 1
#SBATCH -c 1
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
mpirun -np 1 --cpu-set ${core_counter} --bind-to core ./TMIF

EOL

        # Make the script executable
        chmod +x ${case_dir}/TMIF.sh

        # Increment core counter for the next case
        core_counter=$((core_counter + 1))

    done
done
