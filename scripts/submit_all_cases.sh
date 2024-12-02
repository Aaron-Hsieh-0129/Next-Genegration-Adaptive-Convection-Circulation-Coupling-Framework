#!/usr/bin/bash

time_values=(5400)
seed_values=$(seq 0 10 30)
expname="1202_2speed_eva"

# Iterate through the seeds first
for seed in $seed_values; do
    echo "Processing jobs for seed: $seed"
    
    # Then iterate through each time for the current seed
    for time in "${time_values[@]}"; do
        # Construct the expected folder pattern using the time and seed
        for case in /home/Aaron/TMIF_VVM_CSSWM/RUN/${expname}/200_"$time"_*"_seed$seed"_*; do
            if [ -d "$case" ]; then  # Check if the path is a directory
                echo "Running case with time: $time and seed: $seed in folder: $case"
                cd "$case" || { echo "Failed to enter $case"; continue; }  # Change to the case directory

                # Submit the job and capture the job ID
                job_id=$(sbatch TMIF.sh --parsable)
                
                if [ $? -eq 0 ]; then
                    echo "Submitted job with ID: $job_id"
                else
                    echo "Failed to submit job in $case"
                fi

                # Wait for 5 seconds before submitting the next job
                echo "Waiting for 5 seconds before the next submission..."
                # sleep 5
            else
                echo "No directory found for time: $time and seed: $seed"
            fi
        done
    done
done
