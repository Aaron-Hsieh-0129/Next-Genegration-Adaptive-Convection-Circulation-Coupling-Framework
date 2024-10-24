#!/bin/bash

# Loop through each case in the RUN directory
for case in /home/Aaron/TMIF_VVM_CSSWM/RUN/1021_2speed_4interval/200_*; do
    if [ -d "$case" ]; then  # Check if the path is a directory
        echo "Running $case"
        cd "$case" || { echo "Failed to enter $case"; continue; }  # Change to the case directory

        # Submit the job and capture the job ID (optional)
        job_id=$(sbatch TMIF.sh --parsable)
        
        if [ $? -eq 0 ]; then
            echo "Submitted job with ID: $job_id"
        else
            echo "Failed to submit job in $case"
        fi

        # Wait for 30 seconds before submitting the next job
        echo "Waiting for 5 seconds before the next submission..."
        sleep 3
    fi
