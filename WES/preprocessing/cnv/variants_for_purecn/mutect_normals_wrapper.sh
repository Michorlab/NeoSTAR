#!/bin/bash

# Path to sbatch script
SBATCH_SCRIPT="mutect2_normals.sh"
file=normal_sample_paths.txt
mkdir mutect2_normal_slurms

# Loop through each line
while IFS= read -r line; do
    normal_name=$(echo "$line" | awk '{print $1}') 
    normal_bam=$(echo "$line" | awk '{print $2}')

    # Generate unique job name based on sample names
    job_name="mutect2_${normal_name}"
    output="mutect2_normal_slurms/${normal_name}"

    # Submit the sbatch job
    sbatch --job-name="$job_name" --output="$output" "$SBATCH_SCRIPT" "$normal_name" "$normal_bam"

    # Print a message indicating job submission
    echo "Submitted job: $job_name"
done < $file

