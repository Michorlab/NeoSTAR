#!/bin/bash

# Path to script
SBATCH_SCRIPT=/michorlab/jacobg/Ellisen/somatic/mutect2.sh
file="sample_pairs.txt"

# Loop through each line
while IFS= read -r line; do
    # Split the line into tumor and normal sample names
    tumor_name=$(echo "$line" | awk '{print $1}')
    normal_name=$(echo "$line" | awk '{print $2}')
    tumor_bam=$(echo "$line" | awk '{print $3}')
    normal_bam=$(echo "$line" | awk '{print $4}')

    # Generate unique job name based on sample names
    job_name="mutect2_${tumor_name}_${normal_name}"
    output="mutect2_slurms/${tumor_name}"

    # Submit the sbatch job with tumor and normal sample names as arguments
    sbatch --job-name="$job_name" --output="$output" "$SBATCH_SCRIPT" "$tumor_name" "$normal_name" "$tumor_bam" "$normal_bam"

    # Print a message indicating job submission
    echo "Submitted job: $job_name"
done < $file

