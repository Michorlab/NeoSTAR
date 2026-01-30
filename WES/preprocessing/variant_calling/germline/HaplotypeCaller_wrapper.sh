#!/bin/bash
# Path to sbatch script
SBATCH_SCRIPT=/michorlab/jacobg/Ellisen/germline/HaplotypeCaller.sh
file="normals.txt"

# Loop through each line
while IFS= read -r line; do
    # Split the line
    normal_name=$(echo "$line" | awk '{print $1}')
    normal_bam=$(echo "$line" | awk '{print $2}')

    # Generate unique job name based on sample name
    job_name="HaplotypeCaller_${normal_name}"
    output="HaplotypeCaller_slurms/${normal_name}"

    # Submit the sbatch job with normal sample name and bam as arguments
    sbatch --job-name="$job_name" --output="$output" "$SBATCH_SCRIPT" "$normal_name" "$normal_bam"

    # Print a message indicating job submission
    echo "Submitted job: $job_name"
done < $file

