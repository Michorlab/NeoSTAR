#!/bin/bash
# Submit Purecn runs


SBATCH_SCRIPT="run_purecn.sh"
file=tumors.txt 
#file=tumors_test.txt
mkdir purecn_slurms

# Loop through each line
while IFS= read -r line; do
    # Split the line
    tumor_name=$(echo "$line" | awk '{print $1}')
    
    # Generate unique job name based on sample names
    job_name="${tumor_name}_purecn"
    output="purecn_slurms/${tumor_name}"

    # Submit the sbatch job
    sbatch --job-name="$job_name" --output="$output" "$SBATCH_SCRIPT" "$tumor_name"

    # Print a message indicating job submission
    echo "Submitted job: $job_name"
done < $file

