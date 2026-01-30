#!/bin/bash
# Submit Purecn runs

SBATCH_SCRIPT=annotate_and_maf.sh
file=all_tumor_vcf_paths.txt
mkdir annotate_and_maf_slurms

# Loop through each line in sample_list.txt
while IFS= read -r line; do
    # Split the line into variables
    tumor_name=$(echo "$line" | awk '{print $1}')
    vcf_path=$(echo "$line" | awk '{print $2}')

    # Generate unique job name based on sample names
    job_name="${tumor_name}_vcf2maf"
    output="annotate_and_maf_slurms/${tumor_name}"

    # Submit the sbatch job 
    sbatch --job-name="$job_name" --output="$output" "$SBATCH_SCRIPT" "$tumor_name" "$vcf_path"

    # Print a message indicating job submission
    echo "Submitted job: $job_name"
done < $file

