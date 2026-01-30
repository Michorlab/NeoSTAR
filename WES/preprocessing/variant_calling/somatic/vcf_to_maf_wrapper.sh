#!/bin/bash
# Submit Mutect2 jobs for multiple tumor-normal pairs

# Path to your sbatch script
SBATCH_SCRIPT=/michorlab/jacobg/Ellisen/somatic/vcf_to_maf.sh

# Loop through each line
while IFS= read -r line; do
    # Split the line into tumor and normal sample names
    tumor_name=$(echo "$line" | awk '{print $1}')
    normal_name=$(echo "$line" | awk '{print $2}')

    # Generate unique job name based on sample names
    job_name="vcf2maf_${tumor_name}_${normal_name}"
    output="vcf_to_maf_slurms/${tumor_name}"

    # Submit the sbatch job with tumor and normal sample names as arguments
    sbatch --job-name="$job_name" --output="$output" "$SBATCH_SCRIPT" "$tumor_name" "$normal_name"

    # Print a message indicating job submission
    echo "Submitted job: $job_name"
done < "sample_pairs.txt"

