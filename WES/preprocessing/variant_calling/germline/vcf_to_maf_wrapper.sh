#!/bin/bash
# Submit Mutect2 jobs for multiple tumor-normal pairs

# Path to your sbatch script
SBATCH_SCRIPT=/michorlab/jacobg/Ellisen/germline/vcf_to_maf.sh
file="normals.txt"

# Loop through each line in sample_list.txt
while IFS= read -r line; do
    # Split the line into tumor and normal sample names
    normal_name=$(echo "$line" | awk '{print $1}')

    # Generate unique job name based on sample names
    job_name="vcf2maf_${normal_name}"
    output="vcf2maf_slurms/${normal_name}"

    # Submit the sbatch job
    sbatch --job-name="$job_name" --output="$output" "$SBATCH_SCRIPT" "$normal_name"
    echo "Submitted job: $job_name"

done < $file

