#!/bin/bash


SBATCH_SCRIPT=/michorlab/jacobg/Ellisen/run_qc.sh
TABLE='/michorlab/jacobg/Ellisen/qc_prep.tsv'

#OVERWRITES EXISTING 
while IFS= read -r line; do 
    sample=$(echo "$line" | awk '{print $1}')
    bam_file=$(echo "$line" | awk '{print $2}')
    regions=$(echo "$line" | awk '{print $3}')
    job_name="qc_${sample}"
    output="/michorlab/jacobg/Ellisen/qc_slurm/slurm_qc_${sample}.out" 
    sbatch --job-name="$job_name" --output="$output" "$SBATCH_SCRIPT" "$sample" "$bam_file" "$regions" 

    # Print a message indicating job submission
    echo "Submitted job: $job_name"
done < $TABLE
