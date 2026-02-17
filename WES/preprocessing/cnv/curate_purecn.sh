#!/bin/bash
#SBATCH --job-name=purecn
#SBATCH -n 1
#SBATCH --mem=20G       # total memory need
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

# ARG PARSING
if [ "$#" -ne  1 ]; then
    echo "Usage: $0 <tumor_name>"
    exit 1
fi
tumor_name="$1"

echo tumor_name: $tumor_name


module load singularity
cd /liulab/jacobg


singularity exec /liulab/jacobg/Ellisen/purecn/purecn_latest.sif \
Rscript /opt/PureCN/PureCN.R --out /liulab/jacobg/Ellisen/purecn/curation/results/$tumor_name \
  --rds /liulab/jacobg/Ellisen/purecn/curation/${tumor_name}.rds
