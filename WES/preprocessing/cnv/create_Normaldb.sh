#!/bin/bash
#SBATCH --job-name=normaldb
#SBATCH -n 1
#SBATCH --mem=20G       # total memory need
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

module load singularity
cd /liulab/ 

singularity exec /liulab/jacobg/Ellisen/purecn/purecn_latest.sif \
Rscript /opt/PureCN/NormalDB.R --out-dir /liulab/jacobg/Ellisen/purecn/NormalDB --coverage-files /liulab/jacobg/Ellisen/purecn/normal_coverages.list --normal-panel /liulab/jacobg/Ellisen/purecn/normals_combined.vcf.gz --genome hg38 --assay Broad_WES_v6 # purely aesthetic 

