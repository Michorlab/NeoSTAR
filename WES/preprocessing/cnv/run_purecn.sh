#!/bin/bash
#SBATCH --job-name=purecn
#SBATCH -n 16
#SBATCH --mem=50G       # total memory need
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
Rscript /opt/PureCN/PureCN.R --out /liulab/jacobg/Ellisen/purecn/results/$tumor_name \
    --tumor /liulab/jacobg/Ellisen/purecn/coverage/tumors/${tumor_name}_recalibrated_coverage_loess.txt.gz \
    --sampleid $tumor_name \
    --vcf /liulab/jacobg/Ellisen/variant_calling_for_purecn/mutect2/tumors/${tumor_name}/${tumor_name}_unfiltered_modified.MuTect2.vcf \
    --fun-segmentation PSCBS \
    --normaldb /liulab/jacobg/Ellisen/purecn/NormalDB/normalDB_Broad_WES_v6_hg38.rds \
    --mapping-bias-file /liulab/jacobg/Ellisen/purecn/NormalDB/mapping_bias_Broad_WES_v6_hg38.rds \
    --intervals /liulab/jacobg/Ellisen/purecn/Twist_Exome_Targets_hg38_modified_no_chr_intervals.txt \
    --genome hg38 \
    --model betabin \
    --force \
    --post-optimize \
    --seed 507346 \
    --popaf-info-field POP_AF \
    --max-non-clonal 0.3 \
    --max-copy-number 7 \
    --alpha 0.0000000001 \
    --undo-sd 5 \
    --cores 16 

