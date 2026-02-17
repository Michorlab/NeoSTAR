#!/bin/bash
#SBATCH --job-name=intervals
#SBATCH -n 1
#SBATCH --mem=5G       # total memory need
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

# DOWNLOAD BIGWIG FILE
#wget https://s3.amazonaws.com/purecn/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw /michorlab/jacobg/Ellisen/purecn/wgEncodeCrgMapabilityAlign100mer.bigWig

cd /liulab/
module load singularity

singularity exec --bind /liulab/jacobg/Ellisen/iwhale/annotations/:/annotations /liulab/jacobg/Ellisen/purecn/purecn_latest.sif \
Rscript /opt/PureCN/IntervalFile.R --in-file /liulab/jacobg/Ellisen/Twist_Exome_Targets_hg38_modified_no_chr.bed --fasta /annotations/GRCh38.fa --out-file /liulab/jacobg/Ellisen/purecn/Twist_Exome_Targets_hg38_modified_no_chr_target_only_intervals.txt  --genome hg38 --export /liulab/jacobg/Ellisen/purecn/Twist_Exome_Targets_hg38_modified_no_chr_target_only_optimized.bed --mappability /liulab/jacobg/Ellisen/purecn/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw


#singularity exec --bind /liulab/jacobg/Ellisen/iwhale/annotations/:/annotations /liulab/jacobg/Ellisen/purecn/purecn_latest.sif \
#Rscript /opt/PureCN/IntervalFile.R --in-file /liulab/jacobg/Ellisen/Twist_Exome_Targets_hg38_modified_no_chr.bed --fasta /annotations/GRCh38.fa --out-file /liulab/jacobg/Ellisen/purecn/Twist_Exome_Targets_hg38_modified_no_chr_intervals.txt --off-target --genome hg38 --export /liulab/jacobg/Ellisen/purecn/Twist_Exome_Targets_hg38_modified_no_chr_optimized.bed --mappability /liulab/jacobg/Ellisen/purecn/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw
 




