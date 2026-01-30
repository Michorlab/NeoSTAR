#!/bin/bash
#SBATCH --job-name=vcf2maf
#SBATCH -n 1
#SBATCH --mem=10G       # total memory need
#SBATCH --time=72:00:00
##SBATCH --mail-type=END,FAIL # email notification when job ends/fails
##SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

# ARG PARSING
if [ "$#" -ne 1 ]; then
    echo "Usage: $0  <normal_name>"
    exit 1
fi
normal_name="$1"


# SETUP
eval "$(conda shell.bash hook)"
conda activate bioinformatics

cd /michorlab/jacobg/Ellisen/germline/results/${normal_name}




/michorlab/jacobg/miniforge3/envs/bioinformatics/bin/vcf2maf.pl --input-vcf ${normal_name}_HaplotypeCaller.vcf  --output-maf ${normal_name}_HaplotypeCaller.maf  --custom-enst /michorlab/jacobg/Ellisen/ref_files/hg38/vep/myc_isoform_overrides_uniprot  --ref-fasta /michorlab/jacobg/Ellisen/iwhale/annotations/GRCh38.fa --tumor-id ${normal_name} --ncbi-build GRCh38  --inhibit-vep
