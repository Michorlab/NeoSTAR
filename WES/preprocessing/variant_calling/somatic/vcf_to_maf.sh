#!/bin/bash
#SBATCH --job-name=vcf2maf
#SBATCH -n 1
#SBATCH --mem=10G       # total memory need
#SBATCH --time=72:00:00
##SBATCH --mail-type=END,FAIL # email notification when job ends/fails
##SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

# ARG PARSING
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <tumor_name> <normal_name>"
    exit 1
fi
tumor_name="$1"
normal_name="$2"

source /michorlab/jacobg/miniforge3/etc/profile.d/conda.sh
conda activate bioinformatics
cd /michorlab/jacobg/Ellisen/somatic/mutect2/${tumor_name}

/michorlab/jacobg/miniforge3/envs/bioinformatics/bin/vcf2maf.pl --input-vcf ${tumor_name}_filtered.MuTect2.vcf --output-maf ${tumor_name}_filtered.MuTect2.maf --custom-enst /michorlab/jacobg/Ellisen/ref_files/hg38/vep/myc_isoform_overrides_uniprot  --ref-fasta /michorlab/jacobg/Ellisen/iwhale/annotations/GRCh38.fa --tumor-id ${tumor_name} --normal-id ${normal_name}  --ncbi-build GRCh38  --inhibit-vep
