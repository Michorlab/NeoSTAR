#!/bin/bash
#SBATCH --job-name=vcf2maf
#SBATCH -n 1
#SBATCH --mem=10G       # total memory need
#SBATCH --time=72:00:00
##SBATCH --mail-type=END,FAIL # email notification when job ends/fails
##SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

# ARG PARSING
if [ "$#" -ne  2 ]; then
    echo "Usage: <tumor_name> <vcf_path>"
    exit 1
fi
tumor_name="$1"
vcf_path="$2"

echo tumor_name: $tumor_name
echo vcf_path: $vcf_path


# SETUP
cd /liulab/jacobg
module load singularity/3.4.2
outdir=/liulab/jacobg/Ellisen/variant_calling_for_purecn/mafs/${tumor_name}
mkdir $outdir


#VCF2MAF
vep_fp=/opt/conda/bin/vep
vep_path=$(dirname "$vep_fp")
center=Novogene
build=GRCh38


singularity exec --bind /liulab/jacobg/Ellisen/iwhale/annotations/:/annotations /liulab/jacobg/Ellisen/variant_calling_for_purecn/vep.simg \
vcf2maf.pl --input-vcf $vcf_path --output-maf $outdir/$tumor_name.maf \
           --tumor-id $tumor_name \
           --ref-fasta /annotations/GRCh38.fa \
           --vep-data /liulab/jacobg/Ellisen/variant_calling_for_purecn/ \
           --ncbi-build $build \
           --vep-path $vep_path \
           --maf-center $center
           # normal id and alternate isoforms would go here
