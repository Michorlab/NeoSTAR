#!/bin/bash
#SBATCH --job-name=mutect2
#SBATCH -n 1
#SBATCH --mem=15G       # total memory need
#SBATCH --time=72:00:00
##SBATCH --mail-type=END,FAIL # email notification when job ends/fails
##SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

# ARG PARSING
if [ "$#" -ne  4 ]; then
    echo "Usage: $0 <tumor_name> <normal_name> <tumor_bam> <normal_bam>"
    exit 1
fi
tumor_name="$1"
normal_name="$2"
tumor_bam="$3"
normal_bam="$4"

echo tumor_name: $1
echo normal_name: $2
echo tumor_bam: $3
echo normal_bam: $4

# SETUP
cd /liulab/jacobg
module load singularity/3.4.2
mkdir /liulab/jacobg/Ellisen/variant_calling_for_purecn/mutect2/tumors/${tumor_name}

# CALL VARAINTS
singularity exec --bind /liulab/jacobg/Ellisen/iwhale/annotations/:/annotations /liulab/jacobg/Ellisen/somatic/gatk_4.0.5.2.sif \
java -Xmx20G -jar /gatk/gatk.jar Mutect2 \
--reference /annotations/GRCh38.fa \
--intervals /liulab/jacobg/Ellisen/Twist_Exome_Targets_hg38_modified_no_chr.bed \
--input ${tumor_bam} \
--input ${normal_bam} \
--normal-sample ${normal_name} \
--tumor-sample ${tumor_name} \
--output /liulab/jacobg/Ellisen/variant_calling_for_purecn/mutect2/tumors/${tumor_name}/${tumor_name}_unfiltered.MuTect2.vcf \
--genotype-germline-sites true \
--genotype-pon-sites true \
--interval-padding 50 \
--germline-resource /liulab/jacobg/Ellisen/variant_calling_for_purecn/af-only-gnomad.hg38.nochr.vcf.gz
#--germline-resource /michorlab/jacobg/Michalina/af-only-gnomad.hg38.vcf.gz 
#--panel-of-normals /michorlab/jacobg/Michalina/1000g_pon.hg38.vcf.gz 

# FILTER VARIANTS
singularity exec --bind /liulab/jacobg/Ellisen/iwhale/annotations/:/annotations /liulab/jacobg/Ellisen/somatic/gatk_4.0.5.2.sif \
java -Xmx20G -jar /gatk/gatk.jar FilterMutectCalls \
--variant /liulab/jacobg/Ellisen/variant_calling_for_purecn/mutect2/tumors/${tumor_name}/${tumor_name}_unfiltered.MuTect2.vcf \
--reference /annotations/GRCh38.fa \
--output /liulab/jacobg/Ellisen/variant_calling_for_purecn/mutect2/tumors/${tumor_name}/${tumor_name}_filtered.MuTect2.vcf


echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2
