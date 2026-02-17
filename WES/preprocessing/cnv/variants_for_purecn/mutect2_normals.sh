#!/bin/bash
#SBATCH --job-name=mutect2
#SBATCH -n 1
#SBATCH --mem=10G       # total memory need
#SBATCH --time=72:00:00
##SBATCH --mail-type=END,FAIL # email notification when job ends/fails
##SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

# ARG PARSING
if [ "$#" -ne  2 ]; then
    echo "Usage: $0 <normal_name> <normal_bam>"
    exit 1
fi
normal_name="$1"
normal_bam="$2"


echo normal_name: $normal_name
echo normal_bam: $normal_bam


# SETUP
cd /liulab/jacobg
module load singularity/3.4.2
outdir=/liulab/jacobg/Ellisen/variant_calling_for_purecn/mutect2/normals/${normal_name}
mkdir $outdir

# CALL VARAINTS
singularity exec --bind /liulab/jacobg/Ellisen/iwhale/annotations/:/annotations /liulab/jacobg/Ellisen/somatic/gatk_4.0.5.2.sif \
java -Xmx20G -jar /gatk/gatk.jar Mutect2 \
--reference /annotations/GRCh38.fa \
--intervals /liulab/jacobg/Ellisen/Twist_Exome_Targets_hg38_modified_no_chr.bed \
--input ${normal_bam} \
--tumor-sample ${normal_name} \
--output $outdir/${normal_name}_unfiltered.MuTect2.vcf \
--genotype-germline-sites true \
--genotype-pon-sites true \
--interval-padding 50 

# FILTER VARIANTS
singularity exec --bind /liulab/jacobg/Ellisen/iwhale/annotations/:/annotations /liulab/jacobg/Ellisen/somatic/gatk_4.0.5.2.sif \
java -Xmx20G -jar /gatk/gatk.jar FilterMutectCalls \
--variant $outdir/${normal_name}_unfiltered.MuTect2.vcf \
--reference /annotations/GRCh38.fa \
--output $outdir/${normal_name}_filtered.MuTect2.vcf

echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2
