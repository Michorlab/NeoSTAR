#!/bin/bash
#SBATCH --job-name=haplotyper
#SBATCH -n 1
#SBATCH --mem=15G # 10G       # total memory need
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

# ARG PARSING
if [ "$#" -ne  2 ]; then
    echo "Usage: $0 <normal_name> <normal_bam>"
    exit 1
fi
normal_name="$1"
normal_bam="$2"

echo normal_name: $1
echo normal_bam: $2

# SETUP
cd /michorlab/jacobg
module load singularity/3.4.2
mkdir /michorlab/jacobg/Ellisen/germline/results/${normal_name}

# CALL VARAINTS
singularity exec --bind /michorlab/jacobg/Ellisen/iwhale/annotations/:/annotations /michorlab/jacobg/Ellisen/somatic/gatk_4.0.5.2.sif \
java -Xmx20G -jar /gatk/gatk.jar HaplotypeCaller \
--reference /annotations/GRCh38.fa \
--intervals /michorlab/jacobg/Ellisen/Twist_Exome_Targets_hg38_modified_no_chr.bed \
--interval-padding 100 \
--input ${normal_bam} \
--output /michorlab/jacobg/Ellisen/germline/results/${normal_name}/${normal_name}_HaplotypeCaller.vcf


echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2
