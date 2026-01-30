#!/bin/bash
#SBATCH --job-name=GISTIC2
#SBATCH -n 1
#SBATCH --mem=32G       # total memory need
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

#SET PARAMS HERE!!!!!
segfile=/liulab/jacobg/Ellisen/gistic/input_tsv/purecn/usable_samples.tsv
OUTDIR=/liulab/jacobg/Ellisen/gistic/purecn/usable_samples

module load singularity/3.4.2


# RUN CONATINER (https://github.com/ShixiangWang/install_GISTIC/blob/master/run_singularity.sh)
GISTIC_LOC=/opt/GISTIC          # DONT change this position, it refers to docker path


echo --- creating output directory ---
basedir=$OUTDIR
mkdir -p $basedir
echo output directory is $basedir

echo --- running GISTIC ---
echo Set input file paths...

## reference file
refgenefile=$GISTIC_LOC/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat


echo Starting...
DOCKER_OUTDIR="$GISTIC_LOC"/run_result

singularity run --bind $basedir:$DOCKER_OUTDIR,$segfile docker://shixiangwang/gistic \
  -b $DOCKER_OUTDIR -seg $segfile -refgene $refgenefile \
  -rx 1 -genegistic 1 -smallmem 1 -broad 1 -twosize 1 \
  -armpeel 1 -savegene 1 -qvt 0.01 


