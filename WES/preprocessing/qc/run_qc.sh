#!/bin/bash
#SBATCH --job-name=qc
#SBATCH -n 1
#SBATCH --mem=1G       # total memory need
#SBATCH --time=168:00:00
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify

# Script parameters
patient=$1
bam_file=$2
regions=$3

# Ensure necessary directories exist
mkdir -p "qc/${patient}"
mkdir -p "qc/fastqc"

# qc environment
export CONDA_PREFIX=/michorlab/jacobg/miniconda3
export CONDA_ROOT=/michorlab/jacobg/miniconda3
export PATH=/michorlab/jacobg/miniconda3/bin:$PATH
source activate qc

# Task: Coverage
echo "Running coverage task..."
#mosdepth --by "default_value" -t 4 "qc/${patient}" "$bam_file"  # Replace 'default_value' with your actual value
mosdepth --by "$regions" -t 4 "qc/${patient}" "$bam_file"

# Task: FastQC
echo "Running FastQC task..."
tmpdir="qc/fastqc/${patient}.tmp"
mkdir -p "$tmpdir"
fastqc --outdir "$tmpdir" "$bam_file"
mv "$tmpdir/07_recalibrated_fastqc.html" "qc/fastqc/${patient}_fastqc.html"
mv "$tmpdir/07_recalibrated_fastqc.zip" "qc/fastqc/${patient}_fastqc.zip"
rm -r "$tmpdir"

# Task: Stats
echo "Running stats task..."
samtools flagstat "$bam_file" > "qc/${patient}.flagstat"
echo "All tasks completed."
