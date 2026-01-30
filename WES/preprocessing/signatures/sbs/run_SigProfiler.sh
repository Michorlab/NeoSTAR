#!/bin/bash
#SBATCH --job-name=signatures
#SBATCH -n 1
#SBATCH --mem=40G       # total memory need
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify


cd /michorlab/jacobg/Ellisen/signatures/SigProfiler
eval "$(conda shell.bash hook)"
conda activate SigProfiler

python run_SigProfiler.py
