import os
import pandas as pd
import argparse


#READ AND COMBINE MAFS
outdir = "/michorlab/jacobg/Ellisen/signatures/SigProfiler/input_mafs"
rootdir = '/michorlab/jacobg/Ellisen/somatic/mutect2/'

for root, subdirs, files in os.walk(rootdir):
    for file in files:
        if file.endswith(".maf"):
            path = os.path.join(root,file)
            print(path)
            d = pd.read_csv(path, header=1, sep='\t')
            d = d[d['FILTER'] == 'PASS']
            d = d[d["Chromosome"] != 'MT'].dropna(axis=0, how='all')
            dest_path = os.path.join(outdir,file).replace("_filtered.MuTect2", '_somatic')
            d.to_csv(dest_path, index=False, sep='\t')

