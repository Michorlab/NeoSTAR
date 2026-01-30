import os
import pandas as pd

folder = '/michorlab/jacobg/Ellisen/purecn/final_loh_files_with_curation'
out_path = '/michorlab/jacobg/Ellisen/signatures/SigProfiler/all_cnvs.tsv' 
cnv_files = [x for x in sorted(os.listdir(folder)) if x.endswith('loh.csv')]

dfs = []
for file in cnv_files:
    d = pd.read_csv(os.path.join(folder, file))
    d['Normal TCN'] = 2
    d['Normal BCN'] = 1
    d = d.loc[:,['Sampleid', 'chr', 'start', 'end','Normal TCN', 'Normal BCN', 'C', 'M']]
    d.columns = ['sample', 'Chromosome', 'Start Position', 'End Position', 'Normal TCN', 'Normal BCN', 'Tumour TCN', 'Tumour BCN']
    dfs.append(d)

df = pd.concat(dfs)
df = df.dropna(axis=0)
df.to_csv(out_path, sep='\t', index=False)
