import os
import pandas as pd
import numpy as np

folder = '/liulab/jacobg/Ellisen/purecn/final_loh_files_with_curation'
out_path = '/liulab/jacobg/Ellisen/gistic/input_tsv/purecn/usable_samples.tsv' 
cnv_files = [x for x in sorted(os.listdir(folder)) if x.endswith('loh.csv')]

dfs = []
for file in cnv_files:
    d = pd.read_csv(os.path.join(folder, file))
    d['num.markers'] = 1 # we think this doesn't mattter and is ignored
    d['gistic_cn'] = d['C'].apply(lambda x: (np.log2(x) -1) if x > 0 else -5)
    d = d.loc[:,['Sampleid', 'chr', 'start', 'end', 'num.markers', 'gistic_cn']]
    dfs.append(d)

df = pd.concat(dfs)
bad_samples = ['BIDMC2_pre', 'BIDMC3_post', 'BIDMC4_pre', 'BIDMC4_post', 'BIDMC5_post', 'BIDMC6_pre', 
               'BIDMC6_post', 'BIDMC7_post', 'DFC7_pre_FFT', 'DFC9_pre', 'DFC12_post', 'DFC19_pre', 
               'DFC19_pre', 'DFC20_post', 'DFC23_post', 'MGH3_post', 'MGH5_post', 'MGNS4_pre', 
               'NWH2_pre_FFT', 'MGNS5_post', 'DFC1_pre', 'DFC19_post', 'MGNS4_post', 'BIDMC3_post', 
               'DFC1_post']
df = df[df['Sampleid'].isin(bad_samples) == False].reset_index(drop=True)

# handle 2 cases where purecn segmentation has overlapping segments (persists from dnacopy to loh.csv)
for i in df.index[1:]:
    if df.iloc[i]['start'] < df.iloc[i - 1]['end']:
        if df.iloc[i]['chr'] == df.iloc[i - 1]['chr']:
            if df.iloc[i]['gistic_cn'] == df.iloc[i - 1]['gistic_cn']:
                print("overlap!")
                print(df.iloc[i-1:i+1])
                df.at[i,'start'] = df.at[i - 1,'end']
            else: 
                # remove areas of overlap is probably smarter than the mean if the areas are small
                print(df.iloc[i-1:i+1])
                df.at[i-1,'end'], df.at[i,'start']  = df.at[i,'start'], df.at[i-1,'end']
                print(df.iloc[i-1:i+1])                   


print(df)
df.to_csv(out_path, sep='\t', index=False, header=False)
