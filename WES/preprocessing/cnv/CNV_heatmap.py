"""Heatmap used in supp 12B"""
import collections
import logging

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from skgenome.rangelabel import unpack_range

from cnvlib import plots
import os
import pandas as pd

# LOAD FILES
segments_path = '/michorlab/jacobg/Ellisen/purecn/final_loh_files_with_curation'
dfs = []
for f in os.listdir(segments_path):
    df = pd.read_csv(f"{segments_path}/{f}")
    dfs.append(df)
cnvs = pd.concat(dfs)

print(cnvs)
cnvs['C'] = cnvs['C'].apply(lambda x: max(6,x))
cnvs = cnvs.rename(columns={'C': 'chromosome'})
print(cnvs)

def copy_number_to_color(cn):
    """Convert copy number to RGB color"""
    colors = {
        0: (0, 0, 0.5),      # dark blue
        1: (0.4, 0.4, 0.8),  # light blue
        2: (1, 1, 1),        # white
        3: (0.8, 0.4, 0.4),  # light red
        4: (0.7, 0, 0),      # medium red
        5: (0.5, 0, 0),      # dark red
        6: (0.3, 0, 0)       # darkest red
    }
    # Ensure copy number is in valid range
    cn = min(max(int(cn), 0), 6)
    return colors[cn]

def set_colorbar_discrete(axis):
    """Create a discrete colorbar for absolute copy numbers 0-6"""
    # Define colors for each copy number
    colors = {
        0: (0, 0, 0.5),      # dark blue
        1: (0.4, 0.4, 0.8),  # light blue
        2: (1, 1, 1),        # white
        3: (0.8, 0.4, 0.4),  # light red
        4: (0.7, 0, 0),      # medium red
        5: (0.5, 0, 0),      # dark red
        6: (0.3, 0, 0)       # darkest red
    }
    
    # Create a custom colormap with discrete steps
    cmap = mpl.colors.ListedColormap([colors[i] for i in range(7)])
    
    # Set boundaries for discrete values
    bounds = np.arange(-0.5, 6.51, 1.0)  # -0.5 to 6.5 to center the colors
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    # Create colorbar
    mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    cbar = plt.colorbar(mappable, ax=axis, orientation='vertical',
                       fraction=0.04, pad=0.03, shrink=0.6,
                       boundaries=bounds, ticks=range(7))
    cbar.set_label("Copy Number", labelpad=0)

def do_heatmap(df, show_range=None, by_bin=False, ax=None):
    """Plot copy number for multiple samples as a heatmap."""
    if ax is None:
        _fig, axis = plt.subplots()
    else:
        axis = ax
    set_colorbar_discrete(axis)

    # Get unique samples
    samples = df['Sampleid'].unique()
    
    # List sample names on the y-axis
    axis.set_yticks([i + 0.5 for i in range(len(samples))])
    axis.set_yticklabels(samples)
    axis.set_ylim(0, len(samples))
    axis.invert_yaxis()
    axis.set_ylabel("Samples")
    axis.set_facecolor('#DDDDDD')

    # Get chromosome info
    if show_range:
        r_chrom = show_range
        r_start = None
        r_end = None
        if ':' in show_range:
            r_chrom, positions = show_range.split(':')
            if positions:
                r_start, r_end = map(lambda x: int(x) if x else None, 
                                   positions.split('-'))
        # Filter data for range
        plot_df = df[df['chromosome'] == r_chrom].copy()
        if r_start is not None:
            plot_df = plot_df[plot_df['start'] >= r_start]
        if r_end is not None:
            plot_df = plot_df[plot_df['end'] <= r_end]
    else:
        plot_df = df.copy()
        r_chrom = None
        
    def plot_sample_chrom(i, sample_data):
        """Draw the given coordinates and colors as a horizontal series."""
        verts = []
        colors = []
        for _, row in sample_data.iterrows():
            verts.append([
                (row['start'], i),      # bottom left
                (row['start'], i + 1),  # top left
                (row['end'], i + 1),    # top right
                (row['end'], i)         # bottom right
            ])
            colors.append(copy_number_to_color(row['chromosome']))
            
        bars = mpl.collections.PolyCollection(verts,
                                            facecolors=colors,
                                            edgecolors="none")
        axis.add_collection(bars)

    if show_range:
        # Plot just the selected chromosome
        MB = plots.MB
        axis.set_xlabel("Position (Mb)")
        max_pos = plot_df['end'].max()
        axis.set_xlim((r_start or 0) * MB, (r_end or max_pos) * MB)
        axis.set_title(show_range)
        
        # Plot each sample
        for i, sample in enumerate(samples):
            sample_data = plot_df[plot_df['Sampleid'] == sample]
            if len(sample_data) == 0:
                logging.warning(f"Sample {sample} has no datapoints in selection {show_range}")
            sample_data['start'] *= MB
            sample_data['end'] *= MB
            plot_sample_chrom(i, sample_data)
            
    else:
        # Get chromosome sizes and offsets
        chrom_sizes = df.groupby('chromosome')['end'].max()
        curr_offset = 0
        chrom_offsets = {}
        for chrom in sorted(df['chromosome'].unique()):
            chrom_offsets[chrom] = curr_offset
            curr_offset += chrom_sizes[chrom] + 1
            
        # Plot chromosomes with offsets
        for i, sample in enumerate(samples):
            sample_data = plot_df[plot_df['Sampleid'] == sample]
            if len(sample_data) == 0:
                logging.warning(f"Sample {sample} has no datapoints")
                continue
                
            # Add offsets to each chromosome
            for chrom in chrom_offsets:
                chrom_data = sample_data[sample_data['chromosome'] == chrom].copy()
                if len(chrom_data) > 0:
                    chrom_data['start'] += chrom_offsets[chrom]
                    chrom_data['end'] += chrom_offsets[chrom]
                    plot_sample_chrom(i, chrom_data)
                    
        # Add chromosome labels
        for chrom, offset in chrom_offsets.items():
            axis.axvline(offset, color='k', linewidth=0.5, zorder=0)
            axis.text(offset + chrom_sizes[chrom]/2, -0.5, 
                     str(chrom), ha='center', va='top')
            
    return axis

# Use it like this:
plt.figure(figsize=(15, 10))
do_heatmap(cnvs)
plt.tight_layout()
plt.show()

