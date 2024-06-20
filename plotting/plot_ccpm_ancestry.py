import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestry', help='ccpm ancestry labels', required=True)
    parser.add_argument('--ancestry_dir', help='ccpm ancestry results dir', required=True)
    parser.add_argument('--png', help='output png file', required=True)

    return parser.parse_args()

def read_ccpm_ancestry(ccpm_ancestry_file):
    '''
    Read the ccpm ancestry file and return a dictionary with ccpm_id as key and ancestry as value
    @param ccpm_ancestry_file: path to the ccpm ancestry file
    @return: dictionary with ccpm_id as key and ancestry as value
    '''
    ccpm_ancestry = dict()

    with open(ccpm_ancestry_file, 'r') as f:
        # Skip header
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            ccpm_id = line[0]
            ancestry = line[1]
            ccpm_ancestry[ccpm_id] = ancestry

    return ccpm_ancestry

def read_ccpm_hits(chrm_ancestry_file,
                   ccpm_data):
    '''
    Read the ccpm ancestry file and return a dictionary with ccpm id as key and match ids as values
    @param ccpm_ancestry_file: path to the ccpm single chrm ancestry file
    @return: dictionary with ccpm_id, as key and match ids as values
    '''

    with open(chrm_ancestry_file, 'r') as f:
        # Skip header
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            query_id = line[0].strip().split(',')[0]
            for match_id_anc in line[1:]:
                match_id = match_id_anc.split(',')[0]
                try:
                    ccpm_data[query_id][match_id] += 1
                except KeyError:
                    try:
                        ccpm_data[query_id][match_id] = 1
                    except KeyError:
                        ccpm_data[query_id] = {match_id: 1}

    return ccpm_data

def organize_data(ccpm_ancestry_data,
                  ccpm_ancestry_dict,
                  ccpm_ancestry_labels):
    '''
    Create a dictionary of dictionaries with ancestry as key and count as value s
    @param ccpm_ancestry_data: dictionary with ccpm id as key and match id as value
    @param ccpm_ancestry_dict: dictionary with ccpm id as key and ancestry as value
    @return: dictionary of dictionaries with ancestry as key and list of cohort counts as values
    '''
    ancestry_counts = {ancestry: {ancestry: [] for ancestry in ccpm_ancestry_labels} for ancestry in ccpm_ancestry_labels}
    for query_id in ccpm_ancestry_data:
        query_counts = {ancestry: 0 for ancestry in ccpm_ancestry_labels}
        query_ancestry = ccpm_ancestry_dict[query_id]
        for match_id in ccpm_ancestry_data[query_id]:
            match_ancestry = ccpm_ancestry_dict[match_id]
            try:
                query_counts[match_ancestry] += 1
            except KeyError:
                query_counts[match_ancestry] = 1
        for a in query_counts:
            ancestry_counts[query_ancestry][a].append(query_counts[a])

    return ancestry_counts


def plot_data(ancestry_counts, png_file, num_chrom):
    '''
    Plot the ancestry data
    @param ancestry_counts: dictionary of dictionaries with ancestry as key and list of cohort counts as values
    @param png_file: path to the output png file
    '''
    ordered_ancestry = ['Africa', 'America', 'East Asian', 'Europe', 'Middle East', 'Central South Asian']
    ordered_ancestry_labels = ['Africa', 'America', 'East\nAsian', 'Europe', 'Middle\nEast', 'Central\nSouth Asian']
    # Plot with heatmap ancestry labels as x and y in order
    fig, ax = plt.subplots(figsize=(18, 15), dpi=300)
    # get average counts for each ancestry in order
    # avg_counts = {a: [np.mean(ancestry_counts[a][b])/num_chrom for b in ancestry_counts[a]] for a in ancestry_counts}
    avg_counts = {a: [np.mean(ancestry_counts[a][b])/ num_chrom for b in ordered_ancestry] for a in ordered_ancestry}
    df = pd.DataFrame(avg_counts)
    sns.set(font_scale=1.5)
    # Plot the heatmap transposed with range 0 to 20
    sns.heatmap(df.T, cmap='Greys', ax=ax,
                square=True,
                annot=False, fmt='.2f', annot_kws={'size': 20},
                vmin=0, vmax=20)
    cbar = ax.collections[0].colorbar
    cbar.set_label('Average number of hits in cohort', fontsize=20, labelpad=20)
    ax.set_title('Average cohort counts by ancestry\nCCPM (chrm 1-11,20-22)', fontsize=40, pad=20)
    ax.set_xlabel('Cohort Population', fontsize=35, labelpad=20)
    ax.set_ylabel('Query Population', fontsize=35, labelpad=20)
    ax.set_xticklabels(ordered_ancestry_labels, fontsize=20)
    ax.set_yticklabels(ordered_ancestry_labels, fontsize=20)


    plt.tight_layout()
    plt.savefig(png_file)


def main():
    # Parse command line arguments
    args = parse_args()
    ccpm_ancestry_file = args.ancestry
    ancestry_dir = args.ancestry_dir
    png_file = args.png

    num_chrom = 13

    print('Reading ccpm ancestry file')
    ccpm_ancestry = read_ccpm_ancestry(ccpm_ancestry_file)
    ccpm_ancestry_labels = list(set(ccpm_ancestry.values()))

    # {queryID: {matchID: count}}
    ccpm_ancestry_data = defaultdict(dict)
    for chrm_file in os.listdir(ancestry_dir):
        if chrm_file.endswith('.txt'):
            print('Reading ccpm hits for {}'.format(chrm_file))
            chrm_ancestry_file = os.path.join(ancestry_dir, chrm_file)
            ccpm_ancestry_data = read_ccpm_hits(chrm_ancestry_file,
                                                ccpm_ancestry_data)

    ancestry_counts = organize_data(ccpm_ancestry_data,
                                    ccpm_ancestry,
                                    ccpm_ancestry_labels)

    # Plot the data
    plot_data(ancestry_counts, png_file, num_chrom)



if __name__ == '__main__':
    main()