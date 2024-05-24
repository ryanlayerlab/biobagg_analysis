import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ks_2samp
from scipy.stats import norm

import plotting.ancestry_helpers as ancestry_helpers

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestry', help='1KG ancestry labels', required=True)
    parser.add_argument('--k', help='k used for knn', required=True)
    parser.add_argument('--genosis', type=str, required=True, help='Path to genosis cohort file')
    parser.add_argument('--png', help='Output png file with top k percents', required=True)

    return parser.parse_args()

def get_ancestry_colors():
    return {
        'AFR': 'darkorange',
        'AMR': 'mediumpurple',
        'EAS': 'deeppink',
        'EUR': 'steelblue',
        'SAS': 'goldenrod'
    }

def read_subpop_counts(subpop_counts_file):
    '''
    Read subpopulation counts from file
    @param subpop_counts_file: path to subpopulation counts file
    @return: dictionary of subpopulation counts
    '''
    ## query_subpop	match_subpop,count...
    subpop_counts = defaultdict(dict)
    f = open(subpop_counts_file, 'r')
    f.readline()
    for line in f:
        line = line.strip().split()
        query_subpop = line[0]
        match_subpop_counts = line[1:]
        for match_subpop_count in match_subpop_counts:
            match_subpop, count = match_subpop_count.split(',')
            subpop_counts[query_subpop][match_subpop] = int(count)
    return subpop_counts

def plot_subpop_counts_heatmatp(subpop_counts,
                                output_file):
    SUPER_SUBPOPULATIONS = ancestry_helpers.SUPER_SUBPOPULATIONS
    num_subpops = len(subpop_counts)
    ordered_subpops = []

    # plot should be ordered and colored by superpopulation
    # make a matrix of counts
    index_offset = 0
    counts = np.zeros((num_subpops, num_subpops))
    for superpop in SUPER_SUBPOPULATIONS:
        subpops = sorted([subpop for subpop in subpop_counts if ancestry_helpers.SUB_SUPERPOPULATIONS[subpop] == superpop])
        ordered_subpops.extend(subpops)
        for i, query_subpop in enumerate(subpops):
            for j, match_subpop in enumerate(subpops):
                counts[i+index_offset, j+index_offset] = subpop_counts[query_subpop][match_subpop]
        index_offset += len(subpops)


    # make a dataframe
    df = pd.DataFrame(counts, index=ordered_subpops, columns=ordered_subpops)

    # # plot
    # plt.figure(figsize=(10, 10))

    # plot
    plt.figure(figsize=(10, 10), dpi=300)
    sns.set(font_scale=1.5)
    sns.heatmap(df, cmap='Greys', annot=False, fmt='g', cbar=True, square=True, xticklabels=True, yticklabels=True)

    # draw lines between superpopulations
    line_offset = 0
    for i in SUPER_SUBPOPULATIONS:
        num_subpops = len(SUPER_SUBPOPULATIONS[i])
        plt.axhline(num_subpops + line_offset, color='black', linewidth=2)
        plt.axvline(num_subpops + line_offset, color='black', linewidth=2)
        line_offset += num_subpops


    plt.xlabel('Match Subpopulation')
    plt.ylabel('Query Subpopulation')
    plt.title('Subpopulation Counts')
    plt.tight_layout()
    plt.savefig(output_file)

def main():
    args = parse_args()

    genosis_subpop_counts = read_subpop_counts(args.genosis)
    plot_subpop_counts_heatmatp(genosis_subpop_counts, args.png)


if __name__ == '__main__':
    main()