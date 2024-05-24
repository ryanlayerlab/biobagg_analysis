import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm
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
        'EUR': 'dodgerblue',
        'SAS': 'goldenrod'
    }
    # return {
    #     'AFR': 'copper',
    #     'AMR': 'Purples',
    #     'EAS': 'Reds',
    #     'EUR': 'Blues',
    #     'SAS': 'Wistia'
    # }

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
    all_subpops = sorted(subpop_counts.keys())
    ordered_subpops = []
    for superpop in SUPER_SUBPOPULATIONS:
        subpops = sorted([subpop for subpop in all_subpops if ancestry_helpers.SUB_SUPERPOPULATIONS[subpop] == superpop])
        ordered_subpops.extend(subpops)
    colormap = 'Greys'

    # make a matrix of counts
    index_offset = 0
    counts = np.zeros((num_subpops, num_subpops))
    for superpop in SUPER_SUBPOPULATIONS:
        subpops = sorted([subpop for subpop in subpop_counts if ancestry_helpers.SUB_SUPERPOPULATIONS[subpop] == superpop])
        for i, query_subpop in enumerate(subpops):
            for j, match_subpop in enumerate(ordered_subpops):
                counts[i+index_offset, j] = subpop_counts[query_subpop][match_subpop]
                ## handle zeros for log
                if counts[i+index_offset, j] == 0:
                    counts[i+index_offset, j] = 1
        index_offset += len(subpops)

    # make a dataframe
    df = pd.DataFrame(counts, index=ordered_subpops, columns=ordered_subpops)
    # plot
    plt.figure(figsize=(18, 15), dpi=300)
    sns.set(font_scale=1.5)
    row_colors = [get_ancestry_colors()[ancestry_helpers.SUB_SUPERPOPULATIONS[subpop]] for subpop in ordered_subpops]

    sns.heatmap(df, cmap=colormap, annot=False, cbar=True, square=True,
                xticklabels=False, yticklabels=False,
                norm=LogNorm(), cbar_kws={'label': 'Cohort Counts'})

    # add patches for superpopulations
    line_offset = 0
    for i in SUPER_SUBPOPULATIONS:
        num_subpops = len(SUPER_SUBPOPULATIONS[i])
        plt.gca().add_patch(Rectangle((line_offset, line_offset), num_subpops, num_subpops,
                                      edgecolor=get_ancestry_colors()[i],
                                      facecolor=get_ancestry_colors()[i],
                                      alpha=0.3, lw=5))
        line_offset += num_subpops

    # add xtick labels and color by superpopulation
    xtick_labels = []
    ytick_labels = []
    xtick_positions = []
    ytick_positions = []
    tick_offset = 0
    x_i = 0
    y_i = 0
    for i in SUPER_SUBPOPULATIONS:
        num_subpops = len(SUPER_SUBPOPULATIONS[i])
        xtick_labels.extend([subpop for subpop in ordered_subpops[tick_offset:tick_offset+num_subpops]])
        ytick_labels.append(i)
        xtick_positions.extend(range(x_i, x_i + num_subpops))
        ytick_positions.append(tick_offset + num_subpops//2)
        x_i += num_subpops
        tick_offset += num_subpops
    plt.xticks(xtick_positions, xtick_labels,
               rotation=90, fontsize=25)
    for tick in plt.gca().get_xticklabels():
        tick.set_color(get_ancestry_colors()[ancestry_helpers.SUB_SUPERPOPULATIONS[tick.get_text()]])


    plt.yticks(ytick_positions, ytick_labels,
               rotation=90, fontsize=30, fontweight='bold',
               va='center')
    for tick in plt.gca().get_yticklabels():
        tick.set_color(get_ancestry_colors()[tick.get_text()])


    # draw lines between superpopulations
    line_offset = 0
    for i in SUPER_SUBPOPULATIONS:
        num_subpops = len(SUPER_SUBPOPULATIONS[i])
        plt.axhline(num_subpops + line_offset, color='black', linewidth=2)
        plt.axvline(num_subpops + line_offset, color='black', linewidth=2)
        line_offset += num_subpops


    plt.xlabel('Match Population', fontsize=40)
    plt.ylabel('Query Population', fontsize=40)
    plt.title('Total Counts for\nSubpopulation Cohorts', fontsize=50)
    # move plot to make room for title
    plt.subplots_adjust(top=0.9, bottom=0.2, right=0.9, left=0.1)
    plt.tight_layout()
    plt.savefig(output_file)

def main():
    args = parse_args()

    genosis_subpop_counts = read_subpop_counts(args.genosis)
    plot_subpop_counts_heatmatp(genosis_subpop_counts, args.png)


if __name__ == '__main__':
    main()