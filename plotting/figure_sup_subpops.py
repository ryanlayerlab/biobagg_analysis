import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys

sys.path.append(os.path.abspath('plotting/'))
import ancestry_helpers

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestry', help='1KG ancestry labels', required=True)
    parser.add_argument('--k', help='k used for knn', required=True)
    parser.add_argument('--colors', help='file with color codes', required=True)
    # INPUT FILES
    parser.add_argument('--genosis', type=str, required=True, help='Path to genosis cohort file')
    parser.add_argument('--dst', type=str, required=True, help='Path to plink dst top 20 file')
    parser.add_argument('--pihat', type=str, required=True, help='Path to plink pi-hat top 20 file')
    parser.add_argument('--kinship', type=str, required=True, help='Path to plink kinship top 20 file')
    # OUTPUT
    parser.add_argument('--png', help='Output png file with top k percents', required=True)

    return parser.parse_args()

def get_subpop_number_samples(ancestry_file):
    '''
    Get the number of samples for each subpopulation
    @param ancestry_file: path to ancestry file
    @return: dictionary of number of samples for each subpopulation
    '''
    sample_subpopulations = ancestry_helpers.get_subpopulations(ancestry_file)
    subpop_number_samples = defaultdict(int)
    for sample in sample_subpopulations:
        subpop = sample_subpopulations[sample]
        try:
            subpop_number_samples[subpop] += 1
        except KeyError:
            subpop_number_samples[subpop] = 1
    return subpop_number_samples


def read_subpop_counts_old(subpop_counts_file):
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

def read_subpop_counts(subpop_counts_file,
                       header=True):
    '''
    Read subpopulation counts from file
    @param subpop_counts_file: path to subpopulation counts file
    @return: dictionary of subpopulation counts
    '''
    # query_subpop	match_subpop,count...
    subpop_counts = defaultdict(dict)
    f = open(subpop_counts_file, 'r')
    if header:
        f.readline()
    for line in f:
        line = line.strip().split()
        query_subpop = line[0]
        match_subpop = line[1]
        match_subpop_counts = line[2].strip().split(',')
        subpop_counts[query_subpop][match_subpop] = match_subpop_counts

    return subpop_counts

def plot_subpop_counts_heatmap(subpop_counts,
                                subpopulation_number_samples,
                                colors,
                                output_file,
                                title):
    SUPER_SUBPOPULATIONS = ancestry_helpers.SUPER_SUBPOPULATIONS
    num_subpops = len(subpop_counts)
    all_subpops = sorted(subpop_counts.keys())
    ordered_subpops = []
    for superpop in SUPER_SUBPOPULATIONS:
        subpops = sorted([subpop for subpop in all_subpops if ancestry_helpers.SUB_SUPERPOPULATIONS[subpop] == superpop])
        ordered_subpops.extend(subpops)
    colormap = 'Greys'

    # [row_idx, col_idx]

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

                # # normalize by size of subpopulation
                # counts[i+index_offset, j] /= subpopulation_number_samples[query_subpop]
            # print(query_subpop, subpopulation_number_samples[query_subpop])

        index_offset += len(subpops)

    # make a dataframe
    df = pd.DataFrame(counts, index=ordered_subpops, columns=ordered_subpops)
    # plot
    plt.figure(figsize=(18, 15), dpi=300)
    sns.set(font_scale=1.5)
    row_colors = [colors[ancestry_helpers.SUB_SUPERPOPULATIONS[subpop]] for subpop in ordered_subpops]

    # log scale
    sns.heatmap(df, cmap=colormap, annot=False, cbar=True, square=True,
            xticklabels=False, yticklabels=False,
            norm=LogNorm(),
            cbar_kws={'label': 'Cohort Counts'})
    # with log, add this to heatmatp plot above
    # norm = LogNorm()

    # add patches for superpopulations
    line_offset = 0
    for i in SUPER_SUBPOPULATIONS:
        num_subpops = len(SUPER_SUBPOPULATIONS[i])
        plt.gca().add_patch(Rectangle((line_offset, line_offset), num_subpops, num_subpops,
                                      edgecolor=colors[i],
                                      facecolor='none',
                                      alpha=1, lw=15))
        line_offset += num_subpops

    # add xtick labels and color by superpopulation
    xtick_labels = []
    ytick_labels = []
    xtick_positions = []
    ytick_positions = []
    tick_offset = 0
    x_i = 0
    for i in SUPER_SUBPOPULATIONS:
        num_subpops = len(SUPER_SUBPOPULATIONS[i])
        xtick_labels.extend([subpop for subpop in ordered_subpops[tick_offset:tick_offset+num_subpops]])
        ytick_labels.append(i)
        xtick_positions.extend(range(x_i + 1, x_i + num_subpops + 1))
        ytick_positions.append(tick_offset + num_subpops//2)
        x_i += num_subpops
        tick_offset += num_subpops
    plt.xticks(xtick_positions, xtick_labels,
               rotation=90, fontsize=25, ha='right')
    for tick in plt.gca().get_xticklabels():
        tick.set_color(colors[ancestry_helpers.SUB_SUPERPOPULATIONS[tick.get_text()]])


    plt.yticks(ytick_positions, ytick_labels,
               rotation=90, fontsize=30, fontweight='bold',
               va='center')
    for tick in plt.gca().get_yticklabels():
        tick.set_color(colors[tick.get_text()])


    # draw lines between superpopulations
    line_offset = 0
    for i in SUPER_SUBPOPULATIONS:
        num_subpops = len(SUPER_SUBPOPULATIONS[i])
        plt.axhline(num_subpops + line_offset, color='black', linewidth=2)
        plt.axvline(num_subpops + line_offset, color='black', linewidth=2)
        line_offset += num_subpops


    plt.xlabel('Cohort Population', fontsize=40)
    plt.ylabel('Query Population', fontsize=40)
    title += '\nTotal Counts for Subpopulation Cohorts'
    plt.title(title, fontsize=50)
    # move plot to make room for title
    plt.subplots_adjust(top=0.9, bottom=0.2, right=0.9, left=0.1)
    plt.tight_layout()
    plt.savefig(output_file)

def plot_heatmap_counts(subpop_counts,
                        colors,
                        title,
                        output_file):

    SUPER_SUBPOPULATIONS = ancestry_helpers.SUPER_SUBPOPULATIONS

    ordered_subpop = []
    for superpop in SUPER_SUBPOPULATIONS:
        for subpop in SUPER_SUBPOPULATIONS[superpop]:
            ordered_subpop.append(subpop)


    # get average counts for each subpopulation
    subpop_counts_avg = {}
    for subpop in ordered_subpop:
        for match_subpop in ordered_subpop:
            avg = sum([int(count) for count in subpop_counts[subpop][match_subpop]]) / len(subpop_counts[subpop][match_subpop])
            try:
                subpop_counts_avg[subpop][match_subpop] = avg
            except KeyError:
                subpop_counts_avg[subpop] = {match_subpop: avg}

    plt.figure(figsize=(18, 15), dpi=300)
    df = pd.DataFrame(subpop_counts_avg)
    sns.set(font_scale=1.5)
    sns.heatmap(df.T, cmap='Greys', annot=False, cbar=True, square=True,
                vmin=0, vmax=20,
                cbar_kws={'label': 'Average number of hits in cohort'})
    # add patches for superpopulations
    line_offset = 0
    for i in SUPER_SUBPOPULATIONS:
        num_subpops = len(SUPER_SUBPOPULATIONS[i])
        plt.gca().add_patch(Rectangle((line_offset, line_offset), num_subpops, num_subpops,
                                      edgecolor=colors[i],
                                      facecolor='none',
                                      alpha=1, lw=15))
        line_offset += num_subpops

    # add xtick labels and color by superpopulation
    xtick_labels = []
    ytick_labels = []
    xtick_positions = []
    ytick_positions = []
    tick_offset = 0
    x_i = 0
    for i in SUPER_SUBPOPULATIONS:
        num_subpops = len(SUPER_SUBPOPULATIONS[i])
        xtick_labels.extend([subpop for subpop in df.columns[tick_offset:tick_offset+num_subpops]])
        ytick_labels.append(i)
        xtick_positions.extend(range(x_i + 1, x_i + num_subpops + 1))
        ytick_positions.append(tick_offset + num_subpops//2)
        x_i += num_subpops
        tick_offset += num_subpops
    plt.xticks(xtick_positions, xtick_labels,
               rotation=90, fontsize=25, ha='right')
    for tick in plt.gca().get_xticklabels():
        tick.set_color(colors[ancestry_helpers.SUB_SUPERPOPULATIONS[tick.get_text()]])

    plt.yticks(ytick_positions, ytick_labels,
                rotation=90, fontsize=30, fontweight='bold',
                va='center')
    for tick in plt.gca().get_yticklabels():
        tick.set_color(colors[tick.get_text()])

    # draw lines between superpopulations except last one
    line_offset = 0
    for i in SUPER_SUBPOPULATIONS:
        if i == 'SAS':
            break
        num_subpops = len(SUPER_SUBPOPULATIONS[i])
        plt.axhline(num_subpops + line_offset, color='black', linewidth=2)
        plt.axvline(num_subpops + line_offset, color='black', linewidth=2)
        line_offset += num_subpops


    plt.xlabel('Cohort Population', fontsize=40)
    plt.ylabel('Query Population', fontsize=40)
    figure_title = ('Average cohort counts by Subpopulation' +
                    '\n(1KG, ' + title + ')' )

    # set colorbar font and label
    cbar = plt.gca().collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    cbar.set_label('Average number of hits in cohort', fontsize=20, labelpad=20)

    plt.title(figure_title, fontsize=40, pad=20)
    # move plot to make room for title
    plt.tight_layout()
    plt.savefig(output_file)


def main():
    args = parse_args()

    subpopulation_number_samples = get_subpop_number_samples(args.ancestry)
    colors = ancestry_helpers.get_colors(args.colors)

    genosis_subpop_counts = read_subpop_counts(args.genosis,
                                               header=True)
    plot_heatmap_counts(genosis_subpop_counts,
                        colors,
                        'Genosis',
                        args.png + 'genosis_counts.png')

    dst_subpop_counts = read_subpop_counts(args.dst,
                                             header=True)
    plot_heatmap_counts(dst_subpop_counts,
                        colors,
                        'plink DST',
                        args.png + 'dst_counts.png')

    pihat_subpop_counts = read_subpop_counts(args.pihat,
                                                header=True)
    plot_heatmap_counts(pihat_subpop_counts,
                        colors,
                        'plink pi-hat',
                        args.png + 'pihat_counts.png')

    kinship_subpop_counts = read_subpop_counts(args.kinship,
                                                header=True)
    plot_heatmap_counts(kinship_subpop_counts,
                        colors,
                        'plink2 kinship',
                        args.png + 'kinship_counts.png')



if __name__ == '__main__':
    main()