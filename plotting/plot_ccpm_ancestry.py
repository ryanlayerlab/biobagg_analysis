import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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

def read_ccpm_ancestries(chrm_ancestry_file,
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

def organize_ancestry_data(ccpm_ancestry_data,
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

def read_ccpm_genosis_scores(chrm_genosis_file,
                             ccpm_genosis_scores):
    '''
    Read the ccpm genosis scores file and return a dictionary with query id as key and match id and genosis score as values
    @param chrm_genosis_file: path to the ccpm genosis scores file
    @param ccpm_genosis_scores: dictionary with query id as key and match id and genosis score as values
    @return: dictionary with query id as key and match id and genosis score as values
    '''
    with open(chrm_genosis_file, 'r') as f:
        # Skip header
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            query_id = line[0].strip().split(',')[0]
            for match_id_score in line[1:]:
                match_id = match_id_score.split(',')[0]
                genosis_score = float(match_id_score.split(',')[1])
                try:
                    ccpm_genosis_scores[query_id][match_id] += genosis_score
                except KeyError:
                    try:
                        ccpm_genosis_scores[query_id][match_id] = genosis_score
                    except KeyError:
                        ccpm_genosis_scores[query_id] = {match_id: genosis_score}

    return ccpm_genosis_scores

def organize_genosis_data(ccpm_genosis_scores,
                            ccpm_ancestry_dict,
                            ccpm_ancestry_labels):
    '''
    Create a dictionary of dictionaries with ancestry as key and genosis scores as values

    @param ccpm_genosis_scores: dictionary with query id as key and match id and genosis score as values
    @param ccpm_ancestry_dict: dictionary with ccpm id as key and ancestry as value
    @param ccpm_ancestry_labels: list of ancestry labels
    @return:
    '''
    ancestry_scores = {ancestry: {ancestry: [] for ancestry in ccpm_ancestry_labels} for ancestry in ccpm_ancestry_labels}
    for query_id in ccpm_genosis_scores:
        query_scores = {ancestry: [] for ancestry in ccpm_ancestry_labels}
        query_ancestry = ccpm_ancestry_dict[query_id]
        for match_id in ccpm_genosis_scores[query_id]:
            if query_id == match_id:
                continue
            match_ancestry = ccpm_ancestry_dict[match_id]
            try:
                query_scores[match_ancestry].append(ccpm_genosis_scores[query_id][match_id])
            except KeyError:
                query_scores[match_ancestry] = [ccpm_genosis_scores[query_id][match_id]]
        for a in query_scores:
            ancestry_scores[query_ancestry][a].extend(query_scores[a])

    return ancestry_scores


def plot_data(ancestry_counts, png_file, num_chrom):
    '''
    Plot the ancestry data
    @param ancestry_counts: dictionary of dictionaries with ancestry as key and list of cohort counts as values
    @param png_file: path to the output png file
    '''
    ordered_ancestry = ['Africa', 'America', 'East Asian', 'Europe', 'Middle East', 'Central South Asian']
    # ordered_ancestry_labels = ['Africa', 'America', 'East\nAsian', 'Europe', 'Middle\nEast', 'Central\nSouth Asian']
    ordered_ancestry_labels = ['TGP+HGP-\nAFR-like',
                               'TGP+HGP-\nAMR-like',
                               'TGP+HGP-\nEAS-like',
                               'TGP+HGP-\nEUR-like',
                               'TGP+HGP-\nMLE-like',
                               'TGP+HGP-\nSAS-like']

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
                annot=True, fmt='.3f', annot_kws = {'size': 35, 'weight': 'bold'},
                vmin=0, vmax=20)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=25, pad=10)
    cbar.set_ticks([0, 5, 10, 15, 20])
    cbar.set_label('Average counts in cohort', fontsize=25, labelpad=20)
    # ax.set_title('Average cohort counts by ancestry\nCCPM (chrm 1-22)', fontsize=40, pad=20)
    ax.set_xlabel('Cohort Population', fontsize=30, labelpad=20, fontweight='bold')
    ax.set_ylabel('Query Population', fontsize=30, labelpad=20, fontweight='bold')
    ax.set_xticklabels(ordered_ancestry_labels, fontsize=25)
    ax.set_yticklabels(ordered_ancestry_labels, fontsize=25)

    # add rectangles around the diagonal colors
    color_CCPM = {'Africa': 'deepskyblue',
                  'America': 'goldenrod',
                  'East Asian': 'crimson',
                  'Europe': 'yellowgreen',
                  'Middle East': 'darkorange',
                  'Central South Asian': 'mediumpurple'}

    for i, a in enumerate(ordered_ancestry):
        ax.add_patch(Rectangle((i, i), 1, 1,
                               fill=False, edgecolor=color_CCPM[a], lw=8))

    plt.tight_layout()
    plt.savefig(png_file)

def plot_scores_new(ancestry_scores, png_file):
    '''
    plot distribution of genosis scores, one column for each ancestry
    @param ancestry_scores:
    @param png_file:
    @return:
    '''
    ordered_ancestry = ['Africa', 'America', 'East Asian', 'Europe', 'Middle East', 'Central South Asian']
    # ordered_ancestry_labels = ['Africa', 'America', 'East\nAsian', 'Europe', 'Middle\nEast', 'Central\nSouth Asian']
    ordered_ancestry_labels = ['TGP+HGP-\nAFR-like',
                               'TGP+HGP-\nAMR-like',
                               'TGP+HGP-\nEAS-like',
                               'TGP+HGP-\nEUR-like',
                               'TGP+HGP-\nMLE-like',
                               'TGP+HGP-\nSAS-like']

    color_CCPM = {'Africa': 'deepskyblue',
                  'America': 'goldenrod',
                  'East Asian': 'crimson',
                  'Europe': 'yellowgreen',
                  'Middle East': 'darkorange',
                  'Central South Asian': 'mediumpurple'}

    # 6 by 1
    fig, ax = plt.subplots(1, 6, figsize=(15, 3), dpi=300, sharex=True, sharey=True)

    # for i, a in enumerate(ordered_ancestry):
    #     all_scores = []
    #     for b in ordered_ancestry:
    #         all_scores.extend(ancestry_scores[a][b])
    #     sns.histplot(all_scores, ax=ax[i//3, i%3], bins=20, color=color_CCPM[a])
    #     # sns.histplot(ancestry_scores[a][a], ax=ax[i], bins=20, color=color_CCPM[a])
    #     # kde plot
    #     # sns.kdeplot(ancestry_scores[a][a], ax=ax[i], color='black')
    #
    #     # remove background grid and color
    #     # ax[j, i].grid(False)
    #     ax[i//3, i%3].set_facecolor('white')
    #     # ax[j, i].set_ylabel('Count', fontsize=10)
    #     # log scale
    #     ax[i//3, i%3].set_yscale('log')
    #
    #     # remove spines
    #     ax[i//3, i%3].spines['top'].set_visible(False)
    #     ax[i//3, i%3].spines['right'].set_visible(False)
    #
    #     # add title to be the ancestry
    #     ax[i//3, i%3].set_title(ordered_ancestry_labels[i],
    #                             fontsize=15, fontweight='bold', color=color_CCPM[a])
    #     # only label the x-axis for the bottom row
    #     if i//3 == 1:
    #         ax[i//3, i%3].set_xlabel('GenoSiS Score', fontsize=10)
    #     # only label the y-axis for the left column
    #     if i%3 == 0:
    #         ax[i//3, i%3].set_ylabel('Count', fontsize=10)

    for i, a in enumerate(ordered_ancestry):
        all_scores = []
        for b in ordered_ancestry:
            all_scores.extend(ancestry_scores[a][b])
        sns.histplot(all_scores, ax=ax[i], bins=20, color=color_CCPM[a])
        # sns.histplot(ancestry_scores[a][a], ax=ax[i], bins=20, color=color_CCPM[a])
        # kde plot
        # sns.kdeplot(ancestry_scores[a][a], ax=ax[i], color='black')

        # remove background grid and color
        ax[i].grid(False)
        ax[i].set_facecolor('white')
        ax[i].set_ylabel('Count', fontsize=10)
        # log scale
        ax[i].set_yscale('log')

        # remove spines
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)

        # add labels along the left side and bottom
        # add title to be the ancestry
        ax[i].set_title(ordered_ancestry_labels[i], fontsize=10, fontweight='bold', color=color_CCPM[a])
        # only label y-axis for the left column
        if i == 0:
            ax[i].set_ylabel('Count', fontsize=10)
        # label x-axis
        ax[i].set_xlabel('GenoSiS Score', fontsize=10)


        # ax[i].set_ylabel('Count', fontsize=10)
        # ax[i].set_xlabel('GenoSiS Score', fontsize=10)
        # ax[i].set_xticklabels(ordered_ancestry_labels, fontsize=20)
        # ax[i].set_yticklabels(ordered_ancestry_labels, fontsize=20)

    # add text box at bottom "Cohort Population"
    # fig.text(0.5, 0.01, 'Cohort Population', ha='center', fontsize=30, fontweight='bold')
    # add text box at left "Query Population"
    # fig.text(0.03, 0.5, 'Query Population', va='center', rotation='vertical', fontsize=30, fontweight='bold')

    # fig.suptitle('Genosis Score Distribution\nCCPM (chrm 1-22)', fontsize=40)

    # add custom legend

    plt.tight_layout()
    plt.savefig(png_file)



def plot_scores(ancestry_scores, png_file):
    '''
    plot distribution of genosis scores
    @param ancestry_scores: dictionary of dictionaries with ancestry as key and list of genosis scores as values
    @param png_file: path to the output png file
    @return: None
    '''
    ordered_ancestry = ['Africa', 'America', 'East Asian', 'Europe', 'Middle East', 'Central South Asian']
    # ordered_ancestry_labels = ['Africa', 'America', 'East\nAsian', 'Europe', 'Middle\nEast', 'Central\nSouth Asian']
    ordered_ancestry_labels = ['TGP+HGP-\nAFR-like',
                               'TGP+HGP-\nAMR-like',
                               'TGP+HGP-\nEAS-like',
                               'TGP+HGP-\nEUR-like',
                               'TGP+HGP-\nMLE-like',
                               'TGP+HGP-\nSAS-like']

    color_CCPM = {'Africa': 'deepskyblue',
                  'America': 'goldenrod',
                  'East Asian': 'crimson',
                  'Europe': 'yellowgreen',
                  'Middle East': 'darkorange',
                  'Central South Asian': 'mediumpurple'}

    # 6 by 6
    fig, ax = plt.subplots(6, 6, figsize=(18, 15), dpi=300, sharex=True, sharey=True)
    for i, a in enumerate(ordered_ancestry):
        for j, b in enumerate(ordered_ancestry):
            sns.histplot(ancestry_scores[a][b], ax=ax[i, j], bins=20, color=color_CCPM[a])
            # kde plot
            # sns.kdeplot(ancestry_scores[a][b], ax=ax[i, j], color='black')

            # remove background grid and color
            # ax[j, i].grid(False)
            ax[j, i].set_facecolor('white')
            # ax[j, i].set_ylabel('Count', fontsize=10)
            # log scale
            ax[j, i].set_yscale('log')

            # remove spines
            ax[j, i].spines['top'].set_visible(False)
            ax[j, i].spines['right'].set_visible(False)


    # add labels along the left side and bottom
    for i, a in enumerate(ordered_ancestry_labels):
        # ax[i, 0].set_ylabel(a, fontsize=20)
        # ax[5, i].set_xlabel(a, fontsize=20, labelpad=20)
        ax[5, i].set_xlabel(a, fontsize=20, labelpad=10)
        ax[i, 0].set_ylabel(a, fontsize=20, labelpad=10)
        # ax[i, 0].set_xticklabels(ordered_ancestry_labels, fontsize=20)
        # ax.set_yticklabels(ordered_ancestry_labels, fontsize=20)
        # ax[5, i].set_xlabel(a, fontsize=10)

    # shift whole plot up
    plt.subplots_adjust(top=0.99)

    # add text box at bottom "Cohort Population"
    fig.text(0.5, 0.01, 'Cohort Population', ha='center', fontsize=30, fontweight='bold')
    # add text box at left "Query Population"
    fig.text(0.03, 0.5, 'Query Population', va='center', rotation='vertical', fontsize=30, fontweight='bold')

    # fig.suptitle('Genosis Score Distribution\nCCPM (chrm 1-22)', fontsize=40)

    # add custom legend


    # plt.tight_layout()
    plt.savefig(png_file)

def plot_with_rank(ccpm_ancestry_data,
                   ccpm_genosis_scores,
                   ccpm_ancestry,
                   png_file):
    '''
    plot a scatter plot where x-axis is the cohort rank, y-axis is the genosis score,
     and the color is "in pop" or "out pop"
    @param ccpm_ancestry_data: cohort ancestry labels
    @param ccpm_genosis_scores: cohort genosis scores
    @param ccpm_ancestry: dict of samples and their ancestry labels
    @return: None
    '''

    self_color = 'black'
    in_pop_color = 'yellowgreen'
    out_pop_color = 'tomato'


    ancestry_options = ['Africa', 'America', 'East Asian', 'Europe', 'Middle East', 'Central South Asian']

    x_africa, y_africa, color_data_africa = [], [], []
    x_america, y_america, color_data_america = [], [], []
    x_east_asian, y_east_asian, color_data_east_asian = [], [], []
    x_europe, y_europe, color_data_europe = [], [], []
    x_middle_east, y_middle_east, color_data_middle_east = [], [], []
    x_central_south_asian, y_central_south_asian, color_data_central_south_asian = [], [], []

    for query_id in ccpm_ancestry_data:
        query_ancestry = ccpm_ancestry[query_id]
        for match_id in ccpm_ancestry_data[query_id]:
            match_ancestry = ccpm_ancestry[match_id]
            # rank is the index of the match in query's top k hits
            rank = list(ccpm_ancestry_data[query_id]).index(match_id) + 1

            if query_ancestry == 'Africa':
                x_africa.append(rank)
                y_africa.append(ccpm_genosis_scores[query_id][match_id])
                if query_id == match_id:
                    color_data_africa.append(self_color)
                elif query_ancestry == match_ancestry:
                    color_data_africa.append(in_pop_color)
                else:
                    color_data_africa.append(out_pop_color)
            elif query_ancestry == 'America':
                x_america.append(rank)
                y_america.append(ccpm_genosis_scores[query_id][match_id])
                if query_id == match_id:
                    color_data_america.append(self_color)
                elif query_ancestry == match_ancestry:
                    color_data_america.append(in_pop_color)
                else:
                    color_data_america.append(out_pop_color)
            elif query_ancestry == 'East Asian':
                x_east_asian.append(rank)
                y_east_asian.append(ccpm_genosis_scores[query_id][match_id])
                if query_id == match_id:
                    color_data_east_asian.append(self_color)
                elif query_ancestry == match_ancestry:
                    color_data_east_asian.append(in_pop_color)
                else:
                    color_data_east_asian.append(out_pop_color)
            elif query_ancestry == 'Europe':
                x_europe.append(rank)
                y_europe.append(ccpm_genosis_scores[query_id][match_id])
                if query_id == match_id:
                    color_data_europe.append(self_color)
                elif query_ancestry == match_ancestry:
                    color_data_europe.append(in_pop_color)
                else:
                    color_data_europe.append(out_pop_color)
            elif query_ancestry == 'Middle East':
                x_middle_east.append(rank)
                y_middle_east.append(ccpm_genosis_scores[query_id][match_id])
                if query_id == match_id:
                    color_data_middle_east.append(self_color)
                elif query_ancestry == match_ancestry:
                    color_data_middle_east.append(in_pop_color)
                else:
                    color_data_middle_east.append(out_pop_color)
            elif query_ancestry == 'Central South Asian':
                x_central_south_asian.append(rank)
                y_central_south_asian.append(ccpm_genosis_scores[query_id][match_id])
                if query_id == match_id:
                    color_data_central_south_asian.append(self_color)
                elif query_ancestry == match_ancestry:
                    color_data_central_south_asian.append(in_pop_color)
                else:
                    color_data_central_south_asian.append(out_pop_color)
            else:
                print('Ancestry not found...' + query_ancestry)
                continue
        else:
            continue

    # one plot for each of the 6 ancestries
    fig, ax = plt.subplots(2, 3, figsize=(18, 15), dpi=300, sharex=True, sharey=True)
    for i, a in enumerate(ancestry_options):
        print('Plotting for {}'.format(a))
        if a == 'Africa':
            x_data = x_africa
            y_data = y_africa
            color_data = color_data_africa
        elif a == 'America':
            x_data = x_america
            y_data = y_america
            color_data = color_data_america
        elif a == 'East Asian':
            x_data = x_east_asian
            y_data = y_east_asian
            color_data = color_data_east_asian
        elif a == 'Europe':
            x_data = x_europe
            y_data = y_europe
            color_data = color_data_europe
        elif a == 'Middle East':
            x_data = x_middle_east
            y_data = y_middle_east
            color_data = color_data_middle_east
        elif a == 'Central South Asian':
            x_data = x_central_south_asian
            y_data = y_central_south_asian
            color_data = color_data_central_south_asian
        else:
            print('Ancestry not found...' + a)
            continue

        # jitter the x values
        x_data = [x + np.random.normal(0, 0.1, 1)[0] for x in x_data]

        ax[i//3, i%3].scatter(x_data, y_data, c=color_data, s=10, alpha=0.5)
        ax[i//3, i%3].set_title(a, fontsize=20)
        ax[i//3, i%3].set_xlabel('Cohort Rank', fontsize=15)
        ax[i//3, i%3].set_ylabel('Genosis Score', fontsize=15)
        ax[i//3, i%3].set_xticks(range(1, 21, 1))
        ax[i//3, i%3].grid(False)
        # log scale
        ax[i//3, i%3].set_yscale('log')

        # add custom legend
        ax[i//3, i%3].scatter([], [], c=self_color, s=10, label='Self')
        ax[i//3, i%3].scatter([], [], c=in_pop_color, s=10, label='In Population')
        ax[i//3, i%3].scatter([], [], c=out_pop_color, s=10, label='Out Population')
        ax[i//3, i%3].legend(loc='upper right', fontsize=10)

        # remove spines
        ax[i//3, i%3].spines['top'].set_visible(False)
        ax[i//3, i%3].spines['right'].set_visible(False)

        # stop after 1 plot for testing
        # break

    fig.suptitle('Genosis Score vs Cohort Rank\nCCPM (chrm 1-22)', fontsize=40)
    plt.savefig(png_file)


def main():
    # Parse command line arguments
    args = parse_args()
    ccpm_ancestry_file = args.ancestry
    ancestry_dir = args.ancestry_dir
    png_dir = args.png

    num_chrom = 18

    print('Reading ccpm ancestry file')
    ccpm_ancestry = read_ccpm_ancestry(ccpm_ancestry_file)
    ccpm_ancestry_labels = list(set(ccpm_ancestry.values()))

    # {queryID: {matchID: count}}
    ccpm_ancestry_data = defaultdict(dict)
    ccpm_genosis_scores = defaultdict(dict)

    # print how many samples are in each ancestry
    for a in ccpm_ancestry_labels:
        print(f'{a}: {list(ccpm_ancestry.values()).count(a)}')

    for chrm_file in os.listdir(ancestry_dir):
        if chrm_file.endswith('_anc.txt'):
            print('Reading ccpm ancestries for {}'.format(chrm_file))
            chrm_ancestry_file = os.path.join(ancestry_dir, chrm_file)
            ccpm_ancestry_data = read_ccpm_ancestries(chrm_ancestry_file,
                                                ccpm_ancestry_data)
        elif chrm_file.endswith('_scores.txt'):
            print('Reading ccpm scores for {}'.format(chrm_file))
            chrm_genosis_file = os.path.join(ancestry_dir, chrm_file)
            ccpm_genosis_scores = read_ccpm_genosis_scores(chrm_genosis_file,
                                                ccpm_genosis_scores)

    ancestry_counts = organize_ancestry_data(ccpm_ancestry_data,
                                             ccpm_ancestry,
                                             ccpm_ancestry_labels)

    ancestry_scores = organize_genosis_data(ccpm_genosis_scores,
                                             ccpm_ancestry,
                                             ccpm_ancestry_labels)

    # Plot the data
    heatmap_png = png_dir + 'ccpm_ancestry.png'
    distribution_png = png_dir + 'ccpm_genosis_scores.png'
    new_png = png_dir + 'ccpm_genosis_scores_new.png'
    # scatter_png = png_dir + 'ccpm_rank_scores_scatter.png'

    # plot_data(ancestry_counts, heatmap_png, num_chrom)
    # plot_scores(ancestry_scores, distribution_png)
    # plot_with_rank(ccpm_ancestry_data, ccpm_genosis_scores, ccpm_ancestry, scatter_png)
    plot_scores_new(ancestry_scores, new_png)



if __name__ == '__main__':
    main()