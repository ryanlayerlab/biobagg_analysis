import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
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
    # INPUT DIR
    parser.add_argument('--quality_dir', type=str, required=True,
                        help='Path to directory of quality experiments')
    # OUTPUT
    parser.add_argument('--png', help='Output png prefix', required=True)


    return parser.parse_args()

def read_quality_results(query_pop_file):
    '''
    Read quality results from directory
    @param quality_dir: path to directory of quality results
    @return: dictionary of quality results
    '''
    quality_scores = defaultdict(dict)
    f = open(query_pop_file, 'r')
    header = f.readline()
    query_population = header.strip().split(':')[1].strip()
    for line in f:
        if 'database' in line:
            # new database
            database_population = line.strip().split(':')[1].strip()
            continue
        query_subpopulation = line.strip().split('->')[0].strip()
        match_subpopulation = line.strip().split('->')[1].split(':')[0].strip()
        database_quality_scores = line.strip().split('->')[1].split(':')[1].split(',')
        # remove last element which is empty
        database_quality_scores = database_quality_scores[:-1]
        database_quality_scores = [float(score) for score in database_quality_scores]

        try:
            quality_scores[database_population][query_subpopulation][match_subpopulation] = database_quality_scores
        except KeyError:
            try:
                quality_scores[database_population][query_subpopulation] = {}
                quality_scores[database_population][query_subpopulation][match_subpopulation] = database_quality_scores
            except KeyError:
                quality_scores[database_population] = {}
                quality_scores[database_population][query_subpopulation] = {}
                quality_scores[database_population][query_subpopulation][match_subpopulation] = database_quality_scores

    f.close()
    return quality_scores

def combine_subpopulations(quality_scores):
    '''
    Combine subpopulations into superpopulations, keeping track of subpop matches as own
    @param quality_scores: dictionary of quality scores
    @return: dictionary of quality scores by superpopulation
    '''
    super_quality_scores = dict()
    same_subpop = []
    for database_population in quality_scores.keys():
        for query_subpopulation in quality_scores[database_population].keys():
            for match_subpopulation in quality_scores[database_population][query_subpopulation].keys():
                if query_subpopulation == match_subpopulation:
                    same_subpop.extend(quality_scores[database_population][query_subpopulation][match_subpopulation])
                else:
                    try:
                        super_quality_scores[database_population].extend(
                            quality_scores[database_population][query_subpopulation][match_subpopulation])
                    except KeyError:
                        super_quality_scores[database_population] = (
                            quality_scores)[database_population][query_subpopulation][match_subpopulation]
        super_quality_scores['same_subpop'] = same_subpop
    return super_quality_scores

def combine_subpopulations_old(quality_scores):
    '''
    Combine subpopulations into superpopulations
    @param quality_scores: dictionary of quality scores
    @return: dictionary of quality scores by superpopulation
    '''
    super_quality_scores = dict()
    for database_population in quality_scores.keys():
        database_scores = []
        for query_subpopulation in quality_scores[database_population].keys():
            for match_subpopulation in quality_scores[database_population][query_subpopulation].keys():
                database_scores.extend(quality_scores[database_population][query_subpopulation][match_subpopulation])
        super_quality_scores[database_population] = database_scores
    return super_quality_scores

def organize_subpopulations(quality_scores):
    '''
    Organize quality scores by subpopulation
    @param quality_scores: dictionary of quality scores
    @return: dictionary of quality scores by subpopulation
    '''
    super_quality_scores = defaultdict(dict)
    for database_superpop in quality_scores.keys():
        for database_subpop in quality_scores[database_superpop].keys():
            for scores_subpop in quality_scores[database_superpop][database_subpop].keys():
                super_quality_scores[database_subpop][scores_subpop] = (
                    quality_scores)[database_superpop][database_subpop][scores_subpop]
    return super_quality_scores


def plot_super_histogram(AFR_super_scores,
                         AMR_super_scores,
                         EAS_super_scores,
                         EUR_super_scores,
                         SAS_super_scores,
                         colors,
                         png):
    SUPER_SUBPOPULATIONS = ancestry_helpers.SUPER_SUBPOPULATIONS

    fig, ax = plt.subplots(5, 5,
                           figsize=(20, 15), dpi=300,
                           sharex=True, sharey=True)
    alpha = 1
    suppop_alpha = 1
    subpop_color = 'black'


    # African query is first column
    ax[0, 0].set_title('AFR Query', fontsize=20, color=colors['AFR'], fontweight='bold')
    for i, database_population in enumerate(SUPER_SUBPOPULATIONS.keys()):
        if database_population == 'AFR':
            # plot subpop scores in black, with curve
            # ax[i, 0].hist(AFR_super_scores['same_subpop'],
            #                 bins=50, color=subpop_color, alpha=suppop_alpha)
            sns.histplot(AFR_super_scores[database_population], color=colors[database_population],
                         ax=ax[i, 0], binwidth=60, alpha=alpha, linewidth=0.0)
            sns.histplot(AFR_super_scores['same_subpop'], color=subpop_color,
                         ax=ax[i, 0], binwidth=60, alpha=suppop_alpha, fill=False, linewidth=0.3)

            # ax[i, 0].hist(AFR_super_scores[database_population],
            #               bins=50, color=colors[database_population], alpha=alpha)
            # ax[i, 0].hist(AFR_super_scores['same_subpop'],
            #               bins=50, color=subpop_color, alpha=suppop_alpha, histtype='step')
        else:
            sns.histplot(AFR_super_scores[database_population], color=colors[database_population],
                         ax=ax[i, 0], alpha=alpha, binwidth=60, linewidth=0.0)
            # ax[i, 0].hist(AFR_super_scores[database_population],
            #           bins=50, color=colors[database_population], alpha=alpha)
        # set y-axis label
        ax[i, 0].set_ylabel(database_population + '\nDatabase',
                            fontsize=25, color=colors[database_population],
                            fontweight='bold', labelpad=40)
        # set x-axis label
        ax[4, i].set_xlabel('GenoSiS Score', fontsize=20, labelpad=20)
        # log scale
        ax[i, 1].set_yscale('log')

    # American query is second column
    ax[0, 1].set_title('AMR Query', fontsize=20, color=colors['AMR'], fontweight='bold')
    for i, database_population in enumerate(SUPER_SUBPOPULATIONS.keys()):
        if database_population == 'AMR':
            # plot subpop scores in black
            # ax[i, 1].hist(AMR_super_scores['same_subpop'],
            #                 bins=50, color=subpop_color, alpha=suppop_alpha)
            # ax[i, 1].hist(AMR_super_scores[database_population],
            #               bins=50, color=colors[database_population], alpha=alpha)
            sns.histplot(AMR_super_scores[database_population], color=colors[database_population],
                         ax=ax[i, 1], alpha=alpha, binwidth=60, linewidth=0.0)
            sns.histplot(AMR_super_scores['same_subpop'], color=subpop_color,
                            ax=ax[i, 1], alpha=suppop_alpha, binwidth=60, fill=False, linewidth=0.3)
        else:
            # ax[i, 1].hist(AMR_super_scores[database_population],
            #               bins=50, color=colors[database_population], alpha=alpha)
            sns.histplot(AMR_super_scores[database_population], color=colors[database_population],
                            ax=ax[i, 1], alpha=alpha, binwidth=60, linewidth=0.0)
        # set x-axis label
        ax[4, i].set_xlabel('GenoSiS Score', fontsize=20, labelpad=20)


    # # East Asian query is third column
    ax[0, 2].set_title('EAS Query', fontsize=20, color=colors['EAS'], fontweight='bold')
    for i, database_population in enumerate(SUPER_SUBPOPULATIONS.keys()):
        if database_population == 'EAS':
            # plot subpop scores in black
            # ax[i, 2].hist(EAS_super_scores['same_subpop'],
            #                 bins=50, color=subpop_color, alpha=suppop_alpha)
            # ax[i, 2].hist(EAS_super_scores[database_population],
            #               bins=50, color=colors[database_population], alpha=alpha)
            sns.histplot(EAS_super_scores[database_population], color=colors[database_population],
                         ax=ax[i, 2], alpha=alpha, binwidth=60, linewidth=0.0)
            sns.histplot(EAS_super_scores['same_subpop'], color=subpop_color,
                            ax=ax[i, 2], alpha=suppop_alpha, binwidth=60, fill=False, linewidth=0.3)
        else:
            # ax[i, 2].hist(EAS_super_scores[database_population],
            #               bins=50, color=colors[database_population], alpha=alpha)
            sns.histplot(EAS_super_scores[database_population], color=colors[database_population],
                            ax=ax[i, 2], alpha=alpha, binwidth=60, linewidth=0.0)
        # set x-axis label
        ax[4, i].set_xlabel('GenoSiS Score', fontsize=20, labelpad=20)

    # European query is fourth column
    ax[0, 3].set_title('EUR Query', fontsize=20, color=colors['EUR'], fontweight='bold')
    for i, database_population in enumerate(SUPER_SUBPOPULATIONS.keys()):
        if database_population == 'EUR':
            # plot subpop scores in black
            # ax[i, 3].hist(EUR_super_scores['same_subpop'],
            #                 bins=50, color=subpop_color, alpha=suppop_alpha)
            # ax[i, 3].hist(EUR_super_scores[database_population],
            #               bins=50, color=colors[database_population], alpha=alpha)
            sns.histplot(EUR_super_scores[database_population], color=colors[database_population],
                         ax=ax[i, 3], alpha=alpha, binwidth=60, linewidth=0.0)
            sns.histplot(EUR_super_scores['same_subpop'], color=subpop_color,
                            ax=ax[i, 3], alpha=suppop_alpha, binwidth=60, fill=False, linewidth=0.3)
        else:
            # ax[i, 3].hist(EUR_super_scores[database_population],
            #               bins=50, color=colors[database_population], alpha=alpha)
            sns.histplot(EUR_super_scores[database_population], color=colors[database_population],
                            ax=ax[i, 3], alpha=alpha, binwidth=60, linewidth=0.0)
        # set x-axis label
        ax[4, i].set_xlabel('GenoSiS Score', fontsize=20, labelpad=20)

    # South Asian query is fifth column
    ax[0, 4].set_title('SAS Query', fontsize=20, color=colors['SAS'], fontweight='bold')
    for i, database_population in enumerate(SUPER_SUBPOPULATIONS.keys()):
        if database_population == 'SAS':
            # plot subpop scores in black
            # ax[i, 4].hist(SAS_super_scores['same_subpop'],
            #                 bins=50, color=subpop_color, alpha=suppop_alpha)
            # ax[i, 4].hist(SAS_super_scores[database_population],
            #               bins=50, color=colors[database_population], alpha=alpha)
            sns.histplot(SAS_super_scores[database_population], color=colors[database_population],
                         ax=ax[i, 4], alpha=alpha, binwidth=60, linewidth=0.0)
            sns.histplot(SAS_super_scores['same_subpop'], color=subpop_color,
                            ax=ax[i, 4], alpha=suppop_alpha, binwidth=60, fill=False, linewidth=0.3)
        else:
            # ax[i, 4].hist(SAS_super_scores[database_population],
            #               bins=50, color=colors[database_population], alpha=alpha)
            sns.histplot(SAS_super_scores[database_population], color=colors[database_population],
                            ax=ax[i, 4], alpha=alpha, binwidth=60, linewidth=0.0)
        # set x-axis label
        ax[4, i].set_xlabel('GenoSiS\nScore', fontsize=20, labelpad=20)


    # formatting for all subplots
    for i in range(5):
        for j in range(5):
            # legend
            if i == j:
                ax[i, j].legend(['SUPER', 'SUB'],
                                loc='upper right', frameon=False)
            # remove spines
            ax[i, j].spines['top'].set_visible(False)
            ax[i, j].spines['right'].set_visible(False)


    # add text at top for figure title
    fig.suptitle('GenoSiS Scores for Isolated\nQuery and Database Populations',
                 fontsize=40, fontweight='bold')

    # move subplots down to make room for title
    plt.subplots_adjust(top=0.85)

    # plt.tight_layout()
    plt.savefig(png)
    plt.close()

def plot_sub_histogram(AFR_sub_scores,
                       AMR_sub_scores,
                       EAS_sub_scores,
                       EUR_sub_scores,
                       SAS_sub_scores,
                       colors,
                       png):
    SUPER_SUBPOPULATIONS = ancestry_helpers.SUPER_SUBPOPULATIONS

    fig, ax = plt.subplots(len(ancestry_helpers.SUBPOPULATIONS), len(ancestry_helpers.SUBPOPULATIONS),
                           dpi=300, figsize=(40, 35), sharex=True, sharey=True)

    alpha = 1

    col_offsets = [len(ancestry_helpers.SUPER_SUBPOPULATIONS[sup]) for sup in ancestry_helpers.SUPER_SUBPOPULATIONS]

    # [row_idx, col_idx]

    # African queries are the first set of columns
    for col_idx, query_subpopulation in enumerate(ancestry_helpers.SUPER_SUBPOPULATIONS['AFR']):
        # column titles: query populations
        ax[0, col_idx].set_title(query_subpopulation,
                            color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[query_subpopulation]],
                            fontsize=15, fontweight='bold', pad=20)
        for row_idx, score_subpop in enumerate(ancestry_helpers.SUBPOPULATIONS):
            # row titles: database populations
            ax[row_idx, 0].set_ylabel(score_subpop, fontsize=15,
                                      color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[score_subpop]],
                                      fontweight='bold', labelpad=20)
            # plot histograms
            ax[row_idx, col_idx].hist(AFR_sub_scores[query_subpopulation][score_subpop], bins=50,
                                        color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[score_subpop]], alpha=alpha)
            # log scale
            ax[row_idx, col_idx].set_yscale('log')

    # American queries are the second set of columns
    for col_idx, query_subpopulation in enumerate(ancestry_helpers.SUPER_SUBPOPULATIONS['AMR']):
        # column titles: query populations
        ax[0, col_idx + col_offsets[0]].set_title(query_subpopulation,
                            color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[query_subpopulation]],
                            fontsize=15, fontweight='bold', pad=20)
        for row_idx, score_subpop in enumerate(ancestry_helpers.SUBPOPULATIONS):
            # plot histograms
            ax[row_idx, col_idx + col_offsets[0]].hist(AMR_sub_scores[query_subpopulation][score_subpop], bins=50,
                                        color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[score_subpop]], alpha=alpha)
            # log scale
            ax[row_idx, col_idx + col_offsets[0]].set_yscale('log')

    # East Asian queries are the third set of columns
    for col_idx, query_subpopulation in enumerate(ancestry_helpers.SUPER_SUBPOPULATIONS['EAS']):
        # column titles: query populations
        ax[0, col_idx + sum(col_offsets[:2])].set_title(query_subpopulation,
                            color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[query_subpopulation]],
                            fontsize=15, fontweight='bold', pad=20)
        for row_idx, score_subpop in enumerate(ancestry_helpers.SUBPOPULATIONS):
            # plot histograms
            ax[row_idx, col_idx + sum(col_offsets[:2])].hist(EAS_sub_scores[query_subpopulation][score_subpop], bins=50,
                                        color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[score_subpop]], alpha=alpha)
            # log scale
            ax[row_idx, col_idx + sum(col_offsets[:2])].set_yscale('log')

    # European queries are the fourth set of columns
    for col_idx, query_subpopulation in enumerate(ancestry_helpers.SUPER_SUBPOPULATIONS['EUR']):
        # column titles: query populations
        ax[0, col_idx + sum(col_offsets[:3])].set_title(query_subpopulation,
                            color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[query_subpopulation]],
                            fontsize=15, fontweight='bold', pad=20)
        for row_idx, score_subpop in enumerate(ancestry_helpers.SUBPOPULATIONS):
            # plot histograms
            ax[row_idx, col_idx + sum(col_offsets[:3])].hist(EUR_sub_scores[query_subpopulation][score_subpop], bins=50,
                                        color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[score_subpop]], alpha=alpha)
            # log scale
            ax[row_idx, col_idx + sum(col_offsets[:3])].set_yscale('log')

    # South Asian queries are the fifth set of columns
    for col_idx, query_subpopulation in enumerate(ancestry_helpers.SUPER_SUBPOPULATIONS['SAS']):
        # column titles: query populations
        ax[0, col_idx + sum(col_offsets[:4])].set_title(query_subpopulation,
                            color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[query_subpopulation]],
                            fontsize=15, fontweight='bold', pad=20)
        for row_idx, score_subpop in enumerate(ancestry_helpers.SUBPOPULATIONS):
            # plot histograms
            ax[row_idx, col_idx + sum(col_offsets[:4])].hist(SAS_sub_scores[query_subpopulation][score_subpop], bins=50,
                                        color=colors[ancestry_helpers.SUB_SUPERPOPULATIONS[score_subpop]], alpha=alpha)
            # log scale
            ax[row_idx, col_idx + sum(col_offsets[:4])].set_yscale('log')


    # Add textbox at top for "AFR Queries", "AMR Queries", etc.
    for i, superpop in enumerate(ancestry_helpers.SUPER_SUBPOPULATIONS.keys()):
        middle = int(col_offsets[i]/2) + sum(col_offsets[:i])
        ax[0, middle].text(0.5, 2., superpop + ' Query',
                            horizontalalignment='center',
                            verticalalignment='center',
                            transform=ax[0, middle].transAxes,
                            fontsize=30, fontweight='bold',
                            color=colors[superpop])

    # Add textbox at left for "AFR Databases", "AMR Databases", etc.
    for i, superpop in enumerate(ancestry_helpers.SUPER_SUBPOPULATIONS.keys()):
        middle = int(col_offsets[i]/2) + sum(col_offsets[:i])
        ax[middle, 0].text(-1.5, 0.5, superpop + '\nDatabase',
                            horizontalalignment='center',
                            verticalalignment='center',
                            transform=ax[middle, 0].transAxes,
                            fontsize=30, fontweight='bold',
                            color=colors[superpop],
                            rotation=90)

    # add x-axis labels
    for i in range(len(ancestry_helpers.SUBPOPULATIONS)):
        for j in range(len(ancestry_helpers.SUBPOPULATIONS)):
            ax[i, j].set_xlabel('GenoSiS\nScore', fontsize=15, labelpad=20)


    # formatting for all subplots
    for i in range(len(ancestry_helpers.SUBPOPULATIONS)):
        for j in range(len(ancestry_helpers.SUBPOPULATIONS)):
            # remove spines
            ax[i, j].spines['top'].set_visible(False)
            ax[i, j].spines['right'].set_visible(False)

    # add text at top for figure title
    fig.suptitle('GenoSiS Scores for Isolated\nQuery and Database Populations\n(by subpopulation)',
                 fontsize=50, fontweight='bold')

    # move subplots down to make room for title
    plt.subplots_adjust(top=0.85)

    # plt.tight_layout()
    plt.savefig(png)
    plt.close()


def get_size_pop(super_pop):
    '''
    Get size of super population
    @param super_pop: super population
    @return: size of super population
    '''
    return {'AFR': 894, 'AMR': 490, 'EAS': 585, 'EUR': 633, 'SAS': 601}[super_pop]

def plot_quality_results_heatmap(AFR_sub_scores,
                                 AMR_sub_scores,
                                 EAS_sub_scores,
                                 EUR_sub_scores,
                                 SAS_sub_scores,
                                 colors,
                                 png):

    SUPER_SUBPOPULATIONS = ancestry_helpers.SUPER_SUBPOPULATIONS
    ALL_SUBPOPULATIONS = ancestry_helpers.SUBPOPULATIONS
    colormap = 'Greys'

    scores_dict = {'AFR': AFR_sub_scores,
                   'AMR': AMR_sub_scores,
                   'EAS': EAS_sub_scores,
                   'EUR': EUR_sub_scores,
                   'SAS': SAS_sub_scores}

    # one big heatmap
    df = pd.DataFrame()
    for super_pop in scores_dict.keys():
        super_pop_scores = scores_dict[super_pop]
        for query_sub_population in SUPER_SUBPOPULATIONS[super_pop]:
            for match_sub_population in ALL_SUBPOPULATIONS:
                scores = super_pop_scores[query_sub_population][match_sub_population]
                size_pop = get_size_pop(super_pop)
                # single value is max of scores
                single_value = max(scores)
                df.at[query_sub_population, match_sub_population] = single_value

    fig, ax = plt.subplots(figsize=(20, 15), dpi=300)
    sns.heatmap(df, cmap=colormap, annot=False, cbar=True, square=True,
            xticklabels=False, yticklabels=False,
            # norm=LogNorm(),
            cbar_kws={'label': 'Max GenoSiS Scores'})

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
        xtick_labels.extend([subpop for subpop in df.index[tick_offset:tick_offset+num_subpops]])
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

    # add text at top for figure title
    fig.suptitle('Max GenoSiS Scores for Isolated\nQuery and Database Populations',
                 fontsize=40, fontweight='bold')

    plt.savefig(png)
    plt.close()

def main():
    args = parse_args()
    ancestry_file = args.ancestry
    k = args.k
    quality_dir = args.quality_dir
    colors = ancestry_helpers.get_colors(args.colors)

    AFR_quality_scores = read_quality_results(quality_dir + 'AFR.txt')
    AMR_quality_scores = read_quality_results(quality_dir + 'AMR.txt')
    EAS_quality_scores = read_quality_results(quality_dir + 'EAS.txt')
    EUR_quality_scores = read_quality_results(quality_dir + 'EUR.txt')
    SAS_quality_scores = read_quality_results(quality_dir + 'SAS.txt')


    AFR_super_scores = combine_subpopulations(AFR_quality_scores)
    AMR_super_scores = combine_subpopulations(AMR_quality_scores)
    EAS_super_scores = combine_subpopulations(EAS_quality_scores)
    EUR_super_scores = combine_subpopulations(EUR_quality_scores)
    SAS_super_scores = combine_subpopulations(SAS_quality_scores)

    plot_super_histogram(AFR_super_scores,
                         AMR_super_scores,
                         EAS_super_scores,
                         EUR_super_scores,
                         SAS_super_scores,
                         colors,
                         args.png + 'super_histogram.png')

    AFR_sub_scores = organize_subpopulations(AFR_quality_scores)
    AMR_sub_scores = organize_subpopulations(AMR_quality_scores)
    EAS_sub_scores = organize_subpopulations(EAS_quality_scores)
    EUR_sub_scores = organize_subpopulations(EUR_quality_scores)
    SAS_sub_scores = organize_subpopulations(SAS_quality_scores)
    #
    # plot_sub_histogram(AFR_sub_scores,
    #                    AMR_sub_scores,
    #                    EAS_sub_scores,
    #                    EUR_sub_scores,
    #                    SAS_sub_scores,
    #                    colors,
    #                    args.png + 'sub_histogram.png')

    plot_quality_results_heatmap(AFR_sub_scores,
                                 AMR_sub_scores,
                                 EAS_sub_scores,
                                 EUR_sub_scores,
                                 SAS_sub_scores,
                                 colors,
                                 args.png + 'heatmap.png')


if __name__ == '__main__':
    main()