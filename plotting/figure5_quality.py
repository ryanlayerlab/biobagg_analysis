import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import os
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
    parser.add_argument('--png_hist', help='Output png firl for histogram', required=True)

    return parser.parse_args()

def read_quality_results(query_pop_file):
    '''
    Read quality results from directory
    @param quality_dir: path to directory of quality results
    @return: dictionary of quality results
    '''
    superpopulations = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    quality_scores = defaultdict(dict)
    f = open(query_pop_file, 'r')
    header = f.readline()
    query_population = header.strip().split(':')[1].strip()
    for line in f:
        if 'database' in line:
            # new database
            database_population = line.strip().split(':')[1].strip()
            continue
        database_subpopulation = line.strip().split(':')[0].strip()
        database_quality_scores = line.strip().split(':')[1].split(',')
        # remove last element which is empty
        database_quality_scores = database_quality_scores[:-1]
        database_quality_scores = [float(score) for score in database_quality_scores]
        quality_scores[database_population][database_subpopulation] = database_quality_scores
    f.close()
    return quality_scores

def plot_quality_results_heatmap(AFR_quality_scores,
                                   AMR_quality_scores,
                                   EAS_quality,scores,
                                   EUR_quality_scores,
                                   SAS_quality_scores,
                                   png):

    SUPER_SUBPOPULATIONS = ancestry_helpers.SUPER_SUBPOPULATIONS


def plot_quality_results_histogram(AFR_quality_scores,
                                   AMR_quality_scores,
                                   EAS_quality_scores,
                                   EUR_quality_scores,
                                   SAS_quality_scores,
                                   colors,
                                   png):
    SUPER_SUBPOPULATIONS = ancestry_helpers.SUPER_SUBPOPULATIONS

    fig, ax = plt.subplots(5, 5,
                           figsize=(20, 15), dpi=300,
                           sharex=True, sharey=True)

    # African query is first column
    ax[0, 0].set_title('AFR Query', fontsize=20)
    for i, database_population in enumerate(AFR_quality_scores.keys()):
        database_scores = []
        for j, database_subpopulation in enumerate(AFR_quality_scores[database_population].keys()):
            database_scores.extend(AFR_quality_scores[database_population][database_subpopulation])
        ax[i, 0].hist(database_scores, bins=50, color=colors[database_population], alpha=0.5)
        # set y-axis label
        ax[i, 0].set_ylabel(database_population + '\nDatabase',
                            fontsize=25, color=colors[database_population],
                            fontweight='bold', labelpad=40)
        # set x-axis label
        ax[4, i].set_xlabel('GenoSiS Score', fontsize=20, labelpad=20)
        # ignore outliers
        # ax[i, 1].set_xlim(0, 1000)
        # log scale
        ax[i, 1].set_yscale('log')


    # American query is second column
    ax[0, 1].set_title('AMR Query', fontsize=20)
    for i, database_population in enumerate(AMR_quality_scores.keys()):
        database_scores = []
        for j, database_subpopulation in enumerate(AMR_quality_scores[database_population].keys()):
            database_scores.extend(AMR_quality_scores[database_population][database_subpopulation])
        ax[i, 1].hist(database_scores, bins=50, color=colors[database_population], alpha=0.5)
        # set x-axis label
        ax[4, i].set_xlabel('GenoSiS Score', fontsize=20, labelpad=20)


    # # East Asian query is third column
    ax[0, 2].set_title('EAS Query', fontsize=20)
    for i, database_population in enumerate(EAS_quality_scores.keys()):
        database_scores = []
        for j, database_subpopulation in enumerate(EAS_quality_scores[database_population].keys()):
            database_scores.extend(EAS_quality_scores[database_population][database_subpopulation])
        ax[i, 2].hist(database_scores, bins=50, color=colors[database_population], alpha=0.5)
        # set x-axis label
        ax[4, i].set_xlabel('GenoSiS Score', fontsize=20, labelpad=20)

    # European query is fourth column
    ax[0, 3].set_title('EUR Query', fontsize=20)
    for i, database_population in enumerate(EUR_quality_scores.keys()):
        database_scores = []
        for j, database_subpopulation in enumerate(EUR_quality_scores[database_population].keys()):
            database_scores.extend(EUR_quality_scores[database_population][database_subpopulation])
        ax[i, 3].hist(database_scores, bins=50, color=colors[database_population], alpha=0.5)
        # set x-axis label
        ax[4, i].set_xlabel('GenoSiS Score', fontsize=20, labelpad=20)

    # South Asian query is fifth column
    ax[0, 4].set_title('SAS Query', fontsize=20)
    for i, database_population in enumerate(SAS_quality_scores.keys()):
        database_scores = []
        for j, database_subpopulation in enumerate(SAS_quality_scores[database_population].keys()):
            database_scores.extend(SAS_quality_scores[database_population][database_subpopulation])
        ax[i, 4].hist(database_scores, bins=50, color=colors[database_population], alpha=0.5)
        # set x-axis label
        ax[4, i].set_xlabel('GenoSiS\nScore', fontsize=20, labelpad=20)


    # formatting for all subplots
    for i in range(5):
        for j in range(5):
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


    # plot_quality_results_heatmap(AFR_quality_scores,
    #                              AMR_quality_scores,
    #                              EAS_quality_scores,
    #                              EUR_quality_scores,
    #                              SAS_quality_scores,
    #                              colors,
    #                              out_png + '_heatmap.png')


    plot_quality_results_histogram(AFR_quality_scores,
                                   AMR_quality_scores,
                                   EAS_quality_scores,
                                   EUR_quality_scores,
                                   SAS_quality_scores,
                                   colors,
                                   args.png_hist)





if __name__ == '__main__':
    main()