import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import random

from plotting import ancestry_helpers
from src import get_relations

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--trios', type=str, help='top hits file', required=True)
    parser.add_argument('--pop', type=str, help='population query', required=True)
    parser.add_argument('--color', type=str, help='color of plots', required=True, default='olivedrab')

    return parser.parse_args()


def read_trios(trios_file):
    trios_dict = {}

    with open(trios_file, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip().split()
            label = line[0]
            scores = [float(x) for x in line[1:]]
            trios_dict[label] = scores

    return trios_dict

def normalize_pop_scores(scores, pop_min, pop_max):
    # normalize scores: x - pop_min / pop_max(x) - pop_min(x)
    scores = [(x - pop_min) / (pop_max - pop_min) for x in scores]
    return scores

def normalize_scores(scores):
    # normalize scores: x - min() / range()
    scores = [(x - min(scores)) / (max(scores) - min(scores)) for x in scores]
    return scores

def plot_single_figure(trio_scores, color, pop):
    color_dict = {'self': color,
                  'parent/child': color,
                  'subpop': color,
                  'AFR': 'darkorange',
                  'AMR': 'mediumpurple',
                  'EAS': 'hotpink',
                  'EUR': 'steelblue',
                  'SAS': 'goldenrod'}

    # plot a distribution of scores for each category, on a different axes
    categories_labels = {'self': ['self'],
                         'parent/child': ['parent', 'child'],
                         'subpop': ['subpop'],
                         'AFR': ['AFR'],
                         'AMR': ['AMR'],
                         'EAS': ['EAS'],
                         'EUR': ['EUR'],
                         'SAS': ['SAS']}
    categories = ['self', 'parent/child', 'subpop', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']

    # combine scores for all categories
    all_scores = []
    for label in categories:
        for type in categories_labels[label]:
            all_scores.extend(trio_scores[type])
    all_min = min(all_scores)
    all_max = max(all_scores)

    # combined scores for AFR, AMR, EAS, EUR, SAS
    all_pop_names = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    pop_scores = []
    for label in all_pop_names:
        pop_scores.extend(trio_scores[label])
    pop_min = min(pop_scores)
    pop_max = max(pop_scores)

    plt.figure(figsize=(10, 5))
    plt.title('GenoSiS scores for ' + pop + ' trios', fontsize=20)

    for i, category in enumerate(categories):
        scores = []
        for label in categories_labels[category]:
            scores.extend(trio_scores[label])
        if category == 'AFR' or category == 'AMR' or category == 'EAS' or category == 'EUR' or category == 'SAS':
            sns.kdeplot(scores, color=color_dict[category], linewidth=2, fill=True, alpha=0.6)
        else:
            sns.kdeplot(scores, color=color_dict[category], linewidth=2, fill=True, alpha=0.6)

        # remove spines
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.x_label = 'GenoSiS Score'

    # sns.kdeplot(all_scores, color='black', linewidth=2, fill=True, alpha=0.6)
    plt.tight_layout()
    plt.savefig('1kg_trio_plots/' + pop + '_trio_distributions_ONE.png')


def plot_by_category(trio_scores, color, pop):

    color_dict = {'self': color,
                    'parent/child': color,
                    'subpop': color,
                    'AFR': 'darkorange',
                    'AMR': 'mediumpurple',
                    'EAS': 'hotpink',
                    'EUR': 'steelblue',
                    'SAS': 'goldenrod'}

    # plot a distribution of scores for each category, on a different axes
    categories_labels = {'self': ['self'],
                         'parent/child': ['parent', 'child'],
                         'subpop': ['subpop'],
                         'AFR': ['AFR'],
                         'AMR': ['AMR'],
                         'EAS': ['EAS'],
                         'EUR': ['EUR'],
                         'SAS': ['SAS']}
    categories = ['self', 'parent/child', 'subpop', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']

    fig, axs = plt.subplots(len(categories), 1, figsize=(10, 15))
    fig.suptitle('GenoSiS scores for ' + pop + ' trios', fontsize=20)

    # combined scores for AFR, AMR, EAS, EUR, SAS
    all_pop_names = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    pop_scores = []
    for label in all_pop_names:
        pop_scores.extend(trio_scores[label])
    pop_min = min(pop_scores)
    pop_max = max(pop_scores)

    # combine scores for all categories
    all_scores = []
    for label in categories:
        for type in categories_labels[label]:
            all_scores.extend(trio_scores[type])
    all_min = min(all_scores)
    all_max = max(all_scores)

    for i, category in enumerate(categories):
        scores = []
        for label in categories_labels[category]:
            scores.extend(trio_scores[label])
        # plot density plot
        sns.kdeplot(scores, color=color_dict[category], ax=axs[i], linewidth=2, fill=True, alpha=0.6)
        # axs[i].hist(scores, bins=20, color=color)
        # axs[i].hist(scores, bins=100, color=color, density=True, alpha=0.6, histtype='stepfilled', linewidth=1.5)

        axs[i].set_title(category, fontsize=14, loc='left')
        axs[i].set_ylabel('Density')


        if category == 'AFR' or category == 'AMR' or category == 'EAS' or category == 'EUR' or category == 'SAS':
            sns.kdeplot(pop_scores, color='black', ax=axs[i], linewidth=2, fill=True, alpha=0.6)
            # axs[i].set_ylim([0, 1])
        # else:
        #     axs[i].set_ylim([0, 1])

        # remove spines
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)

        # # remove x-axis labels for all but the bottom plot
        # if i != len(categories) - 1:
        #     axs[i].set_xticklabels([])
        #     axs[i].set_xticks([])
        #     axs[i].spines['bottom'].set_visible(False)
        # else:
        #     axs[i].set_xlabel('Score')

    axs[i].set_xlabel('GenoSiS Score')
    # get the maximum x value for all plots
    max_x = max([ax.get_xlim()[1] for ax in axs])
    # make all y-axes and x-axes the same
    for ax in axs:
        # ax.set_ylim([0, 100])
        if max_x > 1: max_x + 50
        ax.set_xlim([0, max_x])


    plt.tight_layout()
    plt.savefig('1kg_trio_plots/' + pop + '_trio_distributions.png')

def plot_simple(trio_scores, color, pop):
    color_dict = {'self': color,
                  'parent/child': color,
                  'inpop': color,
                  'outpop': 'black'}

    outpop_categories = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    # remove pop from outpop_categories
    outpop_categories.remove(pop)
    # plot a distribution of scores for each category, on a different axes
    categories_labels = {'self': ['self'],
                         'parent/child': ['parent', 'child'],
                         'inpop': ['subpop', pop],
                         'outpop': outpop_categories}
    categories = ['self', 'parent/child', 'inpop', 'outpop']

    fig, axs = plt.subplots(len(categories), 1, figsize=(7, 5), sharex=True, sharey=True, dpi=150)
    fig.suptitle('GenoSiS scores for ' + pop + ' trios', fontsize=20)

    for i, category in enumerate(categories):
        scores = []
        for label in categories_labels[category]:
            scores.extend(trio_scores[label])
        # plot density plot
        sns.kdeplot(scores, color=color_dict[category], ax=axs[i], linewidth=2, fill=True, alpha=0.6)
        # axs[i].hist(scores, bins=20, color=color)
        # axs[i].hist(scores, bins=100, color=color, density=True, alpha=0.6, histtype='stepfilled', linewidth=1.5)

        axs[i].set_title(category, fontsize=14, loc='left')
        axs[i].set_ylabel('Density')

        # remove spines
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)

        # # remove x-axis labels for all but the bottom plot
        # if i != len(categories) - 1:
        #     axs[i].set_xticklabels([])
        #     axs[i].set_xticks([])
        #     axs[i].spines['bottom'].set_visible(False)
        # else:
        #     axs[i].set_xlabel('Score')

    axs[i].set_xlabel('GenoSiS Score')
    # get the maximum x value for all plots
    max_x = max([ax.get_xlim()[1] for ax in axs])
    # make all y-axes and x-axes the same
    # for ax in axs:
    #     # ax.set_ylim([0, 100])
    #     if max_x > 1: max_x + 50
    #     ax.set_xlim([0, max_x])

    plt.tight_layout()
    plt.savefig('1kg_trio_plots/' + pop + '_trio_distributions_RANK.png')


def plot_all_populations(trio_scores):
    # plot family members in one category and all other populations in another
    color_dict = {'in_family': 'salmon',
                  'out_family': 'black'}

    # plot a distribution of scores for each category on the same figure
    categories_labels = {'in_family': ['parent', 'child'],
                         'out_family': ['subpop', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']}

    fig, axs = plt.subplots(1, 1, figsize=(7, 5), dpi=150)
    fig.suptitle('GenoSiS scores for trios in all 1KG populations', fontsize=20)

    for i, category in enumerate(categories_labels):
        scores = []
        for label in categories_labels[category]:
            scores.extend(trio_scores[label])
        # plot density plot
        sns.kdeplot(scores, color=color_dict[category], ax=axs, linewidth=2, fill=True, alpha=0.6)

    # remove spines
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.set_ylabel('Density')
    axs.set_xlabel('GenoSiS Score')
    axs.legend(['Parent/Child', 'Unrelated'], loc='upper right')

    plt.tight_layout()
    plt.savefig('1kg_trio_plots/ALL_trio_distributions.png')


def combine_population_scores():
    all_populations = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    data_dir = '1kg_trio_data/'
    file_prefix = 'GenoSiS_'
    file_suffix = '_trio_scores_20.txt'

    pop_scores = {'self': [], 'parent': [], 'child': [], 'subpop': [], 'AFR': [], 'AMR': [], 'EAS': [], 'EUR': [], 'SAS': []}

    for pop in all_populations:
        # open file for population
        with open(data_dir + file_prefix + pop + file_suffix, 'r') as f:
            header = f.readline()
            for line in f:
                line = line.strip().split()
                label = line[0]
                scores = [float(x) for x in line[1:]]
                pop_scores[label].extend(scores)

    return pop_scores

def main():

    args = get_args()
    trios = read_trios(args.trios)
    # plot_by_category(trios, args.color, args.pop)
    # plot_single_figure(trios, args.color, args.pop)
    # plot_simple(trios, args.color, args.pop)

    # get data for all populations
    all_population_trios = combine_population_scores()
    plot_all_populations(all_population_trios)


if __name__ == '__main__':
    main()