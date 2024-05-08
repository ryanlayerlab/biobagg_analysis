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
    parser.add_argument('--hits', type=str, help='top hits file', required=True)
    parser.add_argument('--trios', type=str, help='trios scores', required=True)
    parser.add_argument('--ancestry', type=str, help='population query', required=True)
    parser.add_argument('--ped', type=str, help='pedigree file', required=True)
    parser.add_argument('--out', type=str, help='output file', required=True)

    return parser.parse_args()

def get_trio_scores(trio_data_prefix):
    all_populations = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    file_suffix = '_trio_scores_20.txt'
    pop_scores = {'self': [], 'parent': [], 'child': [], 'subpop': [], 'AFR': [], 'AMR': [], 'EAS': [], 'EUR': [], 'SAS': []}

    for pop in all_populations:
        # open file for population
        with open(trio_data_prefix + pop + file_suffix, 'r') as f:
            header = f.readline()
            for line in f:
                line = line.strip().split()
                label = line[0]
                scores = [float(x) for x in line[1:]]
                pop_scores[label].extend(scores)

    return pop_scores

def plot_trio_all_populations(trio_scores, out_dir):
    # plot family members in one category and all other populations in another
    color_dict = {'in_family': 'red',
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
    plt.savefig(out_dir + 'ALL_trio_distributions.png')


def main():
    args = get_args()
    trio_data_prefix = args.trios
    hits_file = args.hits
    ancestry = args.ancestry
    ped_file = args.ped
    out_dir = args.out

    # plot distribution for trio data, all populations
    trio_scores = get_trio_scores(trio_data_prefix)
    plot_trio_all_populations(trio_scores, out_dir)

    # plot distribution for population data, no trios

if __name__ == '__main__':
    main()