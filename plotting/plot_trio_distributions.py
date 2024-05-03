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

def plot_by_category(trio_scores, color, pop):
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

    fig, axs = plt.subplots(len(categories), 1, figsize=(15, 25))
    fig.suptitle('Distribution of scores for each labeled relationship')
    for i, category in enumerate(categories):
        scores = []
        for label in categories_labels[category]:
            scores.extend(trio_scores[label])
        # plot density plot
        sns.kdeplot(scores, color=color, ax=axs[i], linewidth=2, fill=True, alpha=0.6)
        # kdeplot(np.array(data), bw=0.5)
        # axs[i].hist(scores, bins=20, color=color)
        # axs[i].hist(scores, bins=100, color=color, density=True, alpha=0.6, histtype='stepfilled', linewidth=1.5)


        axs[i].set_title(category)
        axs[i].set_ylabel('Frequency')

        # remove spines
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)

        # remove x-axis labels for all but the bottom plot
        if i != len(categories) - 1:
            axs[i].set_xticklabels([])
            axs[i].set_xticks([])
            axs[i].spines['bottom'].set_visible(False)
        else:
            axs[i].set_xlabel('Score')


    # make all y-axes and x-axes the same
    for ax in axs:
        # ax.set_ylim([0, 100])
        ax.set_xlim([0, 9000])

    plt.tight_layout
    plt.savefig('1kg_trio_plots/' + pop + '_trio_distributions.png')


def main():
    args = get_args()
    trios = read_trios(args.trios)
    plot_by_category(trios, args.color, args.pop)



if __name__ == '__main__':
    main()