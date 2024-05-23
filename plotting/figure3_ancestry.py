import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import ks_2samp
from scipy.stats import norm

def parse_args():
    parser = argparse.ArgumentParser()
    # GenoSiS Scores
    parser.add_argument('--genosis', help='1KG genosis scores', required=True)
    # plink Scores
    parser.add_argument('--dst', help='1KG plink DST scores', required=True)
    parser.add_argument('--pihat', help='1KG plink pi-hat scores', required=True)
    parser.add_argument('--kinship', help='1KG plink kinship scores', required=True)
    # Output
    parser.add_argument('--png', help='Output png file', required=True)

    return parser.parse_args()

def get_ancestry_colors():
    return {
        'AFR': 'darkorange',
        'AMR': 'mediumpurple',
        'EAS': 'deeppink',
        'EUR': 'steelblue',
        'SAS': 'goldenrod'
    }

def get_category_colors(superpop):
    return {
        'superpop': 'gray',
        'subpop': get_ancestry_colors()[superpop],
        'outgroup': 'black'
    }

def read_ancestry_group_scores(ancestry_group_scores_file):
    '''
    Read scores for superpop, subpop, and outpop labels.
    @param ancestry_group_scores_file: path to scores file for 1KG data
    @return: dictionary of scores for each superpop
    '''
    # superpop: category: [scores...]
    scores = defaultdict(dict)
    with open(ancestry_group_scores_file, 'r') as f:
        header = f.readline().strip().split(',')
        for line in f:
            line = line.strip().split(',')
            superpop = line[0]
            category = line[1]
            scores[superpop][category] = [float(x) for x in line[2].strip().split()]
    return scores

def plot_ancestry_group_scores(genosis_scores,
                               dst_scores,
                               pihat_scores,
                               kinship_scores,
                               png_file):

    # Create the combined figure
    combined_figure, axes = plt.subplots(5, 4, figsize=(20, 10), dpi=300)
    alpha_value = 0.5
    genosis_x_max = 1500

    # Plot the GenoSiS scores in first column
    for i, superpop in enumerate(genosis_scores.keys()):
        D = []
        colors = []
        categories = genosis_scores[superpop].keys()
        for j, category in enumerate(genosis_scores[superpop].keys()):
            D.append(genosis_scores[superpop][category])
            colors.append(get_category_colors(superpop)[category])
        sns.kdeplot(D,
                    ax=axes[i, 0],
                    label=superpop,
                    palette=colors,
                    fill=True, alpha=alpha_value)
        axes[i, 0].set_xlabel('GenoSiS Score')
        axes[i, 0].set_ylabel('Density')
        axes[i, 0].set_xlim(0, genosis_x_max)
        axes[i, 0].spines['top'].set_visible(False)
        axes[i, 0].spines['right'].set_visible(False)
        legend = ['Subpopulation', 'Superpopulation', 'Outgroup']
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors[i], ec="k") for i in range(len(legend))]
        axes[i, 0].legend(handles, legend, loc='upper right', frameon=False)
    # Save the figure
    plt.tight_layout()
    combined_figure.savefig(png_file)


def main():
    args = parse_args()


    genosis_scores = read_ancestry_group_scores(args.genosis)
    dst_scores = read_ancestry_group_scores(args.dst)
    pihat_scores = read_ancestry_group_scores(args.pihat)
    kinship_scores = read_ancestry_group_scores(args.kinship)
    png_file = args.png

    plot_ancestry_group_scores(genosis_scores,
                                 dst_scores,
                                 pihat_scores,
                                 kinship_scores,
                                 png_file)


if __name__ == '__main__':
    main()

