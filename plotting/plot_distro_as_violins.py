import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sklearn.decomposition import PCA
import utils
from scipy.stats import gaussian_kde
import argparse
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--topk_file', type=str, required=True)
    parser.add_argument('--pairs_file', type=str, required=True)
    parser.add_argument('--label_file', type=str, required=True)
    parser.add_argument('--out_file', type=str, required=True)
    parser.add_argument('--width', type=float, default=10.0)
    parser.add_argument('--height', type=float, default=5.0)
    parser.add_argument('--y_label', type=str)
    return parser.parse_args()

def add_outlier(O, spop, group, value, i, j, outlier):
    if outlier is not None:
        if value >= outlier:
            O[spop][group].append((i,j))

def summarize_related(O, spop, group, related):
    are_related = 0
    not_related = 0
    for p1, p2 in O[spop][group]:
        if p1 in related and  p2 in related[p1]:
            are_related += 1
        else:
            not_related += 1
    return are_related, not_related

def main():
    args = get_args()

    top_k = utils.get_top_hits(args.topk_file)

    pairs = utils.get_pair_map(args.pairs_file)

    pop_map = utils.get_label_map(args.label_file, 'Population code')
    spop_map = utils.get_label_map(args.label_file, 'Superpopulation code')

    D = {}

    for i in pairs:
        spop = spop_map[i]
        if spop not in D:
            D[spop] = {}
            D[spop]['pop'] = []
            D[spop]['spop'] = []
            D[spop]['topk'] = []

        for j in pairs[i]:
            if i == j : continue
            if pop_map[i] == pop_map[j]:
                D[spop]['pop'].append(pairs[i][j])
            if spop_map[i] == spop_map[j]:
                D[spop]['spop'].append(pairs[i][j])

        for j in top_k[i]:
            if i == j: continue
            if j not in pairs[i]: continue
            D[spop]['topk'].append(pairs[i][j])


    spops = ['EUR', 'EAS', 'AMR', 'SAS', 'AFR']

    super_population_colors = {'AFR': '#fe9d57',
                               'AMR': '#ba75ff',
                               'EAS': '#eb4690',
                               'EUR': '#728eff',
                               'SAS': '#f7d06b'}

    fig = plt.figure(figsize=(args.width, args.height))

    gs = gridspec.GridSpec(1, len(spops), hspace=0.6)

    axs = []
    for i in range(len(spops)):
        ax = fig.add_subplot(gs[i])
        axs.append(ax)
        spop = spops[i]
        ax.set_title(f"{spop}",
                     fontsize=4,
                     loc='left',
                     fontdict={'verticalalignment':'top'})

        violin_data = [D[spop]['spop'],
                       D[spop]['pop'],
                       D[spop]['topk']]

        parts = ax.violinplot(violin_data, showextrema=False, showmeans=True)

        parts['cmeans'].set_color(super_population_colors[spop])
        parts['cmeans'].set_linewidth(0.75)

        for b in parts['bodies']:
            b.set_color(super_population_colors[spop])
            b.set_alpha(0.5)
            b.set_linewidth(0.5)

    max_y = max([ax.get_ylim()[1] for ax in axs])
    min_y = min([ax.get_ylim()[0] for ax in axs])
    for ax in axs:
        ax.set_ylim((min_y,max_y))
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks(range(1, 4))
        ax.set_xticklabels(['Super Pop.', 'Pop.', 'Top K'])

        if ax != axs[0]:
            ax.set_yticks([])
        else:
            ax.set_ylabel(args.y_label, fontsize=4)

        ax.tick_params(axis='both',
                       which='major',
                       labelsize=3,
                       width=0.25,
                       length=1)

    plt.tight_layout()
    plt.savefig(args.out_file, dpi=300)

if __name__ == '__main__':
    main()
