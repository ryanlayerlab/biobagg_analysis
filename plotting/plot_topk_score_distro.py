import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sklearn.decomposition import PCA
import utils
from scipy.stats import gaussian_kde
import argparse
import seaborn as sns

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--topk_file', type=str, required=True)
    parser.add_argument('--label_file', type=str, required=True)
    parser.add_argument('--out_file', type=str, required=True)
    parser.add_argument('--width', type=float, default=10.0)
    parser.add_argument('--height', type=float, default=5.0)
    parser.add_argument('--bins', type=int, default=50)
    parser.add_argument('--log', action='store_true')
    return parser.parse_args()


def main():
    args = get_args()

    top_k = utils.get_top_hits(args.topk_file, get_scores=True)

    pop_map = utils.get_label_map(args.label_file, 'Population code')
    spop_map = utils.get_label_map(args.label_file, 'Superpopulation code')

    D = {}

    for i in top_k:
        spop = spop_map[i]
        if spop not in D:
            D[spop] = []

        for j, value in top_k[i]:
            if i == j: continue
            D[spop].append(value)

    spops = ['EUR', 'EAS', 'AMR', 'SAS', 'AFR']

    fig, axs = plt.subplots(1, len(spops), figsize=(args.width, args.height))

    for i in range(len(spops)):
        spop = spops[i]
        axs[i].set_title(spop, fontsize=8, loc='left')
        axs[i].hist(D[spop], bins=args.bins, density=True)

    max_y = max([ax.get_ylim()[1] for ax in axs])
    min_y = min([ax.get_ylim()[0] for ax in axs])

    max_x = max([ax.get_xlim()[1] for ax in axs])
    min_x = min([ax.get_xlim()[0] for ax in axs])

    for ax in axs:
        #ax.set_ylim((0,30))
        ax.set_xlim((min_x,max_x))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both',
                       which='major',
                       labelsize=4,
                       width=0.25,
                       length=2)
        ax.spines['bottom'].set_linewidth(0.25)
        ax.spines['left'].set_linewidth(0.25)
        if args.log: ax.set_yscale('log')

    axs[0].set_ylabel('Freq.', fontsize=6)

    fig.text(0.5, 0.04, 'Similarity score', ha='center', fontsize=6)


        
    plt.tight_layout()
    plt.savefig(args.out_file, dpi=300)

if __name__ == '__main__':
    main()
