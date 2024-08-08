import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import random

from plotting import ancestry_helpers
from src import get_relations

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ccpm_fst', type=str, help='ccpm fst file', required=True)
    parser.add_argument('--out_dir', type=str, help='output png', required=True)

    return parser.parse_args()

def read_ccpm_fst(ccpm_fst):
    all_fst = defaultdict(dict)
    with open(ccpm_fst, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip().split()
            pop1 = line[0]
            pop2 = line[1]
            fst = float(line[2])
            all_fst[pop1][pop2] = fst
            all_fst[pop2][pop1] = fst

    # add 0s for the diagonal
    for pop in all_fst.keys():
        all_fst[pop][pop] = 0

    return all_fst

def plot_fst(all_fst,
             out_dir,
             color_CCPM):
    '''
    Heatmpat of FST values

    @param all_fst:
    @param out_dir:
    @return:
    '''

    ordered_ancestry = ['Africa', 'America', 'East', 'Europe', 'Middle', 'Central']
    ordered_ancestry_labels = ['TGP+HGP-\nAFR-like',
                               'TGP+HGP-\nAMR-like',
                               'TGP+HGP-\nEAS-like',
                               'TGP+HGP-\nEUR-like',
                               'TGP+HGP-\nMLE-like',
                               'TGP+HGP-\nSAS-like']

    # plot the heatmap
    fig, ax = plt.subplots(figsize=(18, 15), dpi=300)

    fst_values = []
    for pop1 in ordered_ancestry:
        row = []
        for pop2 in ordered_ancestry:
            row.append(all_fst[pop1][pop2])
        fst_values.append(row)

    # plot heatmap with gray
    sns.heatmap(fst_values, cmap='gray', ax=ax,
                square = True,
                annot = True, fmt = '.4f', annot_kws = {'size': 20},
                vmin = 0, vmax = .115)
    sns.set(font_scale=1.5)

    cbar = ax.collections[0].colorbar
    # cbar.set_ticks([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12])
    # set font size of colorbar
    cbar.ax.tick_params(labelsize=15)
    cbar.set_label('FST Value', fontsize=25, labelpad=20)

    ax.set_xticklabels(ordered_ancestry_labels, fontsize=20)
    ax.set_yticklabels(ordered_ancestry_labels, fontsize=20)

    # add rectangles around the diagonal colors
    for i, a in enumerate(ordered_ancestry):
        ax.add_patch(Rectangle((i, i), 1, 1,
                               fill=False, edgecolor=color_CCPM[a], lw=8))


    plt.tight_layout()
    plt.savefig(out_dir + 'ccpm_fst.png')
    plt.close()



def main():
    args = get_args()
    ccpm_fst = args.ccpm_fst
    out_dir = args.out_dir

    ancestry_names = {'Africa': 'Africa',
                      'America': 'America',
                      'East': 'East Asia',
                      'Europe': 'Europe',
                      'Middle': 'Middle East',
                      'Central': 'Central South Asian'}

    # add rectangles around the diagonal colors
    color_CCPM = {'Africa': 'deepskyblue',
                  'America': 'gold',
                  'East': 'crimson',
                  'Europe': 'yellowgreen',
                  'Middle': 'chocolate',
                  'Central': 'mediumpurple'}

    all_fst = read_ccpm_fst(ccpm_fst)

    plot_fst(all_fst,
                out_dir,
                color_CCPM)





if __name__ == '__main__':
    main()