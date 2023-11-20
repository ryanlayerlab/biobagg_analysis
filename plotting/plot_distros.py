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
    parser.add_argument('--ped_file', type=str)
    parser.add_argument('--out_file', type=str, required=True)
    parser.add_argument('--width', type=float, default=10.0)
    parser.add_argument('--height', type=float, default=5.0)
    parser.add_argument('--x_label', type=str)
    parser.add_argument('--alpha', type=float, default=0.25)
    parser.add_argument('--bins', type=int, default=20)
    parser.add_argument('--outliers', type=float)
    return parser.parse_args()

def add_outlier(O, spop, group, value, i, j, outlier):
    if outlier is not None:
        if value >= outlier:
            O[spop][group].append((i,j))

def get_related_map(ped_file):
    if ped_file is None: return None

    related = {}
    
    with open(ped_file) as lines:
        for line in lines:
            sid, fid, mid, sex = line.rstrip().split()
            if fid != '0':
                if sid not in related:
                    related[sid] = {}
                if fid not in related:
                    related[fid] = {}
                related[sid][fid] = 1
                related[fid][sid] = 1
            if mid != '0':
                if sid not in related:
                    related[sid] = {}
                if mid not in related:
                    related[mid] = {}
                related[sid][mid] = 1
                related[mid][sid] = 1
    return related

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

    related = get_related_map(args.ped_file)

    D = {}
    O = {}

    for i in pairs:
        spop = spop_map[i]
        if spop not in D:
            D[spop] = {}
            D[spop]['pop'] = []
            D[spop]['spop'] = []
            D[spop]['topk'] = []

            O[spop] = {}
            O[spop]['pop'] = []
            O[spop]['spop'] = []
            O[spop]['topk'] = []

        for j in pairs[i]:
            if i == j : continue
            if spop_map[i] == spop_map[j]:
                D[spop]['spop'].append(pairs[i][j])
                add_outlier(O, spop, 'spop', pairs[i][j], i, j, args.outliers)
            if pop_map[i] == pop_map[j]:
                D[spop]['pop'].append(pairs[i][j])
                add_outlier(O, spop, 'pop', pairs[i][j], i, j, args.outliers)


        for j in top_k[i]:
            if i == j: continue
            if j not in pairs[i]: continue
            D[spop]['topk'].append(pairs[i][j])
            add_outlier(O, spop, 'topk', pairs[i][j], i, j, args.outliers)


    spops = ['EUR', 'EAS', 'AMR', 'SAS', 'AFR']

    fig = plt.figure(figsize=(args.width, args.height))
    outer_gs = gridspec.GridSpec(len(spops),
                                 1, 
                                 hspace=0.6)

    ax_i = 0
    axs = []
    ab_color = 'black'
    ab_lw = 0.5
    legend_set = False
    for i in range(len(spops)):
        spop = spops[i]

        inner_axs = []

        inner_gs = gridspec.GridSpecFromSubplotSpec(
                3, 
                1,
                subplot_spec=outer_gs[ax_i],
                height_ratios=[1, 1, 1],
                hspace=0.4)

        title_ax = fig.add_subplot(inner_gs[0, :])
        title_ax.set_title(f"{spop}",
                           fontsize=4,
                           loc='left',
                           fontdict={'verticalalignment':'top'})
        title_ax.axis('off')

        ax = fig.add_subplot(inner_gs[0])
        axs.append(ax)
        inner_axs.append(ax)

        ax.hist(D[spop]['spop'],
                density = True,
                bins=args.bins,
                alpha=args.alpha,
                label='Super Pop.',
                color='C0')
        ax.axvline(x=np.mean(D[spop]['spop']),
                   lw=ab_lw,
                   color=ab_color,
                   label='Mean')

        if not legend_set :
            hist_handle0 = mpatches.Patch(color='C0', label='Super pop.')
            hist_handle1 = mpatches.Patch(color='C1', label='Pop.')
            hist_handle2 = mpatches.Patch(color='C2', label='Top K')
            mean_handle = mlines.Line2D([],
                                        [],
                                        lw=ab_lw,
                                        color='black',
                                        label='Mean')
            ax.legend(handles=[hist_handle0,
                               hist_handle1,
                               hist_handle2,
                               mean_handle],
                      loc='right',
                      ncol=4, 
                      fontsize=3,
                      handlelength=0.25,
                      handletextpad=0.2,
                      columnspacing=0.5,
                      frameon=False)
            legend_set = True


        print(spop,
              'spop',
              np.mean(D[spop]['spop']),
              np.median(D[spop]['spop']))
        if args.outliers and related is not None:
            are_related, not_related = summarize_related(O,
                                                         spop,
                                                         'spop',
                                                         related)
            print(spop, 'spop',
                  args.outliers,
                  are_related, 'are related',
                  not_related, 'not related')

        ax = fig.add_subplot(inner_gs[1])
        axs.append(ax)
        inner_axs.append(ax)
    
        ax.hist(D[spop]['pop'],
                density = True,
                bins=args.bins,
                alpha=args.alpha,
                label='Pop.',
                color='C1')
        ax.axvline(x=np.mean(D[spop]['pop']), lw=ab_lw, color=ab_color)
        print(spop,
              'pop',
              np.mean(D[spop]['pop']),
              np.median(D[spop]['pop']))
        if args.outliers:
            are_related, not_related = summarize_related(O,
                                                         spop,
                                                         'pop',
                                                         related)
            print(spop, 'pop',
                  args.outliers,
                  are_related, 'are related',
                  not_related, 'not related')

        ax = fig.add_subplot(inner_gs[2])
        axs.append(ax)
        inner_axs.append(ax)

        ax.hist(D[spop]['topk'],
                density = True,
                bins=args.bins,
                alpha=args.alpha,
                label='Top k',
                color='C2')
        ax.axvline(x=np.mean(D[spop]['topk']), lw=ab_lw, color=ab_color)
        print(spop,
              'topk',
              np.mean(D[spop]['topk']),
              np.median(D[spop]['topk']))
        if args.outliers:
            are_related, not_related = summarize_related(O,
                                                         spop,
                                                         'topk',
                                                         related)
            print(spop, 'topk',
                  args.outliers,
                  are_related, 'are related',
                  not_related, 'not related')


        max_y = max([ax.get_ylim()[1] for ax in inner_axs])
        for ax in inner_axs:
            ax.set_ylim((0,max_y))

        ax_i += 1

    max_x = max([ax.get_xlim()[1] for ax in axs])
    min_x = min([ax.get_xlim()[0] for ax in axs])
    for ax in axs:
        ax.set_xlim((min_x,max_x))
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if ax != axs[-1]:
            ax.set_xticks([])
        else:
            if args.x_label:
                ax.set_xlabel(args.x_label, fontsize=4)

        ax.tick_params(axis='both',
                       which='major',
                       labelsize=3,
                       width=0.25,
                       length=1)

    plt.tight_layout()
    plt.savefig(args.out_file, dpi=300)

if __name__ == '__main__':
    main()
