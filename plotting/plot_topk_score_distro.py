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
    parser.add_argument('--ped_file', type=str)
    parser.add_argument('--width', type=float, default=10.0)
    parser.add_argument('--height', type=float, default=5.0)
    parser.add_argument('--max_x', type=int)
    parser.add_argument('--bins', type=int, default=50)
    parser.add_argument('--log', action='store_true')
    return parser.parse_args()


def main():
    args = get_args()

    top_k = utils.get_top_hits(args.topk_file, get_scores=True)

    pop_map = utils.get_label_map(args.label_file, 'Population code')
    spop_map = utils.get_label_map(args.label_file, 'Superpopulation code')

    related = utils.get_related_map(args.ped_file)

    D = {}

    P = {}
    values = []
    for i in top_k:
        i_spop = spop_map[i]
        i_pop = pop_map[i]

        if ',' in i_spop: continue

        if i_spop not in D:
            D[i_spop] = []


        if i_spop not in P:
            P[i_spop] = {}
            P[i_spop]['out_spop'] = []
            P[i_spop]['in_pop'] = []
            P[i_spop]['in_spop'] = []
            P[i_spop]['in_fam'] = []

        for j, value in top_k[i]:
            if i == j: continue
            j_spop = spop_map[j]
            j_pop = pop_map[j]

            if ','in j_spop: continue

            D[i_spop].append(value)

            values.append(value)

            if j_spop != i_spop:
                #print(i_spop, j_spop, value)
                P[i_spop]['out_spop'].append(value)
            else:
                if i in related and  j in related[i]:
                    P[i_spop]['in_fam'].append(value)
                elif j_pop == i_pop:
                    P[i_spop]['in_pop'].append(value)
                else:
                    P[i_spop]['in_spop'].append(value)



    spops = ['EUR', 'EAS', 'AMR', 'SAS', 'AFR']
    samples = ['out_spop', 'in_spop', 'in_pop','in_fam']
    pretty_names = ['Out Super Pop.', 'In Super Pop.', 'In Pop.', 'In Fam.']

    fig, axs = plt.subplots(2,
                            len(spops),
                            figsize=(args.width, args.height),
                            height_ratios=[2,1])

    fig.subplots_adjust(hspace=0.1)

    ax_i = 0
    for spop in spops:
        data = []
        for sample in samples:
            if len(P[spop][sample]) == 0:
                data.append([0])
            else:
                data.append(P[spop][sample])
        vp = axs[1][ax_i].violinplot(data,
                                     showmeans=True,
                                     vert=False)
        for partname in ('cbars', 'cmins', 'cmaxes', 'cmeans'):
            vp_part = vp[partname]
            vp_part.set_linewidth(0.25)
        ax_i += 1

    for i in range(len(spops)):
        spop = spops[i]
        axs[0][i].set_title(spop, fontsize=8, loc='left')
        axs[0][i].hist(D[spop], bins=args.bins, density=True)

    #max_y = max([axs[1].get_ylim() for ax in axs])
    #min_y = min([axs[1].get_ylim()[0] for ax in axs])

    #max_x = max([axs[1].get_xlim()[1] for ax in axs])
    #min_x = min([axs[1].get_xlim()[0] for ax in axs])

    for ax in axs[0]:
        #ax.set_ylim((0,30))
        if args.max_x is None:
            ax.set_xlim(min(values), max(values))
        else:
            ax.set_xlim(min(values), args.max_x)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both',
                       which='major',
                       labelsize=4,
                       width=0.25,
                       length=2)
        ax.spines['bottom'].set_linewidth(0.25)
        ax.spines['left'].set_linewidth(0.25)
        ax.set_xticks([])
        if args.log: ax.set_yscale('log')

    axs[0][0].set_ylabel('Freq.', fontsize=6)

    for i in range(len(axs[1])):
        ax = axs[1][i]
        if i == 0:
            ax.set_yticks(range(1, len(samples)+1))
            ax.set_yticklabels(pretty_names)
        else:
            ax.set_yticks([])

        if args.max_x is None:
            ax.set_xlim(min(values), max(values))
        else:
            ax.set_xlim(min(values), args.max_x)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both',
                       which='major',
                       labelsize=4,
                       width=0.25,
                       length=2)


#    fig.text(0.5, 0.04, 'Similarity score', ha='center', fontsize=6)
        
    plt.tight_layout()
    plt.savefig(args.out_file, dpi=300)

if __name__ == '__main__':
    main()
