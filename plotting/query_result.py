import argparse
import gzip
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import utils

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--svs_results', type=str, required=True)
    parser.add_argument('--chrm', type=str, required=True)
    parser.add_argument('--target', type=str, required=True)
    parser.add_argument('--label_file', type=str, required=False)
    parser.add_argument('--k', type=int, default=10)
    parser.add_argument('--out_file', type=str, required=True)
    parser.add_argument('--height', type=int, default=3)
    parser.add_argument('--width', type=int, default=10)
    return parser.parse_args()


def get_ranks(L):
    sorted_list = sorted(L, key=lambda x: x[1], reverse=True)
    rank_dict = {i: (rank+1)/len(L) for rank, (i, s) \
            in enumerate(sorted_list)}
    #rank_dict = {i: s for rank, (i, s) \
            #in enumerate(sorted_list)}
    return rank_dict

def get_svs_target_results(svs_file, target, K):
    hits = {}
    seen_dict = {}
    with gzip.open(svs_file, 'rt') as lines:
        for line in lines:
            A = line.rstrip().split()
            chrm = A[0]
            start = int(A[1])
            end = int(A[2])
            s1 = A[3]

            if s1 != target: continue

            s2 = A[4]

            if s1==s2: continue

            sim = float(A[5])


            if (start,end) not in hits:
                hits[(start,end)] = {}
            hits[(start,end)][s2] = sim
            if s2 not in seen_dict:
                seen_dict[s2] = 0
            seen_dict[s2] = seen_dict[s2] + 1

    seen_list = sorted([(s2, seen_dict[s2]) for s2 in seen_dict],
                       key = lambda x: x[1],
                       reverse = True)

    return hits, seen_list[:K]

def main():
    args = get_args()

    svs_results, top_k = get_svs_target_results(args.svs_results,
                                                args.target,
                                                args.k)

    spop_map = None
    if args.label_file:
        spop_map = utils.get_label_map(args.label_file, 'Superpopulation code')

    top_k_samples = [t[0] for t in top_k]

    print(top_k_samples)

    X = []
    Y = []
    X_in = []
    X_out = []
    Y_in = []
    Y_out = []

    for start, end in svs_results:
        for s2 in svs_results[(start, end)]:
            if spop_map is not None:

                if spop_map[args.target] == spop_map[s2]:
                    X_in.append(start)
                    Y_in.append(svs_results[(start,end)][s2])
                else:
                    X_out.append(start)
                    Y_out.append(svs_results[(start,end)][s2])
            else:
                X.append(start)
                Y.append(svs_results[(start,end)][s2])


    fig, ax = plt.subplots( figsize=(args.width, args.height) )


    if spop_map is not None:
        ax.plot(X_out,Y_out,
                'o',
                ms=2,
                markerfacecolor='None',
                markeredgecolor='C1',
                markeredgewidth=0.5,
                alpha=0.5,
                label='Out population')

        ax.plot(X_in,Y_in,
                'o',
                ms=2,
                markerfacecolor='None',
                markeredgecolor='C0',
                markeredgewidth=0.5,
                alpha=0.5,
                label='In population')
        plt.legend(frameon=False,
                   fontsize=6,
                   loc='upper right')
    else:
        ax.plot(X,Y,
                'o',
                ms=2,
                markerfacecolor='None',
                markeredgecolor='grey',
                markeredgewidth=0.5,
                alpha=0.1 if args.k>0 else 0.5)


    if args.k > 0 and args.k <= 5:
        
        colors = ['C0', 'C1', 'C2', 'C3', 'C4']
        labels = ['1st', '2nd', '3rd', '4th', '5th']
        X = []
        Y = []
        C = []
        for start, end in svs_results:
            for s2 in svs_results[(start, end)]:
                if s2 not in top_k_samples: continue
                X.append(start)
                Y.append(svs_results[(start,end)][s2])
                C.append(colors[top_k_samples.index(s2)])



        ax.scatter(X,Y,
                   s=5,
                   lw=1,
                   c=C)

        legend_elements = []
        for i in range(len(colors)):
            legend_elements.append(Line2D([0],
                                          [0],
                                          markerfacecolor=colors[i],
                                          color=colors[i],
                                          marker='o',
                                          lw=0,
                                          markersize=3,
                                          label=labels[i]))
        plt.legend(handles=legend_elements,
                   frameon=False,
                   fontsize=6,
                   loc='upper right')


    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Chromosome ' + args.chrm)
    ax.set_ylabel('Similarity')
    
    plt.tight_layout()
    plt.savefig(args.out_file, dpi=300)

    plt.close()

if __name__ == '__main__':
    main()
