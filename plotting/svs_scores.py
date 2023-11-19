import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
import gzip
import re


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', type=str, required=True)
    parser.add_argument('--out_file', type=str, required=True)
    parser.add_argument('--width', type=float, default=10.0)
    parser.add_argument('--height', type=float, default=5.0)
    return parser.parse_args()

def read_score_file(file):
    scores = []
    with gzip.open(file, 'rt') as lines:
        for line in lines:
            A = line.rstrip().split()
            file_name = A[0].split('/')[1]
            numbers = re.findall(r'\d+', file_name)
            chrm = int(numbers[0])
            segment = int(numbers[1])
            mean = np.mean([float(x) for x in A[1:]])
            scores.append((chrm,segment,mean))

    return scores


def main():
    args = get_args()

    scores = read_score_file(args.in_file) 


    chrms = sorted(list(set([score[0] for score in scores])))

    fig = plt.figure(figsize=(args.width, args.height))

    gs = gridspec.GridSpec(6, 1)

    ax_i = 0
    max_xs = []
    axs = []
    for chrm in chrms:
        ax = fig.add_subplot(gs[ax_i])
        ax.set_title('Chromosome ' + str(chrm), loc='left')
        axs.append(ax)
        chrm_scores = [score for score in scores if score[0] == chrm]
        sorted_scores = sorted(chrm_scores, key=lambda x: x[1])

        max_xs.append(max([score[1] for score in sorted_scores]))

        ax.plot([score[1] for score in sorted_scores],
                    [score[2] for score in sorted_scores])

        sorted_scores = sorted(chrm_scores, key=lambda x: x[2])
        print(chrm, sorted_scores[-5:])
        ax_i += 1


    for ax in axs:
        ax.set_xlim((0,max(max_xs)))
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)


        ax.set_ylabel('Mean score')

        if ax == axs[-1]:
            ax.set_xlabel('Pos.')
        else:
            ax.set_xticklabels([])

    plt.tight_layout()
    plt.savefig(args.out_file, dpi=300)

if __name__ == '__main__':
    main()
