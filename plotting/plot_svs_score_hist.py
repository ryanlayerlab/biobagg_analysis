import argparse
import matplotlib.pyplot as plt
import glob
import sys
import random

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='Path to input file')
    parser.add_argument('-o', type=str, help='Path to output file')
    parser.add_argument('-F', type=int, help='Number of files')
    parser.add_argument('-N', type=int, help='Size of the sample')
    parser.add_argument('--width',
                        type=int,
                        default=6,
                        help='Width of the figure')
    parser.add_argument('--height',
                        type=int,
                        default=3,
                        help='Height of the figure')
    parser.add_argument('--bins',
                        type=int,
                        default=20,
                        help='Number of bins in the histogram')
    return parser.parse_args()

def main():
    args = get_args()

    files = random.sample(glob.glob(args.i), args.F)

    svs_scores = []

    i = 1
    for file in files:
        print(round(i/len(files),2), file, file=sys.stderr)
        with open(file, 'r') as f:
            query = None
            for line in f:
                if len(line) <=1 : continue
                if query is not None and line.startswith(query): continue
                if line.startswith('Query'):
                    query = line.split()[1]
                else:
                    svs_scores.append(float(line.split()[1]))
        i += 1

    svs_scores = random.sample(svs_scores, args.N)

    fig, axs = plt.subplots(1, 2, figsize=(args.width, args.height))

    ax = axs[0]
    ax.hist(svs_scores, bins=args.bins)
    ax.set_xlabel('SVS score')
    ax.set_ylabel('Freq.')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    #svs_scores = [x for x in svs_scores if x < 5.0]

    ax = axs[1]
    ax.hist(svs_scores, bins=args.bins*3)
    ax.set_xlabel('SVS score')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, 5)


    #svs_scores = [x for x in svs_scores if x < 2.0]

    #ax = axs[2]
    #ax.hist(svs_scores, bins=args.bins, density=True)
    #ax.set_xlabel('SVS score')
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)


    fig.tight_layout()
    fig.savefig(args.o)

if __name__ == '__main__':
    main()
