import matplotlib.pyplot as plt
import argparse
import sys
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description='FST for subpopulation queires')
    parser.add_argument('--inputs',
                        type=str,
                        nargs='+',
                        help='Input files')
    parser.add_argument('--labels',
                        type=str,
                        nargs='+',
                        help='Labels')
    parser.add_argument('--output',
                        type=str,
                        default='histogram.png',
                        help='Output file name')
    parser.add_argument('--height',
                        type=float,
                        default=5,
                        help='Height of the figure')
    parser.add_argument('--width',
                        type=float,
                        default=5,
                        help='Width of the figure')
    return parser.parse_args()

def get_file_data(file_name):
    D =  []
    with open(file_name, 'r') as f:
        for line in f:
            src, dst, freq, val = line.rstrip().split()
            #if src == dst: continue
            if int(freq) == 0: continue
            for i in range(int(freq)):
                D.append(float(val))    
    return D

def get_cdf(D, start, stop, step):
    D.sort()
    cdf = []
    for i in np.arange(start, stop, step):
        cdf.append(len([d for d in D if d < i])/len(D))
    return cdf

def main():
    args = get_args()

    cols = 1
    rows = 1
    fig, ax = plt.subplots(rows,
                           cols,
                           figsize=(args.width, args.height))
    for i, input_file in enumerate(args.inputs):
        D = get_file_data(input_file)
        CDF = get_cdf(D, 0.0, 0.15, 0.001)
        ax.plot(np.arange(0.0, 0.15, 0.001), CDF, label=args.labels[i])

    ax.set_xlabel('FST')
    ax.set_ylabel('Perecnt of cohort with at most FST')
    ax.legend(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(args.output)

if __name__ == '__main__':
    main()
