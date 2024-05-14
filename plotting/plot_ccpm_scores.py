import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import seaborn as sns
import scipy.stats as stats
import random

from plotting import ancestry_helpers
from src import get_relations

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ccpm', type=str, help='top hits file', required=True)
    parser.add_argument('--out', type=str, help='output directory', required=True)

    return parser.parse_args()

def read_ccpm_scores(ccpm_chrm_dir):
    all_scores = []
    q_idx = 0
    print(q_idx)
    for query in os.listdir(ccpm_chrm_dir):
        print(q_idx)
        with open(ccpm_chrm_dir + query, 'r') as f:
            query = f.readline().strip().split(': ')[1]
            for line in f:
                line = line.strip().split('\t')
                match = line[0]
                if query == match:
                    continue
                score = float(line[1])
                all_scores.append(score)
        q_idx += 1
    return all_scores

def plot_data(all_scores, out):
    # plot histogram of scores
    plt.figure(figsize=(15, 5), dpi=100)
    sns.histplot(all_scores, bins=100)
    plt.title('CCPM Scores')
    plt.xlabel('Score')
    plt.ylabel('Frequency')
    # log scale
    plt.yscale('log')
    # remove spines
    sns.despine()
    plt.savefig(out + 'ccpm_scores.png')
    plt.close()

def main():
    args = get_args()
    ccpm_dir = args.ccpm
    out = args.out

    # for all chromosomes
    all_scores = []
    for chrm in range(1, 23):
        # if chrm == 22:
        try:
            ccpm_chrm_dir = ccpm_dir + 'top_hits_' + str(chrm) + '/'
            all_scores.extend(read_ccpm_scores(ccpm_chrm_dir))
        except:
            continue

    plot_data(all_scores, out)


if __name__ == '__main__':
    main()