import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import random

from plotting import ancestry_helpers
from src import get_relations

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plink', type=str, help='plink file', required=True)
    parser.add_argument('--dist', type=str, help='ped distances and labels', required=True)
    parser.add_argument('--pop', type=str, help='population query', required=True)
    parser.add_argument('--ancestry', type=str, help='ancestry file', required=True)
    parser.add_argument('--ped', type=str, help='PED file', required=True)
    parser.add_argument('--color', type=str, help='color of plots', required=True, default='olivedrab')

    return parser.parse_args()

SUB_SUPERPOPULATIONS = ancestry_helpers.SUB_SUPERPOPULATIONS
def get_dist(dist_file):
    dist_dict = defaultdict(dict)
    label_dict = defaultdict(dict)

    with open(dist_file, 'r') as f:
        for line in f:
            line = line.strip().split(' ')
            try:
                dist_dict[line[0]].update({line[1]: float(line[2])})
                label_dict[line[0]].update({line[1]: line[3]})
            except KeyError:
                dist_dict[line[0]] = {line[1]: float(line[2])}
                label_dict[line[0]] = {line[1]: line[3]}

    return dist_dict, label_dict

def get_plink_scores(plink_file):
    plink_dict = defaultdict(dict)
    f = open(plink_file, 'r')
    header = f.readline().strip().split('\t')
    for line in f:
        line = line.strip().split()
        sample1 = line[0]
        sample2 = line[2]
        pi_hat = float(line[9])
        dst = float(line[11])
        try:
            plink_dict[sample1].update({sample2: (pi_hat, dst)})
        except KeyError:
            plink_dict[sample1] = {sample2: (pi_hat, dst)}

    return plink_dict

def plot_1KG_trios(plink_dict, dist_dict, label_dict, pop, color, subpopulations, samples):
    #x: [y,y,y,y]
    score_data = defaultdict(list)
    label_data = defaultdict(list)

    plt.figure()
    labels_ordered = {'self': 'self',
                      'parent': 'parent\nchild',
                      'child': 'parent\child',
                      'sibling': 'sibling',
                      'grandparent':'grandparent\ngrandchild',
                      'grandchild':'grandparent\ngrandchild',
                      'subpop':'subpop',
                      'superpop':'superpop',
                      'outpop':'outpop',
                      'AFR':'AFR',
                      'AMR':'AMR',
                      'EAS':'EAS',
                      'EUR':'EUR',
                      'SAS':'SAS',
                      'unknown':'unknown'}
    x_labels_ordered = ['self', 'parent\nchild', 'grandparent\ngrandchild',
                        'subpop', 'superpop',
                        'AFR', 'AMR', 'EAS', 'EUR', 'SAS']

    outpop_color_dict = {'AFR': 'darkorange',
                         'AMR': 'mediumpurple',
                         'EAS': 'hotpink',
                         'EUR': 'steelblue',
                         'SAS': 'goldenrod'}
    outpop_colors = []

    # fill data
    for query in samples:
        for match in plink_dict[query]:
            # dist = dist_dict[query][match]
            try:
                label = label_dict[query][match]
                if label == 'outpop':
                    label = SUB_SUPERPOPULATIONS[subpopulations[match]]
            except KeyError:
                label = SUB_SUPERPOPULATIONS[subpopulations[match]]

            try:
                pi_hat_score = plink_dict[query][match][0]
                dst_score = plink_dict[query][match][1]
                # if label == 'outpop':
                #     if score > 1000:
                #         print(query, match, dist_dict[query][match], hits_dict[query][match])
            except KeyError:
                continue


            try:
                # score_data[dist].append(score)
                label_data[labels_ordered[label]].append(pi_hat_score)
            except KeyError:
                # score_data[dist] = [score]
                label_data[labels_ordered[label]] = [pi_hat_score]
            # print(query, match, dist_dict[query][match], hits_dict[query][match])


    # violin plot for labels
    label_data_plot = []
    for l in x_labels_ordered:
        if len(label_data[l]) == 0:
            label_data_plot.append([0])
        else:
            label_data_plot.append(label_data[l])

        [label_data[l] for l in x_labels_ordered]
    # make two subplots one for violin plot one for scatter plot

    fig, ax = plt.subplots(2, 1, figsize=(28, 20))

    # top plot
    ax[0].violinplot(label_data_plot, showmeans=True)

    # color all parts of the violin plot
    for i in range(0, len(x_labels_ordered)):
        ax[0].get_children()[i].set_facecolor(color)
        ax[0].get_children()[i].set_edgecolor(color)
        for j in range(1, 5):
            ax[0].get_children()[i + j].set_color(color)

    for i in range(0, len(x_labels_ordered)):
        # make x a list of values that is the same length as the y values, jittered
        x = [i + 1 + random.uniform(-0.25, 0.25) for j in range(0, len(label_data_plot[i]))]
        ax[1].scatter(x, label_data_plot[i], color=color, alpha=0.3)


    fig.suptitle('Plink pi_hat Scores for 1KG Trios\n(' + pop + ')', fontsize=40, fontweight='bold', color=color)
    # top plot
    ax[0].set_ylabel('plink pi-hat score', labelpad=10, fontsize=30, fontweight='bold')
    ax[0].tick_params(axis='y', labelsize=20)
    # ax[0].set_yticks(range(0, 10, 1), [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], fontsize=20)
    ax[0].set_xticks([])

    ax[1].tick_params(axis='y', labelsize=20)
    # ax[1].set_yticks(range(0, 10, 1), [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], fontsize=20)
    ax[1].set_xlabel('Relationship', labelpad=10, fontsize=30, fontweight='bold')
    ax[1].set_xticks(range(1, len(x_labels_ordered)+1), x_labels_ordered, fontsize=20)
    ax[1].set_ylabel('plink pi-hat score', labelpad=10, fontsize=30, fontweight='bold')
    # ax[1].set_yticks(fontsize=20)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)


    # plt.xlabel('Relationship', labelpad=10, fontsize=30, fontweight='bold')
    # plt.ylabel('GenoSiS Score', labelpad=10, fontsize=30, fontweight='bold')
    # plt.xticks(range(1, len(x_labels_ordered)+1), x_labels_ordered, fontsize=20)
    # plt.yticks(fontsize=20)
    # plt.gca().spines['right'].set_visible(False)
    # plt.gca().spines['top'].set_visible(False)

    plt.savefig('1kg_trio_plots/1KG_trios_plink_' + pop + '.png')

    # # violin plot
    # plt.violinplot(score_data.values(), showmeans=True)
    # # label x-ticks
    # plt.xticks(range(1, len(score_data.keys())+1), score_data.keys())
    # plt.xlabel('Genetic Distance')
    # plt.ylabel('GenoSiS Score')
    # plt.savefig('1KG_trios.png')



def main():
    args = get_args()
    pop = args.pop
    dist_dict, label_dict = get_dist(args.dist + '_' + pop + '.txt')
    plink_dict = get_plink_scores(args.plink)
    subpopulations = ancestry_helpers.get_subpopulations(args.ancestry)
    samples = get_relations.get_samples(args.ped, pop, subpopulations)

    color = args.color
    plot_1KG_trios(plink_dict, dist_dict, label_dict, pop, color, subpopulations, samples)

    # print(hits_dict)
    # dist_dict = get_hits(args.dist)


if __name__ == '__main__':
    main()