import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import random

from plotting import ancestry_helpers
from src import get_relations

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hits', type=str, help='top hits file', required=True)
    parser.add_argument('--trios', type=str, help='trios scores', required=True)
    parser.add_argument('--ancestry', type=str, help='population query', required=True)
    parser.add_argument('--ped', type=str, help='pedigree file', required=True)
    parser.add_argument('--out', type=str, help='output file', required=True)

    return parser.parse_args()

def get_trio_scores(trio_data_prefix):
    all_populations = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    file_suffix = '_trio_scores_20.txt'
    pop_scores = {'self': [], 'parent': [], 'child': [], 'subpop': [], 'AFR': [], 'AMR': [], 'EAS': [], 'EUR': [], 'SAS': []}

    for pop in all_populations:
        # open file for population
        with open(trio_data_prefix + pop + file_suffix, 'r') as f:
            header = f.readline()
            for line in f:
                line = line.strip().split()
                label = line[0]
                scores = [float(x) for x in line[1:]]
                pop_scores[label].extend(scores)

    return pop_scores

def plot_trio_all_populations(trio_scores, out_dir):
    # plot family members in one category and all other populations in another
    color_dict = {'in_family': 'red',
                  'out_family': 'black'}

    # plot a distribution of scores for each category on the same figure
    categories_labels = {'in_family': ['parent', 'child'],
                         'out_family': ['subpop', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']}

    fig, axs = plt.subplots(1, 1, figsize=(7, 5), dpi=150)
    fig.suptitle('GenoSiS scores for trios in all 1KG populations', fontsize=20)

    for i, category in enumerate(categories_labels):
        scores = []
        for label in categories_labels[category]:
            scores.extend(trio_scores[label])
        # plot density plot
        sns.kdeplot(scores, color=color_dict[category], ax=axs, linewidth=2, fill=True, alpha=0.6)

    # remove spines
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.set_ylabel('Density')
    axs.set_xlabel('GenoSiS Score')
    axs.legend(['Parent/Child', 'Unrelated'], loc='upper right')

    plt.tight_layout()
    plt.savefig(out_dir + 'ALL_trio_distributions.png')

def get_trio_samples(ped_file):
    trio_samples = []
    f = open(ped_file, 'r')
    header = f.readline()
    for line in f:
        line = line.strip().split()
        sampleID = line[0]
        fatherID = line[1]
        motherID = line[2]
        # if sample has a father OR mother; add sample, father, and mother to trio_samples
        if fatherID != '0' or motherID != '0':
            trio_samples.append(sampleID)
            trio_samples.append(fatherID)
            trio_samples.append(motherID)

        # remove duplicates
        trio_samples = list(set(trio_samples))
    f.close()
    return trio_samples

def read_hits_file(hits_file, trio_samples, sample_subpopulations, SUB_SUPERPOPULATIONS):
    # { superpop: {subpop: [scores], superpop: [scores], outgroup: [scores] } }
    hits = {'AFR': {'subpop': [], 'superpop': [], 'outgroup': []},
            'AMR': {'subpop': [], 'superpop': [], 'outgroup': []},
            'EAS': {'subpop': [], 'superpop': [], 'outgroup': []},
            'EUR': {'subpop': [], 'superpop': [], 'outgroup': []},
            'SAS': {'subpop': [], 'superpop': [], 'outgroup': []},}

    # query match,score match,score match,score...
    f = open(hits_file, 'r')
    for line in f:
        line = line.strip().split()
        query = line[0]
        # ignore trios
        if query in trio_samples:
            continue
        else:
            # get subpopulation
            sample_subpop = sample_subpopulations[query]
            # get superpopulation
            sample_superpop = SUB_SUPERPOPULATIONS[sample_subpop]
            # add scores to dictionary
            for match in line[1:]:
                match_ID, match_score = match.split(',')
                match_score = float(match_score)
                # # if match is a trio, ignore
                # if match_ID in trio_samples:
                #     continue
                # if match is self, ignore
                if match_ID == query:
                    continue
                match_subpop = sample_subpopulations[match_ID]
                match_superpop = SUB_SUPERPOPULATIONS[match_subpop]
                # add score to hits dictionary
                # if subpopulation match:
                if match_subpop == sample_subpop:
                    if match_score > 1000:
                        print(sample_superpop, query, match_ID, match_score)
                    hits[sample_superpop]['subpop'].append(float(match_score))
                # if superpopulation match:
                elif match_superpop == sample_superpop:
                    hits[sample_superpop]['superpop'].append(float(match_score))
                # if outgroup match:
                else:
                    hits[sample_superpop]['outgroup'].append(float(match_score))

    f.close()
    return hits

def plot_population_distributions(hits, out_dir):
    pop_color_dict = {'AFR': 'darkorange',
                  'AMR': 'mediumpurple',
                  'EAS': 'hotpink',
                  'EUR': 'steelblue',
                  'SAS': 'goldenrod'}

    fig, axs = plt.subplots(len(pop_color_dict.keys()), 1, figsize=(7, 13), dpi=150)
    fig.suptitle('GenoSiS scores for population in all 1KG', fontsize=20)

    for i, superpop in enumerate(pop_color_dict.keys()):
        for j, category in enumerate(hits[superpop].keys()):
            scores = hits[superpop][category]
            # plot density plot
            category_color_dict = {'subpop': pop_color_dict[superpop],
                                   'superpop': 'gray',
                                   'outgroup': 'black'}
            sns.kdeplot(scores, color=category_color_dict[category], ax=axs[i], linewidth=2, fill=True, alpha=0.6)

        # add population title
        axs[i].set_title(superpop, fontsize=15)
        # same x range
        try:
            max_score = max(max(hits[superpop]['subpop']), max(hits[superpop]['superpop']), max(hits[superpop]['outgroup']))
        except ValueError:
            max_score = max(max(hits[superpop]['subpop']), max(hits[superpop]['superpop']))

        axs[i].set_xlim(0, max_score + 10)

        # format
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)
        axs[i].set_ylabel('Density')
        axs[i].set_xlabel('GenoSiS Score')

        # custom legend with boxes and colors
        legend = ['Subpopulation', 'Superpopulation', 'Outgroup']
        colors = [pop_color_dict[superpop], 'gray', 'black']
        handles = [plt.Rectangle((0,0),1,1, color=colors[i], ec="k") for i in range(len(legend))]
        axs[i].legend(handles, legend, loc='upper right')


    plt.tight_layout()
    plt.savefig(out_dir + '1KG_pop_distributions.png')


def write_pop_hits(hits, out_dir):
    f = open(out_dir + '1KG_pop_hits.txt', 'w')
    for superpop in hits.keys():
        for category in hits[superpop].keys():
            f.write(superpop + ',' + category + ',')
            for score in hits[superpop][category]:
                f.write(str(score) + ' ')
            f.write('\n')
    f.close()

def main():
    args = get_args()
    trio_data_prefix = args.trios
    hits_file = args.hits
    ancestry = args.ancestry
    ped_file = args.ped
    out_dir = args.out

    # plot distribution for trio data, all populations
    trio_scores = get_trio_scores(trio_data_prefix)
    plot_trio_all_populations(trio_scores, out_dir)

    # plot distribution for population data, no trios
    trio_samples = get_trio_samples(ped_file)
    sample_subpop = ancestry_helpers.get_subpopulations(ancestry)
    SUB_SUPERPOPULATIONS = ancestry_helpers.SUB_SUPERPOPULATIONS
    hits = read_hits_file(hits_file, trio_samples, sample_subpop, SUB_SUPERPOPULATIONS)
    write_pop_hits(hits, out_dir)
    plot_population_distributions(hits, out_dir)


if __name__ == '__main__':
    main()