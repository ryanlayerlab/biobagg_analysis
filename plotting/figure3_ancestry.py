import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import sys

sys.path.append(os.path.abspath('plotting/'))
import ancestry_helpers

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestry', help='1KG ancestry labels', required=True)
    parser.add_argument('--k', help='k used for knn', required=True)
    parser.add_argument('--colors', help='file with color codes', required=True)
    # GenoSiS Scores
    parser.add_argument('--genosis_groups', help='1KG genosis scores for groups', required=True)
    parser.add_argument('--genosis_k', help='1KG genosis scores for top K', required=True)
    # plink group scores
    parser.add_argument('--dst_groups', help='1KG plink DST scores for groups', required=True)
    parser.add_argument('--pihat_groups', help='1KG plink pi-hat scores for groups', required=True)
    parser.add_argument('--kinship_groups', help='1KG plink kinship scores for groups', required=True)
    # plink top K scores
    parser.add_argument('--dst_k', help='1KG plink top K DST scores', required=True)
    parser.add_argument('--pihat_k', help='1KG plink top K pi-hat scores', required=True)
    parser.add_argument('--kinship_k', help='1KG plink top K kinship scores', required=True)
    # Output
    parser.add_argument('--png_dist', help='Output png file with density plots', required=True)
    parser.add_argument('--png_k', help='Output png file with top k percents', required=True)

    return parser.parse_args()

def get_category_colors(colors, superpop):
    return {
        'superpop': colors['superpop'],
        'subpop': colors[superpop],
        'outgroup': colors['outgroup']
    }

def read_ancestry_group_scores(ancestry_group_scores_file):
    '''
    Read scores for superpop, subpop, and outpop labels.
    @param ancestry_group_scores_file: path to scores file for 1KG data
    @return: dictionary of scores for each superpop
    '''
    # superpop: category: [scores...]
    scores = defaultdict(dict)
    with open(ancestry_group_scores_file, 'r') as f:
        header = f.readline().strip().split(',')
        for line in f:
            line = line.strip().split(',')
            superpop = line[0]
            category = line[1]
            scores[superpop][category] = [float(x) for x in line[2].strip().split()]
    return scores

def read_top_K(top_k_file,
               subpopulations,
               header=True):
    '''
    Read top K scores and return top K samples, and populations
    @param top_k_file: path top k scores
    @param subpopulations: dictionary of subpopulations for each sample
    @return: dictionary of top K samples and populations for each query
    '''
    top_K_samples = {}
    top_K_subpopulations = {}
    f = open(top_k_file, 'r')
    if header:
        f.readline()
    for line in f:
        line = line.strip().split()
        query = line[0]
        top_K_samples[query] = []
        top_K_subpopulations[query] = []
        matches = line[1:]
        for match_score in matches:
            match = match_score.split(',')[0]
            top_K_samples[query].append(match)
            try:
                top_K_subpopulations[query].append(subpopulations[match])
            except:
                pass # header
    return top_K_samples, top_K_subpopulations

def get_percent_in_top_k(top_K_samples,
                         subpopulations):
    '''

    @param top_K_samples: dictionary with key = sample and value = list of top k samples
    @param top_K_subpopulations: dictionary with key = sample and value = list of top k subpopulations

    @return: dictionary with key = sample and value = percent in top k that are subpop, superpop, outgroup
    '''
    sample_percents = {}
    for query in top_K_samples.keys():
        k = len(top_K_samples[query])
        query_subpop = subpopulations[query]
        query_superpop = ancestry_helpers.SUB_SUPERPOPULATIONS[query_subpop]
        subpop_count = 0
        superpop_count = 0
        outgroup_count = 0
        for match in top_K_samples[query]:
            match_subpop = subpopulations[match]
            match_superpop = ancestry_helpers.SUB_SUPERPOPULATIONS[match_subpop]
            if query_subpop == match_subpop:
                subpop_count += 1
                superpop_count += 1
            elif query_superpop == match_superpop:
                superpop_count += 1
            else:
                outgroup_count += 1
        sample_percents[query] = {'subpop': subpop_count / k,
                                    'superpop': superpop_count / k,
                                    'outgroup': outgroup_count / k}
    return sample_percents


def get_pop_counts(sample_percents,
                   subpopulations):
    '''
    Get the counts of samples in the same subpopulation, superpopulation, and outgroup
    @param sample_percents:
    @return:
    '''
    counts = defaultdict(dict)
    for sample in sample_percents.keys():
        # get superpop of sample
        sample_superpop = ancestry_helpers.SUB_SUPERPOPULATIONS[subpopulations[sample]]
        for category in sample_percents[sample].keys():
            if category not in counts[sample_superpop]:
                counts[sample_superpop][category] = []
            counts[sample_superpop][category].append(sample_percents[sample][category])
    return counts




def get_population_percents(top_K_subpopulations,
                            k,
                            subpopulations):
    '''
    Get the percent of matches in the same subpopulation and superpopulation
    @param top_K_subpopulations:
    @param k: k used for knn
    @return:
    '''
    num_samples = {'AFR': 0,
                     'AMR': 0,
                     'EAS': 0,
                     'EUR': 0,
                     'SAS': 0}

    # start with k as k, and move down to k = 1
    superpop_counts = {'AFR': {k_i:0 for k_i in range(1,21)},
                       'AMR': {k_i:0 for k_i in range(1,21)},
                       'EAS': {k_i:0 for k_i in range(1,21)},
                       'EUR': {k_i:0 for k_i in range(1,21)},
                       'SAS': {k_i:0 for k_i in range(1,21)}}
    subpop_counts = {'AFR': {k_i:0 for k_i in range(1,21)},
                       'AMR': {k_i:0 for k_i in range(1,21)},
                       'EAS': {k_i:0 for k_i in range(1,21)},
                       'EUR': {k_i:0 for k_i in range(1,21)},
                       'SAS': {k_i:0 for k_i in range(1,21)}}
    cohort_size = k
    while cohort_size >= 1:
        for query in top_K_subpopulations.keys():
            query_subpop = subpopulations[query]
            query_superpop = ancestry_helpers.SUB_SUPERPOPULATIONS[query_subpop]
            if cohort_size == k:
                num_samples[query_superpop] += 1
            match_index = 0
            for match_subpop in top_K_subpopulations[query]:
                if match_index >= cohort_size:
                    continue
                match_superpop = ancestry_helpers.SUB_SUPERPOPULATIONS[match_subpop]
                if query_subpop == match_subpop:
                    try:
                        superpop_counts[query_superpop][cohort_size] += 1
                        subpop_counts[query_superpop][cohort_size] += 1
                    except KeyError:
                        superpop_counts[query_superpop]
                        subpop_counts[query_superpop][cohort_size] = 1
                elif query_superpop == match_superpop:
                    try:
                        superpop_counts[query_superpop][cohort_size] += 1
                    except KeyError:
                        superpop_counts[query_superpop][cohort_size] = 1

                match_index += 1

        cohort_size -= 1

    # get percents
    superpop_percents = {'AFR': {p: superpop_counts['AFR'][p] / (p * num_samples['AFR']) for p in superpop_counts['AFR'].keys()},
                         'AMR': {p: superpop_counts['AMR'][p] / (p * num_samples['AMR']) for p in superpop_counts['AMR'].keys()},
                         'EAS': {p: superpop_counts['EAS'][p] / (p * num_samples['EAS']) for p in superpop_counts['EAS'].keys()},
                         'EUR': {p: superpop_counts['EUR'][p] / (p * num_samples['EUR']) for p in superpop_counts['EUR'].keys()},
                         'SAS': {p: superpop_counts['SAS'][p] / (p * num_samples['SAS']) for p in superpop_counts['SAS'].keys()}
                         }
    subpop_percents = {'AFR': {p: subpop_counts['AFR'][p] / (p * num_samples['AFR']) for p in subpop_counts['AFR'].keys()},
                       'AMR': {p: subpop_counts['AMR'][p] / (p * num_samples['AMR']) for p in subpop_counts['AMR'].keys()},
                       'EAS': {p: subpop_counts['EAS'][p] / (p * num_samples['EAS']) for p in subpop_counts['EAS'].keys()},
                       'EUR': {p: subpop_counts['EUR'][p] / (p * num_samples['EUR']) for p in subpop_counts['EUR'].keys()},
                       'SAS': {p: subpop_counts['SAS'][p] / (p * num_samples['SAS']) for p in subpop_counts['SAS'].keys()}
                       }

    return superpop_percents, subpop_percents


def plot_ancestry_group_distributions(genosis_scores,
                                      dst_scores,
                                      pihat_scores,
                                      kinship_scores,
                                      colors,
                                      png_file):

    alpha_value = 0.5

    # Create the combined figure
    combined_figure, axes = plt.subplots(5, 4, figsize=(20, 12), dpi=300)

    col1_xlim = (0, 4000)
    col1_ylim = (0, 0.01)
    col2_xlim = (0.95, 1)
    col2_ylim = (0, 1600)
    col3_xlim = (-0.1, 1)
    col3_ylim = (0, 250)
    col4_xlim = (-0.25, 0.3)
    col4_ylim = (0, 80)

    # set axes limits
    for col_j in range(4):
        for row_i in range(5):
            if col_j == 0:
                axes[row_i, col_j].set_xlim(col1_xlim)
                axes[row_i, col_j].set_ylim(col1_ylim)
            elif col_j == 1:
                axes[row_i, col_j].set_xlim(col2_xlim)
                axes[row_i, col_j].set_ylim(col2_ylim)
            elif col_j == 2:
                axes[row_i, col_j].set_xlim(col3_xlim)
                axes[row_i, col_j].set_ylim(col3_ylim)
            elif col_j == 3:
                axes[row_i, col_j].set_xlim(col4_xlim)
                axes[row_i, col_j].set_ylim(col4_ylim)

    # Plot the GenoSiS scores in first column
    for i, superpop in enumerate(genosis_scores.keys()):
        D = []
        colors_list = []
        for j, category in enumerate(genosis_scores[superpop].keys()):
            D.append(genosis_scores[superpop][category])
            intermediate_colors = get_category_colors(colors, superpop)
            colors_list.append(intermediate_colors[category])
        sns.kdeplot(D,
                    ax=axes[i, 0],
                    label=superpop,
                    palette=colors_list,
                    fill=True, alpha=alpha_value)
        # try:
        #     sns.histplot(D[0], ax=axes[i, 0], color=colors_list[0], kde=False, binwidth=50)
        # except:
        #     pass
        # try:
        #     sns.histplot(D[1], ax=axes[i, 0], color=colors_list[1], kde=False, binwidth=50)
        # except:
        #     pass
        # try:
        #     sns.histplot(D[2], ax=axes[i, 0], color=colors_list[2], kde=False, binwidth=50)
        # except:
        #     pass
        # # log y - axis
        # axes[i, 0].set_yscale('log')
        axes[0, 0].set_title('GenoSiS Scores', fontsize=20, fontweight='bold')
        axes[i, 0].set_ylabel('Density')
        # axes[i, 0].set_xlim(0, genosis_x_max)
        axes[i, 0].spines['top'].set_visible(False)
        axes[i, 0].spines['right'].set_visible(False)
        legend = ['Subpopulation', 'Superpopulation', 'Outgroup']
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors_list[i]) for i in range(len(legend))]
        axes[i, 0].legend(handles, legend, loc='upper right', frameon=False)

    # Plot the plink DST scores in second column
    for i, superpop in enumerate(dst_scores.keys()):
        D = []
        colors_list = []
        for j, category in enumerate(dst_scores[superpop].keys()):
            D.append(dst_scores[superpop][category])
            intermediate_colors = get_category_colors(colors, superpop)
            colors_list.append(intermediate_colors[category])
        sns.kdeplot(D,
                    ax=axes[i, 1],
                    label=superpop,
                    palette=colors_list,
                    fill=True, alpha=alpha_value)
        # try:
        #     sns.histplot(D[0], ax=axes[i, 1], color=colors_list[0], kde=False, binwidth=0.001)
        # except:
        #     pass
        # try:
        #     sns.histplot(D[1], ax=axes[i, 1], color=colors_list[1], kde=False, binwidth=0.001)
        # except:
        #     pass
        # try:
        #     sns.histplot(D[2], ax=axes[i, 1], color=colors_list[2], kde=False, binwidth=0.001)
        # except:
        #     pass
        axes[0, 1].set_title('Plink DST Scores', fontsize=20, fontweight='bold')
        axes[i, 1].set_ylabel('Density')
        axes[i, 1].spines['top'].set_visible(False)
        axes[i, 1].spines['right'].set_visible(False)
        legend = ['Subpopulation', 'Superpopulation', 'Outgroup']
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors_list[i]) for i in range(len(legend))]
        axes[i, 1].legend(handles, legend, loc='upper right', frameon=False)

    # Plot the plink pi-hat scores in third column
    for i, superpop in enumerate(pihat_scores.keys()):
        D = []
        colors_list = []
        for j, category in enumerate(pihat_scores[superpop].keys()):
            D.append(pihat_scores[superpop][category])
            intermediate_colors = get_category_colors(colors, superpop)
            colors_list.append(intermediate_colors[category])
        sns.kdeplot(D,
                    ax=axes[i, 2],
                    label=superpop,
                    palette=colors_list,
                    fill=True, alpha=alpha_value)
        # try:
        #     sns.histplot(D[0], ax=axes[i, 2], color=colors_list[0], kde=False, binwidth=0.01)
        # except:
        #     pass
        # try:
        #     sns.histplot(D[1], ax=axes[i, 2], color=colors_list[1], kde=False, binwidth=0.01)
        # except:
        #     pass
        # try:
        #     sns.histplot(D[2], ax=axes[i, 2], color=colors_list[2], kde=False, binwidth=0.01)
        # except:
        #     pass
        axes[0, 2].set_title('Plink pi-hat Scores', fontsize=20, fontweight='bold')
        axes[i, 2].set_ylabel('Density')
        axes[i, 2].spines['top'].set_visible(False)
        axes[i, 2].spines['right'].set_visible(False)
        legend = ['Subpopulation', 'Superpopulation', 'Outgroup']
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors_list[i]) for i in range(len(legend))]
        axes[i, 2].legend(handles, legend, loc='upper right', frameon=False)

    # Plot the plink kinship scores in fourth column
    for i, superpop in enumerate(kinship_scores.keys()):
        D = []
        colors_list = []
        for j, category in enumerate(kinship_scores[superpop].keys()):
            D.append(kinship_scores[superpop][category])
            intermediate_colors = get_category_colors(colors, superpop)
            colors_list.append(intermediate_colors[category])
        sns.kdeplot(D,
                    ax=axes[i, 3],
                    label=superpop,
                    palette=colors_list,
                    fill=True, alpha=alpha_value)
        # try:
        #     sns.histplot(D[0], ax=axes[i, 3], color=colors_list[0], kde=False, binwidth=0.01)
        # except:
        #     pass
        # try:
        #     sns.histplot(D[1], ax=axes[i, 3], color=colors_list[1], kde=False, binwidth=0.01)
        # except:
        #     pass
        # try:
        #     sns.histplot(D[2], ax=axes[i, 3], color=colors_list[2], kde=False, binwidth=0.001)
        # except:
        #     pass
        axes[0, 3].set_title('Plink kinship Scores', fontsize=20, fontweight='bold')
        axes[i, 3].set_ylabel('Density')
        axes[i, 3].spines['top'].set_visible(False)
        axes[i, 3].spines['right'].set_visible(False)
        legend = ['Subpopulation', 'Superpopulation', 'Outgroup']
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors_list[i]) for i in range(len(legend))]
        axes[i, 3].legend(handles, legend, loc='upper right', frameon=False)

    # shift plots down and to the right
    plt.subplots_adjust(top=0.85, right=1.)

    # create a text label to the left of all plots
    # one label for each row
    for i, superpop in enumerate(genosis_scores.keys()):
        axes[i, 0].text(-0.23, 0.5, superpop, fontsize=20, fontweight='bold',
                        horizontalalignment='center',
                        verticalalignment='center',
                        rotation='vertical',
                        transform=axes[i, 0].transAxes)

    # add text at top
    plt.suptitle('Cohorts for 1KG Data\n(k=20)', fontsize=30, fontweight='bold')

    # Save the figure
    plt.tight_layout()
    combined_figure.savefig(png_file)

def plot_ancestry_top_k(genosis_scores,
                        genosis_superpop_percents, genosis_subpop_percents,
                        dst_superpop_percents, dst_subpop_percents,
                        pihat_superpop_percents, pihat_subpop_percents,
                        kinship_superpop_percents, kinship_subpop_percents,
                        colors,
                        png_file):

    # Create the combined figure
    combined_figure, axes = plt.subplots(5, 3, figsize=(15, 12), dpi=300)
    alpha_value = 0.5
    genosis_x_max = 1500

    # # Plot the GenoSiS scores in first column
    # for i, superpop in enumerate(genosis_scores.keys()):
    #     D = []
    #     colors_list = []
    #     for j, category in enumerate(genosis_scores[superpop].keys()):
    #         D.append(genosis_scores[superpop][category])
    #         intermediate_colors = get_category_colors(colors, superpop)
    #         colors_list.append(intermediate_colors[category])
    #     sns.kdeplot(D,
    #                 ax=axes[i, 0],
    #                 label=superpop,
    #                 palette=colors_list,
    #                 fill=True, alpha=alpha_value)
    #     axes[0, 0].set_title('GenoSiS Scores', fontsize=20, fontweight='bold')
    #     axes[i, 0].set_ylabel('Density')
    #     axes[i, 0].set_xlim(0, genosis_x_max)
    #     axes[i, 0].spines['top'].set_visible(False)
    #     axes[i, 0].spines['right'].set_visible(False)
    #     legend = ['Subpopulation', 'Superpopulation', 'Outgroup']
    #     handles = [plt.Rectangle((0, 0), 1, 1, color=colors_list[i], ec="k") for i in range(len(legend))]
    #     axes[i, 0].legend(handles, legend, loc='upper right', frameon=False)

    # Plot the genoSiS top k percents and dst top K percents in second column
    # line plot for genosis
    # dashed line plot for plink
    for i, superpop in enumerate(genosis_superpop_percents.keys()):
        x = list(genosis_superpop_percents[superpop].keys())
        genosis_y = [list(genosis_superpop_percents[superpop].values()),
                     list(genosis_subpop_percents[superpop].values())]
        D_genosis_superpop = pd.DataFrame({'K': x,
                          'GenoSiS Scores': genosis_y[0]})

        intermediate_colors = get_category_colors(colors, superpop)
        super_color = intermediate_colors['superpop']
        sub_color = intermediate_colors['subpop']

        sns.lineplot(x='K', y='GenoSiS Scores',
                     data=D_genosis_superpop,
                     ax=axes[i, 0],
                     label='Genosis Superpop',
                     color=super_color)
        D_genosis_subpop = pd.DataFrame({'K': x,
                          'GenoSiS Scores': genosis_y[1]})
        sns.lineplot(x='K', y='GenoSiS Scores',
                        data=D_genosis_subpop,
                        ax=axes[i, 0],
                        label='Genosis Subpop',
                        color=sub_color)
        D_dst_superpop = pd.DataFrame({'K': x,
                          'Plink DST Scores': list(dst_superpop_percents[superpop].values())})
        sns.lineplot(x='K', y='Plink DST Scores',
                        data=D_dst_superpop,
                        ax=axes[i, 0],
                        label='Plink DST superpop',
                        color=super_color,
                        linestyle='dashed')
        D_dst_subpop = pd.DataFrame({'K': x,
                          'Plink DST Scores': list(dst_subpop_percents[superpop].values())})
        sns.lineplot(x='K', y='Plink DST Scores',
                        data=D_dst_subpop,
                        ax=axes[i, 0],
                        label='Plink DST subpop',
                        color=sub_color,
                        linestyle='dashed')
        axes[0, 0].set_title('GenoSiS vs.\nPlink DST', fontsize=20, fontweight='bold')
        axes[i, 0].set_ylabel('% in Population')
        axes[i, 0].set_xticks(range(5, 21, 5))
        axes[i, 0].set_ylim(0, 1.1)
        axes[i, 0].spines['top'].set_visible(False)
        axes[i, 0].spines['right'].set_visible(False)

    # Plot the genoSiS top k percents and pihat top K percents in third column
    for i, superpop in enumerate(pihat_superpop_percents.keys()):
        x = list(pihat_superpop_percents[superpop].keys())
        genosis_y = [list(genosis_superpop_percents[superpop].values()),
                     list(genosis_subpop_percents[superpop].values())]
        D_genosis_superpop = pd.DataFrame({'K': x,
                          'GenoSiS Scores': genosis_y[0]})

        intermediate_colors = get_category_colors(colors, superpop)
        super_color = intermediate_colors['superpop']
        sub_color = intermediate_colors['subpop']

        sns.lineplot(x='K', y='GenoSiS Scores',
                     data=D_genosis_superpop,
                     ax=axes[i, 1],
                     label='Genosis Superpop',
                     color=super_color)
        D_genosis_subpop = pd.DataFrame({'K': x,
                          'GenoSiS Scores': genosis_y[1]})
        sns.lineplot(x='K', y='GenoSiS Scores',
                        data=D_genosis_subpop,
                        ax=axes[i, 1],
                        label='Genosis Subpop',
                        color=sub_color)
        D_pihat_superpop = pd.DataFrame({'K': x,
                          'Plink pi-hat Scores': list(pihat_superpop_percents[superpop].values())})
        sns.lineplot(x='K', y='Plink pi-hat Scores',
                        data=D_pihat_superpop,
                        ax=axes[i, 1],
                        label='Plink pi-hat superpop',
                        color=super_color,
                        linestyle='dashed')
        D_pihat_subpop = pd.DataFrame({'K': x,
                          'Plink pi-hat Scores': list(pihat_subpop_percents[superpop].values())})
        sns.lineplot(x='K', y='Plink pi-hat Scores',
                        data=D_pihat_subpop,
                        ax=axes[i, 1],
                        label='Plink pi-hat subpop',
                        color=sub_color,
                        linestyle='dashed')
        axes[0, 1].set_title('GenoSiS vs.\nPlink pi-hat', fontsize=20, fontweight='bold')
        axes[i, 1].set_ylabel('% in Population')
        axes[i, 1].set_xticks(range(5, 21, 5))
        axes[i, 1].set_ylim(0, 1.1)
        axes[i, 1].spines['top'].set_visible(False)
        axes[i, 1].spines['right'].set_visible(False)

    # Plot the genoSiS top k percents and kinship top K percents in fourth column
    for i, superpop in enumerate(kinship_superpop_percents.keys()):
        x = list(kinship_superpop_percents[superpop].keys())
        genosis_y = [list(genosis_superpop_percents[superpop].values()),
                     list(genosis_subpop_percents[superpop].values())]
        D_genosis_superpop = pd.DataFrame({'K': x,
                          'GenoSiS Scores': genosis_y[0]})

        intermediate_colors = get_category_colors(colors, superpop)
        super_color = intermediate_colors['superpop']
        sub_color = intermediate_colors['subpop']

        sns.lineplot(x='K', y='GenoSiS Scores',
                     data=D_genosis_superpop,
                     ax=axes[i, 2],
                     label='Genosis superpop',
                     color=super_color)
        D_genosis_subpop = pd.DataFrame({'K': x,
                          'GenoSiS Scores': genosis_y[1]})
        sns.lineplot(x='K', y='GenoSiS Scores',
                        data=D_genosis_subpop,
                        ax=axes[i, 2],
                        label='Genosis subpop',
                        color=sub_color)
        D_kinship_superpop = pd.DataFrame({'K': x,
                          'Plink kinship Scores': list(kinship_superpop_percents[superpop].values())})
        sns.lineplot(x='K', y='Plink kinship Scores',
                        data=D_kinship_superpop,
                        ax=axes[i, 2],
                        label='Plink kinship superpop',
                        color=super_color,
                        linestyle='dashed')
        D_kinship_subpop = pd.DataFrame({'K': x,
                          'Plink kinship Scores': list(kinship_subpop_percents[superpop].values())})
        sns.lineplot(x='K', y='Plink kinship Scores',
                        data=D_kinship_subpop,
                        ax=axes[i, 2],
                        label='Plink kinship subpop',
                        color=sub_color,
                        linestyle='dashed')
        axes[0, 2].set_title('GenoSiS vs.\nPlink Kinship', fontsize=20, fontweight='bold')
        axes[i, 2].set_ylabel('% in Population')
        axes[i, 2].set_xticks(range(5, 21, 5))
        axes[i, 2].set_ylim(0, 1.1)
        axes[i, 2].spines['top'].set_visible(False)
        axes[i, 2].spines['right'].set_visible(False)

        # shift plots down and to the right
        plt.subplots_adjust(top=0.85, right=1.)

        # create a text label to the left of all plots
        # one label for each row
        for i, superpop in enumerate(genosis_scores.keys()):
            axes[i, 0].text(-0.23, 0.5, superpop, fontsize=20, fontweight='bold',
                            horizontalalignment='center',
                            verticalalignment='center',
                            rotation='vertical',
                            transform=axes[i, 0].transAxes)

        # add text at top
        plt.suptitle('Cohorts for 1KG Data\n(k=20)', fontsize=30, fontweight='bold')

    # Save the figure
    plt.tight_layout()
    combined_figure.savefig(png_file)


def main():
    args = parse_args()

    subpopulations = ancestry_helpers.get_subpopulations(args.ancestry)

    genosis_group_scores = read_ancestry_group_scores(args.genosis_groups)
    dst_group_scores = read_ancestry_group_scores(args.dst_groups)
    pihat_group_scores = read_ancestry_group_scores(args.pihat_groups)
    kinship_group_scores = read_ancestry_group_scores(args.kinship_groups)
    png_dist = args.png_dist

    colors = ancestry_helpers.get_colors(args.colors)

    plot_ancestry_group_distributions(genosis_group_scores,
                                      dst_group_scores,
                                      pihat_group_scores,
                                      kinship_group_scores,
                                      colors,
                                      png_dist)


    # genosis_top_k_samples, genosis_top_k_subpopulations = read_top_K(args.genosis_k, subpopulations, header=False)
    # genosis_superpop_percents, genosis_subpop_percents = get_population_percents(genosis_top_k_subpopulations,
    #                                                                              int(args.k),
    #                                                                              subpopulations)
    # genosis_sample_percents = get_percent_in_top_k(genosis_top_k_samples, subpopulations)
    #
    #
    # dst_top_k_samples, dst_top_k_subpopulations = read_top_K(args.dst_k, subpopulations)
    # dst_superpop_percents, dst_subpop_percents = get_population_percents(dst_top_k_subpopulations,
    #                                                                      int(args.k),
    #                                                                      subpopulations)
    # dst_sample_percents = get_percent_in_top_k(dst_top_k_samples, subpopulations)
    #
    # pihat_top_k_samples, pihat_top_k_subpopulations = read_top_K(args.pihat_k, subpopulations)
    # pihat_superpop_percents, pihat_subpop_percents = get_population_percents(pihat_top_k_subpopulations,
    #                                                                          int(args.k),
    #                                                                          subpopulations)
    # pihat_sample_percents = get_percent_in_top_k(pihat_top_k_samples, subpopulations)
    #
    # kinship_top_k_samples, kinship_top_k_subpopulations = read_top_K(args.kinship_k, subpopulations)
    # kinship_superpop_percents, kinship_subpop_percents = get_population_percents(kinship_top_k_subpopulations,
    #                                                                              int(args.k),
    #                                                                              subpopulations)
    # kinship_sample_percents = get_percent_in_top_k(kinship_top_k_samples, subpopulations)
    #

    # plot_ancestry_top_k(genosis_group_scores,
    #                     genosis_superpop_percents, genosis_subpop_percents,
    #                     dst_superpop_percents, dst_subpop_percents,
    #                     pihat_superpop_percents, pihat_subpop_percents,
    #                     kinship_superpop_percents, kinship_subpop_percents,
    #                     colors,
    #                     args.png_k)

if __name__ == '__main__':
    main()

