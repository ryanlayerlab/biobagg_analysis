import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
from scipy.stats import ks_2samp
from scipy.stats import norm
import sys

sys.path.append(os.path.abspath('plotting/'))
import ancestry_helpers

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--colors', help='file with color codes', required=True)
    # DECODE DATA
    parser.add_argument('--decode_genosis', help='decode genosis scores', required=True)
    parser.add_argument('--decode_ibd', help='decode ibd scores', required=True)
    # 1KG DATA
    parser.add_argument('--AFR_genosis', help='AFR 1kg genosis trio scores', required=True)
    parser.add_argument('--AMR_genosis', help='AMR 1kg genosis trio scores', required=True)
    parser.add_argument('--EAS_genosis', help='EAS 1kg genosis trio scores', required=True)
    parser.add_argument('--EUR_genosis', help='EUR 1kg genosis trio scores', required=True)
    parser.add_argument('--SAS_genosis', help='SAS 1kg genosis trio scores', required=True)
    parser.add_argument('--AFR_dst', help='AFR 1kg plink DST trio scores')
    parser.add_argument('--AMR_dst', help='AMR 1kg plink DST trio scores')
    parser.add_argument('--EAS_dst', help='EAS 1kg plink DST trio scores')
    parser.add_argument('--EUR_dst', help='EUR 1kg plink DST trio scores')
    parser.add_argument('--SAS_dst', help='SAS 1kg plink DST trio scores')
    parser.add_argument('--AFR_pihat', help='AFR 1kg plink pi-hat trio scores')
    parser.add_argument('--AMR_pihat', help='AMR 1kg plink pi-hat trio scores')
    parser.add_argument('--EAS_pihat', help='EAS 1kg plink pi-hat trio scores')
    parser.add_argument('--EUR_pihat', help='EUR 1kg plink pi-hat trio scores')
    parser.add_argument('--SAS_pihat', help='SAS 1kg plink pi-hat trio scores')
    parser.add_argument('--AFR_kin', help='AFR 1kg plink kinship trio scores')
    parser.add_argument('--AMR_kin', help='AMR 1kg plink kinship trio scores')
    parser.add_argument('--EAS_kin', help='EAS 1kg plink kinship trio scores')
    parser.add_argument('--EUR_kin', help='EUR 1kg plink kinship trio scores')
    parser.add_argument('--SAS_kin', help='SAS 1kg plink kinship trio scores')
    # OUTPUT
    parser.add_argument('--png', help='output png file', required=True)
    return parser.parse_args()

def get_num_meiosis(relationship):
    number_meiosis = {'self': 0,
                      'child': 1,
                      'parent': 1,
                      'sibling': 2,
                      'grandchild': 2,
                      'grandparent': 2,
                      'niece-nephew': 3,
                      'aunt-uncle': 3,
                      'great-grandchild': 3,
                      'great-grandparent': 3,
                      '1-cousin': 4,
                      'great-niece-nephew': 4,
                      'great-aunt-uncle': 4,
                      '1-cousin-1-removed': 5,
                      '2-cousin': 6,
                      'unrelated': 10,}
    return number_meiosis[relationship]

def read_decode_data(decode_data_file):
    '''
    Read file with decode data scores
    @param decode_data_file: decode genosis scores or ibd scores
    @return: dictionary where key=number of meoisis, and value=list of scores
    '''
    decode_scores = {}
    with open(decode_data_file, 'r') as f:
        for line in f:
            line = line.strip().split(',')
            relationship_label = line[0]
            # ignore self and undetermined
            if relationship_label == 'self' or relationship_label == 'undetermined':
                continue
            num_meiosis = get_num_meiosis(relationship_label)
            scores = [float(x) for x in line[1:]]
            try:
                decode_scores[num_meiosis].extend(scores)
            except KeyError:
                decode_scores[num_meiosis] = scores

    return decode_scores


def read_1kg_data(scores_file, population):
    '''
    Read file with 1kg data scores
    @param scores_file: 1kg scores file (genosis, dst, pihat, or kin)
    @return: dictionary where key=number of meiosis, and value=list of scores
    '''
    population_scores = {}
    with open(scores_file, 'r') as f:
        # read header
        header = f.readline().strip().split()
        for line in f:
            line = line.strip().split()
            relationship_label = line[0]
            # ignore self, subpop, and pop
            if relationship_label == 'self' or relationship_label == 'subpop' or relationship_label == population:
                continue
            elif relationship_label == 'child' or relationship_label == 'parent':
                num_meiosis = get_num_meiosis(relationship_label)
                population_scores[num_meiosis] = [float(x) for x in line[1:]]
            else:
                relationship_label = 'unrelated'
                num_meiosis = get_num_meiosis(relationship_label)
                try:
                    population_scores[num_meiosis].extend([float(x) for x in line[1:]])
                except KeyError:
                    population_scores[num_meiosis] = [float(x) for x in line[1:]]
    return population_scores

def combine_1kg_pop_scores(AFR_scores, AMR_scores, EAS_scores, EUR_scores, SAS_scores):
    '''
    Combines scores (genosis, dst, pihat, or kin) from 5 populations
    @param AFR_scores: AFR scores file (genosis, dst, pihat, or kin)
    @param AMR_scores: AMR scores file (genosis, dst, pihat, or kin)
    @param EAS_scores: EAS scores file (genosis, dst, pihat, or kin)
    @param EUR_scores: EUR scores file (genosis, dst, pihat, or kin)
    @param SAS_scores: SAS scores file (genosis, dst, pihat, or kin)
    @return: dictionary where key=number of meiosis, and value=list of scores
    '''
    combined_scores = {}
    # only parent-child and unrelated
    for num_meiosis in [1,10]:
        combined_scores[num_meiosis] = (AFR_scores[num_meiosis] +
                                        AMR_scores[num_meiosis] +
                                        EAS_scores[num_meiosis] +
                                        EUR_scores[num_meiosis] +
                                        SAS_scores[num_meiosis])
    return combined_scores



def plot_combined_figures(decode_ibd_data,
                          decode_genosis_data,
                          tg_genosis_scores,
                          tg_dst_scores,
                          plink_label,
                          colors,
                          png_file):
    '''
    Show all plots on one figure
    @param decode_ibd_data: decode ibd dictionary of scores
    @param decode_genosis_data: decode genosis dictionary of scores
    @return: none (saves png)
    '''

    # Create the combined figure
    combined_figure, axes = plt.subplots(2, 2, figsize=(10, 5), dpi=300)
    alpha_value = 0.2

    # row0,col0 = deocde ibd scores
    decode_IBD_plt = axes[0][0]
    for num_meiosis in sorted(decode_ibd_data.keys()):
        if num_meiosis == 10:
            label = 'unrelated'
        else:
            label = str(num_meiosis)
        sns.kdeplot(decode_ibd_data[num_meiosis],
                    color=colors[num_meiosis],
                    ax=decode_IBD_plt,
                    label=label,
                    fill=True, alpha=alpha_value)
    # formatting
    decode_IBD_plt.spines['top'].set_visible(False)
    decode_IBD_plt.spines['right'].set_visible(False)
    decode_IBD_plt.set_yticks(np.arange(0, 0.15, 0.025))
    decode_IBD_plt.set_xlabel('deCODE IBD')
    decode_IBD_plt.set_ylabel('Density')


    # row1,col0 = decode genosis scores
    decode_genosis_plt = axes[1][0]
    for num_meiosis in sorted(decode_genosis_data.keys()):
        if num_meiosis == 10:
            label = 'unrelated'
        else:
            label = str(num_meiosis)
        sns.kdeplot(decode_genosis_data[num_meiosis],
                    color=colors[num_meiosis],
                    ax=decode_genosis_plt,
                    label=label,
                    fill=True, alpha=alpha_value)
    # formatting
    decode_genosis_plt.spines['top'].set_visible(False)
    decode_genosis_plt.spines['right'].set_visible(False)
    decode_genosis_plt.set_yticks(np.arange(0, 0.15, 0.025))
    decode_genosis_plt.set_xlabel('GenoSiS Score')
    decode_genosis_plt.set_ylabel('Density')


    # row0,col1 = plink scores
    tg_plink_plt = axes[0][1]
    for num_meiosis in sorted(tg_dst_scores.keys()):
        if num_meiosis == 10:
            label = 'unrelated'
        else:
            label = str(num_meiosis)
        sns.kdeplot(tg_dst_scores[num_meiosis],
                    color=colors[num_meiosis],
                    ax=tg_plink_plt,
                    label=label,
                    fill=True, alpha=alpha_value)
    # formatting
    tg_plink_plt.spines['top'].set_visible(False)
    tg_plink_plt.spines['right'].set_visible(False)
    tg_plink_plt.set_xlabel(plink_label)
    tg_plink_plt.set_ylabel('Density')

    # row1,col1 = 1kg genosis scores
    tg_genosis_plt = axes[1][1]
    for num_meiosis in sorted(tg_genosis_scores.keys()):
        if num_meiosis == 10:
            label = 'unrelated'
        else:
            label = str(num_meiosis)
        sns.kdeplot(tg_genosis_scores[num_meiosis],
                    color=colors[num_meiosis],
                    ax=tg_genosis_plt,
                    label=label,
                    fill=True, alpha=alpha_value)
    # formatting
    tg_genosis_plt.spines['top'].set_visible(False)
    tg_genosis_plt.spines['right'].set_visible(False)
    tg_genosis_plt.set_xlabel('GenoSiS Score')
    tg_genosis_plt.set_ylabel('Density')


    # title and legend used for left plots
    decode_IBD_plt.set_title('deCODE Family Data', fontsize=14, fontweight='bold')
    decode_IBD_plt.legend(title='Number of Meiosis',
                          fontsize=8,
                          title_fontsize=8,
                          frameon=False,
                          loc='upper right',
                          ncol=2)

    # title and legend used for right plots
    tg_plink_plt.set_title('TGP Trio Data', fontsize=14, fontweight='bold')
    tg_plink_plt.legend(title='Number of Meiosis',
                              fontsize=8,
                              title_fontsize=8,
                              frameon=False,
                              loc='upper center',
                              ncol=2)


    # Save the figure
    plt.tight_layout()
    combined_figure.savefig(png_file)


def main():
    args = parse_args()
    colors = ancestry_helpers.get_colors(args.colors)

    # read deCODE data
    decode_genosis_scores = read_decode_data(args.decode_genosis)
    decode_ibd_scores = read_decode_data(args.decode_ibd)
    # divide IBD by 2
    for num_meiosis in decode_ibd_scores.keys():
        decode_ibd_scores[num_meiosis] = [x/2 for x in decode_ibd_scores[num_meiosis]]

    # 1kg data
    ## genosis scores
    AFR_genosis_scores = read_1kg_data(args.AFR_genosis, 'AFR')
    AMR_genosis_scores = read_1kg_data(args.AMR_genosis, 'AMR')
    EAS_genosis_scores = read_1kg_data(args.EAS_genosis, 'EAS')
    EUR_genosis_scores = read_1kg_data(args.EUR_genosis, 'EUR')
    SAS_genosis_scores = read_1kg_data(args.SAS_genosis, 'SAS')
    tg_genosis_scores = combine_1kg_pop_scores(AFR_genosis_scores,
                                                   AMR_genosis_scores,
                                                   EAS_genosis_scores,
                                                   EUR_genosis_scores,
                                                   SAS_genosis_scores)
    ## plink dst scores
    AFR_dst_scores = read_1kg_data(args.AFR_dst, 'AFR')
    AMR_dst_scores = read_1kg_data(args.AMR_dst, 'AMR')
    EAS_dst_scores = read_1kg_data(args.EAS_dst, 'EAS')
    EUR_dst_scores = read_1kg_data(args.EUR_dst, 'EUR')
    SAS_dst_scores = read_1kg_data(args.SAS_dst, 'SAS')
    tg_dst_scores = combine_1kg_pop_scores(AFR_dst_scores,
                                           AMR_dst_scores,
                                           EAS_dst_scores,
                                           EUR_dst_scores,
                                           SAS_dst_scores)
    # plink pi-hat scores
    AFR_pihat_scores = read_1kg_data(args.AFR_pihat, 'AFR')
    AMR_pihat_scores = read_1kg_data(args.AMR_pihat, 'AMR')
    EAS_pihat_scores = read_1kg_data(args.EAS_pihat, 'EAS')
    EUR_pihat_scores = read_1kg_data(args.EUR_pihat, 'EUR')
    SAS_pihat_scores = read_1kg_data(args.SAS_pihat, 'SAS')
    tg_pihat_scores = combine_1kg_pop_scores(AFR_pihat_scores,
                                             AMR_pihat_scores,
                                             EAS_pihat_scores,
                                             EUR_pihat_scores,
                                             SAS_pihat_scores)
    # plink kinship scores
    AFR_kin_scores = read_1kg_data(args.AFR_kin, 'AFR')
    AMR_kin_scores = read_1kg_data(args.AMR_kin, 'AMR')
    EAS_kin_scores = read_1kg_data(args.EAS_kin, 'EAS')
    EUR_kin_scores = read_1kg_data(args.EUR_kin, 'EUR')
    SAS_kin_scores = read_1kg_data(args.SAS_kin, 'SAS')
    tg_kin_scores = combine_1kg_pop_scores(AFR_kin_scores,
                                           AMR_kin_scores,
                                           EAS_kin_scores,
                                           EUR_kin_scores,
                                           SAS_kin_scores)

    plink_label = 'King-robust coefficient'
    # plot all figures
    plot_combined_figures(decode_ibd_scores,
                          decode_genosis_scores,
                          tg_genosis_scores,
                          tg_kin_scores,
                          plink_label,
                          colors,
                          args.png)



if __name__ == '__main__':
    main()