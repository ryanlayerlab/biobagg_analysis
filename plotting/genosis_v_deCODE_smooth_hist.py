import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from scipy.stats import ks_2samp
from scipy.stats import norm
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--POP', help='POP (svs score) data file')
    parser.add_argument('--IBD', help='IBD (ground truth) data file')
    parser.add_argument('--png', help='output png file')
    parser.add_argument('--height',
                        type=float,
                        default=4,
                        help='height of the violin plot')
    parser.add_argument('--width',
                        type=float,
                        default=6,
                        help='width of the violin plot')
    return parser.parse_args()


def main():
    args = parse_args()
    pop_file = args.POP
    ibd_file = args.IBD
    png_file = args.png

    POP_scores = read_data(pop_file)
    IBD_scores = read_data(ibd_file)
    # plot_violin_by_named_relation(POP_scores, IBD_scores, png_file)
    plot_violin_by_number_meiosis(POP_scores,
                                  IBD_scores,
                                  png_file,
                                  args.height,
                                  args.width)


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
                      'unrelated': 10}
    return number_meiosis[relationship]


def normalize(data, min_val=None, max_val=None):
    # Find the minimum and maximum values
    if min_val is None:
        min_val = min(data)
    if max_val is None:
        max_val = max(data)

    normalized_data = [(x - min_val) / (max_val - min_val) for x in data]

    return normalized_data


def plot_violin_by_number_meiosis(POP_scores,
                                  IBD_scores,
                                  png_file,
                                  height,
                                  width):

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
                      'unrelated': 10}

    meiosis_labels = {0: '0',
                        1: '1',
                        2: '2',
                        3: '3',
                        4: '4',
                        5: '5',
                        6: '6',
                        10: 'unrelated'}

    POP_data = {}
    IBD_data = {}

    for relationship in number_meiosis:
        try:
            POP_data[number_meiosis[relationship]] += POP_scores[relationship]
        except KeyError:
            POP_data[number_meiosis[relationship]] = POP_scores[relationship]
        try:
            IBD_data[number_meiosis[relationship]] += IBD_scores[relationship]
        except KeyError:
            IBD_data[number_meiosis[relationship]] = IBD_scores[relationship]

    min_pop = min([min(POP_data[num_meiosis]) for num_meiosis in POP_data])
    max_pop = max([max(POP_data[num_meiosis]) for num_meiosis in POP_data])
    fig, axs = plt.subplots(2,1,figsize=(width,height) , dpi=200)
    ax = axs[0]
    for num_meiosis in POP_data:
        if num_meiosis == 0: continue
        sns.kdeplot(POP_data[num_meiosis],
                    ax=ax,
                    label=meiosis_labels[num_meiosis],
                    linewidth=1,
                    fill=True,
                    alpha=0.1)
        #ax.hist(POP_data[num_meiosis],
                #bins=20,
                #alpha=0.2,
                #label=meiosis_labels[num_meiosis],
                #linewidth=2,
                #density=True)


    ax.set_xlabel('GenoSiS Score')
    ax.set_ylabel('Frequency')
    ax.legend(fontsize=6,
              title='Number of Meiosis',
              title_fontsize=6,
              frameon=False,
              loc='upper right', 
              ncol=2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yticks([0.0, 0.05])
    ax.set_yticklabels(['0.0', '0.05'])

    min_ibd = min([min(IBD_data[num_meiosis]) for num_meiosis in IBD_data])
    max_ibd = max([max(IBD_data[num_meiosis]) for num_meiosis in IBD_data])
    ax = axs[1]
    for num_meiosis in IBD_data:
        if num_meiosis == 0: continue
        sns.kdeplot(IBD_data[num_meiosis],
                    ax=ax,
                    label=meiosis_labels[num_meiosis],
                    linewidth=1,
                    fill=True,
                    alpha=0.2)
 
    ax.set_xlabel('deCODE IBD')
    ax.set_ylabel('Frequency')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yticks([0.0, 0.05])
    ax.set_yticklabels(['0.0', '0.5'])

    plt.tight_layout()
    plt.savefig(png_file)

    plt.savefig(png_file)

    for num_meiosis in POP_data:

        POP_norm = normalize(POP_data[num_meiosis],
                             min_val=min_pop,
                             max_val=max_pop)
        IBD_norm = normalize(IBD_data[num_meiosis],
                             min_val=min_ibd,
                             max_val=max_ibd)

        ztest = two_sample_z_test(POP_norm, IBD_norm)
        kstest = ks_2samp(POP_norm, IBD_norm)
        print(num_meiosis,
              ztest,
              kstest)

def two_sample_z_test(sample1, sample2):
    # Calculate sample means
    mean1 = np.mean(sample1)
    mean2 = np.mean(sample2)

    std1 = np.std(sample1)
    std2 = np.std(sample2)

    # Calculate sample sizes
    n1 = len(sample1)
    n2 = len(sample2)

    # Calculate the standard error of the difference between the means
    se = np.sqrt((std1**2 / n1) + (std2**2 / n2))

    # Calculate the Z-score
    z_score = (mean1 - mean2) / se

    # Calculate the p-value
    p_value = 2 * norm.sf(np.abs(z_score))  # Two-tailed test

    return z_score, p_value


def read_data(data_file):
    """
	reads results for deCODDE data
	@ input: data_file (str) - path to data file
	@ output: data (dict) - dictionary of data
	"""
    relationship_scores = {}
    with open(data_file, 'r') as f:
        for line in f:
            line = line.strip().split(',')
            relationship = line[0]
            scores = [float(line[i]) for i in range(1, len(line))]
            if relationship not in relationship_scores:
                relationship_scores[relationship] = scores
    return relationship_scores


if __name__ == '__main__':
    main()
