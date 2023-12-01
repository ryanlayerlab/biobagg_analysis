import matplotlib.pyplot as plt
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--POP', help='POP (svs score) data file')
    parser.add_argument('--IBD', help='IBD (ground truth) data file')
    parser.add_argument('--png', help='output png file')
    return parser.parse_args()


def main():
    args = parse_args()
    pop_file = args.POP
    ibd_file = args.IBD
    png_file = args.png

    POP_scores = read_data(pop_file)
    IBD_scores = read_data(ibd_file)
    # plot_violin_by_named_relation(POP_scores, IBD_scores, png_file)
    plot_violin_by_number_meiosis(POP_scores, IBD_scores, png_file)

def plot_violin_by_number_meiosis(POP_scores, IBD_scores, png_file):
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
    fig, ax = plt.subplots(figsize=(25, 15), dpi=200)
    # x-axis is the number of meiosis
    for relationship in number_meiosis:
        try:
            POP_data[number_meiosis[relationship]] += POP_scores[relationship]
        except KeyError:
            POP_data[number_meiosis[relationship]] = POP_scores[relationship]

    # plot POP scores
    ax.violinplot([POP_data[col] for col in POP_data], showmeans=True)
    ax.set_xticks(range(1, len(POP_data) + 1))
    ax.set_xticklabels([label for label in meiosis_labels.values()], rotation=0, size=20)
    ax.set_xlabel('Number of Meiosis', fontsize=20)
    ax.set_ylabel('GeSS Score', fontsize=20, color='olivedrab')
    ax.set_title('deCODE Families', fontsize=30)

    for pc in ax.collections:
        pc.set_facecolor('olivedrab')
        pc.set_edgecolor('olivedrab')
        pc.set_alpha(0.6)

    # plot IBD scores for right y-axis
    ax2 = ax.twinx()
    IBD_data = {}

    for relationship in number_meiosis:
        try:
            IBD_data[number_meiosis[relationship]] += IBD_scores[relationship]
        except KeyError:
            IBD_data[number_meiosis[relationship]] = IBD_scores[relationship]

    ax2.violinplot([IBD_data[col] for col in IBD_data], showmeans=True)
    ax2.set_ylabel('deCODE IBD', fontsize=20, color='black')
    # remove first violin plot
    ax2.collections[0].remove()

    for pc in ax2.collections:
        pc.set_facecolor('grey')
        pc.set_edgecolor('black')
        pc.set_alpha(0.4)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)


    plt.savefig(png_file)



def plot_violin_by_named_relation(POP_scores, IBD_scores, png_file):
    ordered_relations = [['self'],
                         ['child', 'parent'],
                         ['sibling'],
                         ['grandchild', 'grandparent'],
                         ['niece-nephew', 'aunt-uncle'],
                         ['great-grandchild', 'great-grandparent'],
                         ['1-cousin'],
                         ['great-niece-nephew', 'great-aunt-uncle'],
                         ['1-cousin-1-removed'],
                         ['2-cousin'],
                         ['unrelated']]
    POP_data = {}
    fig, ax = plt.subplots(figsize=(25, 15), dpi=200)

    for generation in ordered_relations:
        generation_data = []
        gen_lbl = "\n".join([rl for rl in generation])
        for relationship in generation:
            try:
                generation_data += POP_scores[relationship]
            except KeyError:
                generation_data = [0]
                print('no data for: ', relationship)

        POP_data[gen_lbl] = generation_data

    # plot POP scores
    ax.violinplot([POP_data[col] for col in POP_data], showmeans=True)
    ax.set_xticks(range(1, len(POP_data) + 1))
    ax.set_xticklabels([col for col in POP_data], rotation=45)
    ax.set_xlabel('Relationship')
    ax.set_ylabel('GeSS Score', fontsize=20)
    ax.set_title('deCODE Families', fontsize=30)

    for pc in ax.collections:
        pc.set_facecolor('olivedrab')
        pc.set_edgecolor('olivedrab')
        pc.set_alpha(0.6)

    # # plot IBD scores for right y-axis
    # ax2 = ax.twinx()
    # IBD_data = {}
    # for generation in ordered_relations:
    #     generation_data = []
    #     gen_lbl = "\n".join([rl for rl in generation])
    #     for relationship in generation:
    #         try:
    #             generation_data += IBD_scores[relationship]
    #         except KeyError:
    #             generation_data = [0]
    #             print('no data for: ', relationship)
    #
    #     IBD_data[gen_lbl] = generation_data
    #
    # ax2.violinplot([IBD_data[col] for col in IBD_data], showmeans=True)
    # ax2.set_xticks(range(1, len(IBD_data) + 1))
    # ax2.set_xticklabels([col for col in IBD_data], rotation=45)
    # ax2.set_ylabel('deCODE IBD', fontsize=20)
    # # remove first violin plot
    # # ax2.collections[0].remove()
    #
    # for pc in ax2.collections:
    #     pc.set_facecolor('black')
    #     pc.set_edgecolor('black')
    #     pc.set_alpha(0.4)


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    # ax2.spines['left'].set_visible(False)

    plt.savefig(png_file)


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
