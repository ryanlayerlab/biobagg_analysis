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
    plot_violin(POP_scores, IBD_scores, png_file)
	
def plot_violin(POP_scores, IBD_scores, png_file):
    
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
    
    plot_data = {}
	
    fig, ax = plt.subplots(figsize=(25,13))

    for generation in ordered_relations:
        generation_data = []
        gen_lbl = "\n".join([rl for rl in generation])
        for relationship in generation:
            try:
                generation_data += POP_scores[relationship]
            except KeyError:
                generation_data = [0]
                print('no data for: ', relationship)

        plot_data[gen_lbl] = generation_data

    print(len(plot_data))
    ax.violinplot([plot_data[col] for col in plot_data])
    ax.set_xticks(range(1,len(plot_data)+1))
    ax.set_xticklabels([col for col in plot_data], rotation=45)
    ax.set_xlabel('Relationship')
    ax.set_ylabel('Aggregate SimScore')
    ax.set_title('deCODE pedigree', fontsize=20)

    for pc in ax.collections:
        pc.set_facecolor('olivedrab')
        pc.set_edgecolor('olivedrab')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

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
			score = [float(line[i]) for i in range(1,len(line))]
			if relationship not in relationship_scores:
				relationship_scores[relationship] = [score]
			relationship_scores[relationship].append(score)
	return relationship_scores

if __name__ == '__main__':
	main()
	