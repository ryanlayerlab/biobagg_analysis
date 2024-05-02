import argparse
from collections import defaultdict

import read_plink as rp
import plotting.ancestry_helpers as ah

def parse_args():
    parser = argparse.ArgumentParser(description="Writes scores for top K by relatedness label")
    parser.add_argument("-i", "--input", help="input top hits file")
    parser.add_argument("-s", "--score", help="type of socre (e.g. genosis, plink, etc.)")
    parser.add_argument("-p", "--pop", help="population")
    parser.add_argument("-k", "--knn", help="value of K for hits")
    parser.add_argument("-a", "--ancestry", help="ancestry file")
    parser.add_argument("-o", "--out", help="output dir")
    return parser.parse_args()

def read_top_hits(top_hits_file):
    top_hits_dict = defaultdict(dict)

    with open(top_hits_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            query = line[0]
            top_hits_list = line[1:]
            for hit in top_hits_list:
                match, score = hit.split(',')
                try:
                    score = float(score)
                except ValueError:
                    continue
                try:
                    top_hits_dict[query].update({match: score})
                except KeyError:
                    top_hits_dict[query] = {match: score}

    return top_hits_dict

def read_relationship_labels(pop_labels_file):
    pop_labels = defaultdict(dict)
    with open(pop_labels_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            sample1 = line[0]
            sample2 = line[1]
            relationship = line[3]
            try:
                pop_labels[sample1][sample2] = relationship
            except KeyError:
                pop_labels[sample1] = {}
                pop_labels[sample1][sample2] = relationship
            try:
                pop_labels[sample2][sample1] = relationship
            except KeyError:
                pop_labels[sample2] = {}
                pop_labels[sample2][sample1] = relationship
    return pop_labels

def get_relatedness_dict(top_hits_dict, pop_labels, subpopulations, pop):
    # get population labels for each hit
    relationship_options = ['self', 'parent', 'child', 'subpop', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']

    pop_labels_scores = {label: [] for label in relationship_options}
    for query in top_hits_dict:
        query_superpop = ah.SUB_SUPERPOPULATIONS[subpopulations[query]]
        if query_superpop == pop:
            for match in top_hits_dict[query]:
                match_superpop = ah.SUB_SUPERPOPULATIONS[subpopulations[match]]
                score = top_hits_dict[query][match]
                try:
                    relationship_label = pop_labels[query][match]
                except KeyError:
                    relationship_label = 'outpop'

                if relationship_label not in relationship_options:
                    relationship_label = match_superpop
                try:
                    pop_labels_scores[relationship_label].append(score)
                except KeyError:
                    pop_labels_scores[relationship_label] = [score]
        else:
            continue
    return pop_labels_scores



def write_relationship_scores(relatedness_dict, score_type, output_dir, pop, k):
    output_file = output_dir + "/" + score_type + "_" + pop + "_trio_scores_" + k + ".txt"
    f = open(output_file, 'w')
    with open(output_file, 'w') as f:
        f.write("relationship score\n")
        for relationship in relatedness_dict:
            f.write(f"{relationship} ")
            for score in relatedness_dict[relationship]:
                f.write(f"{score} ")
            f.write("\n")
    f.close()

def main():
    args = parse_args()

    pop_labels_file = 'data/1KG_trios_' + args.pop + '.txt'
    subpopulations = ah.get_subpopulations(args.ancestry)

    top_hits_dict = read_top_hits(args.input)
    pop_labels = read_relationship_labels(pop_labels_file)

    relatedness_dict = get_relatedness_dict(top_hits_dict, pop_labels, subpopulations, args.pop)

    write_relationship_scores(relatedness_dict, args.score, args.out, args.pop, args.knn)

if __name__ == "__main__":
    main()
