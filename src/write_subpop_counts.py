import argparse
from collections import defaultdict

import plotting.ancestry_helpers as ancestry_helpers


def parse_args():
    parser = argparse.ArgumentParser(description='Write subpopulation counts to file')
    parser.add_argument('--ancestry', help='1KG ancestry labels', required=True)
    parser.add_argument('--genosis', type=str, required=True, help='Path to genosis cohort file')
    parser.add_argument('--output_dir', type=str, required=True, help='Path to output directory')
    return parser.parse_args()

def get_subpop_counts(top_k_file,
                      subpopulations,
                      header=True):
    '''
    Read top K scores and return top K samples, and populations
    @param top_k_file: path top k scores
    @param subpopulations: dictionary of subpopulations for each sample
    @return: dictionary of subpopulation counts
    '''
    subpopulation_counts = defaultdict(dict)

    f = open(top_k_file, 'r')
    if header:
        f.readline()
    for line in f:
        line = line.strip().split()
        query = line[0]
        matches = line[1:]
        query_subpopulation = subpopulations[query]
        query_dict = {}
        for match_score in matches:
            match = match_score.split(',')[0]
            # ignore self match
            if query in match:
                continue
            match_subpopulation = subpopulations[match]
            try:
                subpopulation_counts[query_subpopulation][match_subpopulation] += 1
            except KeyError:
                subpopulation_counts[query_subpopulation][match_subpopulation] = 1
    return subpopulation_counts

def write_subpop_counts(subpopulation_counts, output_file):
    '''
    Write subpopulation counts to file
    @param subpopulation_counts: dictionary of subpopulation counts
    @param output_file: path to output file
    '''
    all_subpopulations = sorted(subpopulation_counts.keys())
    with open(output_file, 'w') as f:
        f.write('query_subpop\tmatch_subpop,count...\n')
        for query_subpopulation in all_subpopulations:
            f.write(f'{query_subpopulation}\t')
            for match_subpopulation in all_subpopulations:
                try:
                    count = subpopulation_counts[query_subpopulation][match_subpopulation]
                except KeyError:
                    count = 0
                f.write(f'{match_subpopulation},{count}\t')
            f.write('\n')
    f.close()

def main():
    args = parse_args()

    subpopulations = ancestry_helpers.get_subpopulations(args.ancestry)
    genosis_subpop_counts = get_subpop_counts(args.genosis, subpopulations, header=False)
    output_file = args.output_dir + 'genosis_subpop_counts.tsv'
    write_subpop_counts(genosis_subpop_counts, output_file)


if __name__ == '__main__':
    main()
