import argparse
from collections import defaultdict

import plotting.ancestry_helpers as ancestry_helpers


def parse_args():
    parser = argparse.ArgumentParser(description='Write subpopulation counts to file')
    parser.add_argument('--ancestry', help='1KG ancestry labels', required=True)
    parser.add_argument('--genosis', type=str, required=True, help='Path to genosis cohort file')
    parser.add_argument('--dst', type=str, required=True, help='Path to plink dst top 20 file')
    parser.add_argument('--pihat', type=str, required=True, help='Path to plink pi-hat top 20 file')
    parser.add_argument('--kinship', type=str, required=True, help='Path to plink kinship top 20 file')
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
    all_subpopulations = ancestry_helpers.SUBPOPULATIONS
    subpopulation_counts = {sub: {sub2: [] for sub2 in all_subpopulations} for sub in all_subpopulations}

    f = open(top_k_file, 'r')
    if header:
        f.readline()
    for line in f:
        line = line.strip().split()
        query = line[0]
        matches = line[1:]
        query_subpopulation = subpopulations[query]
        query_dict = {subpop: 0 for subpop in all_subpopulations}
        for match_score in matches:
            match = match_score.split(',')[0]
            # ignore self match
            if query in match:
                continue
            match_subpopulation = subpopulations[match]
            query_dict[match_subpopulation] += 1
        for subpop in query_dict:
            subpopulation_counts[query_subpopulation][subpop].append(query_dict[subpop])

    f.close()
    return subpopulation_counts

def write_subpop_counts(subpopulation_counts, output_file):
    '''
    Write subpopulation counts to file
    @param subpopulation_counts: dictionary of subpopulation counts
    @param output_file: path to output file
    '''
    all_subpopulations = sorted(subpopulation_counts.keys())
    with open(output_file, 'w') as f:
        f.write('query_subpop match_subpop counts,...\n')
        for query_subpopulation in all_subpopulations:
            # f.write(f'{query_subpopulation}\t')
            for match_subpopulation in all_subpopulations:
                counts = subpopulation_counts[query_subpopulation][match_subpopulation]
                f.write(f'{query_subpopulation}\t{match_subpopulation}\t{",".join(map(str, counts))}\t')
                f.write('\n')
    f.close()

def main():
    args = parse_args()

    subpopulations = ancestry_helpers.get_subpopulations(args.ancestry)

    genosis_subpop_counts = get_subpop_counts(args.genosis, subpopulations, header=False)
    genosis_output = args.output_dir + 'genosis_counts.tsv'
    write_subpop_counts(genosis_subpop_counts, genosis_output)

    dst_subpop_counts = get_subpop_counts(args.dst, subpopulations)
    dst_output = args.output_dir + 'dst_counts.tsv'
    write_subpop_counts(dst_subpop_counts, dst_output)

    pihat_subpop_counts = get_subpop_counts(args.pihat, subpopulations)
    pihat_output = args.output_dir + 'pihat_counts.tsv'
    write_subpop_counts(pihat_subpop_counts, pihat_output)

    kinship_subpop_counts = get_subpop_counts(args.kinship, subpopulations)
    kinship_output = args.output_dir + 'kinship_counts.tsv'
    write_subpop_counts(kinship_subpop_counts, kinship_output)


if __name__ == '__main__':
    main()
