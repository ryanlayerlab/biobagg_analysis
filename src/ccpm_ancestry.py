import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestry', help='ccpm ancestry labels', required=True)
    parser.add_argument('--chrom_hits', help='chromosome hits', required=True)
    parser.add_argument('--png', help='output png file', required=True)

    return parser.parse_args()


def get_ccpm_ancestry_labels(ccpm_ancestry_file):
    '''
    Read ccpm ancestry file and get a list of unique ancestry values
    @param ccpm_ancestry_file: path to the ccpm ancestry file
    @return: list of unique ancestry values
    '''
    ancestry_labels = set()
    with open(ccpm_ancestry_file, 'r') as f:
        # Skip header
        f.readline()
        for line in f:
            ancestry = line.strip().split('\t')[1]
            ancestry_labels.add(ancestry)

    return list(ancestry_labels)

def read_ccpm_ancestry(ccpm_ancestry_file):
    '''
    Read the ccpm ancestry file and return a dictionary with ccpm_id as key and ancestry as value
    @param ccpm_ancestry_file: path to the ccpm ancestry file
    @return: dictionary with ccpm_id as key and ancestry as value
    '''
    ccpm_ancestry = dict()

    with open(ccpm_ancestry_file, 'r') as f:
        # Skip header
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            ccpm_id = line[0]
            ancestry = line[1]
            ccpm_ancestry[ccpm_id] = ancestry

    return ccpm_ancestry

def get_cohorts(chrom_hits):
    '''
    Read all the query's hits and return a dictionary with ccpm_id as key and top k hits as value
    @param chrom_hits: path to the chromosome hits dir
    @return: dictionary with ccpm_id as key and top k hits as value
    '''
    num_files = 0
    ccpm_top_k = defaultdict(list)
    for file in os.listdir(chrom_hits):
        if num_files == 100:
            break
        if file.endswith('.knn'):
            query_id = file.split('.')[0]
            with open(os.path.join(chrom_hits, file), 'r') as f:
                for line in f:
                    if 'QUERY' in line:
                        # check that the query id is the same as the file name
                        assert query_id == line.strip().split(':')[1].strip()
                    else:
                        line = line.strip().split('\t')
                        match_id = line[0]
                        ccpm_top_k[query_id].append(match_id)
            num_files += 1

    return ccpm_top_k

def get_cohort_ancestry(ccpm_top_k, ccpm_ancestry):
    '''
    Get the ancestry of the top k hits for each query
    @param ccpm_top_k: dictionary with ccpm_id as key and top k hits as value
    @param ccpm_ancestry: dictionary with ccpm_id as key and ancestry as value
    @return: dictionary with ccpm_id as key and top k hits ancestry as value
    '''
    ccpm_top_k_ancestry = defaultdict(list)
    for query_id, hits in ccpm_top_k.items():
        for hit in hits:

            try:
                hit_ancestry = ccpm_ancestry[hit]
            except KeyError:
                print(f'Error: {hit} not found in ccpm ancestry file')
                continue
            ccpm_top_k_ancestry[query_id].append(hit_ancestry)

    return ccpm_top_k_ancestry

def write_ccpm_hits_ancestry(ccpm_top_k_ancestry, ccpm_ancestry_labels, output_file):
    '''
    Write the ancestry of the top k hits for each query to a file
    @param ccpm_top_k_ancestry: dictionary with ccpm_id as key and top k hits ancestry as value
    @param ccpm_ancestry_labels: list of ancestry labels
    @param output_file: path to the output file
    '''
    with open(output_file, 'w') as f:
        f.write('ccpm_id\t' + '\t'.join(ccpm_ancestry_labels) + '\n')
        for query_id, hits in ccpm_top_k_ancestry.items():
            f.write(query_id + '\t' + '\t'.join(hits) + '\n')

def plot_ccpm_hits(ccpm_top_k_ancestry,
                   ccpm_ancestry_labels,
                   png_file):
    '''
    Plot the ancestry of the top k hits for each query
    @param ccpm_top_k_ancestry: dictionary with ccpm_id as key and top k hits ancestry as value
    '''
    # heatmap with ancestry labels on the both axes
    sns.set()
    fig, ax = plt.subplots()
    ancestry_matrix = np.zeros((len(ccpm_top_k_ancestry), len(ccpm_ancestry_labels)))
    for i, (query_id, hits) in enumerate(ccpm_top_k_ancestry.items()):
        for hit in hits:
            ancestry_matrix[i][ccpm_ancestry_labels.index(hit)] += 1

    sns.heatmap(ancestry_matrix, annot=True, fmt='g', xticklabels=ccpm_ancestry_labels, yticklabels=ccpm_top_k_ancestry.keys())
    plt.show()
    plt.savefig(png_file)


def main():
    # Parse command line arguments
    args = parse_args()
    ccpm_ancestry_file = args.ancestry
    chrom_hits = args.chrom_hits
    png_file = args.png

    # Check if the files exist
    if not os.path.exists(ccpm_ancestry_file):
        sys.exit(f'Error: {ccpm_ancestry_file} does not exist')
    if not os.path.exists(chrom_hits):
        sys.exit(f'Error: {chrom_hits} does not exist')

    ccpm_ancestry_labels = get_ccpm_ancestry_labels(ccpm_ancestry_file)

    print('Reading ccpm ancestry file')
    ccpm_ancestry = read_ccpm_ancestry(ccpm_ancestry_file)
    print('Reading chromosome hits')
    ccpm_top_k = get_cohorts(chrom_hits)
    print('Getting ancestry of top k hits')
    ccpm_top_k_ancestry = get_cohort_ancestry(ccpm_top_k, ccpm_ancestry)
    print('writing chromosome hits ancestry to file')
    write_ccpm_hits_ancestry(ccpm_top_k_ancestry, ccpm_ancestry_labels, 'ccpm_hits_ancestry.txt')

    # print('Plotting ccpm hits')
    # plot_ccpm_hits(ccpm_top_k_ancestry, ccpm_ancestry_labels, png_file)


if __name__ == '__main__':
    main()