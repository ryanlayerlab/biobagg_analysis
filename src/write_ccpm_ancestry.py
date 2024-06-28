import argparse
from collections import defaultdict
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestry', help='ccpm ancestry labels', required=True)
    parser.add_argument('--chrom_hits', help='chromosome hits', required=True)

    return parser.parse_args()

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
    ccpm_top_k_genosis_scores = defaultdict(dict)

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
                        genosis_score = float(line[1])
                        ccpm_top_k[query_id].append(match_id)
                        ccpm_top_k_genosis_scores[query_id][match_id] = genosis_score
            num_files += 1

    return ccpm_top_k, ccpm_top_k_genosis_scores

def write_ccpm_hits_ancestry(ccpm_ancestry,
                             ccpm_top_k,
                             output_file):
    '''
    Write the ancestry of the top k hits for each query to a file
    @param ccpm_top_k_ancestry: dictionary with ccpm_id as key and top k hits ancestry as value
    @param ccpm_ancestry_labels: list of ancestry labels
    @param output_file: path to the output file
    '''
    with open(output_file, 'w') as f:
        f.write('query_id,ancestry\t' + '\thit_id,ancestry\n')
        for query_id, hits in ccpm_top_k.items():
            query_ancestry = ccpm_ancestry[query_id]
            f.write(f'{query_id},{query_ancestry}\t')
            for hit in hits:
                hit_ancestry = ccpm_ancestry[hit]
                f.write(f'{hit},{hit_ancestry}\t')
            f.write('\n')
    f.close()

def write_ccpm_genosis_scores(ccpm_genosis_scores,
                                output_file):
        '''
        Write the genosis scores of the top k hits for each query to a file
        @param ccpm_genosis_scores: dictionary with ccpm_id as key and top k hits genosis scores as value
        @param output_file: path to the output file
        '''
        with open(output_file, 'w') as f:
            f.write('query_id\t' + '\thit_id,genosis_score\n')
            for query_id, hits in ccpm_genosis_scores.items():
                f.write(f'{query_id}\t')
                for hit, genosis_score in hits.items():
                    f.write(f'{hit},{genosis_score}\t')
                f.write('\n')
        f.close()

def main():
    # Parse command line arguments
    args = parse_args()
    ccpm_ancestry_file = args.ancestry
    chrom_hits = args.chrom_hits

    # Check if the files exist
    if not os.path.exists(ccpm_ancestry_file):
        sys.exit(f'Error: {ccpm_ancestry_file} does not exist')
    if not os.path.exists(chrom_hits):
        sys.exit(f'Error: {chrom_hits} does not exist')

    chrm = int(chrom_hits.split('/')[-2].split('_')[-1])

    print('Reading ccpm ancestry file')
    ccpm_ancestry = read_ccpm_ancestry(ccpm_ancestry_file)

    print('Reading chromosome hits')
    ccpm_top_k, ccpm_genosis_scores = get_cohorts(chrom_hits)

    print('writing chromosome hits ancestry to file')
    out_file = str(chrm) + '_ccpm_anc.txt'
    write_ccpm_hits_ancestry(ccpm_ancestry,
                             ccpm_top_k,
                             out_file)

    print('writing chromosome genosis scores to file')
    out_file = str(chrm) + '_ccpm_genosis_scores.txt'
    write_ccpm_genosis_scores(ccpm_genosis_scores,
                             out_file)



if __name__ == '__main__':
    main()