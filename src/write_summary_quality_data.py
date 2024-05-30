import argparse
import os
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--ancestry', type=str, help='ancestry file', required=True)
    parser.add_argument('-d', '--data', type=str, help='data directory', required=True)
    return parser.parse_args()

SUBPOPULATIONS = ['ASW', 'LWK', 'GWD', 'MSL', 'ESN', 'YRI', 'ACB', 
                  'CLM', 'PEL', 'MXL', 'PUR',
                  'CDX', 'CHB', 'JPT', 'KHV', 'CHS',
                  'CEU', 'TSI', 'FIN', 'GBR', 'IBS',
                  'BEB', 'GIH', 'ITU', 'PJL', 'STU']

SUB_SUPERPOPULATIONS = {'ASW': 'AFR', 'LWK': 'AFR', 'GWD': 'AFR', 'MSL': 'AFR', 'ESN': 'AFR', 'YRI': 'AFR', 'ACB': 'AFR',
                        'CLM': 'AMR', 'PEL': 'AMR', 'MXL': 'AMR', 'PUR': 'AMR',
                        'CDX': 'EAS', 'CHB': 'EAS', 'JPT': 'EAS', 'KHV': 'EAS', 'CHS': 'EAS', 
                        'CEU': 'EUR', 'TSI': 'EUR', 'FIN': 'EUR', 'GBR': 'EUR', 'IBS': 'EUR', 
                        'BEB': 'SAS', 'GIH': 'SAS', 'ITU': 'SAS', 'PJL': 'SAS', 'STU': 'SAS'}

SUPER_SUBPOPULATIONS = {'AFR': ['ASW', 'LWK', 'GWD', 'MSL', 'ESN', 'YRI', 'ACB'],
                        'AMR': ['CLM', 'PEL', 'MXL', 'PUR'],
                        'EAS': ['CDX', 'CHB', 'JPT', 'KHV', 'CHS'],
                        'EUR': ['CEU', 'TSI', 'FIN', 'GBR', 'IBS'],
                        'SAS': ['BEB', 'GIH', 'ITU', 'PJL', 'STU']}

def get_subpopulations(ancestry_file):
    sample_subpopulations = {}
    with open(ancestry_file) as f:
        for line in f:
            for subpop in SUBPOPULATIONS:
                if subpop in line:
                    sample_subpopulations[line.split()[0]] = subpop
    return sample_subpopulations

def assert_sample_pop(sample, pop, sample_subpopulations):
    sample_subpop = sample_subpopulations[sample]
    sample_pop = SUB_SUPERPOPULATIONS[sample_subpop]
    if sample_pop != pop:
        #print('!!incorrect population label!!', sample, pop)
        return 0
    else:
        return 1

def read_top_hits(query_file,
                  query_full_dict,
                  database_pop,
                  sample_subpopulations):
    f = open(query_file, 'r')
    header = f.readline()
    query_subpopulation = sample_subpopulations[header.strip().split(':')[1].strip()]
    query_superpopulation = SUB_SUPERPOPULATIONS[query_subpopulation]
    for line in f:
        L = line.strip().split(',')
        hit_ID = L[0]
        # if self, exit early
        if hit_ID == header.strip().split(':')[1].strip():
            continue
        hit_score = float(L[1])
        # get match of hit
        hit_subpopulation = sample_subpopulations[hit_ID]
        # get superpopulation of hit
        hit_superpopulation = SUB_SUPERPOPULATIONS[hit_subpopulation]
        if hit_superpopulation != database_pop:
            print('!!incorrect population label!!', hit_ID, database_pop)
        # assert hit is in expected superpopulation
        assert_sample_pop(hit_ID, database_pop, sample_subpopulations)

        # add score to subpopulation list
        try:
            query_full_dict[database_pop][hit_subpopulation].append(hit_score)
        except KeyError:
            try:
                query_full_dict[database_pop][hit_subpopulation] = [hit_score]
            except KeyError:
                query_full_dict[database_pop] = {hit_subpopulation : [hit_score]}

    #print(query_scores[db_pop])
    f.close()
    return query_full_dict

def write_query_pop_scores(query_pop,
                           query_full_dict):

    o_file = open(query_pop + '.txt', 'w')
    header = 'QUERY POPULATION: ' + query_pop + '\n'
    o_file.write(header)

    for database_pop in query_full_dict.keys():
        
        o_file.write('database: ' + database_pop + '\n')
        for subpop in query_full_dict[database_pop]:
            line = subpop + ':'
            scores = query_full_dict[database_pop][subpop]
            if len(scores) == 0:
                continue
            else:
                for score in scores:
                    line += str(score) + ','
            o_file.write(line + '\n')
    o_file.close()
    

def main():

    args = get_args()
    ancestry_file = args.ancestry
    data_dir = args.data

    sample_subpopulations = get_subpopulations(ancestry_file)

    super_populations = {
                        'AFR': 'African',
                        'AMR': 'American',
                        #'EAS': 'East Asian',
                        'EUR': 'European',
                        'SAS': 'South Asian'}

    # queries were performed by superpopulation
    for query_pop in super_populations.keys():
        print('query...', query_pop)
        # report scores by subpopulation
        query_full_dict = defaultdict(dict)
        # query_subpop_scores = {subpop : [] for subpop in SUBPOPULATIONS}
        for database_pop in super_populations.keys():
            print('...database...', database_pop)
            top_hits_dir = data_dir + query_pop + '_db/' + database_pop + '_top_hits/'
            # iterate through all top hits for each sample
            for query_ID in os.listdir(top_hits_dir):
                q_ID = query_ID.replace('.knn', '')
                # assert sample is in expected superpopulation
                assert_sample_pop(q_ID, query_pop, sample_subpopulations)
                query_file = top_hits_dir + query_ID
                query_full_dict = read_top_hits(query_file,
                                                    query_full_dict,
                                                    database_pop,
                                                    sample_subpopulations)

        write_query_pop_scores(query_pop,
                               query_full_dict)

if __name__ == '__main__':
    main()
