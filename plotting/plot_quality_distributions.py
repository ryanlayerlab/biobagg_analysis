import argparse
from collections import defaultdict
import os
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import random

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query', type=str, help='top hits file', required=True)
    parser.add_argument('--out', type=str, help='output file', required=True)

    return parser.parse_args()

pop_color_dict = {'AFR': 'darkorange',
                  'AMR': 'mediumpurple',
                  'EAS': 'hotpink',
                  'EUR': 'steelblue',
                  'SAS': 'goldenrod'}

def read_summary_data(query_pop_summary):
    query_pop_results = {'AFR':[], 'AMR':[], 'EAS':[], 'EUR':[], 'SAS':[]}
    with open(query_pop_summary, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip().split(',')
            database_pop = line[0]
            scores = [float(x) for x in line[1:]]
            query_pop_results[database_pop].extend(scores)
    return query_pop_results


def read_query_results(db, query_pop_results, database_pop):
    for query in os.listdir(db):
        with open(db + query, 'r') as f:
            query = f.readline().strip().split(': ')[1]
            # only read top 20 hits
            match_count = 0
            for line in f:
                if match_count == 20:
                    break
                line = line.strip().split(',')
                match = line[0]
                # ignore self
                if query == match:
                    continue
                score = float(line[1])
                if score > 100:
                    print(score, query, match)
                query_pop_results[database_pop].append(score)
                match_count += 1
    return query_pop_results

def plot_data(query_pop_results, query_pop, out):
    # density distribution and histogram for each population
    ax = plt.figure(figsize=(8, 5), dpi=200)
    for database_pop, scores in query_pop_results.items():
        # sns.kdeplot(scores, color=pop_color_dict[database_pop], label=database_pop)
        sns.histplot(scores, color=pop_color_dict[database_pop], label=database_pop, bins=100, kde=True, alpha=0.5)

        # remove bars

    plt.legend()
    plt.xlabel('GenoSiS Score')
    plt.ylabel('Frequency')
    sns.despine()
    plt.legend(title='Database Population', labels=['AFR', 'AMR', 'EAS', 'EUR', 'SAS'], frameon=False)

    # title
    plt.title(query_pop + ' Query Population', fontsize=20, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out + query_pop + '_query.png')
    # plot the mean of each population

# def write_summary_data(query_pop_results, query_pop, out):
#     with open(out + query_pop + '_query_summary.txt', 'w') as f:
#         f.write('Database_Population,scores\n')
#         for database_pop, scores in query_pop_results.items():
#             f.write(database_pop + ',' + ','.join([str(x) for x in scores]) + '\n')

def main():
    args = get_args()


    # USED FOR WRITING SUMMARY DATA ONLY
    # query_dir = 'data/' + args.query + '_query/'
    # query_pop_results = {'AFR':[], 'AMR':[], 'EAS':[], 'EUR':[], 'SAS':[]}
    # for database_pop in query_pop_results.keys():
    #     print(database_pop)
    #     query_file = query_dir + database_pop + '_db/'
    #     query_pop_results = read_query_results(query_file, query_pop_results, database_pop)
    # write_summary_data(query_pop_results, args.query, args.out)

    # read in summary data
    query_pop_summary = 'data/' + args.query + '_query_summary.txt'
    query_pop_results = read_summary_data(query_pop_summary)
    # plot summary data
    plot_data(query_pop_results, args.query, args.out)



if __name__ == '__main__':
    main()