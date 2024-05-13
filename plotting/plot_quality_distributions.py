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
    plt.figure(figsize=(10, 6), dpi=200)
    for pop, scores in query_pop_results.items():
        sns.distplot(scores, hist=False, kde=True, kde_kws={'shade': True, 'linewidth': 3},
                     label=pop, color=pop_color_dict[pop])
        # plt.hist(scores, bins=20, alpha=0.6, color=pop_color_dict[pop])
        # log scale
        # plt.yscale('log')

    # make them all share the same x-axis

    sns.despine()
    plt.xlabel('GenoSiS score')
    plt.ylabel('Density')
    plt.title('GenoSiS scores for ' + query_pop)
    plt.legend(loc='upper right', frameon=False, title='Database Population')
    plt.tight_layout()
    plt.savefig(out + query_pop + '_query.png')
    # plot the mean of each population

def write_summary_data(query_pop_results, query_pop, out):
    with open(out + query_pop + '_query_summary.txt', 'w') as f:
        f.write('Database_Population,scores\n')
        for database_pop, scores in query_pop_results.items():
            f.write(database_pop + ',' + ','.join([str(x) for x in scores]) + '\n')



def main():
    args = get_args()
    query_dir = 'data/' + args.query + '_query/'

    query_pop_results = {'AFR':[], 'AMR':[], 'EAS':[], 'EUR':[], 'SAS':[]}

    for database_pop in query_pop_results.keys():
        print(database_pop)
        query_file = query_dir + database_pop + '_db/'
        query_pop_results = read_query_results(query_file, query_pop_results, database_pop)

    plot_data(query_pop_results, args.query, args.out)
    write_summary_data(query_pop_results, args.query, args.out)


if __name__ == '__main__':
    main()