import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import random

from plotting import ancestry_helpers
from src import get_relations

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--times', type=str, help='times dir', required=True)
    parser.add_argument('--out', type=str, help='output directory', required=True)

    return parser.parse_args()

SYSTEM='eureka'
DATA='CCPM'
NUM_SAMPLES=73346
#146692

def get_segment_search_times(segment_file):
    read_embeddings_times = []
    loading_index_times = []
    searching_index_times = []
    scoring_times = []
    full_times = []
    sorting_time = 0

    f = open(segment_file, 'r')
    for line in f:
        if 'config' in line: # chrm8.segment117.config
            segment_idx = int(line.strip().split('segment')[1].split('.')[0])
            print(segment_idx)
            continue
        if 'reading' in line:
            time = float(line.strip().split(':')[1].split()[0])
            read_embeddings_times.append(time)
            continue
        if 'loading' in line:
            time = float(line.strip().split(':')[1].split()[0])
            loading_index_times.append(time)
            continue
        if 'searching' in line:
            time = float(line.strip().split(':')[1].split()[0])
            searching_index_times.append(time)
            continue
        if 'scoring' in line:
            time = float(line.strip().split(':')[1].split()[0])
            scoring_times.append(time)
            continue
        if 'full' in line:
            time = float(line.strip().split(':')[1].split()[0])
            full_times.append(time)
            continue
        if 'sorting' in line:
            time = float(line.strip().split(':')[1].split()[0])
            sorting_time = time
            continue

    f.close()
    return read_embeddings_times, loading_index_times, searching_index_times, scoring_times, full_times, sorting_time

def plot_segment_search_times(reading_embeddings_times,
                              loading_index_times,
                              searching_index_times,
                              scoring_times,
                              full_times,
                              sorting_time,
                              out,
                              chrm,
                              sample):

    # sns.set(style='whitegrid')
    # 5 subplots
    fig, axs = plt.subplots(2, 1, figsize=(15, 10), dpi=300)

    num_segments = len(reading_embeddings_times)
    num_queries = NUM_SAMPLES * 2
    # title = DATA + '-chrm ' + chrm + '\nTIMING BY SEGMENT\n' + SYSTEM
    # fig.suptitle(title, fontsize=25, fontweight='bold')
    out_file = out + 'chrm' + chrm + '_segment_timing.png'

    if sample:
        # divide all scores by number of queries
        reading_embeddings_times = [x/num_queries for x in reading_embeddings_times]
        loading_index_times = [x/num_queries for x in loading_index_times]
        searching_index_times = [x/num_queries for x in searching_index_times]
        scoring_times = [x/num_queries for x in scoring_times]
        full_times = [x/num_queries for x in full_times]
        num_queries = NUM_SAMPLES * 2
        # title = DATA + '-chrm ' + chrm + '\nTIMING BY QUERY\n' + SYSTEM
        # fig.suptitle(title, fontsize=25, fontweight='bold')
        out_file = out + 'chrm' + chrm + '_sample_timing.png'

    color = 'darkorange'
    # sns.histplot(reading_embeddings_times, kde=True, ax=axs[0])
    # sns.histplot(loading_index_times, kde=True, ax=axs[1])
    sns.histplot(searching_index_times, kde=True, ax=axs[0], color=color)
    # sns.histplot(scoring_times, kde=True, ax=axs[3])
    sns.histplot(full_times, kde=True, ax=axs[1], color=color)


    # axs[0].set_title('Reading Embeddings', fontsize=16, fontweight='bold', pad=10, loc='left')
    # axs[0].set_xlim(0, max(reading_embeddings_times) + 0.1)
    # axs[1].set_title('Loading Index', fontsize=16, fontweight='bold', loc='left')
    axs[0].set_title('SVS Search', fontsize=35, fontweight='bold', loc='center', pad=10)
    # axs[3].set_title('Scoring\n(adding pop count in dictionary)', fontsize=16, fontweight='bold', loc='left')
    axs[1].set_title('GenoSiS Search', fontsize=35, fontweight='bold', loc='center', pad=10)


    for i in range(2):
        axs[i].set_xlabel('Time (seconds)', fontsize=20)
        axs[i].set_ylabel('Frequency', fontsize=20)
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)
        # add text box of number of queries and number of segments
        textstr = 'Number of Queries: ' + str(num_queries) + '\nNumber of Segments: ' + str(num_segments)
        # props = dict(boxstyle='round', alpha=0.5, facecolor='white')
        axs[i].text(0.55, 0.55, textstr, transform=axs[i].transAxes, fontsize=20)

    # if not sample or sample:
    #     # subplots 1-4 should have same x-axis limits
    #     max_x = max(max(reading_embeddings_times), max(loading_index_times), max(searching_index_times), max(scoring_times))
    #     # get standard deviation of each
    #     reading_embeddings_std = round(stats.tstd(reading_embeddings_times), 2)
    #     loading_index_std = round(stats.tstd(loading_index_times), 2)
    #     searching_index_std = round(stats.tstd(searching_index_times), 2)
    #     scoring_std = round(stats.tstd(scoring_times), 2)
    #     avg_std = round((reading_embeddings_std + loading_index_std + searching_index_std + scoring_std) / 4, 2)
    #
    #     for i in range(2):
    #         axs[i].set_xlim(0, max_x + avg_std)

    # if sample:
    #     # subplots 1-4 should have same x-axis limits
    #     max_x = max(max(reading_embeddings_times), max(loading_index_times), max(searching_index_times),
    #                 max(scoring_times))
    #     for i in range(4):
    #         axs[i].set_xlim(0, max_x )


    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()

def plot_search_times(loading_index_times,
                      searching_index_times,
                      scoring_times,
                      out,
                      chrm,
                      sample):
    # combine searching and scoring times
    full_times = [x + y + z for x, y, z in zip(searching_index_times, scoring_times, loading_index_times)]
    # plot histogram of full times
    fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
    # title = DATA + '-chrm ' + chrm + '\nTIMING BY SEGMENT\n' + SYSTEM
    title = 'GenoSiS Search\n' + DATA + ', chrm ' + chrm
    # ax.set_title(title, fontsize=20, fontweight='bold')
    out_file = out + 'chrm' + chrm + '_segment_full_timing.png'

    num_segments = len(full_times)
    num_queries = NUM_SAMPLES * 2

    if sample:
        # full_times = [x / num_queries for x in full_times]
        full_times = [(x / num_queries) * num_segments for x in full_times]

        # title = DATA + '-chrm ' + chrm + '\nTIMING BY QUERY\n' + SYSTEM
        title = 'GenoSiS Search\n' + DATA + ', chrm ' + chrm
        # ax.set_title(title, fontsize=20, fontweight='bold')
        out_file = out + 'chrm' + chrm + '_sample_full_timing.png'

    sns.histplot(full_times, kde=True, ax=ax, color='darkorange')
    ax.set_xlabel('Time\n(seconds)')
    ax.set_ylabel('Frequency')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    textstr = 'Number of Queries: ' + str(num_queries) + '\nNumber of Segments: ' + str(num_segments)
    props = dict(boxstyle='round', alpha=0.5, facecolor='white')
    ax.text(0.55, 0.95, textstr, transform=ax.transAxes, fontsize=14,

            verticalalignment='top', bbox=props)

    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()


def plot_outer_loop_timing(full_times,
                            out,
                            chrm,
                            sample):
    # just plot outer loop (full) times
    fig, ax = plt.subplots(figsize=(8, 5), dpi=200)
    title = DATA + '-chrm ' + chrm + '\nTIMING BY SEGMENT\n' + SYSTEM
    ax.set_title(title, fontsize=25, fontweight='bold')
    out_file = out + 'chrm' + chrm + '_segment_full_timing.png'

    num_segments = len(full_times)
    num_queries = NUM_SAMPLES * 2


    if sample:
        full_times = [x / num_queries for x in full_times]
        # title = DATA + '-chrm ' + chrm + '\nTIMING BY QUERY\n' + SYSTEM
        # ax.set_title(title, fontsize=25, fontweight='bold')
        out_file = out + 'chrm' + chrm + '_sample_full_timing.png'

    sns.histplot(full_times, kde=True, ax=ax)
    ax.set_xlabel('Time\n(seconds)')
    ax.set_ylabel('Frequency')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    textstr = 'Number of Queries: ' + str(num_queries) + '\nNumber of Segments: ' + str(num_segments)
    props = dict(boxstyle='round', alpha=0.5, facecolor='white')
    ax.text(0.55, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)

    plt.tight_layout()
    plt.savefig(out_file)


def main():
    args = get_args()
    times_dir = args.times
    out = args.out

    # all_times = get_segment_search_times_OLD(segment_file)
    # plot_segment_search_times_OLD(all_times, out, chrm)

    read_embeddings_times_list = []
    loading_index_times_list = []
    searching_index_times_list = []
    scoring_times_list = []
    full_times_list = []
    sorting_times_list = []

    included_chroms = []

    # for all chromosomes, make long list of times
    for chrom in range(1, 23):
        try:
            segment_file = times_dir + 'chrm' + str(chrom) + '_segments.log'
            (reading_embeddings_times,
             loading_index_times,
             searching_index_times,
             scoring_times,
             full_times,
             sorting_time) = get_segment_search_times(segment_file)

            # append to lists
            read_embeddings_times_list.extend(reading_embeddings_times)
            loading_index_times_list.extend(loading_index_times)
            searching_index_times_list.extend(searching_index_times)
            scoring_times_list.extend(scoring_times)
            full_times_list.extend(full_times)
            sorting_times_list.append(sorting_time)

            included_chroms.append(chrom)

        except FileNotFoundError:
            continue


    chrm = str(min(included_chroms)) + '-' + str(max(included_chroms))


    # plot_segment_search_times(read_embeddings_times_list,
    #                           loading_index_times_list,
    #                           searching_index_times_list,
    #                           scoring_times_list,
    #                           full_times_list,
    #                           sorting_times_list,
    #                           out,
    #                           chrm, False)
    # plot_segment_search_times(read_embeddings_times_list,
    #                           loading_index_times_list,
    #                           searching_index_times_list,
    #                           scoring_times_list,
    #                           full_times_list,
    #                           sorting_times_list,
    #                           out,
    #                           chrm, True)

    # plot_outer_loop_timing(full_times_list, out, chrm, False)
    # plot_outer_loop_timing(full_times_list, out, chrm, True)

    plot_search_times(loading_index_times_list, searching_index_times_list, scoring_times_list, out, chrm, False)
    plot_search_times(loading_index_times_list, searching_index_times_list, scoring_times_list, out, chrm, True)



if __name__ == '__main__':
    main()