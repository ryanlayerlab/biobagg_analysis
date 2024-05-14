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
    parser.add_argument('--segments', type=str, help='top hits file', required=True)
    parser.add_argument('--times', type=str, help='trios scores', required=True)
    parser.add_argument('--chrm', type=str, help='output directory', required=True)
    parser.add_argument('--out', type=str, help='output directory', required=True)

    return parser.parse_args()

SYSTEM='eureka'
DATA='CCPM'
NUM_SAMPLES=100

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


def get_segment_search_times_OLD(segment_file):
    all_times = []
    f = open(segment_file, 'r')
    header = f.readline()
    for line in f:
        if 'config' in line:
            continue
        line = line.strip().split(':')
        time, unit = line[1].split()
        time = float(time)
        all_times.append(time)

    f.close()
    return all_times

def plot_segment_search_times_OLD(all_times, out, chrm):
    # sns.set(style='whitegrid')
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(all_times, kde=True, ax=ax)
    num_segments = len(all_times)
    num_queries = NUM_SAMPLES * 2
    title = DATA + '-chrm ' + chrm + '\nTIMING BY SEGMENT\n' + SYSTEM
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlabel('Time\n(seconds)')
    ax.set_ylabel('Frequency')

    # add text box with data information
    textstr = '\n'.join((
        r'$\mathrm{num\ segments}=%.0f$' % (num_segments, ),
        r'$\mathrm{num\ samples}=%.0f$' % (NUM_SAMPLES, ),
        r'$\mathrm{num\ queries}=%.0f$' % (num_queries, ),
        r'$\mathrm{mean}=%.2f$' % (sum(all_times)/len(all_times), ),
        r'$\mathrm{median}=%.2f$' % (sorted(all_times)[len(all_times)//2], )))

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(out + 'segment_search_times.png')
    plt.close()

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
    fig, axs = plt.subplots(5, 1, figsize=(10, 15), dpi=200)

    num_segments = len(reading_embeddings_times)
    num_queries = NUM_SAMPLES * 2
    title = DATA + '-chrm ' + chrm + '\nTIMING BY SEGMENT\n' + SYSTEM
    fig.suptitle(title, fontsize=25, fontweight='bold')
    out_file = out + 'chrm' + chrm + '_segment_timing.png'

    if sample:
        # divide all scores by number of queries
        reading_embeddings_times = [x/num_queries for x in reading_embeddings_times]
        loading_index_times = [x/num_queries for x in loading_index_times]
        searching_index_times = [x/num_queries for x in searching_index_times]
        scoring_times = [x/num_queries for x in scoring_times]
        full_times = [x/num_queries for x in full_times]
        num_queries = NUM_SAMPLES * 2
        title = DATA + '-chrm ' + chrm + '\nTIMING BY SAMPLE\n' + SYSTEM
        fig.suptitle(title, fontsize=25, fontweight='bold')
        out_file = out + 'chrm' + chrm + '_sample_timing.png'

    color = 'darkblue'
    sns.histplot(reading_embeddings_times, kde=True, ax=axs[0])
    sns.histplot(loading_index_times, kde=True, ax=axs[1])
    sns.histplot(searching_index_times, kde=True, ax=axs[2])
    sns.histplot(scoring_times, kde=True, ax=axs[3])
    sns.histplot(full_times, kde=True, ax=axs[4])


    axs[0].set_title('Reading Embeddings', fontsize=16, fontweight='bold', pad=10, loc='left')
    # axs[0].set_xlim(0, max(reading_embeddings_times) + 0.1)
    axs[1].set_title('Loading Index', fontsize=16, fontweight='bold', loc='left')
    axs[2].set_title('Searching Index', fontsize=16, fontweight='bold', loc='left')
    axs[3].set_title('Scoring\n(adding pop count in dictionary)', fontsize=16, fontweight='bold', loc='left')
    axs[4].set_title('Full\n(outer loop)', fontsize=16, fontweight='bold', loc='left')


    for i in range(5):
        axs[i].set_xlabel('Time\n(seconds)')
        axs[i].set_ylabel('Frequency')
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)

    if not sample or sample:
        # subplots 1-4 should have same x-axis limits
        max_x = max(max(reading_embeddings_times), max(loading_index_times), max(searching_index_times), max(scoring_times))
        # get standard deviation of each
        reading_embeddings_std = round(stats.tstd(reading_embeddings_times), 2)
        loading_index_std = round(stats.tstd(loading_index_times), 2)
        searching_index_std = round(stats.tstd(searching_index_times), 2)
        scoring_std = round(stats.tstd(scoring_times), 2)
        avg_std = round((reading_embeddings_std + loading_index_std + searching_index_std + scoring_std) / 4, 2)

        for i in range(4):
            axs[i].set_xlim(0, max_x + avg_std)

    # if sample:
    #     # subplots 1-4 should have same x-axis limits
    #     max_x = max(max(reading_embeddings_times), max(loading_index_times), max(searching_index_times),
    #                 max(scoring_times))
    #     for i in range(4):
    #         axs[i].set_xlim(0, max_x )


    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()


def main():
    args = get_args()
    segment_file = args.segments
    times_file = args.times
    chrm = args.chrm
    out = args.out

    # all_times = get_segment_search_times_OLD(segment_file)
    # plot_segment_search_times_OLD(all_times, out, chrm)

    (reading_embeddings_times,
     loading_index_times,
     searching_index_times,
     scoring_times,
     full_times,
     sorting_time) = get_segment_search_times(segment_file)

    sample = False
    plot_segment_search_times(reading_embeddings_times,
                              loading_index_times,
                              searching_index_times,
                              scoring_times,
                              full_times,
                              sorting_time,
                              out,
                              chrm, sample)


if __name__ == '__main__':
    main()