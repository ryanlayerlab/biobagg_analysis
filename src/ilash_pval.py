import pysam
import glob
import re
import argparse

from scipy.stats import fisher_exact

def perform_fishers_exact_test(ibd_segments, in_segments, all_segments):
    in_and_ibd = in_segments.intersection(ibd_segments)
    in_and_not_ibd = in_segments.difference(ibd_segments)
    not_in_and_ibd = ibd_segments.difference(in_segments)
    not_in_and_not_ibd = all_segments.difference(in_segments.union(ibd_segments))
    table = [[len(in_and_ibd), len(in_and_not_ibd)], [len(not_in_and_ibd), len(not_in_and_not_ibd)]]
    odds_ratio, p_value = fisher_exact(table)
    return odds_ratio, p_value

def get_args():
    #parser = argparse.ArgumentParser()
    #parser.add_argument('-S', nargs=2, help='The two individuals')
    parser = argparse.ArgumentParser()
    parser.add_argument('--ilash', type=str, help='The iLASH file', required=True)
    parser.add_argument('--hapsis', type=str, help='The hapsis path pattern', required=True)
    parser.add_argument('--segment_bed', type=str, help='The segment bed file', required=True)
    return parser.parse_args()

def get_hapsis_segments(hapsis_path_pattern):
    file_list = glob.glob(hapsis_path_pattern)

    all_segments = []

    hapsis_segments = {}

    pattern = re.compile(r'.*segment(\d+)\.knn$')
    for file_name in file_list:
        match = pattern.search(file_name)
        segment = int(match.group(1))
        all_segments.append(segment)

        src = None
        dst = None
        with open(file_name) as f:
            for line in f:
                if len(line) <= 1:
                    continue
                elif line.startswith('Query:'):
                    src = line.rstrip().split()[1]
                else:
                    dst = line.rstrip().split()[0]
                    if (src, dst) not in hapsis_segments:
                        hapsis_segments[(src, dst)] = []
                    hapsis_segments[(src, dst)].append(int(segment))
    return hapsis_segments, all_segments

def get_segments(A, B, segment_path_pattern):
    #segment_path_pattern = 'svs_results_chrm15-20/chrm15.segment*'
    file_list = glob.glob(segment_path_pattern)

    in_segments = []
    all_segments = []

    pattern = re.compile(r'.*segment(\d+)\.knn$')
    for file_name in file_list:
        match = pattern.search(file_name)
        segment = match.group(1)
        all_segments.append(int(segment))
        with open(file_name) as f:
            in_A = False
            for line in f:
                if line.startswith('Query:'):
                    if A in line:
                        in_A = True
                    elif in_A:
                        break
                elif in_A:
                    if B in line:
                        in_segments.append(int(segment))
                        break

    return in_segments, all_segments

def get_ibd_intervals(ilash_file_name):
    ibd_segments = {}

    with open(ilash_file_name) as f:
        for line in f:
            line = line.strip().split()
            A = line[1]
            B = line[3]
            if (A, B) not in ibd_segments:
                ibd_segments[(A, B)] = []

            ibd_segments[(A,B)].append((int(line[4]), int(line[5]), int(line[6])))
    return ibd_segments

def get_ibd_segments(ibd_intervals, segment_bed_file):
    tabix_file = pysam.TabixFile(segment_bed_file)

    ibd_segments = {}

    cache = {}

    for pair in ibd_intervals:
        ibd_segments[pair] = []
        for ibd_interval in ibd_intervals[pair]:
            if ibd_interval not in cache:
                cache[ibd_interval] = []
                for row in tabix_file.fetch(ibd_interval[0], ibd_interval[1], ibd_interval[2]):
                    cache[ibd_interval].append(int(row.rstrip().split()[3]))
            ibd_segments[pair].extend(cache[ibd_interval])
    return ibd_segments

def main():
    args = get_args()
    ibd_intervals = get_ibd_intervals(args.ilash)

    ibd_segments = get_ibd_segments(ibd_intervals, args.segment_bed)

    hapsis_segments, all_segments  = get_hapsis_segments(args.hapsis)

    all_segments = set(all_segments)

    for pair in ibd_segments:
        if pair not in hapsis_segments:
            continue
        hapsis = set(hapsis_segments[pair])
        ilash = set(ibd_segments[pair])
        odds_ratio, p_value = perform_fishers_exact_test(ilash, hapsis, all_segments)
        print(f'{pair}\t{len(ilash)}\t{len(hapsis)}\t{len(ilash.intersection(hapsis))}\t{odds_ratio}\t{p_value}')

if __name__ == '__main__':
    main()

