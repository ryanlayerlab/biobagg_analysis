import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', type=str, required=True)
    parser.add_argument('--num_samples', type=int, required=True)
    parser.add_argument('--out_file', type=str, required=True)
    parser.add_argument('--width', type=float, default=10.0)
    parser.add_argument('--height', type=float, default=5.0)
    return parser.parse_args()

def read_time_file(file):
    times = []
    with open(file) as lines:
        for line in lines:
            if line[:4] == 'real':
                time_str = line.rstrip().split()[1]
                parts = time_str.replace('s', '').split('m')
                minutes = float(parts[0])
                seconds = float(parts[1])
                total_seconds = minutes * 60 + seconds
                times.append(total_seconds)

    return times


def main():
    args = get_args()
    ms = [ time/args.num_samples * 1000 \
            for time in read_time_file(args.in_file)]

    fig, ax = plt.subplots(figsize=(args.width, args.height))

    ax.hist(ms, bins=50)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylabel('Freq.')
    ax.set_xlabel('Time (ms)')

    print('per sample per segment median (ms)', np.median(ms))
    print('per sample per segment mean (ms)', np.mean(ms))
    print('per sample per segment stdev (ms)', np.std(ms))
    # len(ms) is number of segments, so to get a file search time for running
    # all segments for a sample we multiply then conver to seconds
    print('per sample mean run time (s)', (np.mean(ms) * len(ms))/1000 )

    plt.tight_layout()
    plt.savefig(args.out_file, dpi=300)

if __name__ == '__main__':
    main()
