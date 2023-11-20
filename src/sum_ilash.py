import argparse
import glob

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, required=True)
    parser.add_argument('--out_file', type=str, required=True)
    return parser.parse_args()

def main():
    args = get_args()

    pairs = {}

    for file in glob.glob(args.data_dir + '*.match'):
        with open(file) as lines:
            for line in lines:
                A = line.rstrip().split('\t')
                a = A[1][:-2]
                b = A[3][:-2]
                l = float(A[9])

                pair = tuple(sorted([a,b]))

                if pair not in pairs:
                    pairs[pair] = 0
                pairs[pair] = pairs[pair] + l

    with open(args.out_file, 'w') as f:
        for pair in pairs:
            f.write('\t'.join([pair[0], pair[1] ,str(pairs[pair])]) + '\n')

if __name__ == '__main__':
    main()
