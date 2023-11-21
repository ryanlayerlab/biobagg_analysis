import argparse
import sys
sys.path.insert(0, 'plotting')
import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs_file', type=str, required=True)
    parser.add_argument('--out_file', type=str, required=True)
    parser.add_argument('--N', type=int, required=True)
    return parser.parse_args()

def main():
    args = get_args()

    pairs = utils.get_pair_map(args.pairs_file)

    with open(args.out_file, 'w') as file:
        for i in pairs:
            hits = []
            for j in pairs[i]:
                hits.append((j,pairs[i][j]))
            hits = sorted(hits, key = lambda x: x[1], reverse=True)

            str_hits = [x[0]+','+str(x[1]) for x in hits[:args.N]]

            file.write('\t'.join([i] + str_hits) + '\n')

if __name__ == '__main__':
    main()
