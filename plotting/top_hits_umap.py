import argparse
import numpy as np
import umap
import matplotlib.pyplot as plt
from matplotlib import cm

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', type=str, required=True)
    parser.add_argument('--label_file', type=str, required=True)
    parser.add_argument('--out_file', type=str, required=True)
    return parser.parse_args()

def get_top_hits(file):
    str_hits = {}
    ids = {}

    with open(file) as lines:
        for line in lines:
            A = line.rstrip().split()
            q = A[0]
            ids[q] = len(ids)
            str_hits[q] = [a.split(',')[0] for a in A[1:]]

    hits = {}

    for str_h in str_hits:
        hits[ str_h ]  = [ids[i] for i in str_hits[str_h]]

    return hits

def get_label_map(file_name):
    labels = {}
    header = None
    with open(file_name) as lines:
        for line in lines:
            A = line.rstrip().split('\t')

            if header is None:
                header = A
                continue

            sample = A[0]
            label = A[5]
            labels[sample] = label
    return labels

def main():
    args = get_args()

    top_hits = get_top_hits(args.in_file)
    ids = sorted(top_hits.keys())

    D = []
    for i in ids:
        D.append(top_hits[i])

    label_map = get_label_map(args.label_file)
   
    labels = [label_map[i] for i in ids]

    unique_labels = list(set(labels))

    label_id_map = {value: index for index, value in enumerate(unique_labels)}

    label_idxs = [label_id_map[label] for label in labels]

    reducer = umap.UMAP()
    embedding = reducer.fit_transform(D)

    colormap = cm.get_cmap('tab10', len(unique_labels))

    scatter = plt.scatter(embedding[:, 0],
                          embedding[:, 1],
                          c=label_idxs,
                          cmap=colormap,
                          alpha=0.5)

    plt.title('UMAP projection of the dataset')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')

    cbar = plt.colorbar(scatter, ticks=range(len(unique_labels)))
    cbar.set_ticklabels(unique_labels)

    plt.savefig(args.out_file, dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()
