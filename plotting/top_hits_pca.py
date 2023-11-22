import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sklearn.decomposition import PCA
import utils

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', type=str, required=True)
    parser.add_argument('--label_file', type=str, required=True)
    parser.add_argument('--out_file', type=str, required=True)
    return parser.parse_args()

def main():
    args = get_args()

    top_hits = utils.get_top_hits(args.in_file, integerize=True)
    ids = sorted(top_hits.keys())

    D = []
    for i in ids:
        D.append(top_hits[i])

    label_map = utils.get_label_map(args.label_file, 'Superpopulation code')
   
    labels = [label_map[i] for i in ids]

    unique_labels = list(set(labels))

    label_id_map = {value: index for index, value in enumerate(unique_labels)}

    label_idxs = [label_id_map[label] for label in labels]

    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(D)

    colormap = cm.get_cmap('tab10', len(unique_labels))

    scatter = plt.scatter(principal_components[:, 0],
                          principal_components[:, 1],
                          c=label_idxs,
                          cmap=colormap)

    plt.title('PCA projection of the dataset')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')

    cbar = plt.colorbar(scatter, ticks=range(len(unique_labels)))
    cbar.set_ticklabels(unique_labels)

    plt.savefig(args.out_file, dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()
