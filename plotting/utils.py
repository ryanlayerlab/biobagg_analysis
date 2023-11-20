from scipy.stats import gaussian_kde
import numpy as np

def get_top_hits(file, integerize=False, get_scores=False):
    str_hits = {}
    ids = {}

    with open(file) as lines:
        for line in lines:
            A = line.rstrip().split()
            q = A[0]
            ids[q] = len(ids)
            if get_scores:
                str_hits[q] = \
                    [(a.split(',')[0], float(a.split(',')[1])) for a in A[1:]]
            else:
                str_hits[q] = [a.split(',')[0] for a in A[1:]]

    if integerize:
        hits = {}

        for str_h in str_hits:
            hits[ str_h ]  = [ids[i] for i in str_hits[str_h]]

        return hits
    else:
        return str_hits

def get_label_map(file_name, col_name):
    labels = {}
    header = None
    col_idx = None
    with open(file_name) as lines:
        for line in lines:
            A = line.rstrip().split('\t')

            if header is None:
                header = A
                col_idx = header.index(col_name)
                continue

            sample = A[0]
            label = A[col_idx]
            labels[sample] = label
    return labels

def get_pair_map(file):
    pairs = {}

    with open(file) as lines:
        for line in lines:
            A = line.rstrip().split()
            a = A[0]
            b = A[1]
            dist = float(A[2])

            if a not in pairs:
                pairs[a] = {}
            if b not in pairs:
                pairs[b] = {}

            pairs[a][b] = dist
            pairs[b][a] = dist

    return pairs

def draw_smooth_histo(data, ax, color, label, lw):
    kde = gaussian_kde(data)
    kde.set_bandwidth(bw_method=kde.factor / 3.)
    x_range = np.linspace(min(data), max(data), 500)
    kde_values = kde(x_range)
    ax.plot(x_range, kde_values, color=color, lw=lw, alpha=0.75, label=label)

