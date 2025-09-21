import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

ccpm_probabilities = sys.argv[1]

name_map = {'EUR': 'European',
            'AFR': 'African',
            'AMR': 'Admixed American',
            'SAS': 'South-Central Asian',
            'EAS': 'East Asian',
            'MLE': 'Middle Eastern',
            'OCE': 'Oceania',
            'European': 'EUR',
            'African': 'AFR',
            'Admixed American': 'AMR',
            'South-Central Asian': 'SAS',
            'East Asian': 'EAS',
            'Middle Eastern': 'MLE',
            'Oceania': 'OCE'}

pretty_name_map = {
            'European': 'TGP+HGP\nEUR-like',
            'African': 'TGP+HGP\nAFR-like',
            'Admixed American': 'TGP+HGP\nAMR-like',
            'South-Central Asian': 'TGP+HGP\nSAS-like',
            'East Asian': 'TGP+HGP\nEAS-like',
            'Middle Eastern': 'TGP+HGP\nMLE-like',
            'Oceania': 'OCE-like'}

ordered_ancestry_labels = {1:'TGP+HGP-\nAFR-like',
                               2:'TGP+HGP-\nAMR-like',
                               3:'TGP+HGP-\nEAS-like',
                               4:'TGP+HGP-\nEUR-like',
                               5:'TGP+HGP-\nMLE-like',
                               6:'TGP+HGP-\nSAS-like'}


header = None
D = []
f = open(ccpm_probabilities, 'r')
for l in f:
    A = l.strip().split(',')
    if header is None:
        header = A
        continue
    d = dict(zip(header, A))
    D.append(d)

R = {}
for d in D:
    long_name = d['Inferred_ancestry']
    short_name = name_map[d['Inferred_ancestry']]
    if long_name not in R:
        R[long_name] = []
    R[long_name].append(float(d[short_name]))

rows = int(len(R) / 2)
cols = 2
fig, axs = plt.subplots(1, 6, figsize=(15, 3), dpi=300, sharex=True, sharey=False)

color_CCPM = {'African': 'deepskyblue',
                  'Admixed American': 'goldenrod',
                  'East Asian': 'crimson',
                  'European': 'yellowgreen',
                  'Middle Eastern': 'darkorange',
                  'South-Central Asian': 'mediumpurple'}


for i, (long_name, values) in enumerate(sorted(R.items())):
    ax = axs[i]
    # ax.hist(values, bins=20, color=color_CCPM[long_name], alpha=0.7, edgecolor='black', linewidth=1.2)
    sns.histplot(values, bins=20, color=color_CCPM[long_name], ax=ax)
    ax.set_xlabel('Ancestry inference probability', fontsize=10)
    # remove y-axis labels
    ax.set_title(pretty_name_map[long_name], fontsize=10, fontweight='bold', color=color_CCPM[long_name])

    # spine removal
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # only label y-axis for the left column
    if i == 0:
        ax.set_ylabel('Count', fontsize=10)
    else:
        ax.set_ylabel('', fontsize=10)
        # ax.set_yticklabels([])


# for i, (long_name, values) in enumerate(sorted(R.items())):
#     ax = axs[i // cols, i % cols]
#     # ax.hist(values, bins=20, color='slategrey', alpha=0.7, edgecolor='black', linewidth=1.2)
#     ax.hist(values, bins=20, color=color_CCPM[long_name], alpha=0.7, edgecolor='black', linewidth=1.2)
#     ax.set_title(pretty_name_map[long_name], fontsize=20, fontweight='bold')
#     ax.set_ylabel('Frequency', fontsize=18)
#
#     # spine removal
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)

# add x-axis label and tick labels for only the bottom row
# for ax in axs.flat:
#     if ax in axs[rows - 1]:
#         ax.set_xlabel('Ancestry inference probability', fontsize=18)


plt.tight_layout()
plt.savefig('pub_figures/ccpm_prob.png')