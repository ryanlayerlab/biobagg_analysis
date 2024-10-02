import sys
import numpy as np
import matplotlib.pyplot as plt

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
            'European': 'TGP+HGP EUR-like',
            'African': 'TGP+HGP AFR-like',
            'Admixed American': 'TGP+HGP AMR-like',
            'South-Central Asian': 'TGP+HGP SAS-like',
            'East Asian': 'TGP+HGP EAS-like',
            'Middle Eastern': 'TGP+HGP MLE-like',
            'Oceania': 'OCE-like'}


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
fig, axs = plt.subplots(rows, cols, figsize=(18, 15), sharex=True, sharey=False)

for i, (long_name, values) in enumerate(R.items()):
    ax = axs[i // cols, i % cols]
    ax.hist(values, bins=20, color='gray')
    ax.set_title(pretty_name_map[long_name], fontsize=20, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12)

    # spine removal
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# add x-axis label and tick labels for only the bottom row
for ax in axs.flat:
    if ax in axs[rows - 1]:
        ax.set_xlabel('Ancestry inference probability', fontsize=12)


plt.tight_layout()
plt.savefig('pub_figures/ccpm_prob.png')