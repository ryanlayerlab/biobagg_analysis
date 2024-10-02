import argparse
import sys
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colormaps
from matplotlib.ticker import FormatStrFormatter
from colorspacious import cspace_converter

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath('plotting/'))
import ancestry_helpers


def get_args():
    parser = argparse.ArgumentParser(description="FST for subpopulation queires")
    parser.add_argument("--inputs", type=str, nargs="+", help="Input files")
    parser.add_argument("--labels", type=str, nargs="+", help="Labels")
    parser.add_argument(
        "--output", type=str, default="histogram.png", help="Output file name"
    )
    parser.add_argument("--height", type=float, default=5, help="Height of the figure")
    parser.add_argument("--width", type=float, default=5, help="Width of the figure")
    parser.add_argument('--colors', help='file with color codes', required=True)
    return parser.parse_args()


def get_file_data(file_name):
    D = []
    with open(file_name, "r") as f:
        for line in f:
            src, dst, freq, val = line.rstrip().split()
            # if src == dst: continue
            if int(freq) == 0:
                continue
            for i in range(int(freq)):
                D.append(float(val))
    return D


def get_cdf(D, start, stop, step):
    D.sort()
    cdf = []
    for i in np.arange(start, stop, step):
        cdf.append(len([d for d in D if d < i]) / len(D))
    return cdf


def plot_cdf(inputs: list[str], labels: list[str], ax: plt.Axes):

    colors = ['black',
                'darkorange',
                'limegreen',
                'steelblue']

    line_styles = ['dashed', 'dashed', 'dashed', 'dashed']
    line_styles = ['solid', 'solid', 'solid', 'solid']


    for i, input_file in enumerate(inputs):
        D = get_file_data(input_file)
        CDF = get_cdf(D, 0.0, 0.15, 0.001)
        ax.plot(np.arange(0.0, 0.15, 0.001), CDF, label=labels[i],
                color=colors[i], lw=2, alpha=0.9, linestyle=line_styles[i])

    ax.set_xlabel("FST")
    ax.set_ylabel("Perecnt of cohort with at most FST", )
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))

def plot_fst_heatmap(input: str, ax: plt.Axes, colors):
    df = pd.read_csv(input, sep=" ", names=["A", "B", "hits", "fst"])

    # hard coded subpopulations to group subpops by their superpop
    subpopulations = {
        "ASW": "AFR",
        "LWK": "AFR",
        "GWD": "AFR",
        "MSL": "AFR",
        "ESN": "AFR",
        "YRI": "AFR",
        "ACB": "AFR",
        "CLM": "AMR",
        "PEL": "AMR",
        "MXL": "AMR",
        "PUR": "AMR",
        "CDX": "EAS",
        "CHB": "EAS",
        "JPT": "EAS",
        "KHV": "EAS",
        "CHS": "EAS",
        "CEU": "EUR",
        "TSI": "EUR",
        "FIN": "EUR",
        "GBR": "EUR",
        "IBS": "EUR",
        "BEB": "SAS",
        "GIH": "SAS",
        "ITU": "SAS",
        "PJL": "SAS",
        "STU": "SAS",
    }

    N = len(subpopulations)
    FST = np.zeros((N, N))

    for i, A in enumerate(subpopulations):
        for j, B in enumerate(subpopulations):
            FST[i, j] = df[(df["A"] == A) & (df["B"] == B)]["fst"].values[0]

    cmap_name = "Blues_r"
    im = ax.imshow(FST, cmap=cmap_name, norm=mpl.colors.Normalize(vmin=0, vmax=0.15))

    # Create colorbar only show from 0 to 0.15
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.05, pad=0.04, )
    # plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04) # weird but its perfect

    ax.set(
        xticks=np.arange(N),
        yticks=np.arange(N),
        xticklabels=subpopulations,
        yticklabels=subpopulations,
    )
    # color the x and y tick labels by superpop, and make bold
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_color(colors[subpopulations[label.get_text()]])
        label.set_fontweight("bold")
    for i, label in enumerate(ax.get_yticklabels()):
        label.set_color(colors[subpopulations[label.get_text()]])
        label.set_fontweight("bold")

    # add lines to separate superpopulations
    superpop_indices = [0, 7, 11, 16, 21, 26]
    superpop_labels = ["AFR", "AMR", "EAS", "EUR", "SAS"]
    for i in superpop_indices:
        ax.axhline(i - 0.5, color="black", lw=1)
        ax.axvline(i - 0.5, color="black", lw=1)

    # # add text labels along bottom below heatmap
    # for i, label in enumerate(superpop_labels):
    #     ax.text(
    #         -4,
    #         superpop_indices[i],
    #         label,
    #         color="black",
    #         fontsize=6,
    #         fontweight="bold",
    #         rotation=90
    #     )

    # ax.set_xticks(np.arange(N + 1) - 0.5, minor=True)
    # ax.set_yticks(np.arange(N + 1) - 0.5, minor=True)

    ax.tick_params(axis="both", labelsize=6)
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right", rotation_mode="anchor")


def main():
    args = get_args()

    colors = ancestry_helpers.get_colors(args.colors)

    cols = 2
    rows = 1
    fig, ax = plt.subplots(
        rows,
        cols,
        figsize=(args.width, args.height),
        # gridspec_kw={"width_ratios": [1, 1], "height_ratios": [1]},
    )

    plot_fst_heatmap(args.inputs[0], ax[0], colors)
    plot_cdf(args.inputs, args.labels, ax[1])
    ax[0].set_title(r"$\mathbf{A}$", loc="left")
    ax[1].set_title(r"$\mathbf{B}$", loc="left")
    ax[1].set_aspect(0.15, adjustable="box")

    plt.tight_layout()
    # plt.show()
    plt.savefig(args.output, dpi=600)


if __name__ == "__main__":
    main()
