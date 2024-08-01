import argparse
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description="FST for subpopulation queires")
    parser.add_argument("--inputs", type=str, nargs="+", help="Input files")
    parser.add_argument("--labels", type=str, nargs="+", help="Labels")
    parser.add_argument(
        "--output", type=str, default="histogram.png", help="Output file name"
    )
    parser.add_argument("--height", type=float, default=5, help="Height of the figure")
    parser.add_argument("--width", type=float, default=5, help="Width of the figure")
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
    for i, input_file in enumerate(inputs):
        D = get_file_data(input_file)
        CDF = get_cdf(D, 0.0, 0.15, 0.001)
        ax.plot(np.arange(0.0, 0.15, 0.001), CDF, label=labels[i])

    ax.set_xlabel("FST")
    ax.set_ylabel("Perecnt of cohort with at most FST")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))


def plot_fst_heatmap(input: str, ax: plt.Axes):
    df = pd.read_csv(input, sep=" ", names=["A", "B", "hits", "fst"])

    # hard coded subpopulations to group subpops by their superpop
    subpopulations = [
        "ASW",
        "LWK",
        "GWD",
        "MSL",
        "ESN",
        "YRI",
        "ACB",
        "CLM",
        "PEL",
        "MXL",
        "PUR",
        "CDX",
        "CHB",
        "JPT",
        "KHV",
        "CHS",
        "CEU",
        "TSI",
        "FIN",
        "GBR",
        "IBS",
        "BEB",
        "GIH",
        "ITU",
        "PJL",
        "STU",
    ]
    N = len(subpopulations)
    FST = np.zeros((N, N))
    for i, A in enumerate(subpopulations):
        for j, B in enumerate(subpopulations):
            FST[i, j] = df[(df["A"] == A) & (df["B"] == B)]["fst"].values[0]

    im = ax.imshow(FST, cmap="Blues_r")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04) # weird but its perfect
    ax.set(
        xticks=np.arange(N),
        yticks=np.arange(N),
        xticklabels=subpopulations,
        yticklabels=subpopulations,
    )
    ax.tick_params(axis="both", labelsize=6)
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right", rotation_mode="anchor")


def main():
    args = get_args()

    cols = 2
    rows = 1
    fig, ax = plt.subplots(
        rows,
        cols,
        figsize=(args.width, args.height),
        # gridspec_kw={"width_ratios": [1, 1], "height_ratios": [1]},
    )

    plot_fst_heatmap(args.inputs[0], ax[0])
    plot_cdf(args.inputs, args.labels, ax[1])
    ax[0].set_title(r"$\mathbf{A}$", loc="left")
    ax[1].set_title(r"$\mathbf{B}$", loc="left")
    ax[1].set_aspect(0.15, adjustable="box")

    plt.tight_layout()
    # plt.show()
    plt.savefig(args.output, dpi=600)


if __name__ == "__main__":
    main()
