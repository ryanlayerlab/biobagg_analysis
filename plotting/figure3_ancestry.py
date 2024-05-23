import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import ks_2samp
from scipy.stats import norm

def parse_args():
    parser = argparse.ArgumentParser()
    # GenoSiS Scores
    parser.add_argument('--genosis', help='1KG genosis scores', required=True)
    # plink Scores
    parser.add_argument('--dst', help='1KG plink DST scores', required=True)
    parser.add_argument('--pihat', help='1KG plink pi-hat scores', required=True)
    parser.add_argument('--kinship', help='1KG plink kinship scores', required=True)
    # Output
    parser.add_argument('--png', help='Output png file', required=True)

    return parser.parse_args()

def main():
    args = parse_args()

if __name__ == '__main__':
    main()

