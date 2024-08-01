import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import os

import numpy as np
import seaborn as sns
import scipy.stats as stats
import statistics
import pandas as pd

from plotting import ancestry_helpers

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestry_file', type=str, help='file with ancestry labels', required=True)
    parser.add_argument('--density_dir', type=str, help='directory with density files', required=True)
    parser.add_argument('--distance_dir', type=str, help='directory with distance files', required=True)
    parser.add_argument('--colors', type=str, help='file with colors for each pop', required=True)
    parser.add_argument('--out', type=str, help='output directory', required=True)

    return parser.parse_args()

def get_sample_densities(chrm_dir, density_file, sample_densities):
    segment = density_file.split('.')[1].replace('segment', '')
    seg_idx = int(segment)
    with open(chrm_dir + density_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            sample = line[0]
            density = float(line[1])
            # if density < 100:
            #     print(sample, density, segment)
            sample_densities[sample][seg_idx] = density
    return sample_densities

def get_segment_densities(chrm_dir, density_file, segment_densities):
    segment_idx = density_file.split('.')[1].replace('segment', '')
    with open(chrm_dir + density_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            sample = line[0]
            density = float(line[1])
            segment_densities[segment_idx].append(density)
    return segment_densities

def get_sample_distances(distance_file):
    enc_distances = defaultdict(dict)
    emb_distances = defaultdict(dict)
    with open(distance_file, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip().split()
            sampleA = line[0]
            sampleB = line[1]
            enc_dist = float(line[2])
            emb_dist = float(line[3])
            enc_distances[sampleA][sampleB] = enc_dist
            emb_distances[sampleA][sampleB] = emb_dist
    return enc_distances, emb_distances

def get_single_r2(enc_distances, emb_distances):
    '''
    Calculate the r^2 value for a single sample
    @param enc_distances:
    @param emb_distances:
    @return:
    '''

    sample_r2 = {}
    # compare one sample to all others and get the single samples r^2 value
    for sampleA in enc_distances:
        encoding_distances = []
        embedding_distances = []
        for sampleB in enc_distances[sampleA]:
            encoding_distances.append(enc_distances[sampleA][sampleB])
            embedding_distances.append(emb_distances[sampleA][sampleB])

        r2 = stats.pearsonr(encoding_distances, embedding_distances)[0] ** 2
        sample_r2[sampleA] = r2

    return sample_r2

def plot_density_by_ancestry(sample_densities,
                             sample_subpopulations,
                             sub_to_super,
                             colors,
                             chrm,
                             out):

    # one subplot per superpopulation
    superpops = sorted(set(sub_to_super.values()))
    num_superpops = len(superpops)

    # title = ('Density by Ancestry'
    #          '\nChromosome ' + chrm + ', segments (0,' + str(len(sample_densities['HG00096_0'])) + ')')
    # title = ('Chromosome ' + chrm + ' segments 1-' + str(len(sample_densities['HG00096_0'])+1))

    fig, axes = plt.subplots(1, num_superpops,
                             figsize=(27, 4),
                             sharex=True, sharey=True,
                             dpi=300)

    for superpop, ax in zip(superpops, axes):
        # get sample IDs for this superpopulation
        color = colors[superpop]
        pop_densities = []
        samples = [sample for sample, subpop in sample_subpopulations.items() if sub_to_super[subpop] == superpop]
        for sample in samples:
            # haplotypes
            sample_0 = sample + '_0'
            sample_1 = sample + '_1'
            pop_densities.extend(sample_densities[sample_0])
            pop_densities.extend(sample_densities[sample_1])
            # if 0 in sample_densities[sample_0] or 0 in sample_densities[sample_1]:
            #     print(superpop, sample_0, sample_1)

        # plot densities with histograms
        sns.histplot(pop_densities,
                     ax=ax,
                     color=color,
                     kde=False,
                     bins=20)
        ax.set_title(superpop, fontsize=20, color='black')
        ax.set_xlabel('Number of non-reference alleles')
        ax.set_ylabel('Frequency')

        # add a textbox that shows mean, median, and mode densities
        mean_density = statistics.mean(pop_densities)
        median_density = statistics.median(pop_densities)
        mode_density = statistics.mode(pop_densities)

        ax.text(0.1, 0.95, 'Mean: ' + str(round(mean_density, 2)) + '\n'
                            'Median: ' + str(round(median_density, 2)) + '\n'
                            'Mode: ' + str(round(mode_density, 2)),
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax.transAxes,
                fontsize=12,
                color='black')


        # remove spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # fig.suptitle(title, fontsize=20, fontweight='bold')

    plt.tight_layout()
    fig_name = out + 'density_by_ancestry_chrm' + chrm + '.png'
    plt.savefig(fig_name)

def plot_density_by_segment(segment_densities,
                            sample_densities,
                            sample_subpopulations,
                            sub_to_super,
                            colors,
                            chrm,
                            out):


    superpops = sorted(set(sub_to_super.values()))
    num_superpops = len(superpops)

    sorted_segments = sorted(segment_densities.keys(), key=lambda x: int(x))

    # heatmap where x is segment idx and y is superpop
    print('Plotting heatmap of median density by segment...', chrm)
    fig, ax = plt.subplots(figsize=(8, 20), dpi=300)
    df = {superpop: [] for superpop in superpops}
    for segment in sorted_segments:
        densities = segment_densities[segment]
        for superpop in superpops:
            pop_densities = []
            samples = [sample for sample, subpop in sample_subpopulations.items() if sub_to_super[subpop] == superpop]
            for sample in samples:
                # haplotypes
                sample_0 = sample + '_0'
                sample_1 = sample + '_1'
                pop_densities.append(sample_densities[sample_0][int(segment)])
                pop_densities.append(sample_densities[sample_1][int(segment)])

            # get median density for this segment
            median_density = statistics.median(pop_densities)
            df[superpop].append(median_density)
            x = 0

    densities = pd.DataFrame(df, index=sorted_segments)
    sns.heatmap(densities,
                cmap='Grays',
                ax=ax,
                annot=False, fmt='.2f')

    # label colorbar
    cbar = ax.collections[0].colorbar
    cbar.set_label('Median Density', fontsize=15, fontweight='bold')

    ax.set_xticks(range(num_superpops))
    ax.set_xticklabels(superpops, rotation=45)
    ax.set_yticks(range(len(sorted_segments)))
    ax.set_yticklabels(sorted_segments)
    ax.set_xlabel('Superpopulation')
    ax.set_ylabel('Segment')

    title = ('Median Density by Segment'
             '\nChromosome ') + chrm + ', segments 0-' + str(len(sorted_segments))
    fig.suptitle(title, fontsize=20, fontweight='bold')
    fig_name = out + 'density_by_segment_chrm' + chrm + '_heatmap.png'
    plt.tight_layout()
    plt.savefig(fig_name)

    # histograms for each segment
    print('Plotting histograms of density by segment...', chrm)
    num_cols = 10
    fig, axes = plt.subplots(len(sorted_segments) // num_cols + 1, num_cols,
                             figsize=(20, 20),
                             sharex=True, sharey=True,
                             dpi=300)

    for segment, ax in zip(sorted_segments, axes.flatten()):
        seg_dens = segment_densities[segment]

        # plot densities with histograms
        sns.histplot(seg_dens,
                        ax=ax,
                        color='black',
                        kde=False,
                        bins=50)
        ax.set_title('Segment ' + segment, fontsize=15, fontweight='bold')
        ax.set_xlabel('Density')
        ax.set_ylabel('Frequency')
        # log y-axis
        ax.set_yscale('log')

        # remove spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # remove empty subplots
    for i in range(len(sorted_segments), len(axes.flatten())):
        fig.delaxes(axes.flatten()[
            i])


    title = ('Density by Segment'
             '\nChromosome ') + chrm + ', segments 0-' + str(len(sorted_segments))
    fig.suptitle(title, fontsize=20, fontweight='bold')
    fig_name = out + 'density_by_segment_chrm' + chrm + '_hist.png'
    plt.tight_layout()
    plt.savefig(fig_name)

    # plot median density by mode density per segment, colored by superpopulation
    print('Plotting scatterplot of median density vs. mode density by segment...', chrm)
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    for superpop in superpops:
        color = colors[superpop]
        other_densities = []
        median_densities = []
        for segment in sorted_segments:
            pop_densities = []
            samples = [sample for sample, subpop in sample_subpopulations.items() if sub_to_super[subpop] == superpop]
            for sample in samples:
                # haplotypes
                sample_0 = sample + '_0'
                sample_1 = sample + '_1'
                pop_densities.append(sample_densities[sample_0][int(segment)])
                pop_densities.append(sample_densities[sample_1][int(segment)])

            # get mode density for this segment
            other_density = statistics.mode(pop_densities)
            other_densities.append(other_density)
            # get median density for this segment
            median_density = statistics.median(pop_densities)
            median_densities.append(median_density)

        ax.scatter(other_densities, median_densities, color=color, label=superpop, alpha=0.5)
    ax.set_xlabel('Mode Density')
    ax.set_ylabel('Median Density')
    ax.legend()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    title = (('Median Density vs. Mode Density by Segment'
              '\nChromosome ') + chrm + ', segments 0-' + str(len(sorted_segments)))
    fig.suptitle(title, fontsize=20, fontweight='bold')
    fig_name = out + 'density_by_segment_chrm' + chrm + '_scatter.png'
    plt.tight_layout()
    plt.savefig(fig_name)


def plot_encoding_embedding(good_segment_enc, good_segment_emb,
                            good_segments,
                            sample_subpopulations,
                            sub_to_super,
                            colors,
                            chrm,
                            out):
    # plot encoding vs. embedding distances by ancestry
    print('Plotting encoding vs. embedding distances by ancestry...', chrm)
    superpops = sorted(set(sub_to_super.values()))
    num_superpops = len(superpops)

    # separate plot per segment
    for seg_idx in good_segments:
        fig, axes = plt.subplots(num_superpops, 1,
                                 figsize=(6, 10),
                                 sharex=True, sharey=True,
                                 dpi=300)
        for superpop in superpops:
            # get sample IDs for this superpopulation
            samples = [sample for sample, subpop in sample_subpopulations.items() if sub_to_super[subpop] == superpop]
            encoding_distances = []
            embedding_distances = []
            for sampleA in samples:
                for sampleB in samples:
                    # hap_0 = sampleA + '_0'
                    # hap_1 = sampleA + '_1'
                    sampleA_0 = sampleA + '_0'
                    sampleA_1 = sampleA + '_1'
                    sampleB_0 = sampleB + '_0'
                    sampleB_1 = sampleB + '_1'

                    if sampleA_0 in good_segment_enc[seg_idx] and sampleB_0 in good_segment_enc[seg_idx][sampleA_0]:
                        encoding_distances.append(good_segment_enc[seg_idx][sampleA_0][sampleB_0])
                        embedding_distances.append(good_segment_emb[seg_idx][sampleA_0][sampleB_0])
                    if sampleA_0 in good_segment_enc[seg_idx] and sampleB_1 in good_segment_enc[seg_idx][sampleA_0]:
                        encoding_distances.append(good_segment_enc[seg_idx][sampleA_0][sampleB_1])
                        embedding_distances.append(good_segment_emb[seg_idx][sampleA_0][sampleB_1])
                    if sampleA_1 in good_segment_enc[seg_idx] and sampleB_0 in good_segment_enc[seg_idx][sampleA_1]:
                        encoding_distances.append(good_segment_enc[seg_idx][sampleA_1][sampleB_0])
                        embedding_distances.append(good_segment_emb[seg_idx][sampleA_1][sampleB_0])
                    if sampleA_1 in good_segment_enc[seg_idx] and sampleB_1 in good_segment_enc[seg_idx][sampleA_1]:
                        encoding_distances.append(good_segment_enc[seg_idx][sampleA_1][sampleB_1])
                        embedding_distances.append(good_segment_emb[seg_idx][sampleA_1][sampleB_1])

            if len(encoding_distances) == 0:
                continue
            if len(embedding_distances) == 0:
                continue

            ax = axes[superpops.index(superpop)]
            ax.scatter(encoding_distances, embedding_distances, alpha=0.5, color=colors[superpop])
            ax.set_title(superpop, fontsize=15, color=colors[superpop], fontweight='bold')
            ax.set_xlabel('Encoding Distance')
            ax.set_ylabel('Embedding Distance')

            # get r^2 and correlation coefficient
            r2 = stats.pearsonr(encoding_distances, embedding_distances)[0] ** 2
            cc = np.corrcoef(encoding_distances, embedding_distances)[0, 1]
            ax.text(0.2, 0.9,
                    'R^2: ' + str(round(r2, 2)) + '\n'
                    'Corr. Coeff.: ' + str(round(cc, 2)),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes,
                    fontsize=10,
                    color='black')

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        fig.suptitle('Encoding vs. Embedding Distances\nChromosome ' + chrm + ', segment ' + str(seg_idx),
                     fontsize=20, fontweight='bold')
        fig_name = out + 'encoding_embedding_chrm' + chrm + '_seg' + str(seg_idx) + '.png'
        plt.tight_layout()
        plt.savefig(fig_name)
        plt.close()


def plot_eur_afr_density(sample_densities,
                         sample_subpopulations,
                         sub_to_super,
                         chrm,
                         out):

    # x-axis: EUR density
    # y-axis: AFR density
    # plot scatter plot for all segments

    # get sample IDs for EUR and AFR
    eur_samples = [sample for sample, subpop in sample_subpopulations.items() if sub_to_super[subpop] == 'EUR']
    afr_samples = [sample for sample, subpop in sample_subpopulations.items() if sub_to_super[subpop] == 'AFR']

    eur_segment_densities = {seg: [] for seg in range(len(sample_densities['HG00096_0']))}
    for eur_sample in eur_samples:
        for seg in range(len(sample_densities[eur_sample + '_0'])):
            eur_segment_densities[seg].append(sample_densities[eur_sample + '_0'][seg])
            eur_segment_densities[seg].append(sample_densities[eur_sample + '_1'][seg])

    afr_segment_densities = {seg: [] for seg in range(len(sample_densities['HG00096_0']))}
    for afr_sample in afr_samples:
        for seg in range(len(sample_densities[afr_sample + '_0'])):
            afr_segment_densities[seg].append(sample_densities[afr_sample + '_0'][seg])
            afr_segment_densities[seg].append(sample_densities[afr_sample + '_1'][seg])

    # get mean eur and afr densities for each segment
    eur_densities = []
    afr_densities = []
    for seg in eur_segment_densities:
        eur_densities.append(statistics.mean(eur_segment_densities[seg]))
        afr_densities.append(statistics.mean(afr_segment_densities[seg]))

    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    ax.scatter(eur_densities, afr_densities, alpha=0.5)

    # label each point with segment number
    for i, seg in enumerate(eur_segment_densities):
        ax.text(eur_densities[i], afr_densities[i], str(seg),
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=10,
                fontweight='bold',
                color='black',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

    # label segments 81 and 9 for clarity
    # ax.text(eur_densities[81], afr_densities[81], '81',
    #         horizontalalignment='center',
    #         verticalalignment='center',
    #         fontsize=10,
    #         fontweight='bold',
    #         color='white',
    #         bbox=dict(facecolor='black', edgecolor='black', boxstyle='round,pad=0.5'))
    # ax.text(eur_densities[9], afr_densities[9], '9',
    #         horizontalalignment='center',
    #         verticalalignment='center',
    #         fontsize=10,
    #         fontweight='bold',
    #         color='white',
    #         bbox=dict(facecolor='black', edgecolor='black', boxstyle='round,pad=0.5'))

    ax.set_xlabel('EUR Density')
    ax.set_ylabel('AFR Density')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    title = ('EUR vs. AFR Density'
                '\nChromosome ' + chrm)
    fig.suptitle(title, fontsize=20, fontweight='bold')
    fig_name = out + 'eur_afr_density_chrm' + chrm + '.png'
    plt.tight_layout()
    plt.savefig(fig_name)

def plot_r2_density(good_segment_r2,
                    sample_densities,
                    good_segments,
                    sample_subpopulations,
                    sub_to_super,
                    colors,
                    chrm,
                    out):
    # plot r2 vs. density by ancestry
    # x = density
    # y = r2

    superpops = sorted(set(sub_to_super.values()))
    num_superpops = len(superpops)

    for seg_idx in good_segments:
        fig, axes = plt.subplots(1, num_superpops,
                                figsize=(15, 5),
                                sharex=True, sharey=True,
                                dpi=300)

        for superpop in superpops:
            # get sample IDs for this superpopulation
            samples = [sample for sample, subpop in sample_subpopulations.items() if sub_to_super[subpop] == superpop]
            densities = []
            r2_values = []
            for sample in samples:
                # haplotypes
                sample_0 = sample + '_0'
                sample_1 = sample + '_1'
                try:
                    r2_values.append(good_segment_r2[seg_idx][sample_1])
                    r2_values.append(good_segment_r2[seg_idx][sample_0])
                    densities.append(sample_densities[sample_1][seg_idx])
                    densities.append(sample_densities[sample_0][seg_idx])
                except KeyError:
                    continue

            if len(densities) == 0:
                continue
            if len(r2_values) == 0:
                continue

            ax = axes[superpops.index(superpop)]
            ax.scatter(densities, r2_values, alpha=0.5, color=colors[superpop])
            ax.set_title(superpop, fontsize=15, color=colors[superpop], fontweight='bold')
            ax.set_xlabel('Density')
            ax.set_ylabel('r^2')

            # y scale 0-1
            ax.set_ylim(-0.2, 1.0)

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        fig.suptitle('r^2 vs. Density\nChromosome ' + chrm + ', segment ' + str(seg_idx),
                     fontsize=20, fontweight='bold')

        fig_name = out + 'r2_density_chrm' + chrm + '_seg' + str(seg_idx) + '.png'
        plt.tight_layout()
        plt.savefig(fig_name)
        plt.close()


def read_colors(color_file):
    colors = {}
    with open(color_file, 'r') as f:
        for line in f:
            line = line.strip().split(',')
            colors[line[0]] = line[1]
    return colors

def main():
    args = get_args()
    ancestry_file = args.ancestry_file
    density_dir = args.density_dir
    distance_dir = args.distance_dir
    color_files = args.colors
    out = args.out

    # read ancestry labels
    sample_subpopulations, sub_to_super, super_to_sub = ancestry_helpers.get_population_maps(ancestry_file)

    colors = read_colors(color_files)

    # read density files
    chroms = [str(i) for i in range(8, 9)]
    good_segments = [76, 81, 97, 108, 136]

    num_samples = 6406

    for chrm in chroms:
        num_segments = 0
        max_segments = 300
        sample_densities = defaultdict(lambda: [0] * max_segments)
        segment_densities = defaultdict(list)
        good_segment_enc = {gs: defaultdict(dict) for gs in good_segments}
        good_segment_emb = {gs: defaultdict(dict) for gs in good_segments}
        good_segment_r2 = {gs: defaultdict(dict) for gs in good_segments}

        # open chromosome directory for densities
        dens_chrm_dir = density_dir + 'chrm' + chrm + '_dens/'
        for density_file in os.listdir(dens_chrm_dir):
            if 'segment' in density_file:
                sample_densities = get_sample_densities(dens_chrm_dir, density_file, sample_densities)
                segment_densities = get_segment_densities(dens_chrm_dir, density_file, segment_densities)
                num_segments += 1
            else:
                # ignore
                continue

        # open chrm_dir for distances
        dist_chrm_dir = distance_dir + 'chrm' + chrm + '_dist/'
        for distance_file in os.listdir(dist_chrm_dir):
            if 'segment' in distance_file:
                seg_idx = int(distance_file.split('.')[1].replace('segment', ''))
                enc_distances, emb_distances = get_sample_distances(dist_chrm_dir + distance_file)
                sample_r2 = get_single_r2(enc_distances, emb_distances)
                good_segment_enc[seg_idx] = enc_distances
                good_segment_emb[seg_idx] = emb_distances
                good_segment_r2[seg_idx] = sample_r2

            else:
                # ignore
                continue


        # resize sample densities to be length of num_segments
        for sample, densities in sample_densities.items():
            # sample_densities[sample] = densities[:num_segments]
            sample_densities[sample] = [densities[76]]

        # resize segment densities to be length of num_segments
        for segment, densities in segment_densities.items():
            segment_densities[segment] = densities[:num_samples]

        # plot densities

        print('Plotting densities by ancestry...', chrm)
        plot_density_by_ancestry(sample_densities,
                                 sample_subpopulations,
                                 sub_to_super,
                                 colors,
                                 chrm,
                                 out)

        # print('Plotting densities by segment...', chrm)
        # plot_density_by_segment(segment_densities,
        #                         sample_densities,
        #                         sample_subpopulations,
        #                         sub_to_super,
        #                         colors,
        #                         chrm,
        #                         out)

        # print('Plotting encoding by embedding distances...', chrm)
        # plot_encoding_embedding(good_segment_enc, good_segment_emb,
        #                         good_segments,
        #                         sample_subpopulations,
        #                         sub_to_super,
        #                         colors,
        #                         chrm,
        #                         out)

        # print('Plotting EUR vs. AFR...', chrm)
        # plot_eur_afr_density(sample_densities,
        #                      sample_subpopulations,
        #                      sub_to_super,
        #                      chrm,
        #                      out)

        # print('Plotting r2 by density...', chrm)
        # plot_r2_density(good_segment_r2,
        #                 sample_densities,
        #                 good_segments,
        #                 sample_subpopulations,
        #                 sub_to_super,
        #                 colors,
        #                 chrm,
        #                 out)






if __name__ == '__main__':
    main()