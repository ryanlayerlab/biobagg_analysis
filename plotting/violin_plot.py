from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def plot_subpopulation(sample_subpopulations,
                       sub_to_super,
                       super_to_sub,
                       knn_results,
                       png_dir):
    """
    Plot violin plots for each subpopulation

    :param sample_subpopulation: a dictionary mapping sample IDs to subpopulation codes
    :param sub_to_super: a dictionary mapping subpopulation codes to superpopulation codes
    :param knn_results: KNN results file
    :param png_dir: directory to write png files
    """
    subpopulations = list(set(sample_subpopulations.values()))
    # data: [supopulation: [subpopulation: count]]
    subpopulation_counts = defaultdict(dict)
    for sample in sample_subpopulations:
        sample_subpopulation_code = sample_subpopulations[sample]
        # get all match subpopulations
        for match in knn_results[sample]:
            matchID = match[0]
            match_score = match[1]
            match_subpopulation_code = sample_subpopulations[matchID]
            match_score = match_score
            try:
                subpopulation_counts[sample_subpopulation_code][match_subpopulation_code].append(match_score)
            except KeyError:
                subpopulation_counts[sample_subpopulation_code][match_subpopulation_code] = [match_score]

    # plot violin plots for each subpopulation
    for subpop in subpopulation_counts:
        plot_violin_plot(subpopulation_counts[subpop],
                         sub_to_super,
                         super_to_sub,
                         subpop,
                         png_dir)
    pass


def plot_violin_plot(population_data,
                     sub_to_super,
                     super_to_sub,
                     subpop_name,
                     png_dir):
    super_population_colors = {'AFR': '#fe9d57',
                               'AMR': '#ba75ff',
                               'EAS': '#eb4690',
                               'EUR': '#728eff',
                               'SAS': '#f7d06b'}
    subpop_full_names = get_subpop_full_names()

    subpop_codes = []
    subpop_colors = []
    # make violin plot
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 10)
    violin_data = []
    for super_pop in super_to_sub:
        if super_pop is not None:
            # sort in alphabetical order
            subpops = sorted(list(super_to_sub[super_pop]))
            for sub_pop in subpops:
                if sub_pop is not None:
                    subpop_codes.append(sub_pop)
                    subpop_colors.append(super_population_colors[sub_to_super[sub_pop]])
                    try:
                        violin_data.append(population_data[sub_pop])
                    except KeyError:
                        violin_data.append([0])
                else:  # ignore None
                    pass
        else:  # ignore None
            pass

    # plot!
    parts = ax.violinplot(violin_data, showextrema=False, showmeans=True)

    # add labels
    ax.set_xticks(range(1, len(subpop_codes) + 1))
    ax.set_xticklabels(subpop_codes, rotation=45, fontsize=15)
    ax.set_ylabel('KNN popcount scores', fontsize=20)

    # add title
    subpop_title = "Query Subpopulation:\n" + subpop_full_names[subpop_name]
    ax.set_title(subpop_title, fontsize=30, fontweight='bold')

    # color each x violin plot with superpopulation color
    for i in range(len(subpop_codes)):
        ax.get_xticklabels()[i].set_color(subpop_colors[i])
    parts['cmeans'].set_color(subpop_colors)
    c = 0
    for b in parts['bodies']:
        b.set_color(subpop_colors[c])
        b.set_alpha(0.5)
        b.set_linewidth(2)
        c += 1

    # format plot
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_tick_params(width=0)
    # legend for superpopulation colors
    legend_elements = []
    for super_pop in super_population_colors:
        legend_elements.append(Line2D([0], [0], color=super_population_colors[super_pop], lw=20, label=super_pop))
    # ax legend no border
    ax.legend(handles=legend_elements, fontsize=20, frameon=False)

    # save figure
    png_name = png_dir + subpop_name + '.png'
    plt.savefig(png_name)
    plt.close()


def get_subpop_full_names():
    return {'ACB': 'African Caribbean in Barbados',
            'ASW': 'Americans of African Ancestry in SW USA',
            'BEB': 'Bengali in Bangladesh',
            'CDX': 'Chinese Dai in Xishuangbanna, China',
            'CEU': 'Utah Residents (CEPH) with Northern and Western European Ancestry',
            'CHB': 'Han Chinese in Bejing, China',
            'CHS': 'Southern Han Chinese',
            'CLM': 'Colombians from Medellin, Colombia',
            'ESN': 'Esan in Nigeria',
            'FIN': 'Finnish in Finland',
            'GBR': 'British in England and Scotland',
            'GIH': 'Gujarati Indian in Houston, TX',
            'GWD': 'Gambian in Western Divisions in the Gambia',
            'IBS': 'Iberian Population in Spain',
            'ITU': 'Indian Telugu in the UK',
            'JPT': 'Japanese in Tokyo, Japan',
            'KHV': 'Kinh in Ho Chi Minh City, Vietnam',
            'LWK': 'Luhya in Webuye, Kenya',
            'MSL': 'Mende in Sierra Leone',
            'MXL': 'Mexican Ancestry from Los Angeles USA',
            'PEL': 'Peruvians from Lima, Peru',
            'PJL': 'Punjabi in Lahore, Pakistan',
            'PUR': 'Puerto Ricans from Puerto Rico',
            'STU': 'Sri Lankan Tamil in the UK',
            'TSI': 'Toscani in Italia',
            'YRI': 'Yoruba in Ibadan, Nigeria'}
