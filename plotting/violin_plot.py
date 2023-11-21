from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import ancestry_helpers as ah


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
    """
    Plot violin plot for a subpopulation

    :param population_data: counts for knn results
    :param sub_to_super: mapping of subpopulation to superpopulation
    :param super_to_sub: mapping of superpopulation to subpopulations
    :param subpop_name: name of subpopulation that is being plotted
    :param png_dir: directory to write png files
    :return:
    """
    # colors to match current slide deck
    super_population_colors = {'AFR': '#fe9d57',
                               'AMR': '#ba75ff',
                               'EAS': '#eb4690',
                               'EUR': '#728eff',
                               'SAS': '#f7d06b'}
    # get full names
    subpop_full_names = ah.get_subpop_full_names()

    subpop_codes = []
    subpop_colors = []

    # make violin plot
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 13)
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
    subpop_title = "Query Subpopulation:\n" + subpop_full_names[subpop_name] + "\n(" + subpop_name + ")"
    ax.set_title(subpop_title, fontsize=25, fontweight='bold')

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

