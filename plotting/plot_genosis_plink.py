import matplotlib.pyplot as plt
import numpy as np
import plotting.ancestry_helpers

def get_percent_in_pop(top, pop, subpopulations, k_subset):
    subpop_percents = []
    superpop_percents = []
    outpop_percents = []

    for query in top:
        subpop = 0
        superpop = 0
        outpop = 0
        query_subpop = subpopulations[query]
        query_superpop = plotting.ancestry_helpers.SUB_SUPERPOPULATIONS[query_subpop]
        if query_superpop == pop:
            # check matches in top k_subset
            for i in range(k_subset):
                match, score = top[query][i]
                match_subpop = subpopulations[match]
                match_superpop = plotting.ancestry_helpers.SUB_SUPERPOPULATIONS[match_subpop]
                if match_subpop == query_subpop:
                    subpop += 1
                    superpop += 1
                elif match_superpop == query_superpop:
                    superpop += 1
                else:
                    outpop += 1
        else:
            # skip this query
            continue
        subpop_percents.append(subpop/k_subset)
        superpop_percents.append(superpop/k_subset)
        outpop_percents.append(outpop/k_subset)
    return np.mean(subpop_percents), np.mean(superpop_percents), np.mean(outpop_percents)



def get_y_values(top, pop, k, subpopulations):
    subpop_y = []
    superpop_y = []
    outpop_y = []
    for i in range(1, k+1):
        # get percent in population for k = i
        percent_sub, percent_super, percent_out = get_percent_in_pop(top, pop, subpopulations, i)
        subpop_y.append(percent_sub)
        superpop_y.append(percent_super)
        outpop_y.append(percent_out)

    return subpop_y, superpop_y, outpop_y





def plot_plink_genosis_compare(population_file,
                               genosis_top,
                               plink_dst,
                               plink_pihat,
                               plink_kin,
                               pop, k,
                               color):

    subpopulations = plotting.ancestry_helpers.get_subpopulations(population_file)

    linewidth = 2
    subpop_color = 'black'
    superpop_color = color

    # line plot of genosis vs plink
    # x-axis: K
    # y-axis: % in query population
    x = range(1, k+1)
    genosis_y = get_y_values(genosis_top, pop, k, subpopulations)
    plink_dst_y = get_y_values(plink_dst, pop, k, subpopulations)
    plink_pihat_y = get_y_values(plink_pihat, pop, k, subpopulations)
    plink_kin_y = get_y_values(plink_kin, pop, k, subpopulations)

    # 3 figures, 1 for each plink metric
    fig, axs = plt.subplots(3, figsize=(10, 15), dpi=100)
    fig.suptitle('GenoSiS vs Plink\n' + pop, fontsize=20, fontweight='bold')

    axs[0].plot(x, genosis_y[0], label='Genosis Subpop', color=subpop_color, linewidth=linewidth)
    axs[0].plot(x, genosis_y[1], label='Genosis Superpop', color=superpop_color, linewidth=linewidth)
    # axs[0].plot(x, genosis_y[2], label='Genosis Outpop', color='red')
    axs[0].plot(x, plink_dst_y[0], label='plink Subpop', color=subpop_color, linestyle='dashed', linewidth=linewidth)
    axs[0].plot(x, plink_dst_y[1], label='plink Superpop', color=superpop_color, linestyle='dashed', linewidth=linewidth)
    # axs[0].plot(x, plink_dst_y[2], label='Plink Outpop', color='red', linestyle='dashed')
    axs[0].set_title('plink DST', fontsize=15, fontweight='bold')
    axs[0].set_ylabel('% in Population')
    axs[0].set_xticks([])
    axs[0].set_ylim(0, 1)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)

    axs[1].plot(x, genosis_y[0], label='Genosis Subpop', color=subpop_color, linewidth=linewidth)
    axs[1].plot(x, genosis_y[1], label='Genosis Superpop', color=superpop_color, linewidth=linewidth)
    # axs[1].plot(x, genosis_y[2], label='Genosis Outpop', color='red')
    axs[1].plot(x, plink_pihat_y[0], label='plink Subpop', color=subpop_color, linestyle='dashed', linewidth=linewidth)
    axs[1].plot(x, plink_pihat_y[1], label='plink Superpop', color=superpop_color, linestyle='dashed', linewidth=linewidth)
    # axs[1].plot(x, plink_pihat_y[2], label='Plink Outpop', color='red', linestyle='dashed')
    axs[1].set_title('plink pi-hat', fontsize=15, fontweight='bold')
    axs[1].set_ylabel('% in Population')
    axs[1].set_xticks([])
    axs[1].set_ylim(0, 1)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    axs[1].spines['bottom'].set_visible(False)

    axs[2].plot(x, genosis_y[0], label='Genosis Subpop', color=subpop_color, linewidth=linewidth)
    axs[2].plot(x, genosis_y[1], label='Genosis Superpop', color=superpop_color, linewidth=linewidth)
    # axs[2].plot(x, genosis_y[2], label='Genosis Outpop', color='red')
    axs[2].plot(x, plink_kin_y[0], label='plink Subpop', color=subpop_color, linestyle='dashed', linewidth=linewidth)
    axs[2].plot(x, plink_kin_y[1], label='plink Superpop', color=superpop_color, linestyle='dashed', linewidth=linewidth)
    # axs[2].plot(x, plink_kin_y[2], label='Plink Outpop', color='red', linestyle='dashed')
    axs[2].set_title('plink kinship', fontsize=15, fontweight='bold')
    axs[2].set_xlabel('K')
    axs[2].set_ylabel('% in Population')
    # set x ticks from 1-20 with step 5
    axs[2].set_xticks(np.arange(1, k+1, 5))
    axs[2].set_ylim(0, 1)
    axs[2].spines['top'].set_visible(False)
    axs[2].spines['right'].set_visible(False)

    # one legend for all subplots
    # axs[0].legend()
    # tight layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.legend()

    plt.savefig('genosis_plink_compare_' + pop + '.png')





    # plt.figure(figsize=(15, 7), dpi=100)
    # plt.plot(x, genosis_y[0], label='Genosis Subpop', color='green')
    # plt.plot(x, genosis_y[1], label='Genosis Superpop', color='blue')
    # plt.plot(x, genosis_y[2], label='Genosis Outpop', color='red')
    # plt.plot(x, plink_y[0], label='Plink Subpop', color='green', linestyle='dashed')
    # plt.plot(x, plink_y[1], label='Plink Superpop', color='blue', linestyle='dashed')
    # plt.plot(x, plink_y[2], label='Plink Outpop', color='red', linestyle='dashed')
    # plt.xlabel('K')
    # plt.ylabel('% in Population')
    # plt.legend()
    #
    # # remove spines
    # plt.gca().spines['top'].set_visible(False)
    # plt.gca().spines['right'].set_visible(False)
    # plt.title('Genosis vs Plink ' + plink_metric)
    # plt.savefig('genosis_plink_compare_' + plink_metric + '.png')





# def main():
#
#
# if __name__ == "__main__":
#     main()