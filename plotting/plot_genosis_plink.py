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
        query_superpop = plotting.ancestry_helpers.SUB_SUPERPOPULATIONS
        if query_superpop == pop:
            # check matches in top k_subset
            for i in range(k_subset):
                match = top[query][i]
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
    for i in range(k):
        # get percent in population for k = i
        percent_sub, percent_super, percent_out = get_percent_in_pop(top, pop, subpopulations, i)
        subpop_y.append(percent_sub)
        superpop_y.append(percent_super)
        outpop_y.append(percent_out)




def plot_plink_genosis_compare(population_file, genosis_top, plink_top, pop, k):
    subpopulations = plotting.ancestry_helpers.get_subpopulations(population_file)

    # line plot of genosis vs plink
    # x-axis: K
    # y-axis: % in query population
    x = range(1, k+1)
    genosis_y = get_y_values(genosis_top, pop, k, subpopulations)
    plink_y = get_y_values(plink_top, pop, k, subpopulations)


# def main():
#
#
# if __name__ == "__main__":
#     main()