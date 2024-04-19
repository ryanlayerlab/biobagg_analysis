import argparse
import ancestry_helpers as ah
import violin_plot as vp
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    parser=argparse.ArgumentParser(description="evaluate ancestry for 1KG")
    parser.add_argument("--pop", help="path to population file", required=True)
    parser.add_argument("--knn", help="path to knn results file", required=True)
    parser.add_argument("--plink", help="path to plink results file", required=True)
    parser.add_argument("--out_file", help="output file", required=True)
    parser.add_argument("--height", type=int, default=5, help="figure height")
    parser.add_argument("--width",  type=int, default=10, help="figure width")
    parser.add_argument("--y_min",  type=float, help="min y-axis value")
    parser.add_argument("--y_max",  type=float, help="max y-axis value")
    return parser.parse_args()

def get_hit_rates(sample_knn, sample_subpopulations, sub_to_super):
    hit_rates = {}
    for sample in sample_knn:
        s_subpop = sample_subpopulations[sample]
        s_suppop = sub_to_super[sample_subpopulations[sample]]

        hits = sorted(sample_knn[sample], key=lambda x: x[1], reverse=True)

        if s_subpop not in hit_rates:
            hit_rates[s_subpop] = [[] for h in hits]

        if s_suppop not in hit_rates:
            hit_rates[s_suppop] = [[] for h in hits]

        i = 1
        in_subpop = 0
        in_suppop = 0
        subpop_rates = []
        suppop_rates = []
        for hit, score in sample_knn[sample]:
            h_subpop = sample_subpopulations[hit]
            h_suppop = sub_to_super[sample_subpopulations[hit]]
            if h_subpop == s_subpop:
                in_subpop += 1
            if h_suppop in s_suppop:
                in_suppop += 1
            #else:
                #print(sample, s_suppop, hit, h_suppop)
            subpop_rates.append(in_subpop/i)
            suppop_rates.append(in_suppop/i)
            i+=1

        for j in range(len(subpop_rates)):
            hit_rates[s_subpop][j].append(subpop_rates[j])
            hit_rates[s_suppop][j].append(suppop_rates[j])
    return hit_rates

def main():
    args=parse_args()
    population_file = args.pop
    knn_results = args.knn

    sample_subpopulations, sub_to_super, super_to_sub = ah.get_population_maps(population_file)

    sample_knn = ah.get_knn_results(knn_results)

    plink_knn = ah.get_knn_results(args.plink)

    knn_rates = get_hit_rates(sample_knn, sample_subpopulations, sub_to_super)

    plink_rates = get_hit_rates(plink_knn, sample_subpopulations, sub_to_super)

    fig, ax = plt.subplots(1,1, figsize=(args.width, args.height))

    pops = ['ACB', 'ASW', 'BEB', 'CDX', 'CEU',
            'CHB', 'CHS', 'CLM', 'ESN', 'FIN',
            'GBR', 'GIH', 'GWD', 'IBS', 'ITU',
            'JPT', 'KHV', 'LWK', 'MSL', 'MXL',
            'PEL', 'PJL', 'PUR', 'STU', 'TSI',
            'YRI']

    spops = ['EUR', 'AFR', 'EAS', 'SAS', 'AMR']

    spop_colors = {'EUR':'tab:blue',
                   'AFR':'tab:orange',
                   'EAS':'tab:green',
                   'SAS':'tab:red',
                   'AMR':'tab:purple'}


    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
    color_i = 0
    for pop in spops:
        rates = []
        for i in range(len(knn_rates[pop])):
            rates.append(np.mean(knn_rates[pop][i]))
        ax.plot(rates, label= pop, c=spop_colors[pop])
        color_i += 1

    ax.plot([], [], label='HapSiS', c='black')
    ax.plot([], [], label='Plink', c='black', linestyle='--')

    ax.legend(ncol=1, frameon=False)

    color_i = 0
    for pop in spops:
        rates = []
        for i in range(len(plink_rates[pop])):
           rates.append(np.mean(plink_rates[pop][i]))
        ax.plot(rates, label= 'Plink:' + pop,c=spop_colors[pop], linestyle='--')

        color_i += 1

    ax.set_xlabel('Top K')
    ax.set_ylabel('Percent of K in target population')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if args.y_min:
        ax.set_ylim(bottom=args.y_min)
    if args.y_max:
        ax.set_ylim(top=args.y_max)

    tick_positions = range(0, int(ax.get_xlim()[1])+1)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(range(1, int(ax.get_xlim()[1]+2)))

    plt.tight_layout()
    plt.savefig(args.out_file)


if __name__ == '__main__':
    main()
