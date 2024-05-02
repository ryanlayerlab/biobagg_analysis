import argparse

import read_plink
import read_genosis
import plotting.plot_genosis_plink

def parse_args():
    parser = argparse.ArgumentParser(description="Compare genosis and plink output")
    parser.add_argument("-g", "--genosis", help="genosis top hits file")
    parser.add_argument("-pg", "--plinkgenome", help="plink genome output file")
    parser.add_argument("-pk", "--plinkkin", help="plink genome output file")
    parser.add_argument("-k", "--top_k", help="top K hits to return", default=20)
    parser.add_argument("-a", "--ancestry", help="ancestry file")
    parser.add_argument("-p", "--pop", help="population to investigate")
    parser.add_argument("-c", "--color", help="color for plot")
    return parser.parse_args()

def get_plink_top_hits(plink_dict, K):
    # for each query, return top K samples
    top_hits_dict = {}
    for query in plink_dict:
        top_hits_samples = sorted(plink_dict[query].items(), key=lambda x: x[1], reverse=True)[:K]
        top_hits_dict[query] = top_hits_samples
    return top_hits_dict


def main():
    args = parse_args()

    plink_DST_score_dict = read_plink.get_pairwise_DST_score_dict(args.plinkgenome)
    plink_pihat_score_dict = read_plink.get_pairwise_pihat_score_dict(args.plinkgenome)
    plink_kin_score_dict = read_plink.get_pairwise_kin_score_dict(args.plinkkin)

    genosis_K_dict = read_genosis.get_top_hits_dict(args.genosis)
    plink_DST_K_dict = get_plink_top_hits(plink_DST_score_dict, int(args.top_k))
    plink_pihat_K_dict = get_plink_top_hits(plink_pihat_score_dict, int(args.top_k))
    plink_kin_K_dict = get_plink_top_hits(plink_kin_score_dict, int(args.top_k))

    plotting.plot_genosis_plink.plot_plink_genosis_compare(args.ancestry,
                                                           genosis_K_dict,
                                                           plink_DST_K_dict,
                                                           plink_pihat_K_dict,
                                                           plink_kin_K_dict,
                                                           args.pop,
                                                           int(args.top_k), args.color)

    # plotting.plot_genosis_plink.plot_plink_genosis_compare(args.ancestry,
    #                                                        genosis_K_dict,
    #                                                        plink_pihat_K_dict,
    #                                                        args.pop,
    #                                                        int(args.top_k),
    #                                                        "pihat")
    #
    # plotting.plot_genosis_plink.plot_plink_genosis_compare(args.ancestry,
    #                                                        genosis_K_dict,
    #                                                        plink_kin_K_dict,
    #                                                        args.pop,
    #                                                        int(args.top_k),
    #                                                        "kin")



if __name__ == "__main__":
    main()