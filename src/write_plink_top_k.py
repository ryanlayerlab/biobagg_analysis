import argparse
import read_plink as rp

def parse_args():
    parser = argparse.ArgumentParser(description="Writes top K samples and scores from plink files")
    parser.add_argument("-pg", "--plinkgenome", help="plink genome output file")
    parser.add_argument("-pk", "--plinkkin", help="plink genome output file")
    parser.add_argument("-k", "--knn", help="top K hits to return", default=20)
    parser.add_argument("-o", "--out", help="output directory")
    return parser.parse_args()

def get_plink_top_hits(plink_dict, K):
    # for each query, return top K samples
    top_hits_dict = {}
    for query in plink_dict:
        top_hits_samples = sorted(plink_dict[query].items(), key=lambda x: x[1], reverse=True)[:K]
        top_hits_dict[query] = top_hits_samples
    return top_hits_dict

def write_top_k(top_hits_dict, output_file):
    with open(output_file, 'w') as f:
        f.write("query match,score\n")
        for query in top_hits_dict:
            f.write(f"{query} ")
            for sample in top_hits_dict[query]:
                f.write(f"{sample[0]},{sample[1]} ")
            f.write("\n")
    f.close()

def main():
    args = parse_args()

    plink_DST_score_dict = rp.get_pairwise_DST_score_dict(args.plinkgenome)
    plink_pihat_score_dict = rp.get_pairwise_pihat_score_dict(args.plinkgenome)
    plink_kin_score_dict = rp.get_pairwise_kin_score_dict(args.plinkkin)

    plink_DST_K_dict = get_plink_top_hits(plink_DST_score_dict, int(args.knn))
    plink_pihat_K_dict = get_plink_top_hits(plink_pihat_score_dict, int(args.knn))
    plink_kin_K_dict = get_plink_top_hits(plink_kin_score_dict, int(args.knn))

    write_top_k(plink_DST_K_dict, args.out + "/plink_DST_top_"+str(args.knn)+".txt")
    write_top_k(plink_pihat_K_dict, args.out + "/plink_pihat_top_"+str(args.knn)+".txt")
    write_top_k(plink_kin_K_dict, args.out + "/plink_kin_top_"+str(args.knn)+".txt")

if __name__ == "__main__":
    main()
