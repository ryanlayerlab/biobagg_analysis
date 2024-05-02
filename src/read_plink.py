from collections import defaultdict

## All of these functions help read plink files

def get_pairwise_DST_score_dict(plink_file):
    """
    Get the pairwise DST score from a plink file
    """
    pairsie_DST_score_dict = defaultdict(dict)

    with open(plink_file, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip().split()
            sample_A = line[0]
            sample_B = line[2]
            DST_score = float(line[11])
            try:
                pairsie_DST_score_dict[sample_A][sample_B] = DST_score
            except KeyError:
                pairsie_DST_score_dict[sample_A].update({sample_B: DST_score})
    return pairsie_DST_score_dict

def get_pairwise_pihat_score_dict(plink_file):
    """
    Get the pairwise pihat score from a plink file
    """
    pairsie_pihat_score_dict = defaultdict(dict)

    with open(plink_file, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip().split()
            sample_A = line[0]
            sample_B = line[2]
            pihat_score = float(line[9])
            try:
                pairsie_pihat_score_dict[sample_A][sample_B] = pihat_score
            except KeyError:
                pairsie_pihat_score_dict[sample_A].update({sample_B: pihat_score})
    return pairsie_pihat_score_dict

def get_pairwise_kin_score_dict(plink_file):
    """
    Get the pairwise kin score from a plink file
    """
    pairwise_kin_score_dict = defaultdict(dict)

    with open(plink_file, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip().split()
            sample_A = line[0]
            sample_B = line[1]
            kin_score = float(line[5])
            try:
                pairwise_kin_score_dict[sample_A][sample_B] = kin_score
            except KeyError:
                pairwise_kin_score_dict[sample_A].update({sample_B: kin_score})
    return pairwise_kin_score_dict