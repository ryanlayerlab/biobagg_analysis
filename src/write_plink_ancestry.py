import plotting.ancestry_helpers as ah


def read_hits_file(hits_file, trio_samples, sample_subpopulations, SUB_SUPERPOPULATIONS):
    # { superpop: {subpop: [scores], superpop: [scores], outgroup: [scores] } }
    hits = {'AFR': {'subpop': [], 'superpop': [], 'outgroup': []},
            'AMR': {'subpop': [], 'superpop': [], 'outgroup': []},
            'EAS': {'subpop': [], 'superpop': [], 'outgroup': []},
            'EUR': {'subpop': [], 'superpop': [], 'outgroup': []},
            'SAS': {'subpop': [], 'superpop': [], 'outgroup': []},}

    # query match,score match,score match,score...
    f = open(hits_file, 'r')
    header = f.readline()
    for line in f:
        line = line.strip().split()
        query = line[0]
        # ignore trios
        if query in trio_samples:
            continue
        else:
            # get subpopulation
            sample_subpop = sample_subpopulations[query]
            # get superpopulation
            sample_superpop = SUB_SUPERPOPULATIONS[sample_subpop]
            # add scores to dictionary
            for match in line[1:]:
                match_ID, match_score = match.split(',')
                match_score = float(match_score)
                # # if match is a trio, ignore
                if match_ID in trio_samples:
                    continue
                # if match is self, ignore
                if match_ID == query:
                    continue
                match_subpop = sample_subpopulations[match_ID]
                match_superpop = SUB_SUPERPOPULATIONS[match_subpop]
                # add score to hits dictionary
                # if subpopulation match:
                if match_subpop == sample_subpop:
                    # if match_score > 1000:
                    #     print(sample_superpop, query, match_ID, match_score)
                    hits[sample_superpop]['subpop'].append(float(match_score))
                # if superpopulation match:
                elif match_superpop == sample_superpop:
                    hits[sample_superpop]['superpop'].append(float(match_score))
                # if outgroup match:
                else:
                    hits[sample_superpop]['outgroup'].append(float(match_score))

    f.close()
    return hits

def get_trio_samples(ped_file):
    trio_samples = []
    f = open(ped_file, 'r')
    header = f.readline()
    for line in f:
        line = line.strip().split()
        sampleID = line[0]
        fatherID = line[1]
        motherID = line[2]
        # if sample has a father OR mother; add sample, father, and mother to trio_samples
        if fatherID != '0' or motherID != '0':
            trio_samples.append(sampleID)
            # trio_samples.append(fatherID)
            # trio_samples.append(motherID)

        # remove duplicates
        trio_samples = list(set(trio_samples))
    f.close()
    return trio_samples

def write_plink_hits(hits, out_file):
    f = open(out_file, 'w')
    f.write('superpop,category,scores...\n')
    for superpop in hits.keys():
        for category in hits[superpop].keys():
            f.write(superpop + ',' + category + ',')
            for score in hits[superpop][category]:
                f.write(str(score) + ' ')
            f.write('\n')
    f.close()

def main():
    ped_file = 'data/1kg_trios.txt'
    trio_samples = get_trio_samples(ped_file)

    ancestry_file = 'data/igsr-1000 genomes 30x on grch38.tsv'
    sample_subpop = ah.get_subpopulations(ancestry_file)
    SUB_SUPERPOPULATIONS = ah.SUB_SUPERPOPULATIONS

    plink_dst_hits_file = 'plink_top_K_data/plink_DST_top_20.txt'
    plink_pihat_hits_file = 'plink_top_K_data/plink_pihat_top_20.txt'
    plink_kinship_hits_file = 'plink_top_K_data/plink_kin_top_20.txt'
    dst_hits = read_hits_file(plink_dst_hits_file,
                              trio_samples,
                              sample_subpop,
                              SUB_SUPERPOPULATIONS)
    pihat_hits = read_hits_file(plink_pihat_hits_file,
                                trio_samples,
                                sample_subpop,
                                SUB_SUPERPOPULATIONS)
    kinship_hits = read_hits_file(plink_kinship_hits_file,
                                  trio_samples,
                                  sample_subpop,
                                  SUB_SUPERPOPULATIONS)
    write_plink_hits(dst_hits, 'plink_top_K_data/plink_DST_20_groups.txt')
    write_plink_hits(pihat_hits, 'plink_top_K_data/plink_pihat_20_groups.txt')
    write_plink_hits(kinship_hits, 'plink_top_K_data/plink_kin_20_groups.txt')

if __name__ == '__main__':
    main()
