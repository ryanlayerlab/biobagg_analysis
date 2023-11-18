from collections import defaultdict

def get_population_maps(population_file):
    """
    Make four mappings of population codes

    @param population_file: path to a file containing 1KG population labels
    @return: sample_subpopulation: a dictionary mapping sample IDs to subpopulation codes
                sub_to_super: a dictionary mapping subpopulation codes to superpopulation codes
                super_to_sub: a dictionary mapping superpopulation codes to subpopulation codes
    """
    sample_subpopulations = defaultdict()
    sub_to_super = defaultdict()
    super_to_sub = defaultdict()

    with open(population_file, 'r') as f:
        header = f.readline()
        for line in f:
            L = line.strip().split()
            # get necessary data from file
            sampleID = L[0]
            subpop_code = L[3]
            if ',' in subpop_code:
                subpop_code = subpop_code.split(',')[0]
            superpop_code = find_superpop_code(L)

            # make mappings
            # sample ID: subpopulation code
            sample_subpopulations[sampleID] = subpop_code
            # subpopulation code: superpopulation code
            sub_to_super[subpop_code] = superpop_code
            # superpopulation code: list of subpopulation codes
            try:
                super_to_sub[superpop_code].add(subpop_code)
            except KeyError:
                super_to_sub[superpop_code] = {subpop_code}

    f.close()
    return sample_subpopulations, sub_to_super, super_to_sub

def get_knn_results(knn_results):
    """
    Get KNN results

    @param knn_results: path to KNN results file
    @return: a dictionary mapping sample IDs to KNN results
    """
    sample_knn = defaultdict(list)
    with open(knn_results, 'r') as f:
        for line in f:
            line = line.strip().split()
            queryID = line[0]
            knn = line[1:]
            for k in knn:
                matchID = k.split(',')[0]
                matchScore = float(k.split(',')[1])
                if matchID == queryID: # ignore self
                    continue
                try:
                    sample_knn[queryID].append((matchID, matchScore))
                except KeyError:
                    sample_knn[queryID] = [(matchID, matchScore)]
    f.close()
    return sample_knn

def find_superpop_code(sample_list):
    """
    returns superpopulation code from 1KG file

    :param sample list: line for sample data from 1KG file
    :return: superpopulation code
    """

    # superpop code is one of these: AFR, AMR, EAS, EUR, SAS
    superpop_options = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    for superpop in superpop_options:
        if superpop in sample_list:
            return superpop
