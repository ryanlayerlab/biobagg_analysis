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

def get_subpop_full_names():
    """
    Get full subpopulation names
    :return: mapping of subpopulation codes to full names
    """
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

SUPER_SUBPOPULATIONS = {'AFR': ['ASW', 'LWK', 'GWD', 'MSL', 'ESN', 'YRI', 'ACB'],
                        'AMR': ['CLM', 'PEL', 'MXL', 'PUR'],
                        'EAS': ['CDX', 'CHB', 'JPT', 'KHV', 'CHS'],
                        'EUR': ['CEU', 'TSI', 'FIN', 'GBR', 'IBS'],
                        'SAS': ['BEB', 'GIH', 'ITU', 'PJL', 'STU']}

SUBPOPULATIONS = ['ASW', 'LWK', 'GWD', 'MSL', 'ESN', 'YRI', 'ACB', 'CLM', 'PEL', 'MXL', 'PUR', 'CDX', 'CHB', 'JPT', 'KHV', 'CHS', 'CEU', 'TSI', 'FIN', 'GBR', 'IBS', 'BEB', 'GIH', 'ITU', 'PJL', 'STU']

SUB_SUPERPOPULATIONS = {'ASW': 'AFR', 'LWK': 'AFR', 'GWD': 'AFR', 'MSL': 'AFR', 'ESN': 'AFR', 'YRI': 'AFR', 'ACB': 'AFR', 'CLM': 'AMR', 'PEL': 'AMR', 'MXL': 'AMR', 'PUR': 'AMR', 'CDX': 'EAS', 'CHB': 'EAS', 'JPT': 'EAS', 'KHV': 'EAS', 'CHS': 'EAS', 'CEU': 'EUR', 'TSI': 'EUR', 'FIN': 'EUR', 'GBR': 'EUR', 'IBS': 'EUR', 'BEB': 'SAS', 'GIH': 'SAS', 'ITU': 'SAS', 'PJL': 'SAS', 'STU': 'SAS'}

def get_subpopulations(ancestry_file):
    sample_subpopulations = {}
    with open(ancestry_file) as f:
        for line in f:
            for subpop in SUBPOPULATIONS:
                if subpop in line:
                    sample_subpopulations[line.split()[0]] = subpop
    return sample_subpopulations

def get_colors(color_file):
    '''
    Get colors from file
    @param color_file: path to color file
    @return: dictionary of colors
    '''
    colors = {}
    f = open(color_file, 'r')
    for line in f:
        line = line.strip().split(',')
        try:
            colors[int(line[0])] = line[1]
        except ValueError:
            colors[line[0]] = line[1]
    f.close()
    return colors