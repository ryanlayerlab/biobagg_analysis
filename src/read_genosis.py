from collections import defaultdict

# TOP_HITS.txt file format:
# query match_1,match_1_score match_2,match_2_score ...

def get_top_hits_dict(top_hits_file):
    """
    Get the top hits from a top hits file
    """
    top_hits_dict = {}
    # query: [(match_1, match_1_score), (match_2, match_2_score), ...]

    with open(top_hits_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            query = line[0]
            top_hits_list = line[1:]
            for hit in top_hits_list:
                match, score = hit.split(',')
                score = float(score)
                try:
                    top_hits_dict[query].append((match, score))
                except KeyError:
                    top_hits_dict[query] = [(match, score)]

    return top_hits_dict