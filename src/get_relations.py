import argparse
from plotting import ancestry_helpers
from collections import deque

SUB_SUPERPOPULATIONS = ancestry_helpers.SUB_SUPERPOPULATIONS

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ped', type=str, help='PED file', required=True)
    parser.add_argument('--ancestry', type=str, help='ancestry file', required=True)
    parser.add_argument('--no_header',
                        default=False,
                        action='store_true',
                        help='PED file does not have a header')
    parser.add_argument('--pop', type=str, help='population query', required=True)
    return parser.parse_args()

def get_ped(ped_file, pop, subpopulations, no_header):
    ped = []
    with open(ped_file, 'r') as f:
        if not no_header:
            f.readline()
        for line in f:
            line = line.strip().split(' ')
            if pop == 'ALL':
                ped.append(line)
            # if sample in subpopulation is in the superpopulation query
            elif SUB_SUPERPOPULATIONS[subpopulations[line[0]]] == pop:
                ped.append(line)
    return ped

class FamilyNode:
    # has parents and children
    def __init__(self, name):
        self.name = name
        self.parents = []
        self.children = []

def build_family_tree(ped):
    # build a graph of the family tree
    family_tree = {}
    for person in ped:
        sample_id = person[0]
        father_id = person[1]
        mother_id = person[2]

        if sample_id not in family_tree:
            family_tree[sample_id] = FamilyNode(sample_id)
        if father_id != '0':
            if father_id not in family_tree:
                family_tree[father_id] = FamilyNode(father_id)
            family_tree[sample_id].parents.append(father_id)
            family_tree[father_id].children.append(sample_id)
        if mother_id != '0':
            if mother_id not in family_tree:
                family_tree[mother_id] = FamilyNode(mother_id)
            family_tree[sample_id].parents.append(mother_id)
            family_tree[mother_id].children.append(sample_id)
    return family_tree

def min_path_length(family_tree, sample1, sample2):
    # 0 = self, 1 = parent/child, 2 = grandparent/grandchild/sibling, 3 = great-grandparent/great-grandchild/aunt/uncle/niece/nephew
    # -1 = no traversal

    # traverse the graph to find the shortest path between two nodes
    if sample1 not in family_tree or sample2 not in family_tree:
        return -1
    visited = set()
    queue = deque([(sample1, 0)])
    while queue:
        current, dist = queue.popleft()
        visited.add(current)
        if current == sample2:
            return dist
        for parent in family_tree[current].parents:
            if parent not in visited:
                queue.append((parent, dist + 1))
        for child in family_tree[current].children:
            if child not in visited:
                queue.append((child, dist + 1))
    return -1


def label_relationship(dist, sample1, sample2, family_tree, subpopulations):
    label = 'outpop'
    if dist == 0:
        label = 'self'
    elif dist == 1:
        if sample1 in family_tree[sample2].children or sample2 in family_tree[sample1].parents:
            label = 'child'
        else:
            label = 'parent'
    elif dist == 2:
        # siblings
        s1_parents = family_tree[sample1].parents
        s2_parents = family_tree[sample2].parents
        for p in s1_parents:
            if p in s2_parents:
                label = 'sibling'
        # grandparent/grandchild
        s1_children = family_tree[sample1].children
        s2_children = family_tree[sample2].children
        for c in s1_children:
            s1_grandchildren = family_tree[c].children
            if sample2 in s1_grandchildren:
                label = 'grandchild'
        for c in s2_children:
            s2_grandchildren = family_tree[c].children
            if sample1 in s2_grandchildren:
                label = 'grandparent'
    elif dist > 0:
        # connected somehow, but could be related (e.g. aunt/uncle) or married in (unrelated)
        label = 'unknown'
    else:
        # if not connected, check population
        if dist < 0:
            subpop1 = subpopulations[sample1]
            subpop2 = subpopulations[sample2]
            pop1 = ancestry_helpers.SUB_SUPERPOPULATIONS[subpop1]
            pop2 = ancestry_helpers.SUB_SUPERPOPULATIONS[subpop2]
            if pop1 == pop2:
                if subpop1 == subpop2:
                    label = 'subpop'
                else:
                    label = 'superpop'

    return label

def get_samples(ped, pop, subpopulations):
    samples = []
    f = open(ped, 'r')
    f.readline()
    for line in f:
        line = line.strip().split(' ')
        if pop == 'ALL':
            samples.append(line[0])
        # if sample in subpopulation is in the superpopulation query
        elif SUB_SUPERPOPULATIONS[subpopulations[line[0]]] == pop:
            samples.append(line[0])
    f.close()
    return samples


def main():

    args = get_args()
    subpopulations = ancestry_helpers.get_subpopulations(args.ancestry)
    ped = get_ped(args.ped, args.pop, subpopulations, args.no_header)
    graph = build_family_tree(ped)
    samples = get_samples(args.ped, args.pop, subpopulations)


    o_file = open('1KG_trios_' + args.pop + '.txt', 'w')

    for i in samples:
        for j in samples:
            dist = min_path_length(graph, i, j)
            label = label_relationship(dist, i, j, graph, subpopulations)

            o_file.write(f'{i} {j} {dist} {label}\n')

    o_file.close()


if __name__ == '__main__':
    main()