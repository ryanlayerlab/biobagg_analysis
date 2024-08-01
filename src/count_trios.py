trio_file='/Users/krsc0813/PycharmProjects/biobagg_analysis/data/1kg_info/1kg_trios.txt'

trio_file_open = open(trio_file, 'r')
header = trio_file_open.readline()
num_trios = 0
for line in trio_file_open:

    L = line.strip().split()
    if L[1] == '0' and L[2] == '0':
        continue
    if L[1] == '0' or L[2] == '0':
        print(line)
    num_trios += 1

trio_file_open.close()

print(num_trios)