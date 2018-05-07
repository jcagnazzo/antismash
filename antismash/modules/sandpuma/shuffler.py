import os, sys
from random import shuffle

cwd = os.getcwd()
queries = []

def shuffler(query_file, splits, shuffles):
    groups = {}
    with open(query_file, 'r') as qfile:
        next(qfile)
        for line in qfile:
            line = line.rstrip()
            pair = list(line.split(','))
            queries.append(pair)
    for i in range(shuffles):
        shuffle(queries)
        groupnum = 1
        groupsize = int(len(queries)/splits)
        for j in range(0,len(queries),groupsize):
            groups[groupnum] = queries[j:j+groupsize]
            groupnum += 1
    return groups

shuffled = shuffler(cwd + "\\" 'fullset_with_specs.csv',10,1)
print(shuffled) 


