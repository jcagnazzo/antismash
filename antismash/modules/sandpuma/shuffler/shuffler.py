import os, sys
from random import shuffle

cwd = os.getcwd()

def fasta_to_list(fasta):
    ''' Takes a fasta file and returns a list of lists. Each interior list contains a
    header and its sequence'''

    query_list = []
    pair = []
    with open(fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                pair.append(line.strip('>'))
            else:
                pair.append(line)
                query_list.append(pair)
                pair = []
    return query_list

def shuffler(fasta, splits, jk_file):
    
    qlist = fasta_to_list(fasta)
    groups = {}
    shuffle(qlist)
    groupnum = 1
    groupsize = int(len(qlist)/splits)
    # takes the shuffled qlist and makes the
    # dictionary groups who's keys are the groupnames
    # matched to a dictionary containing the query
    # names matched to the sequence
    for i in range(0,len(qlist),groupsize):
        seq_dict = {}
        for j in range(groupsize):
            if i + j > len(qlist) - 1:
                continue
            else:
                seq_dict[qlist[i + j][0]] = qlist[i + j][1]
                j += 1
        groups['group' + str(groupnum)] = seq_dict
        groupnum += 1
    for gname, indict in groups.items():
        with open(gname + '.faa', 'w') as group_fasta:
            for qname, qseq in indict.items():
                group_fasta.write(qname + '\n')
                group_fasta.write(qseq + '\n')
            with open(jk_file, 'r') as jk:
                with open(gname + 'jackknife.tsv', 'w') as jk_new:
                    for line in jk:
                        if qname in jk:
                            continue
                        else:
                            jk_new.write(line)











shuffler(cwd + '\\shuffler\\' + 'miniset.faa', 10, cwd + '\\data\\sandpuma1_jackknife.tsv' )

