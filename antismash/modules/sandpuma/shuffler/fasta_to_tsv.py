import os, sys

def fasta_to_tsv(fasta):
    fasta_dict = {}
    pair = []
    with open(fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                pair.append(line.strip('>'))
            else:
                pair.append(line)
                fasta_dict[pair[0]] = pair[1]
                pair = []
    return fasta_dict

fd = fasta_to_tsv(sys.argv[1])
print (fd)

