def split_into_groups(fasta, n_groups):
    """ divides a query set into groups to run over SANDPUMA's parallelized pipleline
    Arguments:
        fasta: dictionary of seq names (str) to seqs (str)
        n_groups: number of groups to split into (you can think of this as the number of threads)

    Returns:
        dictionary of groups to list of fasta headers
    """
    n_seqs = len(fasta)
    seqs_per_group = int(n_seqs / n_groups)
    qnum = 0
    groupnum = 1
    groups = {}
    for qname in fasta:
        if (qnum == 0) or (qnum < seqs_per_group):
            groupname = 'group'+str(groupnum)
            if groupname in groups:
                groups[groupname].append(qname)
            else:
                groups[groupname] = [qname]
            qnum += 1
        else:
            groupnum += 1
            groupname = 'group'+str(groupnum)
            groups[groupname] = [qname]
            qnum = 0
    return groups

fasta = {'a': 'AACC', 'b': 'BBCB', 'c':'LLGH'}
split_into_groups(fasta, 3)