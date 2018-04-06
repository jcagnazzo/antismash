# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" SANDPUMA. """

import logging
import os
import re
from typing import Any, Dict, List, Optional, Tuple

from antismash.common import fasta, module_results, pfamdb, subprocessing, pplacer
from antismash.common.secmet import Record, CDSFeature
from antismash.common.secmet.feature import FeatureLocation, PFAMDomain

from sklearn import tree
import numpy as np
from Bio import Phylo
import multiprocessing


class PredicatResults(module_results.ModuleResults):
    """ Results for prediCAT """
    def __init__(self, record_id: str, monophyly: str, forced: str, nn_dist: float, nn_score: float, snn_score: float) -> None:
        super().__init__(record_id)
        self.monophyly = str(monophyly)
        self.forced = str(forced)
        self.nn_dist = float(nn_dist)
        self.nn_score = float(nn_score)
        self.snn_score = float(snn_score)


class SandpumaResults(module_results.ModuleResults):
    """ Results for SANDPUMA """
    def __init__(self, record_id: str, predicat: str, asm: str, svm: str, phmm: str, pid: float, ensemble: str, sandpuma: str) -> None:
        super().__init__(record_id)
        self.predicat = PredicatResults
        self.asm = str(asm)
        self.svm = str(svm)
        self.phmm = str(phmm)
        self.pid = float(pid)
        self.ensemble = str(ensemble)
        self.sandpuma = str(sandpuma)


def cleancall(call: str)->str:
    ''' Cleans up predictions '''
    call = call.replace("_result", "")
    call = call.replace("no_confident", "nocall")
    call = call.replace("N/A", "nocall")
    call = call.replace("no_call", "nocall")
    return call


def is_trustworthy_path(spresult: SandpumaResults, paths: List[str], pathacc: Dict[str, str, Any], cutoff: float) -> bool:
    ''' Decide if a decision tree path is trustworthy '''
    passed = 1
    path = ''
    for pa in paths:
        decisions = pa.split('&')
        for d in decisions:
            if d == 'LEAF_NODE':
                break
            else:
                decision, threshchoice = d.split('%')
                thresh, choice = threshchoice.split('-')
                if decision == 'pid':
                    if choice == 'T': ## Need greater than the thresh to pass
                        if thresh <= spresult.pid:
                            passed = 0
                            break
                        else: ## Need less or equal to thresh to pass
                            if thresh > spresult.pid:
                                passed = 0
                                break
                else: ## Not pid
                    decision = cleancall(decision)
                    a = decision.split('_')
                    spec = a[-1]
                    method = decision.replace('_'+spec, "")
                    tocheck = ''
                    if method == 'SVM':
                        tocheck = spresult.svm
                    elif method == 'prediCAT_SNN':
                        tocheck = spresult.predicat.forced
                    elif method == 'prediCAT_MP':
                        tocheck = spresult.predicat.monophyly
                    elif method == 'pHMM':
                        tocheck = spresult.phmm
                    elif method == 'ASM':
                        tocheck = spresult.asm
                    if choice == 'T': ## than 0.5, so NOT spec
                        if tocheck == spec:
                            passed = 0
                            break
                        else: ## matches spec
                            if tocheck != spec:
                                passed = 0
                                break
        path = pa
    path = re.sub(r"\S+&(LEAF_NODE-\d+)$", "\g<1>", path)
    if pathacc[path]['pct'] < cutoff:
        return False
    else:
        return True


def get_parent(tree, child_clade) -> Phylo.Node:
    ''' Given a tree and a node, returns the parental node '''
    node_path = tree.get_path(child_clade)
    if len(node_path) < 2:
        return None
    return node_path[-2]


def get_child_leaves(tree, parent) -> List:
    ''' Given a tree and a node, returns all children nodes '''
    child_leaves = []
    for leaf in tree.get_terminals():
        for node in tree.get_path(leaf):
            if(node == parent):
                child_leaves.append(leaf)
    return child_leaves


def calcscore(scaleto, distance) -> float:
    ''' Scale distance to score from 0 to 1 '''
    if(distance >= scaleto):
        return 0
    else:
        return float(scaleto - distance) / scaleto


def getscore(scaleto, nd, dist2q, leaf, o) -> float:
    ''' Calculate the SNN '''
    score = 0
    nnspec = leaf[o[0]]['spec']
    for n in o:
        curspec = leaf[o[0]]['spec']
        if(nnspec == curspec):
            tmpscore = calcscore(scaleto, float(dist2q[n]) / nd)
            if(tmpscore > 0):
                score += tmpscore
            else:
                break
        else:
            break
    return score    


def deeperdive(query: int, tree: Phylo.Tree, nearest1: int, nearest2: int, l: Dict[int, str, Any])-> str, str, str:
    """ deeper substrate prediction triggered for non-monophyletic seqs
    Arguments:
        query: index for the query
        tree: tree
        nearest1: index for the nearest neighbor
        nearest2, index for the second nearest neighbor
        l: dictionary of leaf index to str to any
            includes group (str), id (str), spec (str), node (Phylo.node)

    Returns:
        monophyly specificity (str), hit name (str), forced call specificity (str)
    """
    ## Want query to nearest dist to be less than nearest1 to nearest2 dist
    query_to_nn1 = tree.distance(l[query]['node'], l[nearest1]['node'])
    nn1_to_nn2 = tree.distance(l[nearest1]['node'], l[nearest2]['node'])
    query_to_nn2 = tree.distance(l[query]['node'], l[nearest2]['node'])
    if((query_to_nn1 < nn1_to_nn2) and (l[nearest1]['spec'] == l[nearest2]['spec'])):
        return (l[nearest1]['spec'], l[nearest1]['id'], 'no_force_needed')
    elif((query_to_nn1 == query_to_nn2) and (l[nearest1]['spec'] != l[nearest2]['spec'])):
        return ('no_confident_result', 'NA', 'no_confident_result')
    else:
        parent = get_parent(tree, l[query]['node'])
        if parent is None:
            return ('no_confident_result', 'NA', 'no_confident_result')
        sisterdist = {}
        for sister in get_child_leaves(tree, parent):
            sisterdist[sister.name] = {}
            sisterdist[sister.name]['dist'] = tree.distance(l[query]['node'], sister)
            sisterdist[sister.name]['node'] = sister
            ordered_sisterdist = sorted(sisterdist, key=sisterdist.get)
            for name in ordered_sisterdist:
                if(name != l[query]['id']):
                    forced = re.split("_+", name)
                    return ('no_confident_result', 'NA', forced[-1])
            return ('no_confident_result', 'NA', 'no_confident_result') 


def checkclade(query: int, lo: int, hi: int, wc: str, tree: Phylo.Tree, l: Dict[int, str, Any])-> str, str:
    """ recursive substrate prediction for a query & it's sisters in a tree
    Arguments:
        query: index for the query
        lo: index for the lower sister
        hi: index for the higher sister
        wc: wildcard variable
        tree: tree
        l: dictionary of leaf index to str to any
            includes group (str), id (str), spec (str), node (Phylo.node)

    Returns:
        substrate specificity (str), hit name (str)
    """
    if((lo in l) and (hi in l)): ## Not first or last
        if(l[lo]['spec'] == wc): ## lower bound is wildcard
            checkclade(query, lo-1, hi, wc, l)
        elif(l[hi]['spec'] == wc): ## upper bound is wildcard
            checkclade(query, lo, hi+1, wc, l)
        else:
            ## Get the lca's descendants and check specs
            lca = tree.common_ancestor(l[lo]['node'], l[hi]['node'])
            spec = ''
            iname = ''
            passflag = 1
            for child in get_child_leaves(tree, lca):
                split_id = re.split("_+", child.name)
                if(spec != ''): ## Not yet assigned
                    if((split_id[-1] != spec) and (split_id[-1] != wc)): ## Specs are different, Requires a deeper dive
                        passflag = 0
                    else:
                        spec = split_id[-1]
                        iname = split_id[0]
                else:
                    spec = split_id[-1]
                    iname = split_id[0]
            if(passflag == 0):
                return('deeperdive', 'NA')
            else:
                return(spec, iname)
    else: ## First or last
        return('deeperdive', 'NA')        


def predicat(pplacer_tree: str, masscutoff: float, wild: str, snn_thresh: float)-> PredicatResults:
    """ predicat substrate prediction
    Arguments:
        pplacer_tree: pplacer tree (newick)
        masscutoff: cutoff value for pplacer masses
        wild: wildcard variable
        snn_thresh: SNN threshold for confident prediction (default=0.5)

    Returns:
        PredicatResults
            monophyly -> substrate specificity (str)
            forced -> substrate specificity (str)
            nndist -> distance to nearest neighbor (float)
            nn_score -> nearest neighbor score (float)
            snn_score -> scaled nearest neighbor score (float)
    """
    ## predicat settings
    zero_thresh = 0.005 ## Branch distance less than which to call sequences identical
    nppref = ['Q70AZ7_A3', 'Q8KLL5_A3'] ## Leaves used to normalize the nn score to SNN
    npnode = ['', ''] ## initialize node list for leaves above
    dcut = 2.5 ## normalized distance cutoff for nearest neighbor (emperically derived default=2.5)
    ## read tree
    tree = Phylo.read(pplacer_tree, 'newick')
    query = []
    leaves = {}
    ## Loop through leaves, only keep top placement
    for leaf in tree.get_terminals():
        split_id = re.split("_+", leaf.name) ## Split ID on _, split_id[-1] will be specificity
        if re.match(r"^#[123456789]\d*$", split_id[-2]) is not None: ## matches a non-top pplacer placement
            tree.prune(leaf) ## remove it
    ## Loop through leaves to find groups
    last = '' ## Keep track of the last specificity, initialize on ''
    group = 1 ## group number
    leafnum = 1 ## leaf number
    for leaf in tree.get_terminals():
        if(bool(re.search('^'+nppref[0], leaf.name))): ## if node is nppref[0], store the node
            npnode[0] = leaf
        elif(bool(re.search('^'+nppref[1], leaf.name))): ## if node is nppref[1], store the node
            npnode[1] = leaf
        split_id = re.split("_+", leaf.name) ## Split ID on _, split_id[-1] will be specificity
        if(last != ''): ## Every pass except the first
            if((last != split_id[-1]) or (last != wild)): ## begin new group
                group += 1
        if re.match("^#0$", split_id[-2]) is not None: ## matches pplacer query formatting; #0 is the top placement
            last = wild
            mass = float(re.sub(r"^.+#\d+_M=(\d+?\.?\d*)$", "\g<1>", leaf.name))
            if(mass < masscutoff):
                return PredicatResults('no_confident_result', 'no_confident_result', 0, 0, 0)
        else:
            last = split_id[-1]
        leaves[leafnum] = {}
        leaves[leafnum]['id'] = leaf.name
        leaves[leafnum]['group'] = group
        leaves[leafnum]['spec'] = last
        leaves[leafnum]['node'] = leaf
        ## Record queries
        if(last == wild):
            query.append(leafnum)
        leafnum += 1 
    qnum = next(iter(query))
    ## Get distances to knowns
    distfromquery = {}
    for leafnum in leaves:
        if((qnum != leafnum) and (leaves[leafnum]['spec'] != wild)):
            distfromquery[leafnum] = tree.distance(leaves[qnum]['node'], leaves[leafnum]['node'])
    # Sort distances
    ordered_dist = sorted(distfromquery, key=distfromquery.get)
    ## Get zero distances
    zero_dist = []
    for leafnum in ordered_dist:
        if(distfromquery[leafnum] <= zero_thresh):
            if(distfromquery[leafnum] >= 0):
                zero_dist.append(leafnum)
            else:
                break
    forcedpred = 'no_force_needed'
    pred = 'no_call'
    hit = 'NA'
    if(len(zero_dist) > 0): ## query has zero length known neighbors
        pred = leaves[zero_dist[0]]['spec']
        hit = re.search("^(\S+)_.+$", leaves[zero_dist[0]]['id']).groups()[0]
    else:
        ## predict the clade
        pred, hit = checkclade(qnum, qnum-1, qnum+1, wild, tree, leaves)
        if(pred == 'deeperdive'):
            pred, hit, forcedpred = deeperdive(qnum, tree, ordered_dist[0], ordered_dist[1], leaves)
            if(hit != 'NA'):
                hit = re.search("^(\S+)_.+$", hit).groups()[0]
    if forcedpred == 'no_force_needed':
        forcedpred = pred
    normdist = tree.distance(npnode[0], npnode[1])
    nn_dist = float(distfromq[ordered_dist[0]]) / normdist
    nnscore = 0
    snnscore = 0
    if(nn_dist < dcut):
        snnscore = getscore(dcut, normdist, distfromquery, leaves, ordered_dist)
        nnscore = calcscore(dcut, nn_dist)
    if(snnscore < snn_thresh):
        forcedpred = 'no_confident_result'
    return PredicatResults(pred, forcedpred, nn_dist, nnscore, snnscore)


def run_predicat(reference_aln: str, queryfa: Dict[str, str], wildcard: str, ref_tree: str, ref_pkg: str, masscutoff: float) -> PredicatResults:
    """ pplacer and predicat substrate prediciton
    Arguments:
        reference_aln: filename for reference protein fasta, see sandpuma_multithreaded comments for requirements
        queryfa: seq id to seq dictionary
        wildcard: suffix str identifying query sequence (Default= 'UNK' which means headers end in '_UNK')
        ref_tree: reference tree (newick)
        ref_pkg: pplacer reference package
        masscutoff: cutoff value for pplacer masses

    Returns:                                                                                                                            PredicatResults
            monophyly -> substrate specificity (str)
            forced -> substrate specificity (str)
            nndist -> distance to nearest neighbor (float)
            nn_score -> nearest neighbor score (float)
            snn_score -> scaled nearest neighbor score (float)
    """
    query = next(iter(queryfa))
    ## Align query to a single known sequence
    to_align = queryfa
    ref = read_fasta(reference_aln)
    tname = next(iter(known)) ## Grab any training sequence header
    to_align[tname] = known[tname].replace('-', '')
    aligned = subprocessing.run_mafft_predicat_trim(to_align)
    ## trim overhangs
    head = len(re.sub(r'^(-*).+$', r'\g<1>', aligned[tname]))
    tail = len(re.sub(r'^.+(-*)$', r'\g<1>', aligned[tname]))
    trimmed = aligned[query][head:-tail].replace('-', '') ## Removes head and tail then removes gaps
    trimmedfa = {query: trimmed}
    ## Align trimmed seq to reference
    all_aligned = subprocessing.run_muscle_profile_sandpuma(reference_aln, trimmedfa)
    ## Pplacer (NOTE: this is new to SANDPUMA as of antiSMASH5 and needs to be tested
    pplacer_tree = subprocessing.run_pplacer(ref_tree, reference_aln, ref_pkg, all_aligned)
    ## prediCAT
    return predicat(pplacer_tree, masscutoff, wildcard)

def run_asm(queryfa: Dict[str, str], stachfa: Dict[str, str], seedfa: Dict[str, str]) -> str:
    """ Active site motif (ASM) substrate prediction
    Arguments:
        queryfa: seq id to seq dictionary
        stachfa: seq name to seq for stachelhaus codes
        seedfa: seq name to seq for seed alignment for stachelhaus code extraction

    Returns:                                                                                                                            substrate specificity prediction
    """ 
    ## ASM settings
    gapopen = 3.4
    properlen = 117 ## true length
    grsAcode = {4:1,5:1,8:1,47:1,68:1,70:1,91:1,99:1,100:1} ## positions in grsA for code
    ## Alignment
    toalign = {**queryfa, **seedfa}
    aligned2seed = subprocessing.mafft_sandpuma_asm(toalign, gapopen)
    ## Loop through alignment to find new positions for code
    qname = next(iter(queryfa))
    pos = 0
    newcode = []
    for p, val in enumerate(aligned2seed['phe_grsA']):
        if(val=='-'):
            continue
        else:
            pos += 1
            if(pos in grsAcode):
                newcode.append(p)
    ## Extract codes
    extractedcode = {}
    for seqname in aligned2seed:
        code = ''
        for n in newcode:
            code = code + aligned2seq[seqname][n]
            extractedcode[seqname] = code
    ## Error checking
    truth = {'phe_grsA':'DAWTIAAIC', 'asp_stfA-B2':'DLTKVGHIG','orn_grsB3':'DVGEIGSID','val_cssA9':'DAWMFAAVL'}
    for seqname in extractedcode:
        if seqname == qname:
            continue
        else:
            if extractedcode[seqname] != truth[seqname]:
                return('no_call') ## Issue with the alignment
    ## Score each
    scores = {}
    for sname in stachfa:
        match = 0
        split_id = re.split("_+", sname)
        spec = re.split("|", split_id[-1])
        for p, val in enumerate(stachfa[sname]):
            if val == extractedcode[qname][p]:
                match += 1
        if str(match) in scores:
            for s in spec:
                if s in scores[str(match)]:
                    scores[str(match)][s] += 1
                else:
                    scores[str(match)][s] = 1
        else:
            scores[str(match)] = {}
            for s in spec:
                scores[str(match)][s] = 1
    if '9' in scores:
        return('|'.join(sorted(scores['9'])))
    elif '8' in scores:
        return('|'.join(sorted(scores['8'])))
    elif '7' in scores:
        return('|'.join(sorted(scores['7'])))
    else:
        return('no_call')


def run_svm(queryfa: Dict[str, str]) -> str:
    """ Support vector machine (SVM) substrate prediction
    Arguments:
        queryfa: seq id to seq dictionary

    Returns:                                                                                                                            substrate specificity prediction
    """
    ## Set input files
    ref = "A_domains_muscle.fasta" ## From NRPSPredictor2
    ## Set positions
    startpos = 66
    a34positions = [210, 213, 214, 230, 234,
                    235, 236, 237, 238, 239,
                    240, 243, 278, 279, 299,
                    300, 301, 302, 303, 320,
                    321, 322, 323, 324, 325,
                    326, 327, 328, 329, 330,
                    331, 332, 333, 334]
    positions34 = []
    for p in a34positions:
        positions34.append(p-startpos)
    aligned = subprocessing.run_muscle_profile_sandpuma(ref, queryfa)
    refname = "P0C062_A1"
    ## Get the 34 code for the query
    qname = next(iter(queryfa))
    refseq = aligned[refname]
    allp, nongaps = 0, 0
    poslist = []
    while refseq != '':
        if nongaps in positions34 and refseq[0] != '-':
            poslist.append(b)
        if refseq[0] != '-':
            nongaps += 1
        allp += 1
        refseq[1:]
    seq34 = ''
    for j in poslist:
        aa = queryfa[qname][j]
        k, l = j, j
        if aa == '-':
            k += 1
            l = l - 1
            if l not in poslist:
                aa = queryfa[qname][l]
            elif (j+1) not in poslist:
                aa = queryfa[qname][k]
        seq34 = seq34+aa
    return subprocessing.run_svm_sandpuma(seq34)


def get_feature_matrix(spec: str, i2s: List[str]) -> List:
    """ Generate a feature matrix given a specificity and an ordered index2specificity map
    Arguments:
        spec: substrate specificity
        i2s: ordered list of specificities

    Returns:                                     
        List of features (0= absent, 1=present)
    """
    f = []
    for i in i2s:
        if i == spec:
            f.append(1)
        else:
            f.append(0)
    return f


def sandpuma_multithreaded(group: str, fasta: Dict[str, str], knownfaa: str, wildcard: str, snn_thresh: float, knownasm: str, max_depth: int, min_leaf_sup: int, ref_aln: str, ref_tree: str, ref_pkg: str, masscutoff: float, stachfa: Dict[str, str], seedfa: Dict[str, str], clf: DecisionTreeClassifier, i2s: List[str], paths: List[str], pathacc: Dict[str, str, Any]) -> SandpumaResults:
    """ SANDPUMA
    Order of processing:
        predicat: both monophyly and SNN substrate specificity prediction
        asm: active-site motif substrate specificity prediction
        svm: support vector machine substrate specificity prediction
        phmm: profile hidden markov model substrate specificity prediction
        pid: calculate protein percent identity to the known/training set
        ensemble: ML decision tree substrate specificity prediction based on the results above
        rescore: exclude unreliable ensemble tree paths

    Arguments:
        group: prefix group name
        fasta: dictionary of seq names (str) to seqs (str)
        knownfaa: filename for reference protein fasta; assumes each header ends in '_' followed by the <substrate specificity>
        wildcard: str to append to the end of each query sequence; should be different that all specificities (Default= 'UNK')
        snn_thresh: threshold for SNN score (Default= 0.5) NOTE: may need to be adjusted with new pplacer implementation
        knownasm: filename for reference active site motif protein fasta, similar header formatting as knownfaa
        max_depth: maximum depth for the sklearn decision tree; default= 40
        min_leaf_sup: minimum leaf support required within the decision tree; default= 10
        jk: jackknife benchmarking results dictionary
        ref_aln: reference alignment (fasta) file
        ref_tree: reference tree (newick)
        ref_pkg: pplacer reference package
        masscutoff: cutoff value for pplacer masses
        stachfa: seq name to seq for stachelhaus codes
        seedfa: seq name to seq for seed alignment for stachelhaus code extraction
        clf: trained machine learning decision tree
        i2s: ordered list of specificities
        paths: list of possible decision tree paths
        pathacc: dictionary of path accuracies

    Returns:                                                                                                                
        dictionary of SandpumaResults
    """
    sp_results = {}
    for query in fasta:
        wc_name = query+'_'+wildcard
        ## Store as a dictionary for functions that don't
        queryfa = {wc_name: fasta[query]}
        ## PrediCAT
        predicat_result = run_predicat(ref_aln, queryfa, wildcard, ref_tree, ref_pkg, masscutoff, snn_thresh)
        ## ASM
        asm = run_asm(queryfa, stachfa, seedfa)
        ## SVM
        svm = run_svm(queryfa)
        ## pHMM
        hmmdb = 'nrpsA.hmmdb'
        phmm = subprocessing.run_phmm_sandpuma(queryfa, hmmdb)
        ## PID
        pid = subprocessing.run_pid_sandpuma(queryfa, knownfaa)
        ## Ensemble
        query_features = [pid]
        query_features.extend(get_feature_matrix(predicat_result.monophyly, i2s))
        query_features.extend(get_feature_matrix(predicat_result.forced, i2s))
        query_features.extend(get_feature_matrix(svm, i2s))
        query_features.extend(get_feature_matrix(asm, i2s))
        query_features.extend(get_feature_matrix(phmm, i2s))
        ensemble = clf.predict(query_features)[0]
        ## Rescore paths
        sp = SandpumaResults(predicat_result, asm, svm, phmm, pid, ensemble, 'Unchecked')
        if is_trustworthy_path(sp, paths, pathacc, 0.5):
            sp_results[query] = SandpumaResults(predicat_result, asm, svm, phmm, pid, ensemble, ensemble)
        else:
            sp_results[query] = SandpumaResults(predicat_result, asm, svm, phmm, pid, ensemble, 'no_call')
    return sp_results


def split_into_groups(fasta: Dict[str, str], n_groups: int) -> Dict[str, List[str]]:
    """ divides a query set into groups to run over SANDPUMA's parallelized pipleline
    Arguments:
        fasta: dictionary of seq names (str) to seqs (str)
        n_groups: number of groups to split into (you can think of this as the number of threads)

    Returns:
        dictionary of groups to list of fasta headers
    """
    n_seqs = length(fasta)
    seqs_per_group = int(n_seqs / n_groups)
    qnum = 0
    groupnum = 1
    groups = {}
    for qname in fasta:
        if (qnum == 0) or (i < seqs_per_group):
            groupname = 'group'+str(groupnum)
            if groupname in groups:
                groups[groupname].append(qname)
            else:
                groups[groupname] = [qname]
            i += 1
        else:
            groupnum += 1
            groupname = 'group'+str(groupnum)
            groups[groupname] = [qname]
    return groups


def run_sandpuma(name2seq: Dict[str, str], threads: int, knownfaa: str, wildcard: str, snn_thresh: float, knownasm: str, max_depth: int, min_leaf_sup: int, jackknife_data: str, ref_aln: str, ref_tree: str, ref_pkg: str, masscutoff:float, stach_file: str, seed_file: str):
    """ SANDPUMA parallelized pipleline
    Arguments:
        name2seq: dictionary of seq names (str) to seqs (str)
        threads: number of threads
        knownfaa: filename for reference protein fasta; assumes each header ends in '_' followed by the <substrate specificity>
        wildcard: str to append to the end of each query sequence; should be different that all specificities (Default= 'UNK')
        snn_thresh: threshold for SNN score (Default= 0.5) NOTE: may need to be adjusted with new pplacer implementation
        knownasm: filename for reference active site motif protein fasta, similar header formatting as knownfaa
        max_depth: maximum depth for the sklearn decision tree; default= 40
        min_leaf_sup: minimum leaf support required within the decision tree; default= 10
        jackknife_data: filename for jackknife benchmarking results
        ref_aln: reference alignment (fasta) file
        ref_tree: reference tree (newick)
        ref_pkg: pplacer reference package
        masscutoff: cutoff value for pplacer masses
        stach_file: fasta file of stachelhaus codes
        seed_file: seed fasta file (single entry) used for stachelhaus code extraction

    Returns:                                     

    """
    ## Load jackknife data
    jk = {}
    allspec = {}
    with open(jackknife_data, "r") as j:
        next(j) ## skip header
        for line in j:
            line = line.strip()
            l = line.split("\t")
            jk[l[10]] = {'true': l[4],
                         'pid': l[3],
                         'shuf': l[0],
                         'jk': = l[1],
                         'query': l[2],
                         'bin': l[11]}
            called_spec = l[5]
            if l[7] == 'N':
                called_spec = 'no_call'
            jk[l[10]]['method'] = {}
            jk[l[10]]['method'][l[6]] = called_spec
            allspec[l[4]] = -1
            allspec[l[5]] = -1
    ## Map specificities to integers
    i2s = []
    i = 0
    for spec in sorted(allspec, key=allspec.get):
        allspec[spec] = i
        i2s.append(spec)
        i += 1
    ## Prepare features and labels
    allmethods = ['prediCAT', 'forced_prediCAT_snn50', 'svm', 'stach', 'phmm']
    features = []
    labels = []
    for uname in jk:
        for m in allmethods:
            if m in jk[uname]['method']:
                continue
            else:
                jk[uname]['method'][m] = 'no_call'
        labels.append( allspec[jk[uname]['true']] )
        feature_matrix = [ jk[uname]['pid'] ]
        for m in allmethods:
            feature_matrix.extend(get_feature_matrix(jk[uname]['method'][m], i2s))
        features.append(feature_matrix)
    ## Train the decision tree
    clf = tree.DecisionTreeClassifier(min_samples_leaf=msl, max_depth=md)
    clf = clf.fit(features, labels)
    ## Load the nodemap for decision tree
    nodemap = {}
    with open('nodemap.tsv', "r") as nm:
        for line in nm:
            if line[0] == '#':
                continue
            else:
                line = line.strip()
                l = line.split("\t")
                nodemap[l[0]] = {'parent': l[1],
                                 'parent_call': l[2],
                                 'decision': l[3],
                                 'threshold': l[4]}
    ## Define paths
    paths = []
    for n in sorted(nodemap):
        if nodemap[n]['decision'] == 'LEAF_NODE':
            p = nodemap[n]['parent']
            traceback = nodemap[p]['decision']+'%'+node[p]['thresh']+'-'+node[n]['parent_call']+'&LEAF_NODE-'+n
            while(p != 0):
                n = p
                p = node[p]['parent']
                t = node[p]['decision']+'%'+node[p]['thresh']+'-'+node[n]['parent_call']
                traceback = t+'&'+traceback
            paths.append(traceback)
    ## Load path accuracies
    pathacc = {}
    with open('traceback', "r") as tb:
        for line in tb:
            line = line.strip()
            l = line.split("\t")
            l[2] = re.sub(r"\S+&(LEAF_NODE-\d+)$", "\g<1>", l[2])
            pathacc[l[2]] = {'pct': l[0],
                             'n': l[1]}
    ## Load ASM fastas        
    stach_fa = read_fasta(stach_file)
    seed_fa = read_fasta(seed_file)
    ## Split groups
    groups = split_into_groups(name2seq, threads)
    for group in groups:
        toprocess = {}
        for name in name2seq:
            if name in groups[group]:
                toprocess[name] = name2seq[name]
        p = multiprocessing.Process(target=sandpuma_multithreaded, args=(group, toprocess, knownfaa, wildcard, snn_thresh, knownasm, max_depth, min_leaf_sup, ref_aln, ref_tree, ref_pkg, masscutoff, stach_fa, seed_fa, clf, i2s, paths, pathacc))
        p.start()
