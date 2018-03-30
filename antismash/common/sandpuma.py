# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" SANDPUMA. """

import logging
import os
from typing import Any, Dict, List, Optional, Tuple

from antismash.common import fasta, module_results, pfamdb, subprocessing, pplacer
from antismash.common.secmet import Record, CDSFeature
from antismash.common.secmet.feature import FeatureLocation, PFAMDomain

from ete3 import Tree
import multiprocessing

class SandpumaResults(module_results.ModuleResults):
        """ Results for SANDPUMA """
        
        def __init__(self, record_id: str, predicat_monophyly: str, predicat_snn: str, snn_score: float, asm: str, svm: str, phmm: str, pid: float, ensemble: str, sandpuma: str) -> None:
            super().__init__(record_id)
            self.predicat_monophyly = str(predicat_monophyly)
            self.predicat_snn = str(predicat_snn)
            self.snn_score = float(snn_score)
            self.asm = str(asm)
            self.svm = str(svm)
            self.phmm = str(phmm)
            self.pid = float(pid)
            self.ensemble = str(ensemble)
            self.sandpuma = str(sandpuma)







def sandpuma_multithreaded(group: str, fasta: Dict[str, str], knownfaa: str, wildcard: str, snn_thresh: float, knownasm: str, max_depth: int, min_leaf_sup: int, jackknife_data: str):
    """ SANDPUMA

        Arguments:
            group: prefix group name
            fasta: dictionary of seq names (str) to seqs (str)
            knownfaa: filename for reference protein fasta; assumes each header ends in '_' followed by the <substrate specificity>
            wildcard: str to append to the end of each query sequence; should be different that all specificities (Default= 'UNK')
            snn_thresh: threshold for SNN score (Default= 0.5) NOTE: may need to be adjusted with new pplacer implementation
            knownasm: filename for reference active site motif protein fasta, similar header formatting as knownfaa
            max_depth: maximum depth for the sklearn decision tree; default= 40
            min_leaf_sup: minimum leaf support required within the decision tree; default= 10
            jackknife_data: filename for jackknife benchmarking results

        Returns:                                                                                                                             
    """

    for query in fasta:
        tmpfaa = group+'.tmp.faa'
        write_fasta(query+'_'+wildcard, fasta[query], tmpfaa)

        """
        run_predicat(tmpfaa)
        run_asm(tmpfaa)
        run_svm(tmpfaa)
        run_phmm(tmpfaa)
        run_pid(tmpfaa)
        run_ensemble()
        run_rescore()
        """


def split_into_groups(fasta: Dict[str, str], n_groups: int) -> Dict[str, List[str]]:
    n_seqs = length(fasta)
    seqs_per_group = int(n_seqs / n_groups)
    qnum = 0
    groupnum = 1
    groups = []
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

def run_sandpuma(name2seq: Dict[str, str], threads: int, knownfaa: str, wildcard: str, snn_thresh: float, knownasm: str, max_depth: int, min_leaf_sup: int, jackknife_data: str):
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

        Returns:                                                                                                                             
    """
 
    
    groups = split_into_groups(name2seq, threads)
    for group in groups:
        toprocess = {}
        for name in name2seq:
            if name in groups[group]:
                toprocess[name] = name2seq[name]
        p = multiprocessing.Process(target=sandpuma_multithreaded, args=(group, toprocess, knownfaa, wildcard))
        p.start()
