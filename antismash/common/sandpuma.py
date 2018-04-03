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

from ete3 import Tree
import multiprocessing

class PredicatResults(module_results.ModuleResults):
        """ Results for prediCAT """
        def __init__(self, record_id: str, monophyly: str, forced: str, nndist: float, nn_score: float, snn_score: float) -> None:
            super().__init__(record_id)
            self.monophyly = str(monophyly)
            self.forced = str(forced)
            self.nndist = float(nndist)
            self.nn_score = float(nn_score)
            self.snn_score = float(snn_score)

        
class SandpumaResults(module_results.ModuleResults):
        """ Results for SANDPUMA """
        
        def __init__(self, record_id: str, predicat_monophyly: str, predicat_snn: str, snn_score: float, asm: str, svm: str, phmm: str, pid: float, ensemble: str, sandpuma: str) -> None:
            super().__init__(record_id)
            self.predicat = PredicatResults
            self.predicat_snn = str(predicat_snn)
            self.snn_score = float(snn_score)
            self.asm = str(asm)
            self.svm = str(svm)
            self.phmm = str(phmm)
            self.pid = float(pid)
            self.ensemble = str(ensemble)
            self.sandpuma = str(sandpuma)

def run_predicat(reference_aln: str, queryfa: Dict[str, str], wildcard: str, leaf2clade_fi: str, ref_tree: str, ref_pkg: str, masscutoff: float) -> PredicatResults:
        """ SANDPUMA parallelized pipleline

        Arguments:
            reference_aln: filename for reference protein fasta, see sandpuma_multithreaded comments for requirements
            tmpfaa: filename for the single query faa
            wildcard: suffix str identifying query sequence (Default= 'UNK' which means headers end in '_UNK')
            leaf2clade_fi: leaf to clade mapping file
            ref_tree: reference tree (newick)
            ref_pkg: pplacer reference package
            masscutoff: cutoff value for pplacer masses


        Returns:                                                                                                                             PredicatResults
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
        all_aligned = subprocessing.run_muscle_profile_predicat(reference_aln, trimmedfa)
        
        ## Pplacer (NOTE: this is new to SANDPUMA as of antiSMASH5 and needs to be tested
        clade_assignment = pplacer.pplacer_clade_assignment(leaf2clade_fi, ref_tree, reference_aln, ref_pkg, all_aligned, masscutoff)

        ## MORE HERE SOON!


        
        
def sandpuma_multithreaded(group: str, fasta: Dict[str, str], knownfaa: str, wildcard: str, snn_thresh: float, knownasm: str, max_depth: int, min_leaf_sup: int, jackknife_data: str, ref_aln: str, leaf2clade_fi: str, ref_tree: str, ref_pkg: str, masscutoff: float):
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
            jackknife_data: filename for jackknife benchmarking results
            ref_aln: reference alignment (fasta) file
            leaf2clade_fi: leaf to clade mapping file
            ref_tree: reference tree (newick)
            ref_pkg: pplacer reference package
            masscutoff: cutoff value for pplacer masses


        Returns:                                                                                                                             
    """

    for query in fasta:
        wc_name = query+'_'+wildcard
        ## Write fasta for functions that need files
        tmpfaa = group+'.tmp.faa'
        write_fasta(wc_name, fasta[query], tmpfaa) ## Single entry fasta

        ## Store as a dictionary for functions that don't
        queryfa = {wc_name: fasta[query]}

        ## PrediCAT
        run_predicat(ref_aln, queryfa, wildcard, leaf2clade_fi, ref_tree, ref_pkg, masscutoff)
        



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

def run_sandpuma(name2seq: Dict[str, str], threads: int, knownfaa: str, wildcard: str, snn_thresh: float, knownasm: str, max_depth: int, min_leaf_sup: int, jackknife_data: str, ref_aln: str, leaf2clade_fi: str, ref_tree: str, ref_pkg: str, masscutoff:float):
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
            leaf2clade_fi: leaf to clade mapping file
            ref_tree: reference tree (newick)
            ref_pkg: pplacer reference package
            masscutoff: cutoff value for pplacer masses

        Returns:                                                                                                                             
    """
    
    groups = split_into_groups(name2seq, threads)
    for group in groups:
        toprocess = {}
        for name in name2seq:
            if name in groups[group]:
                toprocess[name] = name2seq[name]
        p = multiprocessing.Process(target=sandpuma_multithreaded, args=(group, toprocess, knownfaa, wildcard, snn_thresh, knownasm, max_depth, min_leaf_sup, jackknife_data, ref_aln, leaf2clade_fi, ref_tree, ref_pkg, masscutoff))
        p.start()
