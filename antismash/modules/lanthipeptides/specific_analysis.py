# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
More detailed lanthipeptide analysis using HMMer-based leader peptide
cleavage site prediction as well as prediction of number of lanthionine
bridges and molcular mass.
"""

from collections import defaultdict
import logging
import os
import re
from typing import Any, Dict, List, Set, Optional

from antismash.common.signature import HmmSignature
from antismash.common import all_orfs, path, subprocessing, secmet, \
                             module_results, serialiser, utils
from antismash.common.fasta import get_fasta_from_features

from .rodeo import run_rodeo

KNOWN_PRECURSOR_DOMAINS = set([
    'Antimicr18',
    'Gallidermin',
    'L_biotic_A',
    'lacticin_l',
    'leader_d',
    'leader_abc',
    'leader_eh',
    'mature_ha',
    'lacticin_mat',
    'mature_b',
    'mature_d',
    'mature_ab',
    'mature_a',
    'TIGR03731',
    'LD_lanti_pre',
    'strep_PEQAXS',
])

THRESH_DICT = {'Class-I': -15,
               'Class-II': -7.3,
               'Class-III': -3.5}

# the maximal number of bases in each direction of a core enzyme that a
# precursor can be defined in
MAX_PRECURSOR_DISTANCE = 10000


class LanthiResults(module_results.ModuleResults):
    """ Holds the results of lanthipeptide analysis for a record

    """
    schema_version = 2

    def __init__(self, record_id, *args):
        super().__init__(record_id, *args)
        # keep new CDS features
        self.new_cds_features = set()
        # keep new CDSMotifs by the gene they match to
        # e.g. self.motifs_by_locus[gene_locus] = [motif1, motif2..]
        self.motifs_by_locus = defaultdict(list)
        # keep clusters and which genes in them had precursor hits
        # e.g. self.clusters[cluster_number] = {gene1_locus, gene2_locus}
        self.clusters = defaultdict(set)

    def to_json(self):
        cds_features = [(serialiser.location_to_json(feature.location),
                         feature.get_name()) for feature in self.new_cds_features]
        motifs = {}
        for locus, locus_motifs in self.motifs_by_locus.items():
            motifs[locus] = [motif.to_json() for motif in locus_motifs]
        return {"record_id": self.record_id,
                "schema_version": LanthiResults.schema_version,
                "motifs": motifs,
                "new_cds_features": cds_features,
                "clusters": {key: list(val) for key, val in self.clusters.items()}}

    @staticmethod
    def from_json(json, record) -> "LanthiResults":
        if json.get("schema_version") != LanthiResults.schema_version:
            logging.warning("Discarding Lanthipeptide results, schema version mismatch")
            return None
        results = LanthiResults(json["record_id"])
        for locus, motifs in json["motifs"].items():
            for motif in motifs:
                results.motifs_by_locus[locus].append(LanthipeptideMotif.from_json(motif))
        results.clusters = {int(key): set(val) for key, val in json["clusters"].items()}
        for location, name in json["new_cds_features"]:
            cds = all_orfs.create_feature_from_location(record, location, label=name)
            results.new_cds_features.add(cds)
        return results

    def add_to_record(self, record):
        for feature in self.new_cds_features:
            record.add_cds_feature(feature)

        for motifs in self.motifs_by_locus.values():
            for motif in motifs:
                record.add_cds_motif(motif)


class PrepeptideBase:
    """ A generic prepeptide class for tracking various typical components """
    def __init__(self, start, end, score, rodeo_score=None):
        self.start = start  # same as CDS
        self.end = end  # same as CDS
        self.score = score  # of cleavage site
        self.rodeo_score = rodeo_score
        self._leader = None
        self._core = ''
        self._tail = None
        self._lan_bridges = -1
        self._weight = -1
        self._monoisotopic_weight = -1
        self._alt_weights = None
        self.core_analysis_monoisotopic = None
        self.core_analysis = None

    @property
    def core(self) -> str:
        """ The sequence of the prepeptide core """
        return self._core

    @core.setter
    def core(self, seq: str):
        self.core_analysis_monoisotopic = utils.RobustProteinAnalysis(seq, monoisotopic=True)
        self.core_analysis = utils.RobustProteinAnalysis(seq, monoisotopic=False)
        self._core = seq
        self._calculate_mw()

    @property
    def leader(self) -> Optional[str]:
        """ The sequence of the prepeptide leader """
        return self._leader

    @leader.setter
    def leader(self, seq: str):
        self._leader = seq

    def __repr__(self):
        return "PrepeptideBase(%s..%s, %s, %r, %r)" % (self.start, self.end,
                                      self.score, self.rodeo_score, self._core)

    @property
    def number_of_lan_bridges(self):
        """
        function determines the number of lanthionine bridges in the core peptide
        """
        raise NotImplementedError()

    def _calculate_mw(self):
        """
        (re)calculate the monoisotopic mass and molecular weight
        """
        assert self._core
        raise NotImplementedError()

    @property
    def monoisotopic_mass(self):
        """ weight of the dehydrated core
        """
        if self._monoisotopic_weight is None:
            raise ValueError("No core to calculate weight of")

        self._calculate_mw()
        return self._monoisotopic_weight

    @property
    def molecular_weight(self):
        """ Weight of the dehydrated core
        """
        if self._weight is None:
            raise ValueError("No core to calculate weight of")

        self._calculate_mw()
        return self._weight

    @property
    def alternative_weights(self):
        """ The possible alternative weights assuming one or more of the Ser/Thr
            residues aren't dehydrated
        """
        if self._alt_weights is None:
            raise ValueError("No core to calculate weight of")
        return self._alt_weights


class Lanthipeptide(PrepeptideBase):
    """ Calculates and stores lanthipeptide information
    """
    def __init__(self, start, end, score, rodeo_score, lantype):
        super().__init__(start, end, score, rodeo_score)
        self.lantype = lantype
        self._aminovinyl = False
        self._chlorinated = False
        self._oxygenated = False
        self._lac = False

    def __repr__(self):
        return "Lanthipeptide(%s..%s, %s, %r, %r, %s, %s(%s))" % (self.start,
                    self.end, self.score, self.lantype, self._core,
                    self._lan_bridges, self._monoisotopic_weight, self._weight)

    @property
    def number_of_lan_bridges(self) -> int:
        """ Determines the number of lanthionine bridges in the core peptide
        """
        if not self._core:
            raise ValueError("No core to calculate bridges from")

        amino_counts = self.core_analysis.count_amino_acids()
        no_cys = amino_counts['C']
        no_thr_ser = amino_counts['T'] + amino_counts['S']
        self._lan_bridges = min(no_cys, no_thr_ser)
        if self._aminovinyl:
            self._lan_bridges -= 1
        return self._lan_bridges

    def _calculate_mw(self) -> None:
        """ (re)calculates the monoisotopic mass and molecular weight
        """
        assert self._core, "calculating weight without a core"

        amino_counts = self.core_analysis.count_amino_acids()
        no_thr_ser = amino_counts['T'] + amino_counts['S']

        mol_mass = self.core_analysis.molecular_weight()
        mods = 18.02 * no_thr_ser
        if self._aminovinyl:
            mods += 46
        if self._chlorinated:
            mods -= 34
        if self._oxygenated:
            mods -= 16
        if self._lac:
            mods -= 2
        self._weight = mol_mass - mods

        # every unbridged Ser or Thr might not be dehydrated
        self._alt_weights = []
        for i in range(1, no_thr_ser - amino_counts['C'] + 1):
            self._alt_weights.append(self._weight + 18.02 * i)

        monoisotopic_mass = self.core_analysis_monoisotopic.molecular_weight()
        mods = 18 * no_thr_ser
        if self._aminovinyl:
            mods += 46
        if self._chlorinated:
            mods -= 34
        if self._oxygenated:
            mods -= 16
        if self._lac:
            mods -= 2
        self._monoisotopic_weight = monoisotopic_mass - mods

    @property
    def aminovinyl_group(self) -> bool:
        """ Returns True if lanthipeptide contains an aminovinyl group
        """
        return self._aminovinyl

    @aminovinyl_group.setter
    def aminovinyl_group(self, value: bool) -> None:
        """ Sets whether lanthipeptide contains an aminovinyl group and triggers
            recalculation of the molecular weight
        """
        self._aminovinyl = value
        if self._core:
            self._calculate_mw()
            # recalculate the number of lan bridges
            self._lan_bridges = -1

    @property
    def chlorinated(self) -> bool:
        """ Returns True if lanthipeptide is chlorinated
        """
        return self._chlorinated

    @chlorinated.setter
    def chlorinated(self, value: bool) -> None:
        """ Sets whether lanthipeptide is chlorinated and triggers
            recalculation of the molecular weight
        """
        self._chlorinated = value
        if self._core:
            self._calculate_mw()

    @property
    def oxygenated(self) -> bool:
        """ Returns True if lanthipeptide is oxygenated
        """
        return self._oxygenated

    @oxygenated.setter
    def oxygenated(self, value: bool) -> None:
        """ Sets whether lanthipeptide is oxygenated and triggers
            recalculation of the molecular weight
        """
        self._oxygenated = value
        if self._core:
            self._calculate_mw()

    @property
    def lactonated(self) -> bool:
        """ Returns True if lanthipeptide starts with a lactone
        """
        return self._lac

    @lactonated.setter
    def lactonated(self, value: bool):
        """ Sets whether lanthipeptide has a lactone and triggers
            recalculation of the molecular weight
        """
        self._lac = value
        if self._core:
            self._calculate_mw()


class CleavageSiteHit(object):
    """ A simple container for storing cleavage site information """
    def __init__(self, start, end, score, lantype):
        self.start = start
        self.end = end
        self.score = score
        self.lantype = lantype

    def __repr__(self):
        return "CleavageSiteHit(start=%s, end=%s, score=%s, lantype='%s')" % (
                    self.start, self.end, self.score, self.lantype)


def get_detected_domains(genes: List[secmet.CDSFeature]) -> List[str]:
    """ Gathers all detected domains in a cluster, including some not detected
        by hmm_detection.

        Arguments:
            genes: a list of genes to check

        Returns:
            a list of strings, each string being the name of a domain in the
            cluster
    """
    found_domains = []  # type: List[str]
    # Gather biosynthetic domains
    for feature in genes:
        if not feature.sec_met:
            continue
        found_domains.extend(feature.sec_met.domain_ids)

    # Gather non-biosynthetic domains
    cluster_fasta = get_fasta_from_features(genes)
    assert cluster_fasta
    non_biosynthetic_hmms_by_id = run_non_biosynthetic_phmms(cluster_fasta)
    non_biosynthetic_hmms_found = []  # type: List[str]
    for hsps_found in non_biosynthetic_hmms_by_id.values():
        for hsp in hsps_found:
            if hsp not in non_biosynthetic_hmms_found:
                non_biosynthetic_hmms_found.append(hsp)
    found_domains += non_biosynthetic_hmms_found

    return found_domains


def run_non_biosynthetic_phmms(fasta: str) -> Dict[str, List[str]]:
    """ Finds lanthipeptide-specific domains in the input fasta

        Arguments:
            fasta: a string containing gene sequences in fasta format

        Returns:
            a dictionary mapping the hit id to a list of matching query ids
    """
    with open(path.get_full_path(__file__, "data", "non_biosyn_hmms", "hmmdetails.txt"), "r") as handle:
        hmmdetails = [line.strip().split("\t") for line in handle if line.count("\t") == 3]
    signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]
    non_biosynthetic_hmms_by_id = defaultdict(list)  # type: Dict[str, List[str]]
    for sig in signature_profiles:
        sig.path = path.get_full_path(__file__, "data", "non_biosyn_hmms", sig.path.rpartition(os.sep)[2])
        runresults = subprocessing.run_hmmsearch(sig.path, fasta)
        for runresult in runresults:
            # Store result if it is above cut-off
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    non_biosynthetic_hmms_by_id[hsp.hit_id].append(hsp.query_id)
    return non_biosynthetic_hmms_by_id


def predict_cleavage_site(query_hmmfile, target_sequence, threshold=-100) -> Optional[CleavageSiteHit]:
    """ Extracts from HMMER the start position, end position and score
        of the HMM alignment for a cleavage site

        Arguments:
            query_hmmfile: the path to a HMM file for the cleavage site profile
            target_sequence: the sequence of a CDS feature
            threshold: a minimum bitscore for a HMMer hit, exclusive

        Returns:
            a CleavageSiteHit instance with the information about the hit, or
            None if no hit was above the threshold
    """
    hmmer_res = subprocessing.run_hmmpfam2(query_hmmfile, target_sequence)

    for res in hmmer_res:
        for hits in res:
            lanthi_type = hits.description
            for hsp in hits:
                if hsp.bitscore > threshold:
                    return CleavageSiteHit(hsp.query_start - 1, hsp.query_end, hsp.bitscore, lanthi_type)
    return None


def predict_class_from_genes(focus: secmet.CDSFeature, genes: List[secmet.CDSFeature]) -> Optional[str]:
    """ Predict the lanthipeptide class from the gene cluster

        Arguments:
            genes: a list of genes to check

        Returns:
            a string representing the class, or None if no class predicted
    """
    found_domains = set()
    for feature in genes + [focus]:
        if not feature.sec_met:
            continue
        found_domains.update(set(feature.sec_met.domain_ids))

    if 'Lant_dehyd_N' in found_domains or 'Lant_dehyd_C' in found_domains:
        return 'Class-I'
    if 'DUF4135' in found_domains:
        return 'Class-II'
    if 'Pkinase' in found_domains:
        # this could be class 3 or class 4, but as nobody has seen class 4
        # in vivo yet, we'll ignore that
        return 'Class-III'

    return None


def run_cleavage_site_phmm(fasta, hmmer_profile, threshold):
    """ Try to identify cleavage site using pHMM """
    profile = path.get_full_path(__file__, hmmer_profile)
    return predict_cleavage_site(profile, fasta, threshold)


def run_cleavage_site_regex(fasta):
    """ Try to identify cleavage site using regular expressions"""
    # Regular expressions; try 1 first, then 2, etc.
    rex1 = re.compile('F?LD')
    rex2 = re.compile('[LF]?LQ')

    # For regular expression, check if there is a match that is <10 AA from the end
    if re.search(rex1, fasta) and len(re.split(rex1, fasta)[-1]) > 10:
        start, end = [m.span() for m in rex1.finditer(fasta)][-1]
        end += 16
    elif re.search(rex2, fasta) and len(re.split(rex2, fasta)[-1]) > 10:
        start, end = [m.span() for m in rex2.finditer(fasta)][-1]
        end += 15
    else:
        return [None, None, None]
    return start, end, 0


def determine_precursor_peptide_candidate(record: secmet.Record, query: secmet.CDSFeature,
                                          query_sequence: str, domains: List[str],
                                          hmmer_profile: str) -> Optional[Lanthipeptide]:
    """ Identify precursor peptide candidates and split into two,
        only valid for Class-I lanthipeptides
    """

    # Skip sequences with >200 AA
    if len(query_sequence) > 200 or len(query_sequence) < 20:
        return None

    # Create FASTA sequence for feature under study
    lan_a_fasta = ">%s\n%s" % (query.get_name(), query_sequence)

    # Run sequence against pHMM; if positive, parse into a vector containing START, END and SCORE
    cleavage_result = run_cleavage_site_phmm(lan_a_fasta, hmmer_profile, THRESH_DICT["Class-I"])

    if cleavage_result is not None and cleavage_result.end <= len(query_sequence) - 8:
        start = cleavage_result.start
        end = cleavage_result.end
        score = cleavage_result.score
        lanthi_type = cleavage_result.lantype
    else:
        # If no pHMM hit, try regular expression
        start, end, score = run_cleavage_site_regex(lan_a_fasta)
        if score is None or end > len(query_sequence) - 8:
            # abort, since RODEO will predict duplicates based only on cluster
            # attributes
            return None
        lanthi_type = "lanthipeptide"

    # if the cleavage results in no core, that's not valid
    if end == len(query_sequence):
        return None

    # Run RODEO to assess whether candidate precursor peptide is judged real
    rodeo_result = run_rodeo(record, query, query_sequence[:end], query_sequence[end:], domains)
    if rodeo_result < 14:
        return None
    lanthipeptide = Lanthipeptide(start, end, score, rodeo_result, lanthi_type)

    # Determine the leader and core peptide
    lanthipeptide.leader = query_sequence[:end]
    lanthipeptide.core = query_sequence[end:]

    return lanthipeptide


def run_lanthipred(record: secmet.Record, query: secmet.CDSFeature, lant_class, domains):
    """ Determines if a CDS is a predicted lanthipeptide based on the class
        and any contained domains.

        Arguments:
            record: the parent Record of the feature
            query: the CDSFeature to analyse
            lant_class: a string representing the class
            domains: a list of domain names in the current cluster
    """
    hmmer_profiles = {'Class-I': 'data/class1.hmm',
                      'Class-II': 'data/class2.hmm',
                      'Class-III': 'data/class3.hmm', }
    query_sequence = query.translation
    lan_a_fasta = ">%s\n%s" % (query.get_name(), query_sequence)

    if lant_class in ("Class-II", "Class-III"):
        profile = path.get_full_path(__file__, hmmer_profiles[lant_class])
        cleavage_result = predict_cleavage_site(profile, lan_a_fasta)

        if cleavage_result is None:
            return None

        if THRESH_DICT[lant_class] > cleavage_result.score:
            return None

        # if the cleavage results in no core, that's not valid
        if cleavage_result.end == len(query_sequence):
            return None

        result = Lanthipeptide(cleavage_result.start, cleavage_result.end,
                               cleavage_result.score, "N/A", lant_class)
        result.leader = query_sequence[:result.end]
        result.core = query_sequence[result.end:]

    else:
        result = determine_precursor_peptide_candidate(record, query, query_sequence,
                                                       domains, hmmer_profiles[lant_class])
        if result is None:
            return None

    # extract now (that class is known and thus the END component) the core peptide
    if result.number_of_lan_bridges == 0:
        return None

    query.gene_functions.add(secmet.GeneFunction.ADDITIONAL, "lanthipeptides",
                             "predicted lanthipeptide")
    return result


def find_lan_a_features(cluster: secmet.Cluster) -> List[secmet.CDSFeature]:
    """ Finds all lanthipeptide candidate features """
    lan_a_features = []
    for feature in cluster.cds_children:
        if not feature.is_contained_by(cluster):
            continue

        if len(feature.translation) < 80:
            lan_a_features.append(feature)
            continue
        if feature.sec_met and set(feature.sec_met.domain_ids).intersection(KNOWN_PRECURSOR_DOMAINS):
            lan_a_features.append(feature)

    return lan_a_features


def contains_feature_with_single_domain(genes: List[secmet.CDSFeature], domains: Set[str]) -> bool:
    """ Checks for the existence of a feature within a group that has a single
        domain and that the domain is within the provided set of domains

        Arguments:
            genes: a list of genes to check
            domains: the set of domain names allowable

        Returns:
            True if a feature matching the conditions was found, otherwise False
    """
    for feature in genes:
        if not feature.sec_met or len(feature.sec_met.domain_ids) > 1:
            continue
        if len(domains.intersection(set(feature.sec_met.domain_ids))) == 1:
            return True
    return False


class LanthipeptideMotif(secmet.Prepeptide):
    """ A lanthipeptide-specific feature """
    def __init__(self, location, core_seq, leader_seq,
                 locus_tag, monoisotopic_mass, molecular_weight, alternative_weights,
                 lan_bridges, lanthi_class, score, rodeo_score, aminovinyl,
                 chlorinated, oxygenated, lactonated):
        super().__init__(location, "lanthipeptide", core_seq, locus_tag, lanthi_class,
                         score=score, monoisotopic_mass=monoisotopic_mass,
                         molecular_weight=molecular_weight,
                         alternative_weights=alternative_weights,
                         leader=leader_seq)
        self.lan_bridges = lan_bridges
        self.rodeo_score = rodeo_score
        self.aminovinyl_group = aminovinyl  # bool
        self.chlorinated = chlorinated  # bool
        self.oxygenated = oxygenated  # bool
        self.lactonated = lactonated  # bool
        self._notes_appended = False

    def get_modifications(self) -> List[str]:
        """ Returns the various modifications of the lanthipeptide, if they exist
        """
        mods = []
        if self.aminovinyl_group:
            mods.append("AviCys")
        if self.chlorinated:
            mods.append("Cl")
        if self.oxygenated:
            mods.append("OH")
        if self.lactonated:
            mods.append("Lac")
        return mods

    def to_biopython(self, qualifiers: Dict[str, List] = None):
        notes = []
        if not qualifiers:
            qualifiers = {}
        notes.append('number of bridges: %s' % self.lan_bridges)
        notes.append('RODEO score: %s' % str(self.rodeo_score))
        if self.aminovinyl_group:
            notes.append('predicted additional modification: AviCys')
        if self.chlorinated:
            notes.append('predicted additional modification: Cl')
        if self.oxygenated:
            notes.append('predicted additional modification: OH')
        if self.lactonated:
            notes.append('predicted additional modification: Lac')
        if "note" not in qualifiers:
            qualifiers["note"] = notes
        else:
            qualifiers["note"].extend(notes)
        return super().to_biopython(qualifiers=qualifiers)

    def to_json(self):
        json = super().to_json()
        json["locus_tag"] = self.locus_tag  # not in vars() due to __slots__
        try:
            assert json["locus_tag"]
        except KeyError:
            logging.critical("bad locus tag on motif %s: %s ... %s", self.location, self.locus_tag, json)
        return json

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "LanthipeptideMotif":
        """ Converts a JSON representation of the motif back into an instance
            of LanthipeptideMotif
        """
        args = []
        args.append(serialiser.location_from_json(data["location"]))
        args.append(data["core"])
        for arg_name in ["leader", "locus_tag", "monoisotopic_mass",
                         "molecular_weight", "alternative_weights", "lan_bridges",
                         "peptide_subclass", "score", "rodeo_score", "aminovinyl_group",
                         "chlorinated", "oxygenated", "lactonated"]:
            args.append(data[arg_name])
        # pylint doesn't do well with the splat op, so don't report errors
        return LanthipeptideMotif(*args)  # pylint: disable=no-value-for-parameter


def result_vec_to_feature(orig_feature: secmet.CDSFeature, res_vec: Lanthipeptide) -> LanthipeptideMotif:
    """ Generates a LanthipeptideMotif feature from a CDSFeature and a Lanthipeptide

        Arguments:
            orig_feature: the CDSFeature the lanthipeptide was found in
            res_vec: the Lanthipeptide instance that was calculated

        Returns:
            a LanthipeptideMotif instance
    """
    feature = LanthipeptideMotif(orig_feature.location, res_vec.core, res_vec.leader,
                                 orig_feature.get_name(), res_vec.monoisotopic_mass,
                                 res_vec.molecular_weight, res_vec.alternative_weights,
                                 res_vec.number_of_lan_bridges, res_vec.lantype,
                                 res_vec.score, res_vec.rodeo_score, res_vec.aminovinyl_group,
                                 res_vec.chlorinated, res_vec.oxygenated, res_vec.lactonated)
    return feature


def find_neighbours_in_range(center: secmet.CDSFeature,
                             candidates: List[secmet.CDSFeature]) -> List[secmet.CDSFeature]:
    """ Restrict a set of genes to those within precursor range of a central
        gene.

        Arguments:
            center: the gene to find the neighbours of
            candidates: the genes to filter by range

        Returns:
            a list of genes within range, with the same ordering as the input
    """
    neighbours = []
    for candidate in candidates:
        if candidate < center:
            if center.location.start - candidate.location.start <= MAX_PRECURSOR_DISTANCE:
                neighbours.append(candidate)
        else:
            if candidate.location.end - center.location.end <= MAX_PRECURSOR_DISTANCE:
                neighbours.append(candidate)
            else:
                # skip looking further to the right if the previous one was too far away
                break
    return neighbours


def run_lanthi_on_genes(record: secmet.Record, focus: secmet.CDSFeature,
                        genes: List[secmet.CDSFeature], results: LanthiResults) -> None:
    """ Runs lanthipeptide around a single focus gene which is a core biosynthetic
        enzyme for lanthipeptides.
        Updates the results object with any precursors found.

        Arguments:
            record: the Record instance containing the genes
            focus: a core lanthipeptide gene
            genes: a list of candidate precursor genes
            results: a LanthiResults object to update

        Returns:
            None
    """
    domains = get_detected_domains(genes)
    non_candidate_neighbours = find_neighbours_in_range(focus, focus.cluster.cds_children)
    flavoprotein_found = contains_feature_with_single_domain(non_candidate_neighbours, {"Flavoprotein"})
    halogenase_found = contains_feature_with_single_domain(non_candidate_neighbours, {"Trp_halogenase"})
    oxygenase_found = contains_feature_with_single_domain(non_candidate_neighbours, {"p450"})
    dehydrogenase_found = contains_feature_with_single_domain(non_candidate_neighbours, {"adh_short", "adh_short_C2"})

    lant_class = predict_class_from_genes(focus, genes)
    if not lant_class:
        return

    for candidate in genes:
        result_vec = run_lanthipred(record, candidate, lant_class, domains)
        if result_vec is None:
            continue
        result_vec.aminovinyl_group = flavoprotein_found
        result_vec.chlorinated = halogenase_found
        result_vec.oxygenated = oxygenase_found
        result_vec.lactonated = dehydrogenase_found and result_vec.core.startswith('S')
        motif = result_vec_to_feature(candidate, result_vec)
        results.motifs_by_locus[focus.get_name()].append(motif)
        results.clusters[focus.cluster.get_cluster_number()].add(focus.get_name())
        # track new CDSFeatures if found with all_orfs
        if candidate.cluster is None:
            results.new_cds_features.add(candidate)


def run_specific_analysis(record: secmet.Record) -> LanthiResults:
    """ Runs the full lanthipeptide analysis over the given record

        Arguments:
            record: the Record instance to analyse

        Returns:
            A populated LanthiResults object
    """
    results = LanthiResults(record.id)
    for cluster in record.get_clusters():
        if 'lanthipeptide' not in cluster.products:
            continue

        # find core biosynthetic enzyme locations
        core_domain_names = {'Lant_dehyd_N', 'Lant_dehyd_C', 'DUF4135', 'Pkinase'}
        core_genes = []
        for gene in cluster.cds_children:
            if not gene.sec_met:
                continue
            if core_domain_names.intersection(set(gene.sec_met.domain_ids)):
                core_genes.append(gene)

        precursor_candidates = find_lan_a_features(cluster)
        # Find candidate ORFs that are not yet annotated
        extra_orfs = all_orfs.find_all_orfs(record, cluster)
        for orf in extra_orfs:
            if len(orf.translation) < 80:
                precursor_candidates.append(orf)

        for gene in core_genes:
            neighbours = find_neighbours_in_range(gene, precursor_candidates)
            run_lanthi_on_genes(record, gene, neighbours, results)

    logging.debug("Lanthipeptide module marked %d motifs", sum(map(len, results.motifs_by_locus)))
    return results
