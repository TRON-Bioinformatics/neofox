#!/usr/bin/env python

import os
import os.path
from typing import List

from logzero import logger

from neofox.helpers import intermediate_files
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.literature_features.differential_binding import DifferentialBinding
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


class NeoantigenFitnessCalculator(BlastpRunner):

    def __init__(self, runner, configuration, iedb):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        super().__init__(runner, configuration)
        self.iedb = iedb
        self.differential_binding = DifferentialBinding()

    def _calc_pathogen_similarity(self, fasta_file):
        """
        This function determines the PATHOGENSIMILARITY of epitopes according to Balachandran et al. using a blast
        search against the IEDB pathogenepitope database
        """
        outfile = self.run_blastp(fasta_file=fasta_file, database=os.path.join(self.iedb, "iedb_blast_db"))
        similarity = self.parse_blastp_output(blastp_output_file=outfile)
        os.remove(outfile)
        return similarity

    def wrap_pathogen_similarity(self, mutation):
        fastafile = intermediate_files.create_temp_fasta(sequences=[mutation], prefix="tmpseq", comment_prefix='M_')
        pathsim = None
        try:
            pathsim = self._calc_pathogen_similarity(fastafile)
        except Exception as ex:
            # TODO: do we need this at all? it should not fail and if it fails we probably want to just stop execution
            logger.exception(ex)
        os.remove(fastafile)
        logger.info("Peptide {} has a pathogen similarity of {}".format(mutation, pathsim))
        return pathsim

    def calculate_amplitude_mhc(self, score_mutation, score_wild_type, apply_correction=False):
        """
        This function calculates the amplitude between mutated and wt epitope according to Balachandran et al.
        when affinity is used, use correction from Luksza et al. *1/(1+0.0003*aff_wt)
        """
        amplitude_mhc = "NA"
        try:
            candidate_amplitude_mhc = float(score_wild_type) / float(score_mutation)
            if apply_correction:  #nine_mer or affinity:
                amplitude_mhc = str(candidate_amplitude_mhc * (self._calculate_correction(score_wild_type)))
            else:
                amplitude_mhc = str(candidate_amplitude_mhc)
        except(ZeroDivisionError, ValueError) as e:
            pass
        return amplitude_mhc

    def _calculate_correction(self, score_wild_type):
        return 1 / (1 + 0.0003 * float(score_wild_type))

    def calculate_recognition_potential(
            self, amplitude, pathogen_similarity, mutation_in_anchor, mhc_affinity_mut=None):
        """
        This function calculates the recognition potential, defined by the product of amplitude and pathogensimiliarity of an epitope according to Balachandran et al.
        F_alpha = - max (A_i x R_i)

        Returns (A_i x R_i) value only for nonanchor mutation and epitopes of length 9; only considered by Balachandran
        """
        recognition_potential = None
        try:
            candidate_recognition_potential = str(float(amplitude) * float(pathogen_similarity))
            if mhc_affinity_mut:
                if mutation_in_anchor == "0" and float(mhc_affinity_mut) < 500.0:
                    recognition_potential = candidate_recognition_potential
            else:
                if mutation_in_anchor == "0":
                    recognition_potential = candidate_recognition_potential
        except ValueError:
            pass
        return recognition_potential

    def get_annotations(self,
                        binding_wild_type, binding_mutation,
                        affinity_wild_type, affinity_mutation, rank_wild_type, rank_mutation,
                        affinity_9mer_wild_type, affinity_9mer_mutation,
                        sequence_9mer_mutation: str, sequence_9mer_wild_type: str,
                        sequence_mutation: str, sequence_wild_type: str,
                        sequence_rank_mutation: str, sequence_rank_wild_type: str) -> List[Annotation]:

        amplitude_affinity = self.calculate_amplitude_mhc(
            score_mutation=affinity_mutation, score_wild_type=affinity_wild_type, apply_correction=True)
        amplitude_nemhcpan_rank = self.calculate_amplitude_mhc(
            score_mutation=rank_mutation, score_wild_type=rank_wild_type)
        amplitude_best_affinity_9mers = self.calculate_amplitude_mhc(
            score_mutation=affinity_9mer_mutation, score_wild_type=affinity_9mer_wild_type, apply_correction=True)
        pathogen_similarity_9mer = self.wrap_pathogen_similarity(mutation=sequence_9mer_mutation)
        pathogen_similarity_rank = self.wrap_pathogen_similarity(mutation=sequence_rank_mutation)
        pathogen_similarity_affinity = self.wrap_pathogen_similarity(mutation=sequence_mutation)

        position = EpitopeHelper.position_of_mutation_epitope(wild_type=sequence_wild_type, mutation=sequence_mutation)
        position_9mer = EpitopeHelper.position_of_mutation_epitope(
            wild_type=sequence_9mer_wild_type, mutation=sequence_9mer_mutation)
        position_rank = EpitopeHelper.position_of_mutation_epitope(
            wild_type=sequence_rank_wild_type, mutation=sequence_rank_mutation)

        bdg_cutoff_classical_mhci = 50
        bdg_cutoff_alternative_mhci = 5000
        amplitude_cutoff_mhci = 10

        return [
            AnnotationFactory.build_annotation(name="CDN_mhcI", value=self.differential_binding.classify_adn_cdn(
                score_mutation=affinity_mutation, amplitude=amplitude_affinity,
                bdg_cutoff_classical=bdg_cutoff_classical_mhci, bdg_cutoff_alternative=bdg_cutoff_alternative_mhci,
                amplitude_cutoff=amplitude_cutoff_mhci, category="CDN")),
            AnnotationFactory.build_annotation(name="ADN_mhcI", value=self.differential_binding.classify_adn_cdn(
                score_mutation=affinity_mutation, amplitude=amplitude_affinity,
                bdg_cutoff_classical=bdg_cutoff_classical_mhci, bdg_cutoff_alternative=bdg_cutoff_alternative_mhci,
                amplitude_cutoff=amplitude_cutoff_mhci, category="ADN")),
            AnnotationFactory.build_annotation(name="Amplitude_mhcI_MB", value=self.calculate_amplitude_mhc(
                score_mutation=binding_mutation, score_wild_type=binding_wild_type)),
            AnnotationFactory.build_annotation(
                name="DAI_mhcI_MB", value=self.differential_binding.dai(
                    score_mutation=binding_mutation, score_wild_type=binding_wild_type)),
            AnnotationFactory.build_annotation(
                name="DAI_affinity_filtered", value=self.differential_binding.dai(
                    score_mutation=affinity_mutation, score_wild_type=affinity_wild_type, affin_filtering=True)),
            AnnotationFactory.build_annotation(
                name="DAI_affinity", value=self.differential_binding.dai(
                    score_mutation=affinity_mutation, score_wild_type=affinity_wild_type)),
            AnnotationFactory.build_annotation(
                name="DAI_rank_netmhcpan4", value=self.differential_binding.dai(
                    score_mutation=rank_mutation, score_wild_type=rank_wild_type)),
            AnnotationFactory.build_annotation(name="Amplitude_mhcI_affinity", value=amplitude_affinity),
            AnnotationFactory.build_annotation(name="Amplitude_mhcI_rank_netmhcpan4", value=amplitude_nemhcpan_rank),
            AnnotationFactory.build_annotation(name="Amplitude_mhcI_affinity_9mer_netmhcpan4",
                                               value=amplitude_best_affinity_9mers),
            AnnotationFactory.build_annotation(name="Pathogensimiliarity_mhcI_9mer", value=pathogen_similarity_9mer),
            AnnotationFactory.build_annotation(name="Pathogensimiliarity_mhcI_rank", value=pathogen_similarity_rank),
            AnnotationFactory.build_annotation(name="Pathogensimiliarity_mhcI_affinity_nmers",
                                               value=pathogen_similarity_affinity),
            AnnotationFactory.build_annotation(name="Recognition_Potential_mhcI_affinity",
                                               value=self.calculate_recognition_potential(
                                                   amplitude=amplitude_affinity,
                                                   pathogen_similarity=pathogen_similarity_affinity,
                                                   mutation_in_anchor=EpitopeHelper.position_in_anchor_position(
                                                       position_mhci=position,
                                                       peptide_length=len(sequence_mutation)))),
            AnnotationFactory.build_annotation(name="Recognition_Potential_mhcI_rank_netmhcpan4",
                                               value=self.calculate_recognition_potential(
                                                   amplitude=amplitude_nemhcpan_rank,
                                                   pathogen_similarity=pathogen_similarity_rank,
                                                   mutation_in_anchor=EpitopeHelper.position_in_anchor_position(
                                                       position_mhci=position_rank,
                                                       peptide_length=len(sequence_rank_mutation)))),
            AnnotationFactory.build_annotation(name="Recognition_Potential_mhcI_9mer_affinity",
                                               value=self.calculate_recognition_potential(
                                                   amplitude=amplitude_best_affinity_9mers,
                                                   pathogen_similarity=pathogen_similarity_9mer,
                                                   mutation_in_anchor=EpitopeHelper.position_in_anchor_position(
                                                       position_mhci=position_9mer,
                                                       peptide_length=len(sequence_9mer_mutation)),
                                                   mhc_affinity_mut=affinity_9mer_mutation))
        ]

    def get_annotations_mch2(self, mut_score, wt_score, aff_mut, aff_wt, sc_mut, sc_wt, epitope_mut_mhcii) -> List[Annotation]:

        pathogen_similarity = self.wrap_pathogen_similarity(mutation=epitope_mut_mhcii)
        amplitude_affinity = self.calculate_amplitude_mhc(
            score_mutation=aff_mut, score_wild_type=aff_wt, apply_correction=True)
        bdg_cutoff_classical_mhcii = 1
        bdg_cutoff_alternative_mhcii = 4
        amplitude_cutoff_mhcii = 4

        return [
            AnnotationFactory.build_annotation(
                value=self.calculate_amplitude_mhc(
                    score_mutation=mut_score, score_wild_type=wt_score), name="Amplitude_mhcII_mb"),
            # dai multiple binding mhc II
            AnnotationFactory.build_annotation(
                value=self.differential_binding.dai(score_mutation=mut_score, score_wild_type=wt_score),
                name="DAI_mhcII_MB"),
            AnnotationFactory.build_annotation(value=amplitude_affinity, name="Amplitude_mhcII_affinity"),
            # amplitude rank score mhc II
            AnnotationFactory.build_annotation(value=self.calculate_amplitude_mhc(
                score_mutation=sc_mut, score_wild_type=sc_wt),
                name="Amplitude_mhcII_rank_netmhcpan4"),
            AnnotationFactory.build_annotation(
                value=self.wrap_pathogen_similarity(mutation=epitope_mut_mhcii),
                name="Pathogensimiliarity_mhcII"),
            AnnotationFactory.build_annotation(value=self.calculate_recognition_potential(
                amplitude=amplitude_affinity, pathogen_similarity=pathogen_similarity, mutation_in_anchor="0"),
                name="Recognition_Potential_mhcII_affinity"),
            AnnotationFactory.build_annotation(
                value=self.differential_binding.classify_adn_cdn(
                    score_mutation=sc_mut, amplitude=amplitude_affinity, bdg_cutoff_classical=bdg_cutoff_classical_mhcii,
                    bdg_cutoff_alternative=bdg_cutoff_alternative_mhcii, amplitude_cutoff=amplitude_cutoff_mhcii,
                    category="CDN"),
                name="CDN_mhcII"),
            AnnotationFactory.build_annotation(
                value=self.differential_binding.classify_adn_cdn(
                    score_mutation=sc_mut, amplitude=amplitude_affinity, bdg_cutoff_classical=bdg_cutoff_classical_mhcii,
                    bdg_cutoff_alternative=bdg_cutoff_alternative_mhcii, amplitude_cutoff=amplitude_cutoff_mhcii,
                    category="ADN"),
                name="ADN_mhcII")
            ]
