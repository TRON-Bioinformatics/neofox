#!/usr/bin/env python

import os
import os.path
from typing import List

from logzero import logger

from neofox.helpers import intermediate_files
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.predictors.netmhcpan4.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.predictors.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.literature_features.differential_binding.amplitude import Amplitude


class NeoantigenFitnessCalculator(BlastpRunner):

    def __init__(self, runner, configuration, iedb):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        super().__init__(runner, configuration)
        self.iedb = iedb


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
        amplitude_mhc = None
        try:
            candidate_amplitude_mhc = score_wild_type / score_mutation
            if apply_correction:  # nine_mer or affinity:
                amplitude_mhc = candidate_amplitude_mhc * self._calculate_correction(score_wild_type)
            else:
                amplitude_mhc = candidate_amplitude_mhc
        except(ZeroDivisionError, ValueError, TypeError):
            pass
        return amplitude_mhc

    def _calculate_correction(self, score_wild_type):
        return 1 / (1 + 0.0003 * float(score_wild_type))

    def calculate_recognition_potential(
            self, amplitude: float, pathogen_similarity: float, mutation_in_anchor: bool,
            mhc_affinity_mut: float = None):
        """
        This function calculates the recognition potential, defined by the product of amplitude and pathogensimiliarity
        of an epitope according to Balachandran et al.
        F_alpha = - max (A_i x R_i)

        Returns (A_i x R_i) value only for nonanchor mutation and epitopes of length 9; only considered by Balachandran
        """
        recognition_potential = None
        try:
            candidate_recognition_potential = amplitude * pathogen_similarity
            if mhc_affinity_mut:
                if not mutation_in_anchor and mhc_affinity_mut < 500.0:
                    recognition_potential = candidate_recognition_potential
            else:
                if not mutation_in_anchor:
                    recognition_potential = candidate_recognition_potential
        except (ValueError, TypeError):
            pass
        return recognition_potential

    def get_annotations(self, netmhcpan: BestAndMultipleBinder, amplitude: Amplitude) -> List[Annotation]:
        pathogen_similarity_9mer = self.wrap_pathogen_similarity(mutation=netmhcpan.mhcI_affinity_epitope_9mer)
        pathogen_similarity_rank = self.wrap_pathogen_similarity(mutation=netmhcpan.best4_mhc_epitope)
        pathogen_similarity_affinity = self.wrap_pathogen_similarity(mutation=netmhcpan.best4_affinity_epitope)

        position = EpitopeHelper.position_of_mutation_epitope(wild_type=netmhcpan.best4_affinity_epitope_WT,
                                                              mutation=netmhcpan.best4_affinity_epitope)
        position_9mer = EpitopeHelper.position_of_mutation_epitope(
            wild_type=netmhcpan.mhcI_affinity_epitope_9mer_WT, mutation=netmhcpan.mhcI_affinity_epitope_9mer)
        position_rank = EpitopeHelper.position_of_mutation_epitope(
            wild_type=netmhcpan.best4_mhc_score_WT, mutation=netmhcpan.best4_mhc_score)

        return [
            AnnotationFactory.build_annotation(name="Pathogensimiliarity_MHCI_affinity_9mer",
                                               value=pathogen_similarity_9mer),
            AnnotationFactory.build_annotation(name="Pathogensimiliarity_MHCI_rank_nmers",
                                               value=pathogen_similarity_rank),
            AnnotationFactory.build_annotation(name="Pathogensimiliarity_MHCI_affinity_nmers",
                                               value=pathogen_similarity_affinity),
            AnnotationFactory.build_annotation(name="Recognition_Potential_MHCI_affinity_9mer",
                                               value=self.calculate_recognition_potential(
                                                   amplitude=amplitude.amplitude_mhci_affinity_9mer,
                                                   pathogen_similarity=pathogen_similarity_9mer,
                                                   mutation_in_anchor=EpitopeHelper.position_in_anchor_position(
                                                       position_mhci=position_9mer,
                                                       peptide_length=len(netmhcpan.mhcI_affinity_epitope_9mer)),
                                                   mhc_affinity_mut=netmhcpan.mhcI_affinity_9mer)),
            AnnotationFactory.build_annotation(name="Recognition_Potential_MHCI_affinity_nmers",
                                               value=self.calculate_recognition_potential(
                                                   amplitude=amplitude.amplitude_mhci_affinity,
                                                   pathogen_similarity=pathogen_similarity_affinity,
                                                   mutation_in_anchor=EpitopeHelper.position_in_anchor_position(
                                                       position_mhci=position,
                                                       peptide_length=len(netmhcpan.best4_affinity_epitope)))),
            AnnotationFactory.build_annotation(name="Recognition_Potential_MHCI_rank4",
                                               value=self.calculate_recognition_potential(
                                                   amplitude=amplitude.amplitude_mhci_rank,
                                                   pathogen_similarity=pathogen_similarity_rank,
                                                   mutation_in_anchor=EpitopeHelper.position_in_anchor_position(
                                                       position_mhci=position_rank,
                                                       peptide_length=len(netmhcpan.best4_mhc_epitope))))
        ]

    def get_annotations_mhc2(self, netmhc2pan: BestAndMultipleBinderMhcII, amplitude: Amplitude) -> List[Annotation]:

        pathogen_similarity = self.wrap_pathogen_similarity(mutation=netmhc2pan.best_mhcII_pan_affinity_epitope)
        return [
             AnnotationFactory.build_annotation(
                value=self.wrap_pathogen_similarity(mutation=netmhc2pan.best_mhcII_pan_affinity_epitope),
                name="Pathogensimiliarity_MHCII_affinity"),
            AnnotationFactory.build_annotation(value=self.calculate_recognition_potential(
                amplitude=amplitude.amplitude_mhcii_affinity, pathogen_similarity=pathogen_similarity,
                mutation_in_anchor=False),
                name="Recognition_Potential_MHCII_affinity"),
            ]
