#!/usr/bin/env python

import os
import os.path
from typing import List

from neofox.helpers import intermediate_files
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.predictors.netmhcpan4.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.predictors.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder


class DissimilarityCalculator(BlastpRunner):

    def __init__(self, runner, configuration, proteome_db):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        super().__init__(runner=runner, configuration=configuration)
        self.proteome_db = proteome_db

    def _dissimilarity(self, fasta_file):
        """
        This function determines the dissimilarity to self-proteome of epitopes as described in Richman et al
        """
        outfile = self.run_blastp(
            fasta_file=fasta_file, database=os.path.join(self.proteome_db, "homo_sapiens.mod"))
        similarity = self.parse_blastp_output(blastp_output_file=outfile, a=32)
        dissimilarity = 1 - similarity
        os.remove(outfile)
        return dissimilarity

    def calculate_dissimilarity(self, mhc_mutation, mhc_affinity, filter_binder=False):
        """
        wrapper for dissimilarity calculation
        """
        dissimilarity_score = None
        if mhc_mutation != "-" and (not filter_binder or not mhc_affinity >= 500):
            fastafile = intermediate_files.create_temp_fasta(sequences=[mhc_mutation], prefix="tmp_dissimilarity_",
                                                             comment_prefix='M_')
            dissimilarity_score = self._dissimilarity(fastafile)
            os.remove(fastafile)
        return dissimilarity_score

    def get_annotations(
            self, netmhcpan: BestAndMultipleBinder, netmhcpan2: BestAndMultipleBinderMhcII) -> List[Annotation]:
        """
        returns dissimilarity for MHC I (affinity) MHC II (affinity)
        """
        return [
            AnnotationFactory.build_annotation(
                value=self.calculate_dissimilarity(
                    mhc_mutation=netmhcpan.best4_affinity_epitope, mhc_affinity=netmhcpan.best4_affinity),
                name="Dissimilarity_MHCI"),
            AnnotationFactory.build_annotation(
                value=self.calculate_dissimilarity(
                    mhc_mutation=netmhcpan.best4_affinity_epitope, mhc_affinity=netmhcpan.best4_affinity,
                    filter_binder=True),
                name="Dissimilarity_MHCI_cutoff500nM"),
            AnnotationFactory.build_annotation(
                value=self.calculate_dissimilarity(
                    mhc_mutation=netmhcpan2.best_mhcII_pan_affinity_epitope,
                    mhc_affinity=netmhcpan2.best_mhcII_pan_affinity),
                name="Dissimilarity_MHCII")
            ]
