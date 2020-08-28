#!/usr/bin/env python

import os
import os.path
from typing import List

from neofox.helpers import intermediate_files
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


class DissimilarityCalculator(BlastpRunner):

    def __init__(self, runner, configuration, proteome_db):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        super().__init__(runner=runner, configuration=configuration)
        self.proteome_db = proteome_db

    def _calc_dissimilarity(self, fasta_file):
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
            fastafile = intermediate_files.create_temp_fasta(sequences=[mhc_mutation], prefix="tmp_dissimilarity_", comment_prefix='M_')
            dissimilarity_score = self._calc_dissimilarity(fastafile)
            os.remove(fastafile)
        return dissimilarity_score

    def get_annotations(self, epitope_mhci, affinity_mhci, epitope_mhcii, affinity_mhcii) -> List[Annotation]:
        """
        returns dissimilarity for MHC I (affinity) MHC II (affinity)
        """
        return [
            AnnotationFactory.build_annotation(
                value=self.calculate_dissimilarity(mhc_mutation=epitope_mhci, mhc_affinity=affinity_mhci),
                name="dissimilarity"),
            AnnotationFactory.build_annotation(
                value=self.calculate_dissimilarity(
                    mhc_mutation=epitope_mhci, mhc_affinity=affinity_mhci, filter_binder=True),
                name="dissimilarity_filter500"),
            AnnotationFactory.build_annotation(
                value=self.calculate_dissimilarity(mhc_mutation=epitope_mhcii, mhc_affinity=affinity_mhcii),
                name="dissimilarity_mhcII")
            ]
