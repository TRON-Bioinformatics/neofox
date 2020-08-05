#!/usr/bin/env python

import os
import os.path

from input.helpers import intermediate_files
from input.helpers.blastp_runner import BlastpRunner


class DissimilarityCalculator(BlastpRunner):

    def __init__(self, runner, configuration, proteome_db):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
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
        fastafile = intermediate_files.create_temp_fasta(sequences=[mhc_mutation], prefix="tmp_dissimilarity_", comment_prefix='M_')
        dissim = self._calc_dissimilarity(fastafile)
        os.remove(fastafile)
        sc = dissim
        if filter_binder and float(mhc_affinity) >= 500:
            sc = "NA"
        return sc
