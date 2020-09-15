#!/usr/bin/env python
#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#

import os
import os.path
from typing import List

from neofox.helpers import intermediate_files
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.netmhcpan.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder


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
            fasta_file=fasta_file, database=os.path.join(self.proteome_db, "homo_sapiens"))
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
            self, netmhcpan: BestAndMultipleBinder) -> List[Annotation]:
        """
        returns dissimilarity for MHC I (affinity) MHC II (affinity)
        """
        return [
            AnnotationFactory.build_annotation(
                value=self.calculate_dissimilarity(
                    mhc_mutation=netmhcpan.best4_affinity_epitope, mhc_affinity=netmhcpan.best4_affinity,
                    filter_binder=True),
                name="Dissimilarity_MHCI_cutoff500nM"),
            ]
