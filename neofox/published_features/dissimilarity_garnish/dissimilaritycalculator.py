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

from neofox.helpers.blastp_runner import BlastpRunner
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)
from neofox import AFFINITY_THRESHOLD_DEFAULT
from neofox.MHC_predictors.netmhcpan.combine_netmhcIIpan_pred_multiple_binders import (
    BestAndMultipleBinderMhcII
)


class DissimilarityCalculator:

    def __init__(self, proteome_blastp_runner: BlastpRunner, affinity_threshold=AFFINITY_THRESHOLD_DEFAULT):
        self.affinity_threshold = affinity_threshold
        self.proteome_blastp_runner = proteome_blastp_runner

    def calculate_dissimilarity(self, mutated_peptide, mhc_affinity):
        """
        wrapper for dissimilarity calculation
        """
        dissimilarity = None
        if mutated_peptide != "-" and not mhc_affinity >= self.affinity_threshold:
            similarity = self.proteome_blastp_runner.calculate_similarity_database(
                peptide=mutated_peptide,
                a=32,
            )
            if similarity is not None:
                dissimilarity = 1 - similarity
        return dissimilarity

    def get_annotations(self, netmhcpan: BestAndMultipleBinder, netmhc2pan: BestAndMultipleBinderMhcII) -> List[Annotation]:
        """
        returns dissimilarity for MHC I (affinity) MHC II (affinity)
        """
        dissimilarity_mhci = None
        dissimilarity_mhcii = None
        if netmhcpan.best_epitope_by_affinity.peptide:
            dissimilarity_mhci = self.calculate_dissimilarity(
                mutated_peptide=netmhcpan.best_epitope_by_affinity.peptide,
                mhc_affinity=netmhcpan.best_epitope_by_affinity.affinity_score )
        if netmhc2pan.best_predicted_epitope_affinity.peptide:
            dissimilarity_mhcii = self.calculate_dissimilarity(
                mutated_peptide=netmhc2pan.best_predicted_epitope_affinity.peptide,
                mhc_affinity=netmhc2pan.best_predicted_epitope_affinity.affinity_score )
        annotations = [
            AnnotationFactory.build_annotation(
                value=dissimilarity_mhci,
                name="Dissimilarity_MHCII",
            ),
            AnnotationFactory.build_annotation(
                value=dissimilarity_mhcii,
                name="Dissimilarity_MHCII",
            )
        ]
        return annotations
