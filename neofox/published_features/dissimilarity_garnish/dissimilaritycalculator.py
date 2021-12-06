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
from typing import List
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import PredictedEpitope
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.model.neoantigen import Annotation
from neofox.model.factories import AnnotationFactory
from neofox import AFFINITY_THRESHOLD_DEFAULT


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

    def get_annotations(
            self, mutated_peptide_mhci: PredictedEpitope, mutated_peptide_mhcii: PredictedEpitope) -> List[Annotation]:
        """
        returns dissimilarity for MHC I (affinity) MHC II (affinity)
        """
        dissimilarity_mhci = None
        dissimilarity_mhcii = None
        if mutated_peptide_mhci and mutated_peptide_mhci.peptide:
            dissimilarity_mhci = self.calculate_dissimilarity(
                mutated_peptide=mutated_peptide_mhci.peptide,
                mhc_affinity=mutated_peptide_mhci.affinity_score )
        if mutated_peptide_mhcii and mutated_peptide_mhcii.peptide:
            dissimilarity_mhcii = self.calculate_dissimilarity(
                mutated_peptide=mutated_peptide_mhcii.peptide,
                mhc_affinity=mutated_peptide_mhcii.affinity_score )
        annotations = [
            AnnotationFactory.build_annotation(
                value=dissimilarity_mhci,
                name="Dissimilarity_MHCI",
            ),
            AnnotationFactory.build_annotation(
                value=dissimilarity_mhcii,
                name="Dissimilarity_MHCII",
            )
        ]
        return annotations
