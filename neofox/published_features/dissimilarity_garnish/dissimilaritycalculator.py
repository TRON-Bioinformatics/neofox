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
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.model.neoantigen import Annotation, PredictedEpitope
from neofox.model.factories import AnnotationFactory


class DissimilarityCalculator:

    def __init__(self, proteome_blastp_runner: BlastpRunner):
        self.proteome_blastp_runner = proteome_blastp_runner

    def calculate_dissimilarity(self, epitope: PredictedEpitope):
        """
        wrapper for dissimilarity calculation
        """
        dissimilarity = None
        if epitope.mutated_peptide != "-":
            similarity = self.proteome_blastp_runner.calculate_similarity_database(
                peptide=epitope.mutated_peptide, a=32)
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
        if mutated_peptide_mhci and mutated_peptide_mhci.mutated_peptide:
            dissimilarity_mhci = self.calculate_dissimilarity(epitope=mutated_peptide_mhci)
        if mutated_peptide_mhcii and mutated_peptide_mhcii.mutated_peptide:
            dissimilarity_mhcii = self.calculate_dissimilarity(epitope=mutated_peptide_mhcii)
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

    def get_annotations_epitope(self, epitope: PredictedEpitope) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.calculate_dissimilarity(epitope=epitope),
                name='dissimilarity_score')
        ]
