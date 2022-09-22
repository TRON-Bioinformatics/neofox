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

from neofox.model.neoantigen import Annotation, PredictedEpitope
from neofox.model.factories import AnnotationFactory
from neofox.published_features.differential_binding.amplitude import Amplitude


class DifferentialBinding:

    def dai(self, score_mutation, score_wild_type):
        """
        Calculates DAI: Returns difference between wt and mut MHC binding score.
        """
        score = None
        try:
            score = score_wild_type - score_mutation
        except TypeError:
            pass
        return score

    def classify_adn_cdn(
        self,
        score_mutation,
        amplitude,
        bdg_cutoff_classical,
        bdg_cutoff_alternative,
        amplitude_cutoff,
        category,
    ):
        """
        returns if an epitope belongs to classically and alternatively defined neoepitopes (CDN vs ADN)
        (indicate which category to examine by category)--> Rech et al, 2018
        grouping is based on affinity and affinitiy foldchange between wt and mut
        """
        group = None
        try:
            if category == "CDN":
                group = score_mutation < bdg_cutoff_classical
            elif category == "ADN":
                group = (
                    score_mutation < bdg_cutoff_alternative
                    and amplitude > amplitude_cutoff
                )
        except (ValueError, TypeError):
            pass
        return group

    def get_annotations_dai(self, epitope: PredictedEpitope) -> List[Annotation]:
        dai = None
        if epitope.mutated_peptide and epitope.wild_type_peptide:
            dai = self.dai(
                        score_mutation=epitope.affinity_mutated,
                        score_wild_type=epitope.affinity_wild_type
                    )
        annotations = [
            AnnotationFactory.build_annotation(
                name="DAI_MHCI_bestAffinity",
                value=dai
            ),
        ]
        return annotations

    def get_annotations(
        self, mutated_peptide_mhci: PredictedEpitope, amplitude: Amplitude
    ) -> List[Annotation]:

        bdg_cutoff_classical_mhci = 50
        bdg_cutoff_alternative_mhci = 5000
        amplitude_cutoff_mhci = 10

        cdn = None
        adn = None
        if mutated_peptide_mhci.mutated_peptide:
            cdn = self.classify_adn_cdn(
                        score_mutation=mutated_peptide_mhci.affinity_mutated,
                        amplitude=amplitude.amplitude_mhci_affinity,
                        bdg_cutoff_classical=bdg_cutoff_classical_mhci,
                        bdg_cutoff_alternative=bdg_cutoff_alternative_mhci,
                        amplitude_cutoff=amplitude_cutoff_mhci,
                        category="CDN",
                    )
            adn = self.classify_adn_cdn(
                        score_mutation=mutated_peptide_mhci.affinity_mutated,
                        amplitude=amplitude.amplitude_mhci_affinity,
                        bdg_cutoff_classical=bdg_cutoff_classical_mhci,
                        bdg_cutoff_alternative=bdg_cutoff_alternative_mhci,
                        amplitude_cutoff=amplitude_cutoff_mhci,
                        category="ADN",
                    )
        annotations = [
            AnnotationFactory.build_annotation(
                name="Classically_defined_neopeptide_MHCI",
                value=cdn
            ),
            AnnotationFactory.build_annotation(
                name="Alternatively_defined_neopeptide_MHCI",
                value=adn
            ),
        ]
        return annotations

    def get_annotations_mhc2(
        self, mutated_peptide_mhcii: PredictedEpitope, amplitude: Amplitude
    ) -> List[Annotation]:

        bdg_cutoff_classical_mhcii = 1
        bdg_cutoff_alternative_mhcii = 4
        amplitude_cutoff_mhcii = 4
        cdn = None
        adn = None
        if mutated_peptide_mhcii.mutated_peptide:
            cdn = self.classify_adn_cdn(
                        score_mutation=mutated_peptide_mhcii.rank_mutated,
                        amplitude=amplitude.amplitude_mhcii_rank,
                        bdg_cutoff_classical=bdg_cutoff_classical_mhcii,
                        bdg_cutoff_alternative=bdg_cutoff_alternative_mhcii,
                        amplitude_cutoff=amplitude_cutoff_mhcii,
                        category="CDN",
                    )
            adn = self.classify_adn_cdn(
                        score_mutation=mutated_peptide_mhcii.rank_mutated,
                        amplitude=amplitude.amplitude_mhcii_rank,
                        bdg_cutoff_classical=bdg_cutoff_classical_mhcii,
                        bdg_cutoff_alternative=bdg_cutoff_alternative_mhcii,
                        amplitude_cutoff=amplitude_cutoff_mhcii,
                        category="ADN",
                    )
        annotations = [
            AnnotationFactory.build_annotation(
                value=cdn,
                name="Classically_defined_neopeptide_MHCII",
            ),
            AnnotationFactory.build_annotation(
                value=adn,
                name="Alternatively_defined_neopeptide_MHCII",
            ),
        ]
        return annotations

    def get_annotations_epitope_mhci(self, epitope: PredictedEpitope) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.dai(score_mutation=epitope.affinity_mutated, score_wild_type=epitope.affinity_wild_type),
                name='DAI')
            ]
