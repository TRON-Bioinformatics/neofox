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

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.netmhcpan.combine_netmhcIIpan_pred_multiple_binders import (
    BestAndMultipleBinderMhcII,
)
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)
from neofox import AFFINITY_THRESHOLD_DEFAULT
from neofox.published_features.differential_binding.amplitude import Amplitude


class DifferentialBinding:

    def __init__(self, affinity_threshold=AFFINITY_THRESHOLD_DEFAULT):
        self.affinity_threshold = affinity_threshold

    def dai(self, score_mutation, score_wild_type, affin_filtering=False):
        """
        Calculates DAI: Returns difference between wt and mut MHC binding score.
        """
        score = None
        try:
            if affin_filtering:
                if score_mutation < self.affinity_threshold:
                    score = score_wild_type - score_mutation
            else:
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

    def get_annotations_dai(self, netmhcpan: BestAndMultipleBinder) -> List[Annotation]:
        dai = None
        if netmhcpan.best_epitope_by_affinity.peptide and netmhcpan.best_wt_epitope_by_affinity.peptide:
            dai = self.dai(
                        score_mutation=netmhcpan.best_epitope_by_affinity.affinity_score,
                        score_wild_type=netmhcpan.best_wt_epitope_by_affinity.affinity_score,
                        affin_filtering=True,
                    )
        annotations = [
            AnnotationFactory.build_annotation(
                name="DAI_MHCI_affinity",
                value=dai
            ),
        ]
        return annotations

    def get_annotations(
        self, netmhcpan: BestAndMultipleBinder, amplitude: Amplitude
    ) -> List[Annotation]:

        bdg_cutoff_classical_mhci = 50
        bdg_cutoff_alternative_mhci = 5000
        amplitude_cutoff_mhci = 10

        cdn = None
        adn = None
        if netmhcpan.best_epitope_by_affinity.peptide:
            cdn = self.classify_adn_cdn(
                        score_mutation=netmhcpan.best_epitope_by_affinity.affinity_score,
                        amplitude=amplitude.amplitude_mhci_affinity,
                        bdg_cutoff_classical=bdg_cutoff_classical_mhci,
                        bdg_cutoff_alternative=bdg_cutoff_alternative_mhci,
                        amplitude_cutoff=amplitude_cutoff_mhci,
                        category="CDN",
                    )
            adn = self.classify_adn_cdn(
                        score_mutation=netmhcpan.best_epitope_by_affinity.affinity_score,
                        amplitude=amplitude.amplitude_mhci_affinity,
                        bdg_cutoff_classical=bdg_cutoff_classical_mhci,
                        bdg_cutoff_alternative=bdg_cutoff_alternative_mhci,
                        amplitude_cutoff=amplitude_cutoff_mhci,
                        category="ADN",
                    )
        annotations = [
            AnnotationFactory.build_annotation(
                name="CDN_MHCI",
                value=cdn
            ),
            AnnotationFactory.build_annotation(
                name="ADN_MHCI",
                value=adn
            ),
        ]
        return annotations

    def get_annotations_mhc2(
        self, netmhc2pan: BestAndMultipleBinderMhcII, amplitude: Amplitude
    ) -> List[Annotation]:

        bdg_cutoff_classical_mhcii = 1
        bdg_cutoff_alternative_mhcii = 4
        amplitude_cutoff_mhcii = 4
        cdn = None
        adn = None
        if netmhc2pan.best_predicted_epitope_rank.peptide:
            cdn = self.classify_adn_cdn(
                        score_mutation=netmhc2pan.best_predicted_epitope_rank.rank,
                        amplitude=amplitude.amplitude_mhcii_rank,
                        bdg_cutoff_classical=bdg_cutoff_classical_mhcii,
                        bdg_cutoff_alternative=bdg_cutoff_alternative_mhcii,
                        amplitude_cutoff=amplitude_cutoff_mhcii,
                        category="CDN",
                    )
            adn = self.classify_adn_cdn(
                        score_mutation=netmhc2pan.best_predicted_epitope_rank.rank,
                        amplitude=amplitude.amplitude_mhcii_rank,
                        bdg_cutoff_classical=bdg_cutoff_classical_mhcii,
                        bdg_cutoff_alternative=bdg_cutoff_alternative_mhcii,
                        amplitude_cutoff=amplitude_cutoff_mhcii,
                        category="ADN",
                    )
        annotations = [
            AnnotationFactory.build_annotation(
                value=cdn,
                name="CDN_MHCII",
            ),
            AnnotationFactory.build_annotation(
                value=adn,
                name="ADN_MHCII",
            ),
        ]
        return annotations
