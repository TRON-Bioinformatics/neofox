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


class Amplitude:
    def __init__(self):
        self.amplitude_mhci_affinity_9mer = None
        self.amplitude_mhci_affinity = None
        self.amplitude_mhcii_rank = None

    def calculate_amplitude_mhc(
        self, score_mutation, score_wild_type, apply_correction=False
    ):
        """
        This function calculates the amplitude between mutated and wt epitope according to Balachandran et al.
        when affinity is used, use correction from Luksza et al. *1/(1+0.0003*aff_wt)
        """
        amplitude_mhc = None
        try:
            candidate_amplitude_mhc = score_wild_type / score_mutation
            if apply_correction:  # nine_mer or affinity:
                amplitude_mhc = candidate_amplitude_mhc * self._calculate_correction(
                    score_wild_type
                )
            else:
                amplitude_mhc = candidate_amplitude_mhc
        except (ZeroDivisionError, ValueError, TypeError):
            pass
        return amplitude_mhc

    def _calculate_correction(self, score_wild_type):
        return 1 / (1 + 0.0003 * float(score_wild_type))

    def run(
        self, netmhcpan: BestAndMultipleBinder, netmhc2pan: BestAndMultipleBinderMhcII
    ):
        # MHC I
        if netmhcpan:
            if netmhcpan.best_epitope_by_affinity.peptide and netmhcpan.best_wt_epitope_by_affinity.peptide:
                self.amplitude_mhci_affinity = self.calculate_amplitude_mhc(
                    score_mutation=netmhcpan.best_epitope_by_affinity.affinity_score,
                    score_wild_type=netmhcpan.best_wt_epitope_by_affinity.affinity_score,
                    apply_correction=True,
                )
            if netmhcpan.best_ninemer_epitope_by_affinity.peptide and netmhcpan.best_ninemer_wt_epitope_by_affinity.peptide:
                self.amplitude_mhci_affinity_9mer = self.calculate_amplitude_mhc(
                    score_mutation=netmhcpan.best_ninemer_epitope_by_affinity.affinity_score,
                    score_wild_type=netmhcpan.best_ninemer_wt_epitope_by_affinity.affinity_score,
                    apply_correction=True,
                )
        # MHC II
        if netmhc2pan:
            if netmhc2pan.best_predicted_epitope_rank.peptide and netmhc2pan.best_predicted_epitope_rank_wt.peptide:
                self.amplitude_mhcii_rank = self.calculate_amplitude_mhc(
                    score_mutation=netmhc2pan.best_predicted_epitope_rank.rank,
                    score_wild_type=netmhc2pan.best_predicted_epitope_rank_wt.rank,
                )

    def get_annotations(self) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.amplitude_mhci_affinity_9mer,
                name="Amplitude_MHCI_affinity_9mer",
            ),
            AnnotationFactory.build_annotation(
                value=self.amplitude_mhci_affinity, name="Amplitude_MHCI_affinity"
            ),
        ]

    def get_annotations_mhc2(self) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.amplitude_mhcii_rank, name="Amplitude_MHCII_rank"
            )
        ]
