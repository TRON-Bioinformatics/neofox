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
from typing import List, Set
import scipy.stats as stats

from neofox.MHC_predictors.netmhcpan.netmhcpan_prediction import NetMhcPanPredictor
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import (
    PredictedEpitope,
)
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.neoantigen import Annotation, Mhc1, Zygosity, Mutation
from neofox.model.wrappers import AnnotationFactory


class BestAndMultipleBinder:
    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self._initialise()

    def _initialise(self):
        self.phbr_i = None
        self.epitope_affinities = None
        self.generator_rate = None
        self.mutation_in_anchor_9mer = None

    def calculate_phbr_i(
        self, predictions: List[PredictedEpitope], mhc1_alleles: List[Mhc1]
    ):
        """returns list of multiple binding scores for mhcII considering best epitope per allele, applying different types of means (harmonic ==> PHRB-II, Marty et al).
        2 copies of DRA - DRB1 --> consider this gene 2x when averaging mhcii binding scores
        """
        best_epitopes_per_allele = (
            BestAndMultipleBinder.extract_best_epitope_per_alelle(
                predictions, mhc1_alleles
            )
        )
        phbr_i = None
        if len(best_epitopes_per_allele) == 6:
            phbr_i = stats.hmean(list(map(lambda e: e.rank, best_epitopes_per_allele)))
        return phbr_i

    @staticmethod
    def extract_best_epitope_per_alelle(
        epitopes: List[PredictedEpitope], mhc_isoforms: List[Mhc1]
    ) -> List[PredictedEpitope]:
        """
        This function returns the predicted epitope with the lowest binding score for each patient allele,
        considering homozyogosity
        """
        homozygous_alleles = BestAndMultipleBinder._get_homozygous_mhc1_alleles(
            mhc_isoforms
        )
        hetero_hemizygous_alleles = (
            BestAndMultipleBinder._get_heterozygous_or_hemizygous_mhc1_alleles(
                mhc_isoforms
            )
        )
        return BestAndMultipleBinder._get_sorted_epitopes(
            hetero_hemizygous_alleles, homozygous_alleles, epitopes
        )

    @staticmethod
    def _get_homozygous_mhc1_alleles(mhc_isoforms: List[Mhc1]) -> List[str]:
        """
        Returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles.
        Otherwise retunrs empty list
        """
        return [
            a.name
            for m in mhc_isoforms
            for a in m.alleles
            if m.zygosity == Zygosity.HOMOZYGOUS
        ]

    @staticmethod
    def _get_heterozygous_or_hemizygous_mhc1_alleles(
        mhc_isoforms: List[Mhc1],
    ) -> List[str]:
        """
        Returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles.
        Otherwise retunrs empty list
        """
        return [
            a.name
            for m in mhc_isoforms
            for a in m.alleles
            if m.zygosity in [Zygosity.HETEROZYGOUS, Zygosity.HEMIZYGOUS]
        ]

    @staticmethod
    def _get_sorted_epitopes(
        hetero_hemizygous_alleles,
        homozygous_alleles,
        predictions: List[PredictedEpitope],
    ) -> List[PredictedEpitope]:

        # groups epitopes by allele
        epitopes_by_allele = {}
        for p in predictions:
            epitopes_by_allele.setdefault(p.hla, []).append(p)
        # chooses the best epitope per allele while considering zygosity
        best_epis_per_allele = []
        for list_alleles in epitopes_by_allele.values():
            # sort by rank to choose the best epitope, ties are solved choosing the first peptide in alphabetcial order
            list_alleles.sort(key=lambda x: (x.rank, x.peptide))
            best_epitope = list_alleles[0]
            if best_epitope.hla in hetero_hemizygous_alleles:
                best_epis_per_allele.append(best_epitope)  # adds the epitope once
            if best_epitope.hla in homozygous_alleles:
                best_epis_per_allele.append(best_epitope)
                best_epis_per_allele.append(best_epitope)  # adds the epitope twice
        return best_epis_per_allele

    @staticmethod
    def determine_number_of_binders(predictions: List[PredictedEpitope], threshold=2):
        """
        Determines the number of HLA binders per mutation based on a threshold. Default is set to 2, which is threshold for weak binding using netmhcpan4.
        """
        scores = [epitope.affinity_score for epitope in predictions]
        number_binders = 0
        for score in scores:
            if score < threshold:
                number_binders += 1
        return number_binders if not len(scores) == 0 else None

    def run(
        self,
        mutation: Mutation,
        mhc1_alleles_patient: List[Mhc1],
        mhc1_alleles_available: Set,
    ):
        """
        predicts MHC epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed
        """
        self._initialise()
        netmhcpan = NetMhcPanPredictor(
            runner=self.runner, configuration=self.configuration
        )
        # print alleles
        predictions = netmhcpan.mhc_prediction(
            mhc1_alleles_patient, mhc1_alleles_available, mutation.mutated_xmer
        )
        filtered_predictions = netmhcpan.filter_binding_predictions(
            position_of_mutation=mutation.position, predictions=predictions
        )

        # multiple binding
        self.epitope_affinities = "/".join(
            [str(epitope.affinity_score) for epitope in filtered_predictions]
        )

        # best prediction
        self.best_epitope_by_rank = netmhcpan.select_best_by_rank(filtered_predictions)
        self.best_epitope_by_affinity = netmhcpan.select_best_by_affinity(
            filtered_predictions
        )

        # best predicted epitope of length 9
        ninemer_predictions = netmhcpan.filter_for_9mers(filtered_predictions)
        self.best_ninemer_epitope_by_rank = netmhcpan.select_best_by_rank(
            ninemer_predictions
        )
        self.best_ninemer_epitope_by_affinity = netmhcpan.select_best_by_affinity(
            ninemer_predictions
        )

        # multiple binding based on affinity
        self.generator_rate = self.determine_number_of_binders(
            predictions=filtered_predictions, threshold=50
        )

        # PHBR-I
        self.phbr_i = self.calculate_phbr_i(
            predictions=filtered_predictions, mhc1_alleles=mhc1_alleles_patient
        )

        # wt
        predictions_wt = netmhcpan.mhc_prediction(
            mhc1_alleles_patient, mhc1_alleles_available, mutation.wild_type_xmer
        )
        filtered_predictions_wt = netmhcpan.filter_binding_predictions(
            position_of_mutation=mutation.position, predictions=predictions_wt
        )

        # best prediction
        self.best_wt_epitope_by_rank = None
        if self.best_epitope_by_rank:
            self.best_wt_epitope_by_rank = netmhcpan.select_best_by_rank(
                netmhcpan.filter_wt_predictions_from_best_mutated(
                    filtered_predictions_wt, self.best_epitope_by_rank
                )
            )
        self.best_wt_epitope_by_affinity = None
        if self.best_epitope_by_affinity:
            self.best_wt_epitope_by_affinity = netmhcpan.select_best_by_affinity(
                netmhcpan.filter_wt_predictions_from_best_mutated(
                    filtered_predictions_wt, self.best_epitope_by_affinity
                )
            )

        # best predicted epitope of length 9
        ninemer_predictions_wt = netmhcpan.filter_for_9mers(filtered_predictions_wt)
        self.best_ninemer_wt_epitope_by_rank = None
        if self.best_ninemer_epitope_by_rank:
            self.best_ninemer_wt_epitope_by_rank = netmhcpan.select_best_by_rank(
                netmhcpan.filter_wt_predictions_from_best_mutated(
                    ninemer_predictions_wt, self.best_ninemer_epitope_by_rank
                )
            )
        self.best_ninemer_wt_epitope_by_affinity = None
        if self.best_ninemer_epitope_by_affinity:
            self.best_ninemer_wt_epitope_by_affinity = netmhcpan.select_best_by_affinity(
                netmhcpan.filter_wt_predictions_from_best_mutated(
                    ninemer_predictions_wt, self.best_ninemer_epitope_by_affinity
                )
            )

    def get_annotations(self) -> List[Annotation]:
        annotations = []
        if self.best_epitope_by_rank:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.rank, name="Best_rank_MHCI_score"
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.peptide,
                    name="Best_rank_MHCI_score_epitope",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.hla, name="Best_rank_MHCI_score_allele"
                )])
        if self.best_epitope_by_affinity:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.affinity_score,
                    name="Best_affinity_MHCI_score",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.peptide,
                    name="Best_affinity_MHCI_epitope",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.hla,
                    name="Best_affinity_MHCI_allele",
                )])
        if self.best_ninemer_epitope_by_rank:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.rank,
                    name="Best_rank_MHCI_9mer_score",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.peptide,
                    name="Best_rank_MHCI_9mer_epitope",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.hla,
                    name="Best_rank_MHCI_9mer_allele",
                )])
        if self.best_ninemer_epitope_by_affinity:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.affinity_score,
                    name="Best_affinity_MHCI_9mer_score",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.hla,
                    name="Best_affinity_MHCI_9mer_allele",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.peptide,
                    name="Best_affinity_MHCI_9mer_epitope",
                )])
            # wt
        if self.best_wt_epitope_by_affinity:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_wt_epitope_by_affinity.affinity_score,
                    name="Best_affinity_MHCI_score_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_wt_epitope_by_affinity.peptide,
                    name="Best_affinity_MHCI_epitope_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_wt_epitope_by_affinity.hla,
                    name="Best_affinity_MHCI_allele_WT",
                )])
        if self.best_wt_epitope_by_rank:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_wt_epitope_by_rank.rank, name="Best_rank_MHCI_score_WT"
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_wt_epitope_by_rank.peptide,
                    name="Best_rank_MHCI_score_epitope_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_wt_epitope_by_rank.hla,
                    name="Best_rank_MHCI_score_allele_WT",
                )])
        if self.best_ninemer_wt_epitope_by_rank:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_wt_epitope_by_rank.rank,
                    name="Best_rank_MHCI_9mer_score_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_wt_epitope_by_rank.peptide,
                    name="Best_rank_MHCI_9mer_epitope_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_wt_epitope_by_rank.hla,
                    name="Best_rank_MHCI_9mer_allele_WT",
                )])
        if self.best_ninemer_wt_epitope_by_affinity:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_wt_epitope_by_affinity.affinity_score,
                    name="Best_affinity_MHCI_9mer_score_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_wt_epitope_by_affinity.hla,
                    name="Best_affinity_MHCI_9mer_allele_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_wt_epitope_by_affinity.peptide,
                    name="Best_affinity_MHCI_9mer_epitope_WT",
                )])
        annotations.extend([
            # generator rate
            AnnotationFactory.build_annotation(
                value=self.generator_rate, name="Generator_rate"
            ),
            AnnotationFactory.build_annotation(value=self.phbr_i, name="PHBR-I")
        ])
        annotations.extend(self._get_positions_and_mutation_in_anchor())
        return annotations

    def _get_positions_and_mutation_in_anchor(self):
        """
        returns if mutation is in anchor position for best affinity epitope over all lengths and best 9mer affinity
        """
        annotations = []
        if self.best_ninemer_wt_epitope_by_affinity:
            position_9mer = EpitopeHelper.position_of_mutation_epitope(
                wild_type=self.best_ninemer_wt_epitope_by_affinity.peptide,
                mutation=self.best_ninemer_epitope_by_affinity.peptide,
            )
            mutation_in_anchor_9mer = EpitopeHelper.position_in_anchor_position(
                position_mhci=position_9mer,
                peptide_length=len(self.best_ninemer_epitope_by_affinity.peptide),
            )
            annotations = [
                AnnotationFactory.build_annotation(
                    value=position_9mer, name="Best_affinity_MHCI_9mer_position_mutation"
                ),
                AnnotationFactory.build_annotation(
                    value=mutation_in_anchor_9mer,
                    name="Best_affinity_MHCI_9mer_anchor_mutated",
                ),
            ]
        return annotations
