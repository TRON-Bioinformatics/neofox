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
from logzero import logger
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import (
    PredictedEpitope,
)
from neofox.MHC_predictors.netmhcpan.netmhcIIpan_prediction import NetMhcIIPanPredictor
from neofox.model.neoantigen import Annotation, Mhc2, Zygosity, Mhc2Isoform, Mutation
from neofox.model.wrappers import AnnotationFactory


class BestAndMultipleBinderMhcII:
    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self._initialise()

    def _initialise(self):
        self.phbr_ii = None
        self.best_predicted_epitope_rank = PredictedEpitope(
            peptide="-",
            pos=None,
            hla=Mhc2Isoform(name=None),
            affinity_score=None,
            rank=None,
        )
        self.best_predicted_epitope_affinity = PredictedEpitope(
            peptide="-",
            pos=None,
            hla=Mhc2Isoform(name=None),
            affinity_score=None,
            rank=None,
        )
        self.best_predicted_epitope_rank_wt = PredictedEpitope(
            peptide="-",
            pos=None,
            hla=Mhc2Isoform(name=None),
            affinity_score=None,
            rank=None,
        )
        self.best_predicted_epitope_affinity_wt = PredictedEpitope(
            peptide="-",
            pos=None,
            hla=Mhc2Isoform(name=None),
            affinity_score=None,
            rank=None,
        )

    def calculate_phbr_ii(self, best_epitope_per_allele_mhc2: List[PredictedEpitope]):
        """
        harmonic mean of best MHC II binding scores per MHC II allele
        :param best_epitope_per_allele_mhc2: list of best MHC II epitopes per allele
        :return: PHBR-II score, Marty et al
        """
        best_epitope_per_allele_mhc2_new = list(best_epitope_per_allele_mhc2)
        phbr_ii = None
        for allele_with_score in best_epitope_per_allele_mhc2:
            # add DRB1
            if "DRB1" in allele_with_score.hla:
                best_epitope_per_allele_mhc2_new.append(allele_with_score)
        if len(best_epitope_per_allele_mhc2_new) == 12:
            # 12 genes gene copies should be included into PHBR_II
            best_mhc_ii_scores_per_allele = [
                epitope.rank for epitope in best_epitope_per_allele_mhc2_new
            ]
            phbr_ii = stats.hmean(best_mhc_ii_scores_per_allele)
        return phbr_ii

    def run(
        self,
        mutation: Mutation,
        mhc2_alleles_patient: List[Mhc2],
        mhc2_alleles_available: Set,
    ):
        """predicts MHC II epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed"""
        # mutation
        self._initialise()
        netmhc2pan = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration
        )
        allele_combinations = netmhc2pan.generate_mhc2_alelle_combinations(
            mhc2_alleles_patient
        )
        # TODO: migrate the available alleles into the model for alleles
        patient_mhc2_isoforms = self._get_only_available_combinations(
            allele_combinations, mhc2_alleles_available
        )

        predictions = netmhc2pan.mhcII_prediction(
            patient_mhc2_isoforms, mutation.mutated_xmer
        )
        if len(mutation.mutated_xmer) >= 15:
            filtered_predictions = netmhc2pan.filter_binding_predictions(
                mutation.position, predictions
            )
            # multiple binding
            best_predicted_epitopes_per_alelle = (
                self.extract_best_epitope_per_mhc2_alelle(
                    filtered_predictions, mhc2_alleles_patient
                )
            )
            self.phbr_ii = self.calculate_phbr_ii(best_predicted_epitopes_per_alelle)

            # best prediction
            self.best_predicted_epitope_rank = netmhc2pan.select_best_by_rank(
                filtered_predictions
            )
            self.best_predicted_epitope_affinity = netmhc2pan.select_best_by_affinity(
                filtered_predictions
            )

        # wt
        predictions = netmhc2pan.mhcII_prediction(
            patient_mhc2_isoforms, mutation.wild_type_xmer
        )
        if len(mutation.wild_type_xmer) >= 15:
            filtered_predictions_wt = netmhc2pan.filter_binding_predictions(
                mutation.position, predictions
            )

            # best prediction
            self.best_predicted_epitope_rank_wt = None
            if self.best_predicted_epitope_rank:
                self.best_predicted_epitope_rank_wt = netmhc2pan.select_best_by_rank(
                    netmhc2pan.filter_wt_predictions_from_best_mutated(
                        filtered_predictions_wt, self.best_predicted_epitope_rank
                    )
                )
            self.best_predicted_epitope_affinity_wt = None
            if self.best_predicted_epitope_affinity:
                self.best_predicted_epitope_affinity_wt = (
                    netmhc2pan.select_best_by_affinity(
                        netmhc2pan.filter_wt_predictions_from_best_mutated(
                            filtered_predictions_wt, self.best_predicted_epitope_affinity
                        )
                    )
                )

    @staticmethod
    def _get_only_available_combinations(allele_combinations, set_available_mhc):
        patients_available_alleles = list(
            set(allele_combinations).intersection(set(set_available_mhc))
        )
        patients_not_available_alleles = list(
            set(allele_combinations).difference(set(set_available_mhc))
        )
        if len(patients_not_available_alleles) > 0:
            logger.warning(
                "MHC II alleles {} are not supported by NetMHC2pan and no binding or derived features will "
                "include it".format(",".join(patients_not_available_alleles))
            )
        return patients_available_alleles

    def get_annotations(self) -> List[Annotation]:
        annotations = []
        if self.best_predicted_epitope_rank:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_rank.rank,
                    name="Best_rank_MHCII_score",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_rank.peptide,
                    name="Best_rank_MHCII_score_epitope",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_rank.hla,
                    name="Best_rank_MHCII_score_allele",
                )])
        if self.best_predicted_epitope_affinity:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity.affinity_score,
                    name="Best_affinity_MHCII_score",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity.peptide,
                    name="Best_affinity_MHCII_epitope",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity.hla,
                    name="Best_affinity_MHCII_allele",
                )])
        if self.best_predicted_epitope_rank_wt:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_rank_wt.rank,
                    name="Best_rank_MHCII_score_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_rank_wt.peptide,
                    name="Best_rank_MHCII_score_epitope_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_rank_wt.hla,
                    name="Best_rank_MHCII_score_allele_WT",
                )])
        if self.best_predicted_epitope_affinity_wt:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity_wt.affinity_score,
                    name="Best_affinity_MHCII_score_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity_wt.peptide,
                    name="Best_affinity_MHCII_epitope_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity_wt.hla,
                    name="Best_affinity_MHCII_allele_WT",
                )])
        annotations.append(AnnotationFactory.build_annotation(value=self.phbr_ii, name="PHBR-II"))
        return annotations

    @staticmethod
    def _get_homozygous_mhc2_alleles(mhc_isoforms: List[Mhc2]) -> List[Mhc2]:
        """
        Returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles.
        Otherwise returns empty list
        """
        return [
            m
            for m in mhc_isoforms
            for g in m.genes
            if g.zygosity == Zygosity.HOMOZYGOUS
        ]

    @staticmethod
    def _get_heterozygous_or_hemizygous_mhc2_alleles(
        mhc_isoforms: List[Mhc2],
    ) -> List[Mhc2]:
        """
        Returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles.
        Otherwise returns empty list
        """
        return [
            m
            for m in mhc_isoforms
            for g in m.genes
            if g.zygosity in [Zygosity.HETEROZYGOUS, Zygosity.HEMIZYGOUS]
        ]

    @staticmethod
    def _get_sorted_epitopes_mhc2(
        hetero_hemizygous_alleles,
        homozygous_alleles,
        predictions: List[PredictedEpitope],
    ) -> List[PredictedEpitope]:

        # groups epitopes by allele
        epitopes_by_allele = {}
        for p in predictions:
            allele = p.hla
            epitopes_by_allele.setdefault(allele, []).append(p)

        # chooses the best epitope per allele and considers zygosity
        best_epitopes_per_allele = []
        for allele, epitopes in epitopes_by_allele.items():
            # sort by rank to choose the best epitope, fixes ties with peptide by alphabetical order
            epitopes.sort(key=lambda e: (e.rank, e.peptide))
            best_epitope = epitopes[0]
            num_repetitions = 0
            if (
                best_epitope.hla in hetero_hemizygous_alleles
                or best_epitope.hla in hetero_hemizygous_alleles
            ):
                # adds the epitope once if alleles heterozygous
                num_repetitions = 1
            if (
                best_epitope.hla in homozygous_alleles
            ):
                # adds the epitope twice if one allele is homozygous
                num_repetitions = 2
            best_epitopes_per_allele.extend(
                [best_epitope for _ in range(num_repetitions)]
            )
        return best_epitopes_per_allele

    def extract_best_epitope_per_mhc2_alelle(
            self, predictions: List[PredictedEpitope], mhc_isoforms: List[Mhc2]
    ) -> List[PredictedEpitope]:
        """
        This function returns the predicted epitope with the lowest binding score for each patient allele,
        considering homozygosity
        """
        netmhc2pan = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration
        )
        homozygous_alleles = BestAndMultipleBinderMhcII._get_homozygous_mhc2_alleles(
            mhc_isoforms
        )
        homozygous_alleles = netmhc2pan.generate_mhc2_alelle_combinations(homozygous_alleles)
        homozygous_alleles = list(set(homozygous_alleles))
        hetero_hemizygous_alleles = (
            BestAndMultipleBinderMhcII._get_heterozygous_or_hemizygous_mhc2_alleles(
                mhc_isoforms
            )
        )
        hetero_hemizygous_alleles = netmhc2pan.generate_mhc2_alelle_combinations(hetero_hemizygous_alleles)
        return BestAndMultipleBinderMhcII._get_sorted_epitopes_mhc2(
            hetero_hemizygous_alleles, homozygous_alleles, predictions
        )
