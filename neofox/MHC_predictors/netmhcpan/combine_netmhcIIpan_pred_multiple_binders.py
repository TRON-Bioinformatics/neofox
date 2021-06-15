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
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import AbstractNetMhcPanPredictor
from neofox.MHC_predictors.netmhcpan.netmhcIIpan_prediction import NetMhcIIPanPredictor
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.runner import Runner
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Annotation, Mhc2, Zygosity, Mhc2Isoform, Mutation, Mhc2GeneName
from neofox.model.wrappers import AnnotationFactory
from neofox.references.references import DependenciesConfiguration

LENGTH_MHC2_EPITOPE = 15


class BestAndMultipleBinderMhcII:
    def __init__(self, runner: Runner, configuration: DependenciesConfiguration, mhc_parser: MhcParser,
                 blastp_runner: BlastpRunner):
        self.runner = runner
        self.configuration = configuration
        self.mhc_parser = mhc_parser
        self.proteome_blastp_runner = blastp_runner
        self._initialise()

    def _initialise(self):
        self.phbr_ii = None
        self.generator_rate = None
        self.generator_rate_adn = None
        self.generator_rate_cdn = None
        self.best_predicted_epitope_rank = PredictedEpitope(
            peptide=None,
            pos=None,
            hla=Mhc2Isoform(name=None),
            affinity_score=None,
            rank=None,
        )
        self.best_predicted_epitope_affinity = PredictedEpitope(
            peptide=None,
            pos=None,
            hla=Mhc2Isoform(name=None),
            affinity_score=None,
            rank=None,
        )
        self.best_predicted_epitope_rank_wt = PredictedEpitope(
            peptide=None,
            pos=None,
            hla=Mhc2Isoform(name=None),
            affinity_score=None,
            rank=None,
        )
        self.best_predicted_epitope_affinity_wt = PredictedEpitope(
            peptide=None,
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
            if Mhc2GeneName.DRB1.name in allele_with_score.hla.name:
                best_epitope_per_allele_mhc2_new.append(allele_with_score)
        if len(best_epitope_per_allele_mhc2_new) == 12:
            # 12 genes gene copies should be included into PHBR_II
            best_mhc_ii_scores_per_allele = [
                epitope.rank for epitope in best_epitope_per_allele_mhc2_new
            ]
            phbr_ii = stats.hmean(best_mhc_ii_scores_per_allele)
        return phbr_ii

    @staticmethod
    def determine_number_of_binders(predictions: List[PredictedEpitope], threshold=1):
        """
        Determines the number of HLA II binders per mutation based on a rank threshold. Default is set to 1, which is threshold used in generator rate.
        """
        scores = [epitope.rank for epitope in predictions]
        number_binders = 0
        for score in scores:
            if score < threshold:
                number_binders += 1
        return number_binders if not len(scores) == 0 else None

    @staticmethod
    def determine_number_of_alternative_binders(predictions: List[PredictedEpitope],
                                                predictions_wt: List[PredictedEpitope], threshold=4):
        """
        Determines the number of HLA II neoepitope candidates that bind stronger (4:1) to HLA in comparison to corresponding WT.
        With the netMHCIIpan4.0 the rank score can get a value of 0.0.  If this is the case, the next smaller possible
        value is used as a measure for best possible binding.
        This is 0.01 in case of netMHCIIpan4.0
        TODO: this is not optimal. is there a better long-term solution?
        """
        number_binders = 0
        values = []
        for epitope in predictions:
            values.append(epitope.rank)
            if epitope.rank < 4:
                wt_peptide = AbstractNetMhcPanPredictor.select_best_by_rank(
                    predictions=AbstractNetMhcPanPredictor.filter_wt_predictions_from_best_mutated(
                        predictions_wt, epitope
                    )
                )
                rank_mutation = epitope.rank
                if rank_mutation == 0:
                    rank_mutation = 0.01
                dai = wt_peptide.rank / rank_mutation
                if dai > threshold:
                    number_binders += 1
        return number_binders if not len(values) == 0 else None

    @staticmethod
    def determine_number_of_alternative_binders_alternative(predictions: List[PredictedEpitope],
                                                            predictions_wt: List[PredictedEpitope], threshold=4):
        """
        Determines the number of HLA I neoepitope candidates that bind stronger (10:1) to HLA in comparison to corresponding WT
        """
        number_binders = 0
        dai_values = []
        for mut, wt in zip(predictions, predictions_wt):
            dai_values.append(mut.rank)
            if mut.rank < 4:
                dai = mut.affinity_score / wt.affinity_score
                if dai > threshold:
                    number_binders += 1
        return number_binders if not len(dai_values) == 0 else None

    def run(
        self,
        mutation: Mutation,
        mhc2_alleles_patient: List[Mhc2],
        mhc2_alleles_available: Set,
        uniprot
    ):
        """predicts MHC II epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed"""
        # mutation
        self._initialise()
        netmhc2pan = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        allele_combinations = netmhc2pan.generate_mhc2_alelle_combinations(
            mhc2_alleles_patient
        )
        # TODO: migrate the available alleles into the model for alleles
        patient_mhc2_isoforms = self._get_only_available_combinations(
            allele_combinations, mhc2_alleles_available
        )

        predictions = netmhc2pan.mhc2_prediction(
            patient_mhc2_isoforms, mutation.mutated_xmer
        )
        if len(mutation.mutated_xmer) >= LENGTH_MHC2_EPITOPE:
            if mutation.wild_type_xmer:
                # make sure that predicted epitopes cover mutation in case of SNVs
                predictions = netmhc2pan.filter_peptides_covering_snv(
                    position_of_mutation=mutation.position, predictions=predictions
                )
            filtered_predictions = netmhc2pan.remove_peptides_in_proteome(
                predictions, uniprot
            )
            if len(filtered_predictions) > 0:
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
                self.generator_rate_cdn = self.determine_number_of_binders(
                    predictions=filtered_predictions
                )
                # MHC binding predictions for WT pepti
                if mutation.wild_type_xmer:
                    predictions = netmhc2pan.mhc2_prediction(
                        patient_mhc2_isoforms, mutation.wild_type_xmer
                    )
                    if len(mutation.wild_type_xmer) >= LENGTH_MHC2_EPITOPE:
                        filtered_predictions_wt = netmhc2pan.filter_peptides_covering_snv(
                            mutation.position, predictions
                        )
                        # best prediction
                        if self.best_predicted_epitope_rank:
                            self.best_predicted_epitope_rank_wt = netmhc2pan.select_best_by_rank(
                                netmhc2pan.filter_wt_predictions_from_best_mutated(
                                    filtered_predictions_wt, self.best_predicted_epitope_rank
                                )
                            )
                        if self.best_predicted_epitope_affinity:
                            self.best_predicted_epitope_affinity_wt = (
                                netmhc2pan.select_best_by_affinity(
                                    netmhc2pan.filter_wt_predictions_from_best_mutated(
                                        filtered_predictions_wt, self.best_predicted_epitope_affinity
                                    )
                                )
                            )
                        if len(mutation.mutated_xmer) >= LENGTH_MHC2_EPITOPE:
                            self.generator_rate_adn = self.determine_number_of_alternative_binders(
                                predictions=filtered_predictions, predictions_wt=filtered_predictions_wt
                            )
                            if self.generator_rate_adn is not None:
                                if self.generator_rate_cdn is not None:
                                    self.generator_rate = self.generator_rate_adn + self.generator_rate_cdn
                else:
                    # alternative mutation classes
                    # do BLAST search for all predicted epitopes  covering mutation to identify WT peptide and
                    # predict MHC binding for the identified peptide sequence
                    peptides_wt = netmhc2pan.find_wt_epitope_for_alternative_mutated_epitope(filtered_predictions)
                    filtered_predictions_wt = []
                    for wt_peptide, mut_peptide in zip(peptides_wt, filtered_predictions):
                        if wt_peptide is not None:
                            filtered_predictions_wt.extend(netmhc2pan.mhc2_prediction_peptide(mut_peptide.hla, wt_peptide))
                    if self.best_predicted_epitope_rank:
                        self.best_predicted_epitope_rank_wt = netmhc2pan.filter_wt_predictions_from_best_mutated_alernative(
                            mut_predictions=filtered_predictions, wt_predictions=filtered_predictions_wt,
                            best_mutated_epitope=self.best_predicted_epitope_rank)
                        if self.best_predicted_epitope_affinity:
                            self.best_predicted_epitope_affinity_wt = \
                                netmhc2pan.filter_wt_predictions_from_best_mutated_alernative(
                                    mut_predictions=filtered_predictions, wt_predictions=filtered_predictions_wt,
                                    best_mutated_epitope=self.best_predicted_epitope_affinity
                                )
                        if len(mutation.mutated_xmer) >= LENGTH_MHC2_EPITOPE:
                            # generator rate for MHC II
                            self.generator_rate_cdn = self.determine_number_of_binders(
                                predictions=filtered_predictions
                            )
                            self.generator_rate_adn = self.determine_number_of_alternative_binders_alternative(
                                predictions=filtered_predictions, predictions_wt=filtered_predictions_wt
                            )

                    if self.generator_rate_adn is not None:
                        if self.generator_rate_cdn is not None:
                            self.generator_rate = self.generator_rate_adn + self.generator_rate_cdn

    @staticmethod
    def _get_only_available_combinations(allele_combinations: List[Mhc2Isoform], set_available_mhc: List[str]) -> List[str]:

        # parses isoforms into internal representation
        parsed_allele_combinations = NetMhcIIPanPredictor.represent_mhc2_isoforms(allele_combinations)
        patients_available_alleles = list(
            set(parsed_allele_combinations).intersection(set(set_available_mhc))
        )
        patients_not_available_alleles = list(
            set(parsed_allele_combinations).difference(set(set_available_mhc))
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
                    value=self.best_predicted_epitope_rank.hla.name,
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
                    value=self.best_predicted_epitope_affinity.hla.name,
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
                    value=self.best_predicted_epitope_rank_wt.hla.name,
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
                    value=self.best_predicted_epitope_affinity_wt.hla.name,
                    name="Best_affinity_MHCII_allele_WT",
                )])
        annotations.extend(
            [
                AnnotationFactory.build_annotation(value=self.phbr_ii, name="PHBR_II"),
                # generator rate
                AnnotationFactory.build_annotation(value=self.generator_rate, name="Generator_rate_MHCII"),
                AnnotationFactory.build_annotation(value=self.generator_rate_cdn, name="Generator_rate_CDN_MHCII"),
                AnnotationFactory.build_annotation(value=self.generator_rate_adn, name="Generator_rate_ADN_MHCII"),
            ]

        )
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
        hetero_hemizygous_alleles: List[Mhc2Isoform],
        homozygous_alleles: List[Mhc2Isoform],
        predictions: List[PredictedEpitope],
    ) -> List[PredictedEpitope]:

        hetero_hemizygous_allele_names = [a.name for a in hetero_hemizygous_alleles]
        homozygous_allele_names = [a.name for a in homozygous_alleles]

        # groups epitopes by allele
        epitopes_by_allele = {}
        for p in predictions:
            allele = p.hla.name
            epitopes_by_allele.setdefault(allele, []).append(p)

        # chooses the best epitope per allele and considers zygosity
        best_epitopes_per_allele = []
        for allele, epitopes in epitopes_by_allele.items():
            # sort by rank to choose the best epitope, fixes ties with peptide by alphabetical order
            epitopes.sort(key=lambda e: (e.rank, e.peptide))
            best_epitope = epitopes[0]
            num_repetitions = 0
            if (
                best_epitope.hla.name in hetero_hemizygous_allele_names
                or best_epitope.hla.name in hetero_hemizygous_allele_names
            ):
                # adds the epitope once if alleles heterozygous
                num_repetitions = 1
            if (
                best_epitope.hla in homozygous_allele_names
            ):
                # adds the epitope twice if one allele is homozygous
                num_repetitions = 2
            best_epitopes_per_allele.extend(
                [best_epitope for _ in range(num_repetitions)]
            )
        return best_epitopes_per_allele

    @staticmethod
    def extract_best_epitope_per_mhc2_alelle(
            predictions: List[PredictedEpitope], mhc_isoforms: List[Mhc2]
    ) -> List[PredictedEpitope]:
        """
        This function returns the predicted epitope with the lowest binding score for each patient allele,
        considering homozygosity
        """
        homozygous_alleles = BestAndMultipleBinderMhcII._get_homozygous_mhc2_alleles(
            mhc_isoforms
        )
        homozygous_alleles = NetMhcIIPanPredictor.generate_mhc2_alelle_combinations(homozygous_alleles)
        # removes duplicated homozygous combinations
        homozygous_alleles = list({a.name: a for a in homozygous_alleles}.values())
        hetero_hemizygous_alleles = (
            BestAndMultipleBinderMhcII._get_heterozygous_or_hemizygous_mhc2_alleles(
                mhc_isoforms
            )
        )
        hetero_hemizygous_alleles = NetMhcIIPanPredictor.generate_mhc2_alelle_combinations(hetero_hemizygous_alleles)
        return BestAndMultipleBinderMhcII._get_sorted_epitopes_mhc2(
            hetero_hemizygous_alleles, homozygous_alleles, predictions
        )

    @staticmethod
    def transform_mhc2allele(mhc2_allele):
        if "DR" in mhc2_allele:
            allele_out = mhc2_allele.rstrip("HLA-").replace("*","_").replace(":","")
        else:
            allele_out = mhc2_allele.replace("*", "").replace(":", "")
        return allele_out

