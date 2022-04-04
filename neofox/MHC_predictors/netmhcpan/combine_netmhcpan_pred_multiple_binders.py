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
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Annotation, Mhc1, Zygosity, Mutation, MhcAllele, PredictedEpitope
from neofox.model.factories import AnnotationFactory
from neofox.references.references import DependenciesConfiguration, ORGANISM_HOMO_SAPIENS


class BestAndMultipleBinder:
    def __init__(self, runner: Runner, configuration: DependenciesConfiguration, mhc_parser: MhcParser,
                 blastp_runner: BlastpRunner):
        self.runner = runner
        self.configuration = configuration
        self.mhc_parser = mhc_parser
        self.organism = mhc_parser.mhc_database.organism
        self.blastp_runner = blastp_runner
        self._initialise()
        self.netmhcpan = NetMhcPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.blastp_runner
        )

    def _initialise(self):
        self.phbr_i = None
        self.predictions = None
        self.generator_rate = None
        self.mutation_in_anchor_9mer = None
        self.generator_rate = None
        self.generator_rate_adn = None
        self.generator_rate_cdn = None
        self.best_epitope_by_rank = self._get_empty_epitope()
        self.best_epitope_by_affinity = self._get_empty_epitope()
        self.best_ninemer_epitope_by_affinity = self._get_empty_epitope()
        self.best_ninemer_epitope_by_rank = self._get_empty_epitope()

    @staticmethod
    def _get_empty_epitope():
        return PredictedEpitope(
            peptide=None,
            position=None,
            hla=MhcAllele(name=None),
            affinity_score=None,
            rank=None,
        )

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
        homozygous_alleles = BestAndMultipleBinder._get_homozygous_mhc1_alleles(mhc_isoforms)
        hetero_hemizygous_alleles = (
            BestAndMultipleBinder._get_heterozygous_or_hemizygous_mhc1_alleles(mhc_isoforms))
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
            epitopes_by_allele.setdefault(p.allele.name, []).append(p)
        # chooses the best epitope per allele while considering zygosity
        best_epis_per_allele = []
        for list_alleles in epitopes_by_allele.values():
            # sort by rank to choose the best epitope, ties are solved choosing the first peptide in alphabetcial order
            list_alleles.sort(key=lambda x: (x.rank, x.mutated_peptide))
            best_epitope = list_alleles[0]
            if best_epitope.allele.name in hetero_hemizygous_alleles:
                best_epis_per_allele.append(best_epitope)  # adds the epitope once
            if best_epitope.allele.name in homozygous_alleles:
                best_epis_per_allele.append(best_epitope)
                best_epis_per_allele.append(best_epitope)  # adds the epitope twice
        return best_epis_per_allele

    @staticmethod
    def determine_number_of_binders(predictions: List[PredictedEpitope], threshold=50):
        """
        Determines the number of HLA I binders per mutation based on an affinity threshold. Default is set to 50, which is threshold used in generator rate.
        """
        scores = [epitope.affinity_score for epitope in predictions]
        number_binders = 0
        for score in scores:
            if score < threshold:
                number_binders += 1
        return number_binders if not len(scores) == 0 else None

    @staticmethod
    def determine_number_of_alternative_binders(predictions: List[PredictedEpitope], threshold=10):
        """
        Determines the number of HLA I neoepitope candidates that bind stronger (10:1) to HLA in comparison to corresponding WT
        """
        number_binders = 0
        dai_values = []
        for epitope in predictions:
            dai_values.append(epitope.affinity_score)
            if epitope.affinity_score < 5000:
                if epitope.wild_type_peptide is not None and epitope.affinity_score_wild_type is not None:
                    dai = epitope.affinity_score_wild_type / epitope.affinity_score
                    if dai > threshold:
                        number_binders += 1

        if len(predictions) == 0:
            number_binders = None
        return number_binders

    @staticmethod
    def determine_number_of_alternative_binders_alternative(predictions: List[PredictedEpitope],
                                                predictions_wt: List[PredictedEpitope], threshold=10):
        """
        Determines the number of HLA I neoepitope candidates that bind stronger (10:1) to HLA in comparison to corresponding WT
        """
        number_binders = 0
        dai_values = []
        for mut, wt in zip(predictions, predictions_wt):
            dai_values.append(mut.affinity_score)
            if mut.affinity_score < 5000 and wt.affinity_score:
                dai = mut.affinity_score / wt.affinity_score
                if dai > threshold:
                    number_binders += 1
        return number_binders if not len(dai_values) == 0 else None

    def run(
        self,
        mutation: Mutation,
        mhc1_alleles_patient: List[Mhc1],
        mhc1_alleles_available: Set,
        uniprot,
    ):
        """
        predicts MHC epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed
        """
        self._initialise()

        # gets all predictions overlapping the mutation and not present in the WT proteome
        predictions = self._get_predictions(mhc1_alleles_available, mhc1_alleles_patient, mutation, uniprot)
        if mutation.wild_type_xmer:
            # SNVs with available WT
            # runs the netMHCpan WT predictions and then pair them with previous predictions
            # based on length, position within neoepitope and HLA allele
            predictions_wt = self._get_wt_predictions(mhc1_alleles_available, mhc1_alleles_patient, mutation)
            predictions = self._pair_predictions(predictions=predictions, predictions_wt=predictions_wt)
        else:
            # alternative mutation classes or missing WT
            # do BLAST search for all predicted epitopes to identify the closest WT peptide and
            # predict MHC binding for the identified peptide sequence
            predictions = self.netmhcpan.set_wt_epitope_by_homology(predictions)
            predictions = self._set_wt_netmhcpan_scores(mhc1_alleles_available, predictions)

        self.predictions = predictions

        if len(predictions) > 0:
            # best prediction
            self.best_epitope_by_rank = self.netmhcpan.select_best_by_rank(
                predictions, none_value=self._get_empty_epitope())
            self.best_epitope_by_affinity = self.netmhcpan.select_best_by_affinity(
                predictions, none_value=self._get_empty_epitope())

            # best predicted epitope of length 9
            ninemer_predictions = self.netmhcpan.filter_for_9mers(predictions)
            self.best_ninemer_epitope_by_rank = self.netmhcpan.select_best_by_rank(
                ninemer_predictions, none_value=self._get_empty_epitope())
            self.best_ninemer_epitope_by_affinity = self.netmhcpan.select_best_by_affinity(
                ninemer_predictions, none_value=self._get_empty_epitope())

            # multiple binding based on affinity
            self.generator_rate_cdn = self.determine_number_of_binders(predictions=predictions, threshold=50)
            self.generator_rate_adn = self.determine_number_of_alternative_binders(predictions=predictions)
            if self.generator_rate_adn is not None and self.generator_rate_cdn is not None:
                    self.generator_rate = self.generator_rate_adn + self.generator_rate_cdn

            # PHBR-I
            self.phbr_i = self.calculate_phbr_i(predictions=predictions, mhc1_alleles=mhc1_alleles_patient)

    def _pair_predictions(self, predictions, predictions_wt) -> List[PredictedEpitope]:
        for prediction in predictions:
            for prediction_wt in predictions_wt:
                if len(prediction_wt.peptide) == len(prediction.peptide) and \
                        prediction.position == prediction_wt.position and \
                        prediction.hla.name == prediction_wt.hla.name:
                    prediction.wild_type_peptide = prediction_wt.peptide
                    prediction.rank_wild_type = prediction_wt.rank
                    prediction.affinity_score_wild_type = prediction_wt.affinity_score
                    break
        return predictions

    def _set_wt_netmhcpan_scores(self, mhc1_alleles_available, predictions) -> List[PredictedEpitope]:
        for p in predictions:
            if p.wild_type_peptide is not None:
                wt_predictions = self.netmhcpan.mhc_prediction(
                    mhc_alleles=[Mhc1(zygosity=Zygosity.HOMOZYGOUS, alleles=[p.hla])],
                    set_available_mhc=mhc1_alleles_available,
                    sequence=p.wild_type_peptide, peptide_mode=True)
                if len(wt_predictions) >= 1:
                    # NOTE: netmhcpan in peptide mode should return only one epitope
                    p.rank_wild_type = wt_predictions[0].rank
                    p.affinity_score_wild_type = wt_predictions[0].affinity_score
        return predictions

    def _get_predictions(self, mhc1_alleles_available, mhc1_alleles_patient, mutation, uniprot) -> List[PredictedEpitope]:
        predictions = self.netmhcpan.mhc_prediction(
            mhc1_alleles_patient, mhc1_alleles_available, mutation.mutated_xmer
        )
        if mutation.wild_type_xmer:
            # make sure that predicted epitopes cover mutation in case of SNVs
            predictions = self.netmhcpan.filter_peptides_covering_snv(
                position_of_mutation=mutation.position, predictions=predictions
            )
        # make sure that predicted neoepitopes are not part of the WT proteome
        filtered_predictions = self.netmhcpan.remove_peptides_in_proteome(
            predictions=predictions, uniprot=uniprot
        )
        return filtered_predictions

    def _get_wt_predictions(self, mhc1_alleles_available, mhc1_alleles_patient, mutation) -> List[PredictedEpitope]:
        predictions = []
        if mutation.wild_type_xmer:
            predictions = self.netmhcpan.mhc_prediction(
                mhc1_alleles_patient, mhc1_alleles_available, mutation.wild_type_xmer
            )
            # make sure that predicted epitopes cover mutation in case of SNVs
            predictions = self.netmhcpan.filter_peptides_covering_snv(
                position_of_mutation=mutation.position, predictions=predictions
            )
        return predictions

    def get_annotations(self, mutation) -> List[Annotation]:
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
                    value=self.best_epitope_by_rank.hla.name, name="Best_rank_MHCI_score_allele"
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.rank_wild_type, name="Best_rank_MHCI_score_WT"
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.wild_type_peptide,
                    name="Best_rank_MHCI_score_epitope_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.hla.name,
                    name="Best_rank_MHCI_score_allele_WT",
                )
            ])
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
                    value=self.best_epitope_by_affinity.hla.name,
                    name="Best_affinity_MHCI_allele",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.affinity_score_wild_type,
                    name="Best_affinity_MHCI_score_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.wild_type_peptide,
                    name="Best_affinity_MHCI_epitope_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.hla.name,
                    name="Best_affinity_MHCI_allele_WT",
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
                    value=self.best_ninemer_epitope_by_rank.hla.name,
                    name="Best_rank_MHCI_9mer_allele",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.rank_wild_type,
                    name="Best_rank_MHCI_9mer_score_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.wild_type_peptide,
                    name="Best_rank_MHCI_9mer_epitope_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.hla.name,
                    name="Best_rank_MHCI_9mer_allele_WT",
                )
            ])
        if self.best_ninemer_epitope_by_affinity:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.affinity_score,
                    name="Best_affinity_MHCI_9mer_score",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.hla.name,
                    name="Best_affinity_MHCI_9mer_allele",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.peptide,
                    name="Best_affinity_MHCI_9mer_epitope",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.affinity_score_wild_type,
                    name="Best_affinity_MHCI_9mer_score_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.hla.name,
                    name="Best_affinity_MHCI_9mer_allele_WT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.wild_type_peptide,
                    name="Best_affinity_MHCI_9mer_epitope_WT",
                )
            ])

        if self.organism == ORGANISM_HOMO_SAPIENS:
            annotations.extend([AnnotationFactory.build_annotation(value=self.phbr_i, name="PHBR_I")])

        annotations.extend([
            # generator rate
            AnnotationFactory.build_annotation(value=self.generator_rate, name="Generator_rate_MHCI"),
            AnnotationFactory.build_annotation(value=self.generator_rate_cdn, name="Generator_rate_CDN_MHCI"),
            AnnotationFactory.build_annotation(value=self.generator_rate_adn, name="Generator_rate_ADN_MHCI")
        ])
        annotations.extend(self._get_positions_and_mutation_in_anchor())
        return annotations

    def _get_positions_and_mutation_in_anchor(self):
        """
        returns if mutation is in anchor position for best affinity epitope over all lengths and best 9mer affinity
        """
        position_9mer = None
        mutation_in_anchor_9mer = None
        if self.best_ninemer_epitope_by_affinity.peptide and self.best_ninemer_epitope_by_affinity.wild_type_peptide:
            position_9mer = EpitopeHelper.position_of_mutation_epitope(epitope=self.best_ninemer_epitope_by_affinity)
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
