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
from neofox.helpers.mhc_helper import MhcHelper
from neofox.helpers.runner import Runner
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Annotation, Mhc1, PredictedEpitope, Neoantigen
from neofox.model.factories import AnnotationFactory
from neofox.references.references import DependenciesConfiguration, ORGANISM_HOMO_SAPIENS


class BestAndMultipleBinder:
    def __init__(
            self, runner: Runner, configuration: DependenciesConfiguration, mhc_parser: MhcParser,
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
        self.generator_rate = None
        self.mutation_in_anchor_9mer = None
        self.generator_rate = None
        self.generator_rate_adn = None
        self.generator_rate_cdn = None
        self.best_epitope_by_rank = EpitopeHelper.get_empty_epitope()
        self.best_epitope_by_affinity = EpitopeHelper.get_empty_epitope()
        self.best_ninemer_epitope_by_affinity = EpitopeHelper.get_empty_epitope()
        self.best_ninemer_epitope_by_rank = EpitopeHelper.get_empty_epitope()
        self.predictions = []

    def calculate_phbr_i(
        self, predictions: List[PredictedEpitope], mhc1_alleles: List[Mhc1]):
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
            phbr_i = stats.hmean(list(map(lambda e: e.rank_mutated, best_epitopes_per_allele)))
        return phbr_i

    @staticmethod
    def extract_best_epitope_per_alelle(
            epitopes: List[PredictedEpitope], mhc_isoforms: List[Mhc1]) -> List[PredictedEpitope]:
        """
        This function returns the predicted epitope with the lowest binding score for each patient allele,
        considering homozyogosity
        """
        homozygous_alleles = MhcHelper.get_homozygous_mhc1_alleles(mhc_isoforms)
        hetero_hemizygous_alleles = (MhcHelper.get_heterozygous_or_hemizygous_mhc1_alleles(mhc_isoforms))
        return BestAndMultipleBinder._get_sorted_epitopes(hetero_hemizygous_alleles, homozygous_alleles, epitopes)

    @staticmethod
    def _get_sorted_epitopes(
        hetero_hemizygous_alleles, homozygous_alleles, predictions: List[PredictedEpitope]) -> List[PredictedEpitope]:

        # groups epitopes by allele
        epitopes_by_allele = {}
        for p in predictions:
            epitopes_by_allele.setdefault(p.allele_mhc_i.name, []).append(p)
        # chooses the best epitope per allele while considering zygosity
        best_epis_per_allele = []
        for list_alleles in epitopes_by_allele.values():
            # sort by rank to choose the best epitope, ties are solved choosing the first peptide in alphabetcial order
            list_alleles.sort(key=lambda x: (x.rank_mutated, x.mutated_peptide))
            best_epitope = list_alleles[0]
            if best_epitope.allele_mhc_i.name in hetero_hemizygous_alleles:
                best_epis_per_allele.append(best_epitope)  # adds the epitope once
            if best_epitope.allele_mhc_i.name in homozygous_alleles:
                best_epis_per_allele.append(best_epitope)
                best_epis_per_allele.append(best_epitope)  # adds the epitope twice
        return best_epis_per_allele

    @staticmethod
    def determine_number_of_binders(predictions: List[PredictedEpitope], threshold=50):
        """
        Determines the number of HLA I binders per mutation based on an affinity threshold. Default is set to 50, which is threshold used in generator rate.
        """
        scores = [epitope.affinity_mutated for epitope in predictions]
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
            dai_values.append(epitope.affinity_mutated)
            if epitope.affinity_mutated < 5000:
                if epitope.wild_type_peptide is not None and epitope.affinity_wild_type is not None:
                    dai = epitope.affinity_wild_type / epitope.affinity_mutated
                    if dai > threshold:
                        number_binders += 1

        if len(predictions) == 0:
            number_binders = None
        return number_binders

    def run(
        self,
        neoantigen: Neoantigen,
        mhc1_alleles_patient: List[Mhc1],
        mhc1_alleles_available: Set,
        uniprot,
    ):
        """
        predicts MHC epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed
        """
        self._initialise()

        # gets all predictions overlapping the mutation and not present in the WT proteome
        available_alleles = self.netmhcpan.get_only_available_alleles(mhc1_alleles_patient, mhc1_alleles_available)
        predictions = self.netmhcpan.get_predictions(available_alleles, neoantigen, uniprot)
        if neoantigen.wild_type_xmer:
            # SNVs with available WT
            # runs the netMHCpan WT predictions and then pair them with previous predictions
            # based on length, position within neoepitope and HLA allele
            predictions_wt = self.netmhcpan.get_wt_predictions(available_alleles, neoantigen)
            predictions = EpitopeHelper.pair_predictions(predictions=predictions, predictions_wt=predictions_wt)
        else:
            # alternative mutation classes or missing WT
            # do BLAST search for all predicted epitopes to identify the closest WT peptide and
            # predict MHC binding for the identified peptide sequence
            predictions = EpitopeHelper.set_wt_epitope_by_homology(predictions, self.blastp_runner)
            predictions = self.netmhcpan.set_wt_netmhcpan_scores(predictions)

        self.predictions = predictions

        if len(predictions) > 0:
            # best prediction
            self.best_epitope_by_rank = EpitopeHelper.select_best_by_rank(predictions)
            self.best_epitope_by_affinity = EpitopeHelper.select_best_by_affinity(predictions)

            # best predicted epitope of length 9
            ninemer_predictions = EpitopeHelper.filter_for_9mers(predictions)
            self.best_ninemer_epitope_by_rank = EpitopeHelper.select_best_by_rank(ninemer_predictions)
            self.best_ninemer_epitope_by_affinity = EpitopeHelper.select_best_by_affinity(ninemer_predictions)

            # multiple binding based on affinity
            self.generator_rate_cdn = self.determine_number_of_binders(predictions=predictions, threshold=50)
            self.generator_rate_adn = self.determine_number_of_alternative_binders(predictions=predictions)
            if self.generator_rate_adn is not None and self.generator_rate_cdn is not None:
                    self.generator_rate = self.generator_rate_adn + self.generator_rate_cdn

            # PHBR-I
            self.phbr_i = self.calculate_phbr_i(predictions=predictions, mhc1_alleles=mhc1_alleles_patient)

    def get_annotations(self) -> List[Annotation]:
        annotations = []
        if self.best_epitope_by_rank:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.rank_mutated, name="NetMHCpan_MHCI_bestRank_rank"
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.mutated_peptide,
                    name="NetMHCpan_bestRank_peptide",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.allele_mhc_i.name, name="NetMHCpan_bestRank_allele"
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.rank_wild_type, name="NetMHCpan_bestRank_rankWT"
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_rank.wild_type_peptide,
                    name="NetMHCpan_bestRank_peptideWT",
                )
            ])
        if self.best_epitope_by_affinity:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.affinity_mutated,
                    name="NetMHCpan_bestAffinity_affinity",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.mutated_peptide,
                    name="NetMHCpan_bestAffinity_peptide",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.allele_mhc_i.name,
                    name="NetMHCpan_bestAffinity_allele",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.affinity_wild_type,
                    name="NetMHCpan_bestAffinity_affinityWT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_epitope_by_affinity.wild_type_peptide,
                    name="NetMHCpan_bestAffinity_peptideWT",
                )])
        if self.best_ninemer_epitope_by_rank:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.rank_mutated,
                    name="NetMHCpan_bestRank9mer_rank",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.mutated_peptide,
                    name="NetMHCpan_bestRank9mer_peptide",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.allele_mhc_i.name,
                    name="NetMHCpan_bestRank9mer_allele",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.rank_wild_type,
                    name="NetMHCpan_bestRank9mer_rankWT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_rank.wild_type_peptide,
                    name="NetMHCpan_bestRank9mer_peptideWT",
                )
            ])
        if self.best_ninemer_epitope_by_affinity:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.affinity_mutated,
                    name="NetMHCpan_bestAffinity9mer_affinity",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.allele_mhc_i.name,
                    name="NetMHCpan_bestAffinity9mer_allele",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.mutated_peptide,
                    name="NetMHCpan_bestAffinity9mer_peptide",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.affinity_wild_type,
                    name="NetMHCpan_bestAffinity9mer_affinityWT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_ninemer_epitope_by_affinity.wild_type_peptide,
                    name="NetMHCpan_bestAffinity9mer_peptideWT",
                )
            ])

        if self.organism == ORGANISM_HOMO_SAPIENS:
            annotations.extend([AnnotationFactory.build_annotation(value=self.phbr_i, name="PHBR_I")])

        annotations.extend([
            # generator rate
            AnnotationFactory.build_annotation(value=self.generator_rate, name="GeneratorRate_MHCI"),
            AnnotationFactory.build_annotation(value=self.generator_rate_cdn, name="GeneratorRate_CDN_MHCI"),
            AnnotationFactory.build_annotation(value=self.generator_rate_adn, name="GeneratorRate_ADN_MHCI")
        ])
        annotations.extend(self._get_positions_and_mutation_in_anchor())
        return annotations

    def _get_positions_and_mutation_in_anchor(self):
        """
        returns if mutation is in anchor position for best affinity epitope over all lengths and best 9mer affinity
        """
        position_9mer = None
        mutation_in_anchor_9mer = None
        if self.best_ninemer_epitope_by_affinity.mutated_peptide and self.best_ninemer_epitope_by_affinity.wild_type_peptide:
            position_9mer = EpitopeHelper.position_of_mutation_epitope(epitope=self.best_ninemer_epitope_by_affinity)
            mutation_in_anchor_9mer = EpitopeHelper.position_in_anchor_position(
                position_mhci=position_9mer,
                peptide_length=len(self.best_ninemer_epitope_by_affinity.mutated_peptide),
            )
        annotations = [
            AnnotationFactory.build_annotation(
                value=position_9mer, name="NetMHCpan_bestAffinity9mer_positionMutation"
            ),
            AnnotationFactory.build_annotation(
                value=mutation_in_anchor_9mer,
                name="NetMHCpan_bestAffinity9mer_anchorMutated",
            ),
        ]
        return annotations

    @staticmethod
    def get_annotations_epitope_mhci(epitope: PredictedEpitope) -> List[Annotation]:
        position = EpitopeHelper.position_of_mutation_epitope(epitope=epitope)
        mutation_in_anchor = EpitopeHelper.position_in_anchor_position(
            position_mhci=position, peptide_length=len(epitope.mutated_peptide))
        return [
            AnnotationFactory.build_annotation(
                value=position,
                name='position_mutation'),
            AnnotationFactory.build_annotation(
                value=mutation_in_anchor,
                name='anchor_mutated')
        ]
