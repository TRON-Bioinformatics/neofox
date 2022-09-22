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
from neofox.MHC_predictors.netmhcpan.netmhcIIpan_prediction import NetMhcIIPanPredictor
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Annotation, Mhc2, Zygosity, Mhc2Isoform, Mhc2GeneName, PredictedEpitope, Neoantigen
from neofox.model.factories import AnnotationFactory
from neofox.references.references import DependenciesConfiguration, ORGANISM_HOMO_SAPIENS

MIN_LENGTH_MHC2_EPITOPE = 15


class BestAndMultipleBinderMhcII:
    def __init__(self, runner: Runner, configuration: DependenciesConfiguration, mhc_parser: MhcParser,
                 blastp_runner: BlastpRunner):
        self.runner = runner
        self.configuration = configuration
        self.mhc_parser = mhc_parser
        self.organism = mhc_parser.mhc_database.organism
        self.proteome_blastp_runner = blastp_runner
        self.netmhc2pan = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        self._initialise()

    def _initialise(self):
        self.phbr_ii = None
        self.generator_rate = None
        self.generator_rate_adn = None
        self.generator_rate_cdn = None
        self.best_predicted_epitope_rank = PredictedEpitope(
            mutated_peptide=None,
            position=None,
            isoform_mhc_i_i=Mhc2Isoform(name=None),
            affinity_mutated=None,
            rank_mutated=None,
        )
        self.best_predicted_epitope_affinity = PredictedEpitope(
            mutated_peptide=None,
            position=None,
            isoform_mhc_i_i=Mhc2Isoform(name=None),
            affinity_mutated=None,
            rank_mutated=None,
        )
        self.predictions = []

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
            if Mhc2GeneName.DRB1.name in allele_with_score.isoform_mhc_i_i.name:
                best_epitope_per_allele_mhc2_new.append(allele_with_score)
        if len(best_epitope_per_allele_mhc2_new) == 12:
            # 12 genes gene copies should be included into PHBR_II
            best_mhc_ii_scores_per_allele = [epitope.rank_mutated for epitope in best_epitope_per_allele_mhc2_new]
            phbr_ii = stats.hmean(best_mhc_ii_scores_per_allele)
        return phbr_ii

    @staticmethod
    def determine_number_of_binders(predictions: List[PredictedEpitope], threshold=1):
        """
        Determines the number of HLA II binders per mutation based on a rank threshold. Default is set to 1, which is threshold used in generator rate.
        """
        scores = [epitope.rank_mutated for epitope in predictions]
        number_binders = 0
        for score in scores:
            if score < threshold:
                number_binders += 1
        return number_binders if not len(scores) == 0 else None

    @staticmethod
    def determine_number_of_alternative_binders(predictions: List[PredictedEpitope], threshold=4):
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
            values.append(epitope.rank_mutated)
            if epitope.rank_mutated < 4:
                rank_mutation = epitope.rank_mutated
                if rank_mutation == 0:
                    rank_mutation = 0.01
                if epitope.wild_type_peptide is not None:
                    dai = epitope.rank_wild_type / rank_mutation
                    if dai > threshold:
                        number_binders += 1
        return number_binders if not len(values) == 0 else None

    def run(self, neoantigen: Neoantigen, mhc2_alleles_patient: List[Mhc2], mhc2_alleles_available: Set, uniprot):
        """predicts MHC II epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed"""
        # mutation
        self._initialise()

        allele_combinations = self.netmhc2pan.generate_mhc2_alelle_combinations(mhc2_alleles_patient)
        # TODO: migrate the available alleles into the model for alleles
        patient_mhc2_isoforms = self._get_only_available_combinations(allele_combinations, mhc2_alleles_available)

        # only process neoepitopes with a minimum length
        if len(neoantigen.mutated_xmer) >= MIN_LENGTH_MHC2_EPITOPE:

            predictions = self.netmhc2pan.get_predictions(neoantigen, patient_mhc2_isoforms, uniprot)

            if neoantigen.wild_type_xmer and len(neoantigen.wild_type_xmer) >= MIN_LENGTH_MHC2_EPITOPE:

                # SNVs with available WT
                # runs the netMHCIIpan WT predictions and then pair them with previous predictions
                # based on length, position within neoepitope and HLA allele
                predictions_wt = self.netmhc2pan.get_wt_predictions(neoantigen, patient_mhc2_isoforms)
                predictions = EpitopeHelper.pair_mhcii_predictions(predictions=predictions, predictions_wt=predictions_wt)
            else:

                # alternative mutation classes or missing WT
                # do BLAST search for all predicted epitopes to identify the closest WT peptide and
                # predict MHC binding for the identified peptide sequence
                predictions = EpitopeHelper.set_wt_epitope_by_homology(predictions, blastp_runner=self.proteome_blastp_runner)
                predictions = self.netmhc2pan.set_wt_netmhcpan_scores(predictions)

            self.predictions = predictions

            if len(predictions) > 0:

                # multiple binding
                best_predicted_epitopes_per_alelle = (
                    self.extract_best_epitope_per_mhc2_alelle(predictions, mhc2_alleles_patient))
                self.phbr_ii = self.calculate_phbr_ii(best_predicted_epitopes_per_alelle)

                # best prediction
                self.best_predicted_epitope_rank = EpitopeHelper.select_best_by_rank(predictions)
                self.best_predicted_epitope_affinity = EpitopeHelper.select_best_by_affinity(predictions)
                self.generator_rate_cdn = self.determine_number_of_binders(predictions=predictions)
                self.generator_rate_adn = self.determine_number_of_alternative_binders(predictions=predictions)
                if self.generator_rate_adn is not None and self.generator_rate_cdn is not None:
                    self.generator_rate = self.generator_rate_adn + self.generator_rate_cdn

    def _get_only_available_combinations(self, allele_combinations: List[Mhc2Isoform], set_available_mhc: Set[str]) -> List[str]:

        # parses isoforms into internal representation
        parsed_allele_combinations = self.netmhc2pan.represent_mhc2_isoforms(allele_combinations)
        patients_available_alleles = list(
            set(parsed_allele_combinations).intersection(set_available_mhc)
        )
        patients_not_available_alleles = list(
            set(parsed_allele_combinations).difference(set_available_mhc)
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
                    value=self.best_predicted_epitope_rank.rank_mutated,
                    name="NetMHCIIpan_bestRank_rank",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_rank.mutated_peptide,
                    name="NetMHCIIpan_bestRank_peptide",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_rank.isoform_mhc_i_i.name,
                    name="NetMHCIIpan_bestRank_allele",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_rank.rank_wild_type,
                    name="NetMHCIIpan_bestRank_rankWT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_rank.wild_type_peptide,
                    name="NetMHCIIpan_bestRank_peptideWT",
                ),
            ])
        if self.best_predicted_epitope_affinity:
            annotations.extend([
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity.affinity_mutated,
                    name="NetMHCIIpan_bestAffinity_affinity",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity.mutated_peptide,
                    name="NetMHCIIpan_bestAffinity_peptide",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity.isoform_mhc_i_i.name,
                    name="NetMHCIIpan_bestAffinity_allele",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity.affinity_wild_type,
                    name="NetMHCIIpan_bestAffinity_affinityWT",
                ),
                AnnotationFactory.build_annotation(
                    value=self.best_predicted_epitope_affinity.wild_type_peptide,
                    name="NetMHCIIpan_bestAffinity_peptideWT",
                )
            ])

        if self.organism == ORGANISM_HOMO_SAPIENS:
            annotations.extend([AnnotationFactory.build_annotation(value=self.phbr_ii, name="PHBR_II")])

        annotations.extend([
            # generator rate
            AnnotationFactory.build_annotation(value=self.generator_rate, name="GeneratorRate_MHCII"),
            AnnotationFactory.build_annotation(value=self.generator_rate_cdn, name="GeneratorRate_CDN_MHCII"),
            AnnotationFactory.build_annotation(value=self.generator_rate_adn, name="GeneratorRate_ADN_MHCII"),
        ])
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
            allele = p.isoform_mhc_i_i.name
            epitopes_by_allele.setdefault(allele, []).append(p)

        # chooses the best epitope per allele and considers zygosity
        best_epitopes_per_allele = []
        for allele, epitopes in epitopes_by_allele.items():
            # sort by rank to choose the best epitope, fixes ties with peptide by alphabetical order
            epitopes.sort(key=lambda e: (e.rank_mutated, e.mutated_peptide))
            best_epitope = epitopes[0]
            num_repetitions = 0
            if (
                best_epitope.isoform_mhc_i_i.name in hetero_hemizygous_allele_names
                or best_epitope.isoform_mhc_i_i.name in hetero_hemizygous_allele_names
            ):
                # adds the epitope once if alleles heterozygous
                num_repetitions = 1
            if (
                best_epitope.isoform_mhc_i_i in homozygous_allele_names
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

