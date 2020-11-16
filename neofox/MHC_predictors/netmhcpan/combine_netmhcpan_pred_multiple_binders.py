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

from neofox.MHC_predictors.netmhcpan.multiple_binders import MultipleBinding
from neofox.MHC_predictors.netmhcpan.netmhcpan_prediction import NetMhcPanPredictor
from neofox.helpers import intermediate_files
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.neoantigen import Annotation, Mhc1
from neofox.model.wrappers import AnnotationFactory
import neofox.helpers.casting as casting


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
        self.best4_mhc_score = None
        self.best4_mhc_epitope = None
        self.best4_mhc_allele = None
        self.best4_mhc_position = None
        self.best4_affinity = None
        self.best4_affinity_epitope = "-"
        self.best4_affinity_allele = None
        self.best4_affinity_position = None
        self.epitope_affinities = None
        self.generator_rate = None
        self.mhcI_score_9mer = None
        self.mhcI_score_allele_9mer = None
        self.mhcI_score_position_9mer = None
        self.mhcI_score_epitope_9mer = "-"
        self.mhcI_affinity_9mer = None
        self.mhcI_affinity_allele_9mer = None
        self.mhcI_affinity_position_9mer = None
        self.mhcI_affinity_epitope_9mer = "-"
        # WT features
        self.best4_mhc_score_WT = None
        self.best4_mhc_epitope_WT = None
        self.best4_mhc_allele_WT = None
        self.best4_affinity_WT = None
        self.best4_affinity_epitope_WT = None
        self.best4_affinity_allele_WT = None
        self.mhcI_score_9mer_WT = None
        self.mhcI_score_allele_9mer_WT = None
        self.mhcI_score_epitope_9mer_WT = None
        self.mhcI_affinity_9mer_WT = None
        self.mhcI_affinity_allele_9mer_WT = None
        self.mhcI_affinity_epitope_9mer_WT = None
        self.position_9mer = None
        self.mutation_in_anchor_9mer = None

    def calculate_phbr_i(self, best_mhc_scores_per_allele):
        """returns list of multiple binding scores for mhcII considering best epitope per allele, applying different types of means (harmonic ==> PHRB-II, Marty et al).
        2 copies of DRA - DRB1 --> consider this gene 2x when averaging mhcii binding scores
        """
        list_best_mhc_scores_per_allele = list(best_mhc_scores_per_allele)
        phbr_i = None
        if len(list_best_mhc_scores_per_allele) == 6:
            phbr_i = stats.hmean(list_best_mhc_scores_per_allele)
        return phbr_i

    def run(self, sequence_wt: str, sequence_mut: str, mhc1_alleles_patient: List[Mhc1], mhc1_alleles_available: Set):
        """
        predicts MHC epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed
        """
        # mutation
        self._initialise()
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        netmhcpan = NetMhcPanPredictor(runner=self.runner, configuration=self.configuration)
        tmp_fasta = intermediate_files.create_temp_fasta(sequences=[sequence_mut], prefix="tmp_singleseq_")
        # print alleles
        netmhcpan.mhc_prediction(mhc1_alleles_patient, mhc1_alleles_available, tmp_fasta, tmp_prediction)
        position_of_mutation = netmhcpan.mut_position_xmer_seq(sequence_mut=sequence_mut, sequence_wt=sequence_wt)
        predicted_neoepitopes = netmhcpan.filter_binding_predictions(position_of_mutation=position_of_mutation, tmppred=tmp_prediction)
        # multiple binding
        predicted_neoepitopes_transformed = MultipleBinding.transform_mhc_prediction_output(predicted_neoepitopes)
        self.epitope_affinities = "/".join([str(epitope[1]) for epitope in predicted_neoepitopes_transformed])
        # best prediction
        best_predicted_epitope_rank = netmhcpan.minimal_binding_score(predicted_neoepitopes)
        self.best4_mhc_score = casting.to_float(netmhcpan.add_best_epitope_info(best_predicted_epitope_rank, "%Rank"))
        self.best4_mhc_epitope = netmhcpan.add_best_epitope_info(best_predicted_epitope_rank, "Peptide")
        self.best4_mhc_allele = netmhcpan.add_best_epitope_info(best_predicted_epitope_rank, "HLA")
        self.best4_mhc_position = netmhcpan.add_best_epitope_info(best_predicted_epitope_rank, "Pos")
        best_predicted_epitope_affinity = netmhcpan.minimal_binding_score(predicted_neoepitopes, rank=False)
        self.best4_affinity = casting.to_float(netmhcpan.add_best_epitope_info(best_predicted_epitope_affinity, "Aff(nM)"))
        self.best4_affinity_epitope = netmhcpan.add_best_epitope_info(best_predicted_epitope_affinity, "Peptide")
        self.best4_affinity_allele = netmhcpan.add_best_epitope_info(best_predicted_epitope_affinity, "HLA")
        self.best4_affinity_position = netmhcpan.add_best_epitope_info(best_predicted_epitope_affinity, "Pos")
        # best predicted epitope of length 9
        predicted_epitopes_9mer = netmhcpan.filter_for_9mers(predicted_neoepitopes)
        best_predicted_epitope_9mer_rank = netmhcpan.minimal_binding_score(predicted_epitopes_9mer)
        self.mhcI_score_9mer = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_rank, "%Rank")
        self.mhcI_score_allele_9mer = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_rank, "HLA")
        self.mhcI_score_position_9mer = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_rank, "Pos")
        self.mhcI_score_epitope_9mer = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_rank, "Peptide")
        best_predicted_epitope_9mer_affinity = netmhcpan.minimal_binding_score(predicted_epitopes_9mer, rank=False)
        self.mhcI_affinity_9mer = casting.to_float(netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_affinity, "Aff(nM)"))
        self.mhcI_affinity_allele_9mer = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_affinity, "HLA")
        self.mhcI_affinity_position_9mer = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_affinity, "Pos")
        self.mhcI_affinity_epitope_9mer = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_affinity, "Peptide")
        # multiple binding based on affinity
        all_affinities = [epitope[1] for epitope in predicted_neoepitopes_transformed]
        self.generator_rate = MultipleBinding.determine_number_of_binders(list_scores=all_affinities, threshold=50)
        best_epitopes_per_allele = MultipleBinding.extract_best_epitope_per_alelle(
            predicted_neoepitopes_transformed, mhc1_alleles_patient)
        # PHBR-I
        best_epitopes_per_allele = [epitope[0] for epitope in best_epitopes_per_allele]
        self.phbr_i = self.calculate_phbr_i(best_epitopes_per_allele)


        # wt
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        netmhcpan = NetMhcPanPredictor(runner=self.runner, configuration=self.configuration)
        tmp_fasta = intermediate_files.create_temp_fasta(sequences=[sequence_wt],
                                                         prefix="tmp_singleseq_")
        netmhcpan.mhc_prediction(mhc1_alleles_patient, mhc1_alleles_available, tmp_fasta, tmp_prediction)
        predicted_neoepitopes_wt = netmhcpan.filter_binding_predictions(position_of_mutation=position_of_mutation,
                                                              tmppred=tmp_prediction)
        # best prediction
        best_predicted_epitope_rank_wt = \
            netmhcpan.filter_for_WT_epitope_position(predicted_neoepitopes_wt, self.best4_mhc_epitope,
                                              position_epitope=self.best4_mhc_position)
        self.best4_mhc_score_WT = casting.to_float(netmhcpan.add_best_epitope_info(best_predicted_epitope_rank_wt, "%Rank"))
        self.best4_mhc_epitope_WT = netmhcpan.add_best_epitope_info(best_predicted_epitope_rank_wt, "Peptide")
        self.best4_mhc_allele_WT = netmhcpan.add_best_epitope_info(best_predicted_epitope_rank_wt, "HLA")

        best_predicted_epitope_affinity_wt = \
            netmhcpan.filter_for_WT_epitope_position(predicted_neoepitopes_wt, self.best4_affinity_epitope,
                                              position_epitope=self.best4_affinity_position)
        self.best4_affinity_WT = casting.to_float(netmhcpan.add_best_epitope_info(best_predicted_epitope_affinity_wt, "Aff(nM)"))
        self.best4_affinity_epitope_WT = netmhcpan.add_best_epitope_info(best_predicted_epitope_affinity_wt, "Peptide")
        self.best4_affinity_allele_WT = netmhcpan.add_best_epitope_info(best_predicted_epitope_affinity_wt, "HLA")
        # best predicted epitope of length 9
        predicted_epitopes_9mer_wt = netmhcpan.filter_for_9mers(predicted_neoepitopes_wt)
        best_predicted_epitope_9mer_rank_wt = \
            netmhcpan.filter_for_WT_epitope_position(predicted_epitopes_9mer_wt, self.mhcI_score_epitope_9mer,
                                              position_epitope=self.mhcI_score_position_9mer)
        best_predicted_epitope_9mer_affinity_wt = \
            netmhcpan.filter_for_WT_epitope_position(predicted_epitopes_9mer_wt, sequence_mut=self.mhcI_affinity_epitope_9mer,
                                              position_epitope=self.mhcI_affinity_position_9mer)
        self.mhcI_score_9mer_WT = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_rank_wt, "%Rank")
        self.mhcI_score_allele_9mer_WT = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_rank_wt, "HLA")
        self.mhcI_score_epitope_9mer_WT = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_rank_wt, "Peptide")
        self.mhcI_affinity_9mer_WT = casting.to_float(netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_affinity_wt, "Aff(nM)"))
        self.mhcI_affinity_allele_9mer_WT = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_affinity_wt, "HLA")
        self.mhcI_affinity_epitope_9mer_WT = netmhcpan.add_best_epitope_info(best_predicted_epitope_9mer_affinity_wt, "Peptide")

    def get_annotations(self) -> List[Annotation]:
        annotations = [AnnotationFactory.build_annotation(value=self.best4_mhc_score, name="Best_rank_MHCI_score"),
                       AnnotationFactory.build_annotation(value=self.best4_mhc_epitope,
                                                          name="Best_rank_MHCI_score_epitope"),
                       AnnotationFactory.build_annotation(value=self.best4_mhc_allele,
                                                          name="Best_rank_MHCI_score_allele"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity, name="Best_affinity_MHCI_score"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity_epitope,
                                                          name="Best_affinity_MHCI_epitope"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity_allele,
                                                          name="Best_affinity_MHCI_allele"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_9mer, name="Best_rank_MHCI_9mer_score"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_epitope_9mer,
                                                          name="Best_rank_MHCI_9mer_epitope"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_allele_9mer,
                                                          name="Best_rank_MHCI_9mer_allele"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_9mer,
                                                          name="Best_affinity_MHCI_9mer_score"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_allele_9mer,
                                                          name="Best_affinity_MHCI_9mer_allele"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_epitope_9mer,
                                                          name="Best_affinity_MHCI_9mer_epitope"),
                       # wt
                       AnnotationFactory.build_annotation(value=self.best4_affinity_WT,
                                                          name="Best_affinity_MHCI_score_WT"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity_epitope_WT,
                                                          name="Best_affinity_MHCI_epitope_WT"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity_allele_WT,
                                                          name="Best_affinity_MHCI_allele_WT"),
                       AnnotationFactory.build_annotation(value=self.best4_mhc_score_WT,
                                                          name="Best_rank_MHCI_score_WT"),
                       AnnotationFactory.build_annotation(value=self.best4_mhc_epitope_WT,
                                                          name="Best_rank_MHCI_score_epitope_WT"),
                       AnnotationFactory.build_annotation(value=self.best4_mhc_allele_WT,
                                                          name="Best_rank_MHCI_score_allele_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_9mer_WT,
                                                          name="Best_rank_MHCI_9mer_score_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_epitope_9mer_WT,
                                                          name="Best_rank_MHCI_9mer_epitope_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_allele_9mer_WT,
                                                          name="Best_rank_MHCI_9mer_allele_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_9mer_WT,
                                                          name="Best_affinity_MHCI_9mer_score_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_allele_9mer_WT,
                                                          name="Best_affinity_MHCI_9mer_allele_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_epitope_9mer_WT,
                                                          name="Best_affinity_MHCI_9mer_epitope_WT"),
                       # generator rate
                       AnnotationFactory.build_annotation(value=self.generator_rate, name="Generator_rate"),
                       AnnotationFactory.build_annotation(value=self.phbr_i, name="PHBR-I")
                       ]

        annotations.extend(self._get_positions_and_mutation_in_anchor())
        return annotations

    def _get_positions_and_mutation_in_anchor(self):
        """
        returns if mutation is in anchor position for best affinity epitope over all lengths and best 9mer affinity
        """
        position_9mer = EpitopeHelper.position_of_mutation_epitope(
            wild_type=self.mhcI_affinity_epitope_9mer_WT, mutation=self.mhcI_affinity_epitope_9mer)
        mutation_in_anchor_9mer = EpitopeHelper.position_in_anchor_position(
                                                       position_mhci=position_9mer,
                                                       peptide_length=len(self.mhcI_affinity_epitope_9mer))
        return [
            AnnotationFactory.build_annotation(value=position_9mer, name="Best_affinity_MHCI_9mer_position_mutation"),
            AnnotationFactory.build_annotation(value=mutation_in_anchor_9mer,
                                               name="Best_affinity_MHCI_9mer_anchor_mutated")
            ]
