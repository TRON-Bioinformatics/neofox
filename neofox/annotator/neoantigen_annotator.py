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

from logzero import logger
from datetime import datetime
import neofox
from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMHC2pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.MHC_predictors.prime import Prime
from neofox.annotator.abstract_annotator import AbstractAnnotator
from neofox.annotator.neoantigen_mhc_binding_annotator import NeoantigenMhcBindingAnnotator
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.model.factories import AnnotationFactory
from neofox.model.mhc_parser import MhcParser
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.neoag.neoag_gbm_model import NeoagCalculator
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.published_features.expression import Expression
from neofox.model.neoantigen import Patient, Neoantigen, Annotations, PredictedEpitope
from neofox.published_features.vaxrank.vaxrank import VaxRank
from neofox.references.references import (
    ReferenceFolder,
    DependenciesConfiguration,
    ORGANISM_HOMO_SAPIENS
)


class NeoantigenAnnotator(AbstractAnnotator):
    def __init__(self, references: ReferenceFolder, configuration: DependenciesConfiguration,
                 tcell_predictor: TcellPrediction, self_similarity: SelfSimilarityCalculator,
                 rank_mhci_threshold=neofox.RANK_MHCI_THRESHOLD_DEFAULT,
                 rank_mhcii_threshold=neofox.RANK_MHCII_THRESHOLD_DEFAULT):
        """class to annotate neoantigens"""

        super().__init__(references, configuration, tcell_predictor, self_similarity)
        self.proteome_db = references.proteome_db
        self.available_alleles = references.get_available_alleles()
        self.rank_mhci_threshold = rank_mhci_threshold
        self.rank_mhcii_threshold = rank_mhcii_threshold

        # NOTE: these resources do not read any file thus can be initialised fast
        self.neoag_calculator = NeoagCalculator(runner=self.runner, configuration=configuration)
        self.expression_calculator = Expression()
        self.mhc_database = references.get_mhc_database()
        self.mhc_parser = MhcParser.get_mhc_parser(self.mhc_database)

        self.neoantigen_mhc_binding_annotator = NeoantigenMhcBindingAnnotator(
            references=references, configuration=configuration, proteome_blastp_runner=self.proteome_blastp_runner,
            uniprot=self.uniprot)

        self.resources_versions = references.get_resources_versions()

    def get_annotated_neoantigen(self, neoantigen: Neoantigen, patient: Patient, with_all_neoepitopes=False) -> Neoantigen:
        """Calculate new epitope features and add to dictionary that stores all properties"""
        neoantigen.neofox_annotations = Annotations(
            annotator="NeoFox",
            annotator_version=neofox.VERSION,
            timestamp="{:%Y%m%d%H%M%S%f}".format(datetime.now()),
            resources=self.resources_versions,
            annotations=[]
        )

        # Runs netmhcpan, netmhc2pan, mixmhcpred and mixmhc2prd in parallel
        (
            mixmhc2pred,
            mixmhcpred,
            netmhc2pan,
            netmhcpan,
            prime
        ) = self.neoantigen_mhc_binding_annotator.get_mhc_binding_annotations(neoantigen=neoantigen, patient=patient)

        # HLA I predictions: NetMHCpan
        if netmhcpan:
            neoantigen.neofox_annotations.annotations.extend(netmhcpan.get_annotations())
            neoantigen.neoepitopes_mhc_i = [e for e in netmhcpan.predictions if e.rank_mutated < self.rank_mhci_threshold]

        # HLA II predictions: NetMHCIIpan
        if netmhc2pan:
            neoantigen.neofox_annotations.annotations.extend(netmhc2pan.get_annotations())
            neoantigen.neoepitopes_mhc_i_i = [e for e in netmhc2pan.predictions if e.rank_mutated < self.rank_mhcii_threshold]

        # MixMHCpred
        if mixmhcpred is not None:
            neoantigen.neofox_annotations.annotations.extend(mixmhcpred.get_annotations())
            neoantigen.neoepitopes_mhc_i = AnnotationFactory.annotate_epitopes_with_other_scores(
                epitopes=neoantigen.neoepitopes_mhc_i,
                annotated_epitopes=mixmhcpred.results,
                annotation_name=MixMHCpred.ANNOTATION_PREFIX)

        # PRIME
        if prime is not None:
            neoantigen.neofox_annotations.annotations.extend(prime.get_annotations())
            neoantigen.neoepitopes_mhc_i = AnnotationFactory.annotate_epitopes_with_other_scores(
                epitopes=neoantigen.neoepitopes_mhc_i,
                annotated_epitopes=prime.results,
                annotation_name=Prime.ANNOTATION_PREFIX)

        # MixMHC2pred
        if mixmhc2pred is not None:
            neoantigen.neofox_annotations.annotations.extend(mixmhc2pred.get_annotations())
            neoantigen.neoepitopes_mhc_i_i = AnnotationFactory.annotate_epitopes_with_other_scores(
                epitopes=neoantigen.neoepitopes_mhc_i_i,
                annotated_epitopes=mixmhc2pred.results,
                annotation_name=MixMHC2pred.ANNOTATION_PREFIX)

        # MHC binding independent features
        expression_annotation = self.expression_calculator.get_annotations(neoantigen=neoantigen)
        neoantigen.neofox_annotations.annotations.extend(expression_annotation)

        sequence_not_in_uniprot = self.uniprot.is_sequence_not_in_uniprot(
            neoantigen.mutated_xmer
        )
        neoantigen.neofox_annotations.annotations.extend(
            self.uniprot.get_annotations(sequence_not_in_uniprot)
        )

        # Amplitude
        self.amplitude.run(netmhcpan=netmhcpan, netmhc2pan=netmhc2pan)
        neoantigen.neofox_annotations.annotations.extend(self.amplitude.get_annotations())
        neoantigen.neofox_annotations.annotations.extend(self.amplitude.get_annotations_mhc2())

        # Neoantigen fitness
        neoantigen.neofox_annotations.annotations.extend(
            self.neoantigen_fitness_calculator.get_annotations(
                mutated_peptide_mhci=netmhcpan.best_ninemer_epitope_by_affinity if netmhcpan else None,
                amplitude=self.amplitude.amplitude_mhci_affinity_9mer,
                mutated_peptide_mhcii=netmhc2pan.best_predicted_epitope_affinity if netmhc2pan else None
            )
        )
        neoantigen.neofox_annotations.annotations.extend(
            self.neoantigen_fitness_calculator.get_annotations_extended(
                mutated_peptide_mhci=netmhcpan.best_epitope_by_affinity if netmhcpan else None,
                amplitude=self.amplitude.amplitude_mhci_affinity
            )
        )

        # Differential Binding
        if netmhcpan:
            neoantigen.neofox_annotations.annotations.extend(
                self.differential_binding.get_annotations_dai(epitope=netmhcpan.best_epitope_by_affinity)
            )
            neoantigen.neofox_annotations.annotations.extend(
                self.differential_binding.get_annotations(mutated_peptide_mhci=netmhcpan.best_epitope_by_affinity,
                                                            amplitude=self.amplitude)
            )
        if netmhc2pan:
            neoantigen.neofox_annotations.annotations.extend(
                self.differential_binding.get_annotations_mhc2(mutated_peptide_mhcii=netmhc2pan.best_predicted_epitope_rank,
                                                               amplitude=self.amplitude)
            )

        # T cell predictor
        if netmhcpan:
            neoantigen.neofox_annotations.annotations.extend(
                self.tcell_predictor.get_annotations(
                    neoantigen=neoantigen, netmhcpan=netmhcpan
                )
            )

        # self-similarity
        neoantigen.neofox_annotations.annotations.extend(
            self.self_similarity.get_annnotations(
                epitope_mhci=netmhcpan.best_epitope_by_rank if netmhcpan else None,
                epitope_mhcii=netmhc2pan.best_predicted_epitope_affinity if netmhc2pan else None
            )
        )

        # number of mismatches and priority score
        if netmhcpan:
            neoantigen.neofox_annotations.annotations.extend(
                self.priority_score_calculator.get_annotations(
                    netmhcpan=netmhcpan,
                    neoantigen=neoantigen,
                    mut_not_in_prot=sequence_not_in_uniprot,
                )
            )

        # neoag immunogenicity model
        if netmhcpan and netmhcpan.best_epitope_by_affinity:
            neoantigen.neofox_annotations.annotations.append(
                self.neoag_calculator.get_annotation(
                    epitope_mhci=netmhcpan.best_epitope_by_affinity,
                    neoantigen=neoantigen)
            )

        # IEDB immunogenicity
        if self.organism == ORGANISM_HOMO_SAPIENS:
            neoantigen.neofox_annotations.annotations.extend(
                self.iedb_immunogenicity.get_annotations(
                    mutated_peptide_mhci=netmhcpan.best_epitope_by_affinity if netmhcpan else None,
                    mutated_peptide_mhcii=netmhc2pan.best_predicted_epitope_affinity if netmhc2pan else None
                )
            )

        # dissimilarity to self-proteome
        neoantigen.neofox_annotations.annotations.extend(
            self.dissimilarity_calculator.get_annotations(
                mutated_peptide_mhci=netmhcpan.best_epitope_by_affinity if netmhcpan else None,
                mutated_peptide_mhcii=netmhc2pan.best_predicted_epitope_affinity if netmhc2pan else None)
        )

        # vaxrank
        if netmhcpan and netmhcpan.predictions:
            neoantigen.neofox_annotations.annotations.extend(VaxRank().get_annotations(
                epitope_predictions=netmhcpan.predictions,
                expression_score=expression_annotation[0].value,
            ))

        # hex
        # TODO: hex is failing for mouse with the current IEDB fasta with only 2 entries
        if self.organism == ORGANISM_HOMO_SAPIENS:
            neoantigen.neofox_annotations.annotations.extend(
                self.hex.get_annotation(
                    mutated_peptide_mhci=netmhcpan.best_epitope_by_affinity if netmhcpan else None,
                    mutated_peptide_mhcii=netmhc2pan.best_predicted_epitope_affinity if netmhc2pan else None)
            )

        # annotate neoepitopes
        if with_all_neoepitopes:
            neoantigen.neoepitopes_mhc_i = [
                self.get_additional_annotations_neoepitope_mhci(
                    epitope=e, neoantigen=neoantigen)
                for e in neoantigen.neoepitopes_mhc_i]
            neoantigen.neoepitopes_mhc_i_i = [
                self.get_additional_annotations_neoepitope_mhcii(epitope=e) for e in neoantigen.neoepitopes_mhc_i_i]

        return neoantigen
