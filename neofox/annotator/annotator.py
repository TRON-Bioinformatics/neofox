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

from logzero import logger
from datetime import datetime
import neofox
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.annotator.neoantigen_mhc_binding_annotator import NeoantigenMhcBindingAnnotator
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.model.factories import AnnotationFactory
from neofox.model.mhc_parser import MhcParser
from neofox.published_features.differential_binding.amplitude import Amplitude
from neofox.published_features.differential_binding.differential_binding import DifferentialBinding
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator
from neofox.published_features.neoag.neoag_gbm_model import NeoagCalculator
from neofox.published_features.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.published_features.iedb_immunogenicity.iedb import IEDBimmunogenicity
from neofox.published_features.expression import Expression
from neofox.published_features.priority_score import PriorityScore
from neofox.published_features.hex.hex import Hex
from neofox.model.neoantigen import Patient, Neoantigen, Annotations, PredictedEpitope
from neofox.published_features.vaxrank.vaxrank import VaxRank
from neofox.references.references import (
    ReferenceFolder,
    DependenciesConfiguration,
    ORGANISM_HOMO_SAPIENS
)


class NeoantigenAnnotator:
    def __init__(
        self,
        references: ReferenceFolder,
        configuration: DependenciesConfiguration,
        tcell_predictor: TcellPrediction,
        self_similarity: SelfSimilarityCalculator,
        rank_mhci_threshold=neofox.RANK_MHCI_THRESHOLD_DEFAULT,
        rank_mhcii_threshold=neofox.RANK_MHCII_THRESHOLD_DEFAULT
    ):
        """class to annotate neoantigens"""
        self.runner = Runner()
        self.configuration = configuration
        self.proteome_db = references.proteome_db
        self.available_alleles = references.get_available_alleles()
        self.tcell_predictor = tcell_predictor
        self.self_similarity = self_similarity
        self.organism = references.organism
        self.rank_mhci_threshold = rank_mhci_threshold
        self.rank_mhcii_threshold = rank_mhcii_threshold

        # NOTE: this one loads a big file, but it is faster loading it multiple times than passing it around
        self.uniprot = Uniprot(references.uniprot_pickle)

        # initialise proteome and IEDB BLASTP runners
        self.proteome_blastp_runner = BlastpRunner(
            runner=self.runner, configuration=configuration,
            database=references.get_proteome_database())
        self.iedb_blastp_runner = BlastpRunner(
            runner=self.runner, configuration=configuration,
            database=references.get_iedb_database())

        # NOTE: these resources do not read any file thus can be initialised fast
        self.dissimilarity_calculator = DissimilarityCalculator(proteome_blastp_runner=self.proteome_blastp_runner)
        self.neoantigen_fitness_calculator = NeoantigenFitnessCalculator(iedb_blastp_runner=self.iedb_blastp_runner)
        self.neoag_calculator = NeoagCalculator(runner=self.runner, configuration=configuration)
        self.differential_binding = DifferentialBinding()
        self.priority_score_calculator = PriorityScore()
        self.iedb_immunogenicity = IEDBimmunogenicity()
        self.amplitude = Amplitude()
        self.hex = Hex(runner=self.runner, configuration=configuration, references=references)
        self.mhc_database = references.get_mhc_database()
        self.mhc_parser = MhcParser.get_mhc_parser(self.mhc_database)

        self.neoantigen_mhc_binding_annotator = NeoantigenMhcBindingAnnotator(
            references=references, configuration=configuration, proteome_blastp_runner=self.proteome_blastp_runner,
            uniprot=self.uniprot)

        self.resources_versions = references.get_resources_versions()

    def  _annotate_epitopes_with_other_scores(
            self,
            epitopes: List[PredictedEpitope],
            annotated_epitopes: List[PredictedEpitope],
            annotation_name: str) -> List[PredictedEpitope]:

        merged_epitopes = []
        if annotated_epitopes is not None:
            annotated_epitopes_dict = {EpitopeHelper.get_epitope_id(e): e for e in annotated_epitopes}
            for e in epitopes:

                # intialise annotations for the epitope if not done already
                if e.neofox_annotations is None:
                    e.neofox_annotations = Annotations(annotations=[])

                # adds new annotations if any
                paired_epitope = annotated_epitopes_dict.get(EpitopeHelper.get_epitope_id(e))
                if paired_epitope is not None:
                    if paired_epitope.affinity_mutated is not None:
                        e.neofox_annotations.annotations.append(
                            AnnotationFactory.build_annotation(
                                name=annotation_name + '_affinity_score', value=paired_epitope.affinity_mutated))
                    if paired_epitope.rank_mutated is not None:
                        e.neofox_annotations.annotations.append(
                            AnnotationFactory.build_annotation(
                                name=annotation_name + '_rank', value=paired_epitope.rank_mutated))

                # updates epitope
                merged_epitopes.append(e)
        else:
            # if there are no results to annotate with it returns the input list as is
            merged_epitopes = epitopes

        return merged_epitopes

    def get_additional_annotations_neoepitope_mhci(
            self, epitope: PredictedEpitope, neoantigen: Neoantigen, vaf_rna) -> PredictedEpitope:

        epitope.neofox_annotations.annotations.extend(
            BestAndMultipleBinder.get_annotations_epitope_mhci(epitope=epitope) +
            self.amplitude.get_annotations_epitope_mhci(epitope=epitope)
        )

        # NOTE: this extend() call cannot be joined with the previous as some of the previous annotations are expected
        epitope.neofox_annotations.annotations.extend(
            self.neoantigen_fitness_calculator.get_annotations_epitope_mhci(epitope=epitope) +
            self.differential_binding.get_annotations_epitope_mhci(epitope=epitope) +
            self.self_similarity.get_annotations_epitope_mhci(epitope=epitope) +
            self.uniprot.get_annotations_epitope(epitope=epitope) +
            self.dissimilarity_calculator.get_annotations_epitope(epitope=epitope)
        )

        # restricted to 9-mers
        if len(epitope.mutated_peptide) == 9:
            epitope.neofox_annotations.annotations.extend(self.tcell_predictor.get_annotations_epitope_mhci(
                epitope=epitope, neoantigen=neoantigen))

        num_mismatches = EpitopeHelper.number_of_mismatches(
            epitope_wild_type=epitope.wild_type_peptide, epitope_mutation=epitope.mutated_peptide, )
        epitope.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
            value=num_mismatches,
            name='number_of_mismatches'))

        epitope.neofox_annotations.annotations.extend(
            self.priority_score_calculator.get_annotations_epitope_mhci(
                epitope=epitope, neoantigen=neoantigen, vaf_rna=vaf_rna))

        if self.organism == ORGANISM_HOMO_SAPIENS:
            epitope.neofox_annotations.annotations.extend(
                self.iedb_immunogenicity.get_annotations_epitope_mhci(epitope=epitope) +
                self.hex.get_annotations_epitope(epitope=epitope))

        return epitope

    def get_additional_annotations_neoepitope_mhcii(self, epitope: PredictedEpitope) -> PredictedEpitope:

        epitope.neofox_annotations.annotations.extend(
            self.amplitude.get_annotations_epitope_mhcii(epitope=epitope) +
            self.neoantigen_fitness_calculator.get_annotations_epitope_mhcii(epitope=epitope) +
            self.self_similarity.get_annotations_epitope_mhcii(epitope=epitope) +
            self.uniprot.get_annotations_epitope(epitope=epitope) +
            self.dissimilarity_calculator.get_annotations_epitope(epitope=epitope))

        if self.organism == ORGANISM_HOMO_SAPIENS:

            epitope.neofox_annotations.annotations.extend(
                self.iedb_immunogenicity.get_annotations_epitope_mhcii(epitope=epitope) +
                self.hex.get_annotations_epitope(epitope=epitope))

        return epitope

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
            neoantigen.neoepitopes_mhc_i = self._annotate_epitopes_with_other_scores(
                epitopes=neoantigen.neoepitopes_mhc_i,
                annotated_epitopes=mixmhcpred.results,
                annotation_name='MixMHCpred')

        # PRIME
        if prime is not None:
            neoantigen.neofox_annotations.annotations.extend(prime.get_annotations())
            neoantigen.neoepitopes_mhc_i = self._annotate_epitopes_with_other_scores(
                epitopes=neoantigen.neoepitopes_mhc_i,
                annotated_epitopes=prime.results,
                annotation_name='PRIME')

        # MixMHC2pred
        if mixmhc2pred is not None:
            neoantigen.neofox_annotations.annotations.extend(mixmhc2pred.get_annotations())
            neoantigen.neoepitopes_mhc_i_i = self._annotate_epitopes_with_other_scores(
                epitopes=neoantigen.neoepitopes_mhc_i_i,
                annotated_epitopes=mixmhc2pred.results,
                annotation_name='MixMHC2pred')

        # decides which VAF to use
        vaf_rna = neoantigen.rna_variant_allele_frequency
        if not patient.is_rna_available and neoantigen.dna_variant_allele_frequency is not None:
            logger.warning(
                "Using the DNA VAF to estimate the RNA VAF as the patient does not have RNA available"
            )
            # NOTE: does not overwrite value in the neoantigen object
            vaf_rna = neoantigen.dna_variant_allele_frequency

        # MHC binding independent features
        expression_calculator = Expression(
            transcript_expression=neoantigen.rna_expression, vaf_rna=vaf_rna
        )
        neoantigen.neofox_annotations.annotations.extend(expression_calculator.get_annotations())

        sequence_not_in_uniprot = self.uniprot.is_sequence_not_in_uniprot(
            neoantigen.mutation.mutated_xmer
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
                    vaf_transcr=vaf_rna,
                    vaf_tum=neoantigen.dna_variant_allele_frequency,
                    expr=neoantigen.rna_expression,
                    mut_not_in_prot=sequence_not_in_uniprot,
                )
            )

        # neoag immunogenicity model
        if netmhcpan and netmhcpan.best_epitope_by_affinity:
            neoantigen.neofox_annotations.annotations.append(
                self.neoag_calculator.get_annotation(
                    epitope_mhci=netmhcpan.best_epitope_by_affinity,
                    mutation=neoantigen.mutation)
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
                expression_score=expression_calculator.expression,
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
                self.get_additional_annotations_neoepitope_mhci(epitope=e, neoantigen=neoantigen, vaf_rna=vaf_rna)
                for e in neoantigen.neoepitopes_mhc_i]
            neoantigen.neoepitopes_mhc_i_i = [
                self.get_additional_annotations_neoepitope_mhcii(epitope=e) for e in neoantigen.neoepitopes_mhc_i_i]

        return neoantigen
