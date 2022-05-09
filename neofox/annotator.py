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
from distributed import get_client, secede, rejoin
import neofox
import time
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMhc2Pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.MHC_predictors.netmhcpan.combine_netmhcIIpan_pred_multiple_binders import (
    BestAndMultipleBinderMhcII,
)
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)
from neofox.model.factories import AnnotationFactory
from neofox.model.mhc_parser import MhcParser
from neofox.published_features.differential_binding.amplitude import Amplitude
from neofox.published_features.differential_binding.differential_binding import (
    DifferentialBinding,
)
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import (
    TcellPrediction,
)
from neofox.published_features.dissimilarity_garnish.dissimilaritycalculator import (
    DissimilarityCalculator,
)
from neofox.published_features.neoag.neoag_gbm_model import NeoagCalculator
from neofox.published_features.neoantigen_fitness.neoantigen_fitness import (
    NeoantigenFitnessCalculator,
)
from neofox.published_features.self_similarity.self_similarity import (
    SelfSimilarityCalculator,
)
from neofox.published_features.iedb_immunogenicity.iedb import IEDBimmunogenicity
from neofox.published_features.expression import Expression
from neofox.published_features.priority_score import PriorityScore
from neofox.MHC_predictors.prime import Prime
from neofox.published_features.hex.hex import Hex
from neofox.model.neoantigen import Patient, Neoantigen, Annotations, PredictedEpitope
from neofox.published_features.vaxrank.vaxrank import VaxRank
from neofox.references.references import (
    ReferenceFolder,
    DependenciesConfiguration,
    AvailableAlleles, ORGANISM_HOMO_SAPIENS
)


class NeoantigenAnnotator:
    def __init__(
        self,
        references: ReferenceFolder,
        configuration: DependenciesConfiguration,
        tcell_predictor: TcellPrediction,
        self_similarity: SelfSimilarityCalculator,
        affinity_threshold =neofox.AFFINITY_THRESHOLD_DEFAULT
    ):
        """class to annotate neoantigens"""
        self.runner = Runner()
        self.configuration = configuration
        self.proteome_db = references.proteome_db
        self.available_alleles = references.get_available_alleles()
        self.tcell_predictor = tcell_predictor
        self.self_similarity = self_similarity
        self.organism = references.organism

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
        self.dissimilarity_calculator = DissimilarityCalculator(
            proteome_blastp_runner=self.proteome_blastp_runner, affinity_threshold=affinity_threshold)
        self.neoantigen_fitness_calculator = NeoantigenFitnessCalculator(iedb_blastp_runner=self.iedb_blastp_runner)
        self.neoag_calculator = NeoagCalculator(
            runner=self.runner, configuration=configuration, affinity_threshold=affinity_threshold
        )
        self.differential_binding = DifferentialBinding(affinity_threshold=affinity_threshold)
        self.priority_score_calculator = PriorityScore()
        self.iedb_immunogenicity = IEDBimmunogenicity(affinity_threshold=affinity_threshold)
        self.amplitude = Amplitude()
        self.hex = Hex(runner=self.runner, configuration=configuration, references=references)
        self.mhc_database = references.get_mhc_database()
        self.mhc_parser = MhcParser.get_mhc_parser(self.mhc_database)

        self.resources_versions = references.get_resources_versions()

    def _annotate_epitopes_with_other_scores(
            self,
            epitopes: List[PredictedEpitope],
            annotated_epitopes: List[PredictedEpitope],
            annotation_name: str) -> List[PredictedEpitope]:

        merged_epitopes = []
        annotated_epitopes_dict = {EpitopeHelper.get_epitope_id(e): e for e in annotated_epitopes}
        for e in epitopes:

            # intialise annotations for the epitope if not done already
            if e.neofox_annotations is None:
                e.neofox_annotations = Annotations(annotations=[])

            # adds new annotations if any
            paired_epitope = annotated_epitopes_dict.get(EpitopeHelper.get_epitope_id(e))
            if paired_epitope is not None:
                if paired_epitope.affinity_score is not None:
                    e.neofox_annotations.annotations.append(
                        AnnotationFactory.build_annotation(
                            name=annotation_name + '_affinity_score', value=paired_epitope.affinity_score))
                if paired_epitope.rank is not None:
                    e.neofox_annotations.annotations.append(
                        AnnotationFactory.build_annotation(
                            name=annotation_name + '_rank', value=paired_epitope.rank))

            # updates epitope
            merged_epitopes.append(e)

        return merged_epitopes

    def get_annotation(self, neoantigen: Neoantigen, patient: Patient, with_all_neoepitopes=False) -> Neoantigen:
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
        ) = self._compute_long_running_tasks(neoantigen, patient)

        # HLA I predictions: NetMHCpan
        if netmhcpan:
            neoantigen.neofox_annotations.annotations.extend(netmhcpan.get_annotations())
            neoantigen.neoepitopes_mhc_i = netmhcpan.predictions
            if with_all_neoepitopes:
                for e in neoantigen.neoepitopes_mhc_i:
                    position = EpitopeHelper.position_of_mutation_epitope(epitope=e)
                    mutation_in_anchor = EpitopeHelper.position_in_anchor_position(
                        position_mhci=position, peptide_length=len(e.peptide),)
                    e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                        value=position,
                        name='position_mutation'))
                    e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                        value=mutation_in_anchor,
                        name='anchor_mutated'))

        # HLA II predictions: NetMHCIIpan
        if netmhc2pan:
            neoantigen.neofox_annotations.annotations.extend(netmhc2pan.get_annotations())
            neoantigen.neoepitopes_mhc_i_i = netmhc2pan.predictions

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
            # TODO: overwrite value in the neoantigen object
            vaf_rna = neoantigen.dna_variant_allele_frequency

        # MHC binding independent features
        start = time.time()
        expression_calculator = Expression(
            transcript_expression=neoantigen.rna_expression, vaf_rna=vaf_rna
        )
        neoantigen.neofox_annotations.annotations.extend(expression_calculator.get_annotations())
        end = time.time()
        logger.info(
            "Expression annotation elapsed time {} seconds".format(
                round(end - start, 3)
            )
        )

        start = time.time()
        sequence_not_in_uniprot = self.uniprot.is_sequence_not_in_uniprot(
            neoantigen.mutation.mutated_xmer
        )
        neoantigen.neofox_annotations.annotations.extend(
            self.uniprot.get_annotations(sequence_not_in_uniprot)
        )
        end = time.time()
        logger.info(
            "Uniprot annotation elapsed time {} seconds".format(round(end - start, 3))
        )

        # Amplitude
        start = time.time()
        self.amplitude.run(netmhcpan=netmhcpan, netmhc2pan=netmhc2pan)
        neoantigen.neofox_annotations.annotations.extend(self.amplitude.get_annotations())
        neoantigen.neofox_annotations.annotations.extend(self.amplitude.get_annotations_mhc2())
        if with_all_neoepitopes:
            for e in neoantigen.neoepitopes_mhc_i:
                e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                    value=self.amplitude.calculate_amplitude_mhc(
                        score_mutation=e.affinity_score, score_wild_type=e.affinity_score_wild_type,
                        apply_correction=True),
                    name='amplitude'))
            for e in neoantigen.neoepitopes_mhc_i_i:
                e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                    value=self.amplitude.calculate_amplitude_mhc(
                        score_mutation=e.rank, score_wild_type=e.rank_wild_type),
                    name='amplitude'))
        end = time.time()
        logger.info(
            "Amplitude annotation elapsed time {} seconds".format(round(end - start, 3))
        )

        # Neoantigen fitness
        start = time.time()
        neoantigen.neofox_annotations.annotations.extend(
            self.neoantigen_fitness_calculator.get_annotations(
                mutated_peptide_mhci=netmhcpan.best_ninemer_epitope_by_affinity if netmhcpan else None,
                mutation_in_anchor=netmhcpan.mutation_in_anchor_9mer if netmhcpan else None,
                amplitude=self.amplitude.amplitude_mhci_affinity_9mer,
                mutated_peptide_mhcii=netmhc2pan.best_predicted_epitope_affinity if netmhc2pan else None
            )
        )
        if with_all_neoepitopes:
            for e in neoantigen.neoepitopes_mhc_i:
                pathogen_similarity = self.neoantigen_fitness_calculator.get_pathogen_similarity(peptide=e.peptide)
                e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                    value=pathogen_similarity,
                    name='pathogen_similarity'))
                amplitude = self.amplitude.calculate_amplitude_mhc(
                    score_mutation=e.affinity_score, score_wild_type=e.affinity_score_wild_type,
                    apply_correction=True)
                # TODO: this is computed twice for each epitope, will it be quicker to fetch it from the annotations?
                position = EpitopeHelper.position_of_mutation_epitope(epitope=e)
                mutation_in_anchor = EpitopeHelper.position_in_anchor_position(
                    position_mhci=position, peptide_length=len(e.peptide), )
                recognition_potential = self.neoantigen_fitness_calculator.calculate_recognition_potential(
                    amplitude=amplitude, pathogen_similarity=pathogen_similarity,
                    mutation_in_anchor=mutation_in_anchor, mhc_affinity_mut=e.affinity_score
                )
                e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                    value=recognition_potential,
                    name='recognition_potential'))
            for e in neoantigen.neoepitopes_mhc_i_i:
                pathogen_similarity = self.neoantigen_fitness_calculator.get_pathogen_similarity(peptide=e.peptide)
                e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                    value=pathogen_similarity,
                    name='pathogen_similarity'))

        end = time.time()
        logger.info(
            "Neoantigen annotation elapsed time {} seconds".format(
                round(end - start, 3)
            )
        )

        # Differential Binding
        start = time.time()
        if netmhcpan:
            neoantigen.neofox_annotations.annotations.extend(
                self.differential_binding.get_annotations_dai(epitope=netmhcpan.best_epitope_by_affinity)
            )
            neoantigen.neofox_annotations.annotations.extend(
                self.differential_binding.get_annotations(mutated_peptide_mhci=netmhcpan.best_epitope_by_affinity,
                                                            amplitude=self.amplitude)
            )
            if with_all_neoepitopes:
                for e in neoantigen.neoepitopes_mhc_i:
                    e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                        value=self.differential_binding.dai(
                            score_mutation=e.affinity_score, score_wild_type=e.affinity_score_wild_type),
                        name='DAI'))
        if netmhc2pan:
            neoantigen.neofox_annotations.annotations.extend(
                self.differential_binding.get_annotations_mhc2(mutated_peptide_mhcii=netmhc2pan.best_predicted_epitope_rank,
                                                               amplitude=self.amplitude)
            )
        end = time.time()
        logger.info(
            "Differential binding annotation elapsed time {} seconds".format(
                round(end - start, 3)
            )
        )

        # T cell predictor
        if netmhcpan:
            start = time.time()
            neoantigen.neofox_annotations.annotations.extend(
                self.tcell_predictor.get_annotations(
                    neoantigen=neoantigen, netmhcpan=netmhcpan
                )
            )
            end = time.time()
            logger.info(
                "T-cell predictor annotation elapsed time {} seconds".format(
                    round(end - start, 3)
                )
            )

        # self-similarity
        start = time.time()
        neoantigen.neofox_annotations.annotations.extend(
            self.self_similarity.get_annnotations(
                epitope_mhci=netmhcpan.best_epitope_by_rank if netmhcpan else None,
                epitope_mhcii=netmhc2pan.best_predicted_epitope_affinity if netmhc2pan else None
            )
        )
        if with_all_neoepitopes:
            for e in neoantigen.neoepitopes_mhc_i:
                e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                    value=self.self_similarity.is_improved_binder(
                        score_mutation=e.rank, score_wild_type=e.rank_wild_type),
                    name='Improved_Binder_MHCI'))
                e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                    value=self.self_similarity.get_self_similarity(
                        mutated_peptide=e.peptide, wt_peptide=e.wild_type_peptide),
                    name='Selfsimilarity_MHCI'))

            for e in neoantigen.neoepitopes_mhc_i_i:
                e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                    value=self.self_similarity.get_self_similarity(
                        mutated_peptide=e.peptide, wt_peptide=e.wild_type_peptide),
                    name='Selfsimilarity_MHCII'))

        end = time.time()
        logger.info(
            "Self similarity annotation elapsed time {} seconds".format(
                round(end - start, 3)
            )
        )

        # number of mismatches and priority score
        if netmhcpan:
            start = time.time()
            neoantigen.neofox_annotations.annotations.extend(
                self.priority_score_calculator.get_annotations(
                    netmhcpan=netmhcpan,
                    vaf_transcr=vaf_rna,
                    vaf_tum=neoantigen.dna_variant_allele_frequency,
                    expr=neoantigen.rna_expression,
                    mut_not_in_prot=sequence_not_in_uniprot,
                )
            )
            if with_all_neoepitopes:
                for e in neoantigen.neoepitopes_mhc_i:

                    mutation_not_in_proteome = self.uniprot.is_sequence_not_in_uniprot(e.peptide)
                    e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                        value=mutation_not_in_proteome,
                        name='mutation_not_found_in_proteome'))

                    num_mismatches = EpitopeHelper.number_of_mismatches(
                        epitope_wild_type=e.wild_type_peptide, epitope_mutation=e.peptide, )
                    e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                        value=num_mismatches,
                        name='number_of_mismatches'))

                    e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                        value=self.priority_score_calculator.calc_priority_score(
                            vaf_tumor=neoantigen.dna_variant_allele_frequency,
                            vaf_rna=vaf_rna,
                            transcript_expr=neoantigen.rna_expression,
                            no_mismatch=num_mismatches,
                            score_mut=e.rank,
                            score_wt=e.rank_wild_type,
                            mut_not_in_prot=mutation_not_in_proteome),
                        name='Priority_score'))

                for e in neoantigen.neoepitopes_mhc_i_i:
                    e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                        value=self.uniprot.is_sequence_not_in_uniprot(e.peptide),
                        name='mutation_not_found_in_proteome'))
            end = time.time()
            logger.info(
                "Priotity score annotation elapsed time {} seconds".format(
                    round(end - start, 3)
                )
            )

        # neoag immunogenicity model
        if netmhcpan and netmhcpan.best_epitope_by_affinity:
            start = time.time()
            neoantigen.neofox_annotations.annotations.append(
                self.neoag_calculator.get_annotation(
                    epitope_mhci=netmhcpan.best_epitope_by_affinity,
                    mutation=neoantigen.mutation)
            )
            if with_all_neoepitopes:
                for e in neoantigen.neoepitopes_mhc_i:
                    e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                        value=self.neoag_calculator.calculate_neoag_score(epitope=e),
                        name='neoag_immunogenicity'))
            end = time.time()
            logger.info(
                "Neoag annotation elapsed time {} seconds".format(round(end - start, 3))
            )

        # IEDB immunogenicity
        if self.organism == ORGANISM_HOMO_SAPIENS:
            start = time.time()
            neoantigen.neofox_annotations.annotations.extend(
                self.iedb_immunogenicity.get_annotations(
                    mutated_peptide_mhci=netmhcpan.best_epitope_by_affinity if netmhcpan else None,
                    mutated_peptide_mhcii=netmhc2pan.best_predicted_epitope_affinity if netmhc2pan else None
                )
            )

            if with_all_neoepitopes:
                for e in neoantigen.neoepitopes_mhc_i:
                    iedb = self.iedb_immunogenicity.calculate_iedb_immunogenicity(
                        peptide=e.peptide, mhc_allele=e.hla, mhc_score=e.affinity_score)
                    e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                        value=iedb,
                        name='IEDB_Immunogenicity_MHCI'))

                for e in neoantigen.neoepitopes_mhc_i_i:
                    iedb = self.iedb_immunogenicity.calculate_iedb_immunogenicity(
                        peptide=e.peptide, mhc_allele=e.hla, mhc_score=e.affinity_score)
                    e.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
                        value=iedb,
                        name='IEDB_Immunogenicity_MHCII'))

            end = time.time()
            logger.info(
                "IEDB annotation elapsed time {} seconds".format(round(end - start, 3))
            )

        # dissimilarity to self-proteome
        start = time.time()
        neoantigen.neofox_annotations.annotations.extend(
            self.dissimilarity_calculator.get_annotations(
                mutated_peptide_mhci=netmhcpan.best_epitope_by_affinity if netmhcpan else None,
                mutated_peptide_mhcii=netmhc2pan.best_predicted_epitope_affinity if netmhc2pan else None)
        )
        end = time.time()
        logger.info(
            "Dissimilarity annotation elapsed time {} seconds".format(
                round(end - start, 3)
            )
        )

        # vaxrank
        if netmhcpan and netmhcpan.predictions:
            start = time.time()
            neoantigen.neofox_annotations.annotations.extend(VaxRank().get_annotations(
                epitope_predictions=netmhcpan.predictions,
                expression_score=expression_calculator.expression,
            ))
            end = time.time()
            logger.info(
                "Vaxrank annotation elapsed time {} seconds".format(round(end - start, 3))
            )

        # hex
        # TODO: hex is failing for mouse with the current IEDB fasta with only 2 entries
        if self.organism == ORGANISM_HOMO_SAPIENS:
            start = time.time()
            neoantigen.neofox_annotations.annotations.extend(
                self.hex.get_annotation(
                    mutated_peptide_mhci=netmhcpan.best_epitope_by_affinity if netmhcpan else None,
                    mutated_peptide_mhcii=netmhc2pan.best_predicted_epitope_affinity if netmhc2pan else None)
            )
            end = time.time()
            logger.info(
                "Hex annotation elapsed time {} seconds".format(round(end - start, 3))
            )

        return neoantigen

    def _compute_long_running_tasks(self, neoantigen, patient, sequential=True):

        has_mhc1 = patient.mhc1 is not None and len(patient.mhc1) > 0
        has_mhc2 = patient.mhc2 is not None and len(patient.mhc2) > 0

        netmhcpan = None
        netmhc2pan = None
        mixmhcpred = None
        mixmhc2pred = None
        prime = None

        if sequential:
            if has_mhc1:
                netmhcpan = self.run_netmhcpan(
                    self.runner,
                    self.configuration,
                    self.available_alleles,
                    self.mhc_parser,
                    neoantigen,
                    patient)
            if has_mhc2:
                netmhc2pan = self.run_netmhc2pan(
                    self.runner,
                    self.configuration,
                    self.available_alleles,
                    self.mhc_parser,
                    neoantigen,
                    patient
                )
            # avoids running MixMHCpred and PRIME for non human organisms
            if self.organism == ORGANISM_HOMO_SAPIENS:
                if self.configuration.mix_mhc2_pred is not None and has_mhc2:
                    mixmhc2pred = self.run_mixmhc2pred(
                        self.runner,
                        self.configuration,
                        self.mhc_parser,
                        neoantigen,
                        patient,
                    )
                if self.configuration.mix_mhc_pred is not None and has_mhc1:
                    mixmhcpred = self.run_mixmhcpred(
                        self.runner,
                        self.configuration,
                        self.mhc_parser,
                        neoantigen,
                        patient,
                    )
                if self.configuration.mix_mhc_pred is not None and has_mhc1:
                    prime = self.run_prime(
                        self.runner,
                        self.configuration,
                        self.mhc_parser,
                        neoantigen,
                        patient,
                    )
        else:
            dask_client = get_client()

            netmhcpan_future = None
            if has_mhc1:
                netmhcpan_future = dask_client.submit(
                    self.run_netmhcpan,
                    self.runner,
                    self.references,
                    self.configuration,
                    self.available_alleles,
                    self.mhc_parser,
                    neoantigen,
                    patient,
                )
            netmhc2pan_future = None
            if has_mhc2:
                netmhc2pan_future = dask_client.submit(
                    self.run_netmhc2pan,
                    self.runner,
                    self.configuration,
                    self.available_alleles,
                    self.mhc_parser,
                    neoantigen,
                    patient,
                )
            # avoids running MixMHCpred and PRIME for non human organisms
            mixmhc2pred_future = None
            mixmhcpred_future = None
            prime_future = None
            if self.organism == ORGANISM_HOMO_SAPIENS:
                if self.configuration.mix_mhc2_pred is not None and has_mhc2:
                    mixmhc2pred_future = dask_client.submit(
                        self.run_mixmhc2pred,
                        self.runner,
                        self.configuration,
                        self.mhc_parser,
                        neoantigen,
                        patient,
                    )
                if self.configuration.mix_mhc_pred is not None and has_mhc1:
                    mixmhcpred_future = dask_client.submit(
                        self.run_mixmhcpred,
                        self.runner,
                        self.configuration,
                        self.mhc_parser,
                        neoantigen,
                        patient,
                    )
                if self.configuration.mix_mhc_pred is not None and has_mhc1:
                    prime_future = dask_client.submit(
                        self.run_prime,
                        self.runner,
                        self.configuration,
                        self.mhc_parser,
                        neoantigen,
                        patient,
                    )

            secede()

            if netmhcpan_future:
                netmhcpan = dask_client.gather([netmhcpan_future])[0]
            if netmhc2pan_future:
                netmhc2pan = dask_client.gather([netmhc2pan_future])[0]

            if self.organism == ORGANISM_HOMO_SAPIENS:
                if mixmhcpred_future:
                    mixmhcpred = dask_client.gather([mixmhcpred_future])[0]
                if mixmhc2pred_future:
                    mixmhc2pred = dask_client.gather([mixmhc2pred_future])[0]
                if prime_future:
                    prime = dask_client.gather([prime_future])[0]
            rejoin()

        return mixmhc2pred, mixmhcpred, netmhc2pan, netmhcpan, prime

    def run_netmhcpan(
            self,
            runner: Runner,
            configuration: DependenciesConfiguration,
            available_alleles: AvailableAlleles,
            mhc_parser: MhcParser,
            neoantigen: Neoantigen,
            patient: Patient,
    ):
        netmhcpan = BestAndMultipleBinder(runner=runner, configuration=configuration, mhc_parser=mhc_parser,
                                          blastp_runner=self.proteome_blastp_runner)
        netmhcpan.run(
            mutation=neoantigen.mutation,
            mhc1_alleles_patient=patient.mhc1,
            mhc1_alleles_available=available_alleles.get_available_mhc_i(),
            uniprot=self.uniprot,
        )
        return netmhcpan

    def run_netmhc2pan(
            self,
            runner: Runner,
            configuration: DependenciesConfiguration,
            available_alleles: AvailableAlleles,
            mhc_parser: MhcParser,
            neoantigen: Neoantigen,
            patient: Patient,
    ):
        netmhc2pan = BestAndMultipleBinderMhcII(
            runner=runner, configuration=configuration, mhc_parser=mhc_parser,
            blastp_runner=self.proteome_blastp_runner)
        netmhc2pan.run(
            mutation=neoantigen.mutation,
            mhc2_alleles_patient=patient.mhc2,
            mhc2_alleles_available=available_alleles.get_available_mhc_ii(),
            uniprot=self.uniprot
        )
        return netmhc2pan

    def run_mixmhcpred(
            self,
            runner: Runner,
            configuration: DependenciesConfiguration,
            mhc_parser: MhcParser,
            neoantigen: Neoantigen,
            patient: Patient,
    ):
        mixmhc = MixMHCpred(runner, configuration, mhc_parser)
        mixmhc.run(mutation=neoantigen.mutation, mhc=patient.mhc1, uniprot=self.uniprot)
        return mixmhc

    def run_prime(
            self,
            runner: Runner,
            configuration: DependenciesConfiguration,
            mhc_parser: MhcParser,
            neoantigen: Neoantigen,
            patient: Patient,
    ):
        prime = Prime(runner, configuration, mhc_parser)
        prime.run(mutation=neoantigen.mutation, mhc=patient.mhc1, uniprot=self.uniprot)
        return prime

    def run_mixmhc2pred(
            self,
            runner: Runner,
            configuration: DependenciesConfiguration,
            mhc_parser: MhcParser,
            neoantigen: Neoantigen,
            patient: Patient,
    ):
        mixmhc2 = MixMhc2Pred(runner, configuration, mhc_parser)
        mixmhc2.run(mhc=patient.mhc2, mutation=neoantigen.mutation, uniprot=self.uniprot)
        return mixmhc2
