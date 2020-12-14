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
from distributed import get_client, secede, rejoin
import neofox
import time
from neofox.annotation_resources.uniprot.uniprot import Uniprot
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
from neofox.published_features.vaxrank import vaxrank
from neofox.published_features.iedb_immunogenicity.iedb import IEDBimmunogenicity
from neofox.published_features.expression import Expression
from neofox.published_features.priority_score import PriorityScore
from neofox.model.neoantigen import Patient, Neoantigen, NeoantigenAnnotations
from neofox.references.references import (
    ReferenceFolder,
    DependenciesConfiguration,
    AvailableAlleles,
)


class NeoantigenAnnotator:
    def __init__(
        self,
        references: ReferenceFolder,
        configuration: DependenciesConfiguration,
        tcell_predictor: TcellPrediction,
        self_similarity: SelfSimilarityCalculator,
    ):
        """class to annotate neoantigens"""
        self.runner = Runner()
        self.configuration = configuration
        self.available_alleles = references.get_available_alleles()
        self.tcell_predictor = tcell_predictor
        self.self_similarity = self_similarity

        # NOTE: this one loads a big file, but it is faster loading it multiple times than passing it around
        self.uniprot = Uniprot(references.uniprot)

        # NOTE: these resources do not read any file thus can be initialised fast
        self.dissimilarity_calculator = DissimilarityCalculator(
            runner=self.runner,
            configuration=configuration,
            proteome_db=references.proteome_db,
        )
        self.neoantigen_fitness_calculator = NeoantigenFitnessCalculator(
            runner=self.runner, configuration=configuration, iedb=references.iedb
        )
        self.neoag_calculator = NeoagCalculator(
            runner=self.runner, configuration=configuration
        )
        self.differential_binding = DifferentialBinding()
        self.priority_score_calculator = PriorityScore()
        self.iedb_immunogenicity = IEDBimmunogenicity()
        self.amplitude = Amplitude()

    def get_annotation(
        self, neoantigen: Neoantigen, patient: Patient
    ) -> NeoantigenAnnotations:
        """Calculate new epitope features and add to dictonary that stores all properties"""
        self._initialise_annotations(neoantigen)

        # Runs netmhcpan, netmhc2pan, mixmhcpred and mixmhc2prd in parallel
        (
            mixmhc2pred_annotations,
            mixmhcpred_annotations,
            netmhc2pan,
            netmhcpan,
        ) = self._compute_long_running_tasks(neoantigen, patient)

        # HLA I predictions: NetMHCpan
        if netmhcpan:
            self.annotations.annotations.extend(netmhcpan.get_annotations())

        # HLA II predictions: NetMHCIIpan
        if netmhc2pan:
            self.annotations.annotations.extend(netmhc2pan.get_annotations())

        # MixMHCpred
        if mixmhcpred_annotations is not None:
            self.annotations.annotations.extend(mixmhcpred_annotations)

        # MixMHC2pred
        if mixmhc2pred_annotations is not None:
            self.annotations.annotations.extend(mixmhc2pred_annotations)

        # decides which VAF to use
        vaf_rna = neoantigen.rna_variant_allele_frequency
        if not patient.is_rna_available:
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
        self.annotations.annotations.extend(expression_calculator.get_annotations())
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
        self.annotations.annotations.extend(
            self.uniprot.get_annotations(sequence_not_in_uniprot)
        )
        end = time.time()
        logger.info(
            "Uniprot annotation elapsed time {} seconds".format(round(end - start, 3))
        )

        # Amplitude
        start = time.time()
        self.amplitude.run(netmhcpan=netmhcpan, netmhc2pan=netmhc2pan)
        self.annotations.annotations.extend(self.amplitude.get_annotations())
        self.annotations.annotations.extend(self.amplitude.get_annotations_mhc2())
        end = time.time()
        logger.info(
            "Amplitude annotation elapsed time {} seconds".format(round(end - start, 3))
        )

        # Neoantigen fitness
        if netmhcpan:
            start = time.time()
            self.annotations.annotations.extend(
                self.neoantigen_fitness_calculator.get_annotations(
                    netmhcpan, self.amplitude
                )
            )
            end = time.time()
            logger.info(
                "Neoantigen annotation elapsed time {} seconds".format(
                    round(end - start, 3)
                )
            )

        # Differential Binding
        start = time.time()
        if netmhcpan:
            self.annotations.annotations.extend(
                self.differential_binding.get_annotations_dai(netmhcpan)
            )
            self.annotations.annotations.extend(
                self.differential_binding.get_annotations(netmhcpan, self.amplitude)
            )
        if netmhc2pan:
            self.annotations.annotations.extend(
                self.differential_binding.get_annotations_mhc2(netmhc2pan, self.amplitude)
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
            self.annotations.annotations.extend(
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
        if netmhcpan:
            start = time.time()
            self.annotations.annotations.extend(
                self.self_similarity.get_annnotations(netmhcpan=netmhcpan)
            )
            end = time.time()
            logger.info(
                "Self similarity annotation elapsed time {} seconds".format(
                    round(end - start, 3)
                )
            )

        # number of mismatches and priority score
        if netmhcpan:
            start = time.time()
            self.annotations.annotations.extend(
                self.priority_score_calculator.get_annotations(
                    netmhcpan=netmhcpan,
                    vaf_transcr=vaf_rna,
                    vaf_tum=neoantigen.dna_variant_allele_frequency,
                    expr=neoantigen.rna_expression,
                    mut_not_in_prot=sequence_not_in_uniprot,
                )
            )
            end = time.time()
            logger.info(
                "Priotity score annotation elapsed time {} seconds".format(
                    round(end - start, 3)
                )
            )

        # neoag immunogenicity model
        if netmhcpan:
            start = time.time()
            peptide_variant_position = EpitopeHelper.position_of_mutation_epitope(
                wild_type=netmhcpan.best_wt_epitope_by_affinity.peptide,
                mutation=netmhcpan.best_epitope_by_affinity.peptide,
            )
            self.annotations.annotations.append(
                self.neoag_calculator.get_annotation(
                    sample_id=patient.identifier,
                    netmhcpan=netmhcpan,
                    peptide_variant_position=peptide_variant_position,
                )
            )
            end = time.time()
            logger.info(
                "Neoag annotation elapsed time {} seconds".format(round(end - start, 3))
            )

        # IEDB immunogenicity
        if netmhcpan:
            start = time.time()
            self.annotations.annotations.extend(
                self.iedb_immunogenicity.get_annotations(
                    netmhcpan=netmhcpan, mhci_allele=netmhcpan.best_epitope_by_affinity.hla
                )
            )
            end = time.time()
            logger.info(
                "IEDB annotation elapsed time {} seconds".format(round(end - start, 3))
            )

        # dissimilarity to self-proteome
        if netmhcpan:
            start = time.time()
            self.annotations.annotations.extend(
                self.dissimilarity_calculator.get_annotations(netmhcpan=netmhcpan)
            )
            end = time.time()
            logger.info(
                "Dissimilarity annotation elapsed time {} seconds".format(
                    round(end - start, 3)
                )
            )

        # vaxrank
        if netmhcpan:
            start = time.time()
            vaxrankscore = vaxrank.VaxRank()
            vaxrankscore.run(
                mutation_scores=netmhcpan.epitope_affinities,
                expression_score=expression_calculator.expression,
            )
            self.annotations.annotations.extend(vaxrankscore.get_annotations())
            end = time.time()
            logger.info(
                "Vaxrank annotation elapsed time {} seconds".format(round(end - start, 3))
            )

        return self.annotations

    def _compute_long_running_tasks(self, neoantigen, patient, sequential=True):

        has_mhc1 = patient.mhc1 is not None and len(patient.mhc1) > 0
        has_mhc2 = patient.mhc2 is not None and len(patient.mhc2) > 0

        netmhcpan = None
        netmhc2pan = None
        mixmhcpred_annotations = None
        mixmhc2pred_annotations = None

        if sequential:
            if has_mhc1:
                netmhcpan = self.run_netmhcpan(
                    self.runner,
                    self.configuration,
                    self.available_alleles,
                    neoantigen,
                    patient)
            if has_mhc2:
                netmhc2pan = self.run_netmhc2pan(
                    self.runner,
                    self.configuration,
                    self.available_alleles,
                    neoantigen,
                    patient
                )
            if self.configuration.mix_mhc2_pred is not None and has_mhc2:
                mixmhc2pred_annotations = self.run_mixmhc2pred(
                    self.runner,
                    self.configuration,
                    neoantigen,
                    patient,
                )
            if self.configuration.mix_mhc_pred is not None and has_mhc1:
                mixmhcpred_annotations = self.run_mixmhcpred(
                    self.runner,
                    self.configuration,
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
                    self.configuration,
                    self.available_alleles,
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
                    neoantigen,
                    patient,
                )
            mixmhc2pred_future = None
            if self.configuration.mix_mhc2_pred is not None and has_mhc2:
                mixmhc2pred_future = dask_client.submit(
                    self.run_mixmhc2pred,
                    self.runner,
                    self.configuration,
                    neoantigen,
                    patient,
                )
            mixmhcpred_future = None
            if self.configuration.mix_mhc_pred is not None and has_mhc1:
                mixmhcpred_future = dask_client.submit(
                    self.run_mixmhcpred,
                    self.runner,
                    self.configuration,
                    neoantigen,
                    patient,
                )
            secede()

            if netmhcpan_future:
                netmhcpan = dask_client.gather([netmhcpan_future])[0]
            if netmhc2pan_future:
                netmhc2pan = dask_client.gather([netmhc2pan_future])[0]
            if mixmhcpred_future:
                mixmhcpred_annotations = dask_client.gather([mixmhcpred_future])[0]
            if mixmhc2pred_future:
                mixmhc2pred_annotations = dask_client.gather([mixmhc2pred_future])[0]
            rejoin()

        return mixmhc2pred_annotations, mixmhcpred_annotations, netmhc2pan, netmhcpan

    def _initialise_annotations(self, neoantigen):
        self.annotations = NeoantigenAnnotations()
        self.annotations.neoantigen_identifier = neoantigen.identifier
        self.annotations.annotator = "Neofox"
        self.annotations.annotator_version = neofox.VERSION
        self.annotations.timestamp = "{:%Y%m%d%H%M%S%f}".format(datetime.now())
        # TODO: set the hash fro the resources
        self.annotations.annotations = []

    @staticmethod
    def run_netmhcpan(
        runner: Runner,
        configuration: DependenciesConfiguration,
        available_alleles: AvailableAlleles,
        neoantigen: Neoantigen,
        patient: Patient,
    ):
        netmhcpan = BestAndMultipleBinder(runner=runner, configuration=configuration)
        netmhcpan.run(
            mutation=neoantigen.mutation,
            mhc1_alleles_patient=patient.mhc1,
            mhc1_alleles_available=available_alleles.get_available_mhc_i(),
        )
        return netmhcpan

    @staticmethod
    def run_netmhc2pan(
        runner: Runner,
        configuration: DependenciesConfiguration,
        available_alleles: AvailableAlleles,
        neoantigen: Neoantigen,
        patient: Patient,
    ):
        netmhc2pan = BestAndMultipleBinderMhcII(
            runner=runner, configuration=configuration
        )
        netmhc2pan.run(
            mutation=neoantigen.mutation,
            mhc2_alleles_patient=patient.mhc2,
            mhc2_alleles_available=available_alleles.get_available_mhc_ii(),
        )
        return netmhc2pan

    @staticmethod
    def run_mixmhcpred(
        runner: Runner,
        configuration: DependenciesConfiguration,
        neoantigen: Neoantigen,
        patient: Patient,
    ):
        mixmhc = MixMHCpred(runner, configuration)
        return mixmhc.get_annotations(mutation=neoantigen.mutation, mhc=patient.mhc1)

    @staticmethod
    def run_mixmhc2pred(
        runner: Runner,
        configuration: DependenciesConfiguration,
        neoantigen: Neoantigen,
        patient: Patient,
    ):
        mixmhc2 = MixMhc2Pred(runner, configuration)
        return mixmhc2.get_annotations(mhc=patient.mhc2, mutation=neoantigen.mutation)
