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
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.helpers.available_alleles import AvailableAlleles
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMhc2Pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.MHC_predictors.netmhcpan.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.published_features.differential_binding.amplitude import Amplitude
from neofox.published_features.differential_binding.differential_binding import DifferentialBinding
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator
from neofox.published_features.neoag.neoag_gbm_model import NeoagCalculator
from neofox.published_features.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.published_features.vaxrank import vaxrank
from neofox.published_features.iedb_immunogenicity.iedb import IEDBimmunogenicity
from neofox.published_features.expression import Expression
from neofox.published_features.priority_score import PriorityScore
from neofox.model.neoantigen import Patient, Neoantigen, NeoantigenAnnotations
from neofox.references.references import ReferenceFolder, DependenciesConfiguration


class NeoantigenAnnotator:

    def __init__(self):
        """class to annotate neoantigens"""
        references = ReferenceFolder()
        configuration = DependenciesConfiguration()
        runner = Runner()
        self.dissimilarity_calculator = DissimilarityCalculator(
            runner=runner, configuration=configuration, proteome_db=references.proteome_db)
        self.neoantigen_fitness_calculator = NeoantigenFitnessCalculator(
            runner=runner, configuration=configuration, iedb=references.iedb)
        self.neoag_calculator = NeoagCalculator(runner=runner, configuration=configuration)
        self.netmhc2pan = BestAndMultipleBinderMhcII(runner=runner, configuration=configuration)
        self.mixmhc2 = MixMhc2Pred(runner=runner, configuration=configuration)
        self.netmhcpan = BestAndMultipleBinder(runner=runner, configuration=configuration)
        self.mixmhc = MixMHCpred(runner=runner, configuration=configuration)
        self.available_alleles = AvailableAlleles(references)
        self.uniprot = Uniprot(references.uniprot)
        self.tcell_predictor = TcellPrediction()
        self.self_similarity = SelfSimilarityCalculator()

        # NOTE: these resources do not read any file thus can be initialised fast
        self.differential_binding = DifferentialBinding()
        self.priority_score_calculator = PriorityScore()
        self.iedb_immunogenicity = IEDBimmunogenicity()
        self.differential_binding = DifferentialBinding()
        self.amplitude = Amplitude()

    def get_annotation(self, neoantigen: Neoantigen, patient: Patient) -> NeoantigenAnnotations:
        """Calculate new epitope features and add to dictonary that stores all properties"""

        self._initialise_annotations(neoantigen)

        # decides which VAF to use
        vaf_rna = neoantigen.rna_variant_allele_frequency
        if not patient.is_rna_available:
            logger.warning("Using the DNA VAF to estimate the RNA VAF as the patient does not have RNA available")
            # TODO: overwrite value in the neoantigen object
            vaf_rna = neoantigen.dna_variant_allele_frequency

        # TODO: this is needed by the T cell predictor, move this construction inside by passing the neoantigen
        substitution = "{}{}{}".format(
            neoantigen.mutation.wild_type_aminoacid, neoantigen.mutation.position,
            neoantigen.mutation.mutated_aminoacid)

        # MHC binding independent features
        expression_calculator = Expression(
            transcript_expression=neoantigen.rna_expression, vaf_rna=vaf_rna)
        self.annotations.annotations.extend(expression_calculator.get_annotations())
        sequence_not_in_uniprot = self.uniprot.is_sequence_not_in_uniprot(neoantigen.mutation.mutated_xmer)
        self.annotations.annotations.extend(self.uniprot.get_annotations(sequence_not_in_uniprot))

        # HLA I predictions: NetMHCpan
        self.netmhcpan.run(
            sequence_mut=neoantigen.mutation.mutated_xmer, sequence_wt=neoantigen.mutation.wild_type_xmer,
            alleles=patient.mhc_i_alleles, set_available_mhc=self.available_alleles.get_available_mhc_i())
        self.annotations.annotations.extend(self.netmhcpan.get_annotations())

        # HLA II predictions: NetMHCIIpan
        self.netmhc2pan.run(
            sequence=neoantigen.mutation.mutated_xmer, sequence_reference=neoantigen.mutation.wild_type_xmer,
            alleles=patient.mhc_i_i_alleles, set_available_mhc=self.available_alleles.get_available_mhc_ii())
        self.annotations.annotations.extend(self.netmhc2pan.get_annotations())

        # Amplitude
        self.amplitude.run(netmhcpan=self.netmhcpan, netmhc2pan=self.netmhc2pan)
        self.annotations.annotations.extend(self.amplitude.get_annotations())

        # Neoantigen fitness
        self.annotations.annotations.extend(
            self.neoantigen_fitness_calculator.get_annotations(self.netmhcpan, self.amplitude))

        # Differential Binding
        self.annotations.annotations.extend(self.differential_binding.get_annotations_dai(self.netmhcpan))
        self.annotations.annotations.extend(self.differential_binding.get_annotations(self.netmhcpan, self.amplitude))
        self.annotations.annotations.extend(
            self.differential_binding.get_annotations_mhc2(self.netmhc2pan, self.amplitude))

        # T cell predictor
        self.annotations.annotations.extend(self.tcell_predictor.get_annotations(
            gene=neoantigen.gene.gene, substitution=substitution, netmhcpan=self.netmhcpan))

        # self-similarity
        self.annotations.annotations.extend(self.self_similarity.get_annnotations(
            netmhcpan=self.netmhcpan))

        # number of mismatches and priority score
        self.annotations.annotations.extend(self.priority_score_calculator.get_annotations(
            netmhcpan=self.netmhcpan, vaf_transcr=vaf_rna,
            vaf_tum=neoantigen.dna_variant_allele_frequency,
            expr=neoantigen.rna_expression, mut_not_in_prot=sequence_not_in_uniprot))

        # neoag immunogenicity model
        peptide_variant_position = EpitopeHelper.position_of_mutation_epitope(
            wild_type=self.netmhcpan.best4_affinity_epitope_WT, mutation=self.netmhcpan.best4_affinity_epitope)
        self.annotations.annotations.append(self.neoag_calculator.get_annotation(
            sample_id=patient.identifier, netmhcpan=self.netmhcpan, peptide_variant_position=peptide_variant_position))

        # IEDB immunogenicity
        self.annotations.annotations.extend(self.iedb_immunogenicity.get_annotations(
            netmhcpan=self.netmhcpan,
            mhci_allele=self.netmhcpan.best4_affinity_allele))

        # MixMHCpred
        self.mixmhc.run(
            sequence_wt=neoantigen.mutation.wild_type_xmer, sequence_mut=neoantigen.mutation.mutated_xmer,
            alleles=patient.mhc_i_alleles)
        self.annotations.annotations.extend(self.mixmhc.get_annotations())

        # MixMHC2pred
        self.mixmhc2.run(
            alleles=patient.mhc_i_i_alleles, sequence_wt=neoantigen.mutation.wild_type_xmer,
            sequence_mut=neoantigen.mutation.mutated_xmer)
        self.annotations.annotations.extend(self.mixmhc2.get_annotations())

        # dissimilarity to self-proteome
        self.annotations.annotations.extend(self.dissimilarity_calculator.get_annotations(
            netmhcpan=self.netmhcpan))

        # vaxrank
        vaxrankscore = vaxrank.VaxRank()
        vaxrankscore.run(mutation_scores=self.netmhcpan.epitope_affinities,
                         expression_score=expression_calculator.expression)
        self.annotations.annotations.extend(vaxrankscore.get_annotations())

        return self.annotations

    def _initialise_annotations(self, neoantigen):
        self.annotations = NeoantigenAnnotations()
        self.annotations.neoantigen_identifier = neoantigen.identifier
        self.annotations.annotator = "Neofox"
        self.annotations.annotator_version = neofox.VERSION
        self.annotations.timestamp = "{:%Y%m%d%H%M%S%f}".format(datetime.now())
        # TODO: set the hash fro the resources
        self.annotations.annotations = []
