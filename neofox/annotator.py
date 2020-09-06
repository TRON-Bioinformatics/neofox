#!/usr/bin/env python

from logzero import logger
import re
from datetime import datetime
import neofox
from neofox.aa_index.aa_index import AminoacidIndex
from neofox.annotation_resources.gtex.gtex import GTEx
from neofox.annotation_resources.nmer_frequency.nmer_frequency import AminoacidFrequency, FourmerFrequency
from neofox.annotation_resources.provean.provean import ProveanAnnotator
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.helpers.available_alleles import AvailableAlleles
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.literature_features.differential_expression import DifferentialExpression
from neofox.literature_features.differential_binding.amplitude import Amplitude
from neofox.literature_features.differential_binding.differential_binding import DifferentialBinding
from neofox.model.conversion import ModelConverter
from neofox.model.wrappers import AnnotationFactory
from neofox.predictors.MixMHCpred.mixmhc2pred import MixMhc2Pred
from neofox.predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.predictors.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.predictors.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator
from neofox.predictors.neoag.neoag_gbm_model import NeoagCalculator
from neofox.predictors.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from neofox.predictors.netmhcpan4.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.predictors.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.references.references import ReferenceFolder, DependenciesConfiguration
from neofox.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.vaxrank import vaxrank
from neofox.predictors.iedb.iedb import IEDBimmunogenicity
from neofox.literature_features.expression import Expression
from neofox.literature_features.priority_score import PriorityScore
from neofox.model.neoantigen import Patient, Neoantigen, NeoantigenAnnotations


class NeoantigenAnnotator:

    def __init__(self, uniprot: Uniprot, gtex: GTEx):
        """class to annotate neoantigens"""
        references = ReferenceFolder()
        configuration = DependenciesConfiguration()
        runner = Runner()
        self.provean_annotator = ProveanAnnotator(provean_file=references.prov_scores_mapped3)
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
        self.uniprot = uniprot
        self.gtex = gtex
        self.aa_frequency = AminoacidFrequency()
        self.fourmer_frequency = FourmerFrequency()
        self.aa_index = AminoacidIndex()
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
            transcript_expression=neoantigen.rna_expression, vaf_rna=vaf_rna,
            tumor_content=patient.estimated_tumor_content)
        self.annotations.annotations.extend(expression_calculator.get_annotations())
        self.add_differential_expression_features(
            neoantigen.gene.gene, expression_tumor=neoantigen.rna_expression, tissue=patient.tissue)
        self.add_aminoacid_index_features(
            mutation_aminoacid=neoantigen.mutation.mutated_aminoacid,
            wild_type_aminoacid=neoantigen.mutation.wild_type_aminoacid)
        self.add_provean_score_features(neoantigen)
        sequence_not_in_uniprot = self.uniprot.is_sequence_not_in_uniprot(neoantigen.mutation.mutated_xmer)
        self.annotations.annotations.extend(self.uniprot.get_annotations(sequence_not_in_uniprot))

        # HLA I predictions: NetMHCpan
        self.netmhcpan.run(
            xmer_mut=neoantigen.mutation.mutated_xmer, xmer_wt=neoantigen.mutation.wild_type_xmer,
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
        self.annotations.annotations.extend(self.amplitude.get_annotations_mhc2())

        # Neoantigen fitness metrics
        self.annotations.annotations.extend(
            self.neoantigen_fitness_calculator.get_annotations(self.netmhcpan, self.amplitude))
        self.annotations.annotations.extend(
            self.neoantigen_fitness_calculator.get_annotations_mhc2(self.netmhc2pan, self.amplitude))

        # Differential Binding
        self.annotations.annotations.extend(self.differential_binding.get_annotations(self.netmhcpan, self.amplitude))
        self.annotations.annotations.extend(
            self.differential_binding.get_annotations_mhc2(self.netmhc2pan, self.amplitude))

        # T cell predictor
        self.annotations.annotations.extend(self.tcell_predictor.get_annotations(
            gene=neoantigen.gene.gene, substitution=substitution, netmhcpan=self.netmhcpan))

        # frequencies
        self.annotations.annotations.extend(self.aa_frequency.get_annotations(
            neoantigen.mutation.mutated_aminoacid, self.netmhcpan.best4_mhc_epitope))
        self.annotations.annotations.extend(self.fourmer_frequency.get_annotations(self.netmhcpan.best4_mhc_epitope))

        # self-similarity
        self.annotations.annotations.extend(self.self_similarity.get_annnotations(
            netmhcpan=self.netmhcpan, netmhcpan2=self.netmhc2pan))

        # number of mismatches and priority scores
        self.annotations.annotations.extend(self.priority_score_calculator.get_annotations(
            netmhcpan=self.netmhcpan, netmhcpan2=self.netmhc2pan, vaf_transcr=vaf_rna,
            vaf_tum=neoantigen.dna_variant_allele_frequency,
            expr=neoantigen.rna_expression, mut_not_in_prot=sequence_not_in_uniprot))

        # neoag immunogenicity model
        peptide_variant_position = EpitopeHelper.position_of_mutation_epitope(
            wild_type=self.netmhcpan.best4_affinity_epitope_WT, mutation=self.netmhcpan.best4_affinity_epitope)
        self.annotations.annotations.append(self.neoag_calculator.get_annotation(
            sample_id=patient.identifier, netmhcpan=self.netmhcpan, peptide_variant_position=peptide_variant_position))

        # IEDB immunogenicity
        self.annotations.annotations.extend(self.iedb_immunogenicity.get_annotations(
            netmhcpan=self.netmhcpan, netmhcpan2=self.netmhc2pan,
            mhci_allele=self.netmhcpan.best4_affinity_allele,
            mhcii_allele=self.netmhc2pan.best_mhcII_pan_allele))

        # MixMHCpred
        self.mixmhc.run(
            xmer_wt=neoantigen.mutation.wild_type_xmer, xmer_mut=neoantigen.mutation.mutated_xmer,
            alleles=patient.mhc_i_alleles)
        self.annotations.annotations.extend(self.mixmhc.get_annotations())

        # MixMHC2pred
        self.mixmhc2.run(
            alleles=patient.mhc_i_i_alleles, xmer_wt=neoantigen.mutation.wild_type_xmer,
            xmer_mut=neoantigen.mutation.mutated_xmer)
        self.annotations.annotations.extend(self.mixmhc2.get_annotations())

        # dissimilarity to self-proteome
        self.annotations.annotations.extend(self.dissimilarity_calculator.get_annotations(
            netmhcpan=self.netmhcpan, netmhcpan2=self.netmhc2pan))

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

    def add_provean_score_features(self, neoantigen):
        # PROVEAN score
        transcript_identifier_without_version = re.sub(r'.\d+$', '', neoantigen.gene.transcript_identifier)
        self.annotations.annotations.append(self.provean_annotator.get_provean_annotation(
            mutated_aminoacid=neoantigen.mutation.mutated_aminoacid,
            protein_id=transcript_identifier_without_version,
            position=neoantigen.mutation.position))

    def add_aminoacid_index_features(self, mutation_aminoacid, wild_type_aminoacid):
        # amino acid index
        self.annotations.annotations.extend(self.aa_index.get_annotations(
            wild_type_aminoacid=wild_type_aminoacid, mutation_aminoacid=mutation_aminoacid))

    def add_differential_expression_features(self, gene, expression_tumor, tissue):
        # differential expression
        gtex_mean, gtex_sum, gtex_sd = self.gtex.get_metrics(gene, tissue)
        self.annotations.annotations.extend(self.gtex.get_annotations(gtex_mean, gtex_sd, gtex_sum))
        self.annotations.annotations.extend(DifferentialExpression().get_annotations(
            expression_tumor=expression_tumor, expression_reference=gtex_mean,
            expression_reference_sd=gtex_sd, expression_reference_sum=gtex_sum))
