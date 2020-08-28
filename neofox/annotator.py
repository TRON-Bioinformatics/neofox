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
from neofox.literature_features.differential_expression import DifferentialExpression
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
from neofox.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.vaxrank import vaxrank
from neofox.predictors.iedb.iedb import IEDBimmunogenicity
from neofox.literature_features.differential_binding import DifferentialBinding
from neofox.literature_features.expression import Expression
from neofox.literature_features.priority_score import PriorityScore
from neofox.model.neoantigen import Patient, Neoantigen, NeoantigenAnnotations


class NeoantigenAnnotator:

    def __init__(self,
                 provean_annotator: ProveanAnnotator,
                 uniprot: Uniprot,
                 dissimilarity_calculator: DissimilarityCalculator,
                 neoantigen_fitness_calculator: NeoantigenFitnessCalculator,
                 neoag_calculator: NeoagCalculator,
                 netmhcpan2: BestAndMultipleBinderMhcII,
                 mixmhc2: MixMhc2Pred,
                 netmhcpan: BestAndMultipleBinder,
                 mixmhc: MixMHCpred,
                 available_alleles: AvailableAlleles):
        """class to annotate neoantigens"""
        self.provean_annotator = provean_annotator
        self.dissimilarity_calculator = dissimilarity_calculator
        self.neoantigen_fitness_calculator = neoantigen_fitness_calculator
        self.neoag_calculator = neoag_calculator
        self.netmhcpan2 = netmhcpan2
        self.mixmhc2 = mixmhc2
        self.netmhcpan = netmhcpan
        self.mixmhc = mixmhc
        self.available_alleles = available_alleles
        self.uniprot = uniprot

        # TODO: assess which of this parse files and should be created only once due to performance reasons
        self.differential_binding = DifferentialBinding()
        self.gtex = GTEx()
        self.aa_frequency = AminoacidFrequency()
        self.fourmer_frequency = FourmerFrequency()
        self.aa_index = AminoacidIndex()
        self.priority_score_calculator = PriorityScore()
        self.tcell_predictor = TcellPrediction()
        self.iedb_immunogenicity = IEDBimmunogenicity()
        self.self_similarity = SelfSimilarityCalculator()

    def get_annotation(self, neoantigen: Neoantigen, patient: Patient) -> NeoantigenAnnotations:
        """Calculate new epitope features and add to dictonary that stores all properties"""

        self._initialise_annotations(neoantigen)

        # extract neoantigen data
        xmer_wt = neoantigen.mutation.wild_type_xmer
        xmer_mut = neoantigen.mutation.mutated_xmer
        gene = neoantigen.gene.gene
        vaf_tumor = neoantigen.dna_variant_allele_frequency
        vaf_rna = neoantigen.rna_variant_allele_frequency
        if not patient.is_rna_available:
            logger.warning("Using the DNA VAF to estimate the RNA VAF as the patient does not have RNA available")
            # TODO: overwrite value in the neoantigen object
            vaf_rna = vaf_tumor
        transcript_expr = neoantigen.rna_expression
        mutated_aminoacid = neoantigen.mutation.mutated_aminoacid
        wild_type_aminoacid = neoantigen.mutation.wild_type_aminoacid
        # TODO: this is needed by the T cell predictor, move this construction inside by passing the neoantigen
        substitution = "{}{}{}".format(
            neoantigen.mutation.wild_type_aminoacid, neoantigen.mutation.position, neoantigen.mutation.mutated_aminoacid)
        alleles_mhc_i = patient.mhc_i_alleles
        alleles_mhc_ii = patient.mhc_i_i_alleles
        tumor_content = patient.estimated_tumor_content

        # MHC binding independent features
        expression_calculator = Expression(
            transcript_expression=transcript_expr, vaf_rna=vaf_rna, tumor_content=tumor_content)
        self.annotations.annotations.extend(expression_calculator.get_annotations())
        self.add_differential_expression_features(gene, expression_tumor=transcript_expr, tissue=patient.tissue)
        self.add_aminoacid_index_features(mutation_aminoacid=mutated_aminoacid, wild_type_aminoacid=wild_type_aminoacid)
        self.add_provean_score_features(neoantigen)
        sequence_not_in_uniprot = self.uniprot.is_sequence_not_in_uniprot(xmer_mut)
        self.annotations.annotations.extend(self.uniprot.get_annotations(sequence_not_in_uniprot))

        # HLA I predictions: NetMHCpan
        self.netmhcpan.run(
            xmer_mut=xmer_mut, xmer_wt=xmer_wt, alleles=alleles_mhc_i,
            set_available_mhc=self.available_alleles.get_available_mhc_i())
        self.annotations.annotations.extend(self.netmhcpan.get_annotations())

        # read results from NetMHCpan
        epitope_wt_affinity = self.netmhcpan.best4_affinity_epitope_WT
        epitope_mut_affinity = self.netmhcpan.best4_affinity_epitope
        epitope_mut_affinitiy_9mer = self.netmhcpan.mhcI_affinity_epitope_9mer
        epitope_wt_affinitiy_9mer = self.netmhcpan.mhcI_affinity_epitope_9mer_WT
        epitope_wt_rank = self.netmhcpan.best4_mhc_epitope_WT
        epitope_mut_rank = self.netmhcpan.best4_mhc_epitope

        # MHC affinities/scores
        # TODO: ensure these conversions are not needed
        affinity_wt = self._2float(self.netmhcpan.best4_affinity_WT)
        affinity_mut = self._2float(self.netmhcpan.best4_affinity)
        affinity_wt_9mer = self._2float(self.netmhcpan.mhcI_affinity_9mer_WT)
        affinity_mut_9mer = self._2float(self.netmhcpan.mhcI_affinity_9mer)
        mhc_rank_mut = self._2float(self.netmhcpan.best4_mhc_score)
        mhc_rank_wt = self._2float(self.netmhcpan.best4_mhc_score_WT)
        wild_type_multiple_binding_score = self.netmhcpan.MHC_score_best_per_alelle_WT[1]
        mutation_multiple_binding_score = self.netmhcpan.MHC_score_best_per_alelle[1]

        # Neoantigen fitness metrics
        self.neoantigen_fitness_calculator.get_annotations(
            binding_wild_type=wild_type_multiple_binding_score, binding_mutation=mutation_multiple_binding_score,
            affinity_wild_type=affinity_wt, affinity_mutation=affinity_mut,
            rank_wild_type=mhc_rank_wt, rank_mutation=mhc_rank_mut,
            affinity_9mer_wild_type=affinity_wt_9mer, affinity_9mer_mutation=affinity_mut_9mer,
            sequence_wild_type=epitope_wt_affinity, sequence_mutation=epitope_mut_affinity,
            sequence_9mer_wild_type=epitope_wt_affinitiy_9mer, sequence_9mer_mutation=epitope_mut_affinitiy_9mer,
            sequence_rank_wild_type=epitope_wt_rank, sequence_rank_mutation=epitope_mut_rank
        )

        # T cell predictor
        self.add_tcell_predictor_features(gene, substitution=substitution, affinity=affinity_mut_9mer,
                                          epitope=epitope_mut_affinitiy_9mer)
        self.add_aminoacid_frequency_features(aa_freq=self.aa_frequency, mutation_mhci=epitope_mut_rank,
                                              nmer_freq=self.fourmer_frequency, mutated_aminoacid=mutated_aminoacid)

        # netMHCIIpan predictions
        self.netmhcpan2.run(sequence=xmer_mut, sequence_reference=xmer_wt, alleles=alleles_mhc_ii,
                            set_available_mhc=self.available_alleles.get_available_mhc_ii())
        self.annotations.annotations.extend(self.netmhcpan2.get_annotations())

        # MHC II scores
        affinity_wt_mhcii = self._2float(self.netmhcpan2.best_mhcII_affinity_WT)
        affinity_mut_mhcii = self._2float(self.netmhcpan2.best_mhcII_pan_affinity)
        rank_wt_mhcii = self._2float(self.netmhcpan2.best_mhcII_pan_score_WT)
        rank_mut_mhcii = self._2float(self.netmhcpan2.best_mhcII_pan_score)
        wild_type_multiple_binding_ii = self.netmhcpan2.MHCII_score_top10_WT[1]
        mutation_multiple_binding_ii = self.netmhcpan2.MHCII_score_top10[1]

        # MHC II epitopes
        epitope_wt_rank_mhcii = self.netmhcpan2.best_mhcII_pan_epitope_WT
        epitope_mut_rank_mhcii = self.netmhcpan2.best_mhcII_pan_epitope
        epitope_wt_affinity_mhcii = self.netmhcpan2.best_mhcII_affinity_epitope_WT
        epitope_mut_affinity_mhcii = self.netmhcpan2.best_mhcII_pan_affinity_epitope

        self.annotations.annotations.extend(self.neoantigen_fitness_calculator.get_annotations_mch2(
            mut_score=mutation_multiple_binding_ii,
            wt_score=wild_type_multiple_binding_ii,
            aff_wt=affinity_wt_mhcii, aff_mut=affinity_mut_mhcii,
            sc_wt=rank_wt_mhcii, sc_mut=rank_mut_mhcii, epitope_mut_mhcii=epitope_mut_affinity_mhcii
        ))

        # self-similarity
        self.annotations.annotations.extend(self.self_similarity.get_annnotations(
            epitope_mut_mhci=epitope_mut_rank, epitope_wt_mhci=epitope_wt_rank,
            rank_mut_mhci=mhc_rank_mut, rank_wt_mhci=mhc_rank_wt,
            epitope_mut_mhcii=epitope_mut_rank_mhcii,
            epitope_wt_mhcii=epitope_wt_rank_mhcii,
            rank_mut_mhcii=rank_mut_mhcii, rank_wt_mhcii=rank_wt_mhcii))

        # number of mismatches and priority scores
        self.annotations.annotations.extend(self.priority_score_calculator.get_annotations(
            epi_wt_mhci=epitope_wt_rank, epi_mut_mhci=epitope_mut_rank,
            epi_wt_mhcii=epitope_mut_rank_mhcii, epi_mut_mhcii=epitope_wt_rank_mhcii,
            rank_mut=mhc_rank_mut, rank_wt=mhc_rank_wt,
            mb_mut=mutation_multiple_binding_score,
            mb_wt=wild_type_multiple_binding_score,
            vaf_transcr=vaf_rna, vaf_tum=vaf_tumor, expr=transcript_expr, mut_not_in_prot=sequence_not_in_uniprot))

        # neoag immunogenicity model
        peptide_variant_position = EpitopeHelper.position_of_mutation_epitope(
            wild_type=self.netmhcpan.best4_affinity_epitope_WT, mutation=self.netmhcpan.best4_affinity_epitope)
        self.annotations.annotations.append(self.neoag_calculator.get_annotation(
            sample_id=patient.identifier, mut_peptide=epitope_mut_affinity, score_mut=affinity_mut,
            ref_peptide=epitope_wt_affinity, peptide_variant_position=peptide_variant_position))

        # IEDB immunogenicity
        self.iedb_immunogenicity.get_annotations(
            epitope_mhci=epitope_mut_affinity, affinity_mhci=affinity_mut, epitope_mhcii=epitope_mut_rank_mhcii,
            mhci_allele=self.netmhcpan.best4_affinity_allele, mhcii_allele=self.netmhcpan2.best_mhcII_pan_allele)

        # MixMHCpred
        self.mixmhc.run(xmer_wt=xmer_wt, xmer_mut=xmer_mut, alleles=alleles_mhc_i)
        self.annotations.annotations.extend(self.mixmhc.get_annotations())

        # MixMHC2pred
        self.mixmhc2.run(alleles=alleles_mhc_ii, xmer_wt=xmer_wt, xmer_mut=xmer_mut)
        self.annotations.annotations.extend(self.mixmhc2.get_annotations())

        ###########

        # dissimilarity to self-proteome
        self.annotations.annotations.extend(self.dissimilarity_calculator.get_annotations(
            epitope_mhci=epitope_mut_affinity, affinity_mhci=affinity_mut,
            epitope_mhcii=epitope_mut_affinity_mhcii, affinity_mhcii=affinity_mut_mhcii))

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

    def add_tcell_predictor_features(self, gene, substitution, epitope, affinity):
        # T cell predictor
        self.annotations.annotations.append(self.tcell_predictor.get_annotation(
            gene=gene, substitution=substitution, epitope=epitope, score=affinity))
        self.annotations.annotations.append(self.tcell_predictor.get_annotation(
            gene=gene, substitution=substitution, epitope=epitope, score=affinity, threshold=500))

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

    def add_aminoacid_frequency_features(self, aa_freq, mutation_mhci, nmer_freq, mutated_aminoacid):
        # amino acid frequency
        self.annotations.annotations.extend(aa_freq.get_annotations(mutated_aminoacid, mutation_mhci))
        self.annotations.annotations.extend(nmer_freq.get_annotations(mutation_mhci))

    def add_differential_expression_features(self, gene, expression_tumor, tissue):
        # differential expression
        gtex_mean, gtex_sd, gtex_sum = self.gtex.get_metrics(gene, tissue)
        self.annotations.annotations.extend(self.gtex.get_annotations(gtex_mean, gtex_sd, gtex_sum))
        self.annotations.annotations.extend(DifferentialExpression().get_annotations(
            expression_tumor=expression_tumor, expression_reference=gtex_mean,
            expression_reference_sd=gtex_sd, expression_reference_sum=gtex_sum))

    def _2float(self, value):
        # TODO: this is temporary! aghhh
        return float(value) if value != "NA" else "NA"
