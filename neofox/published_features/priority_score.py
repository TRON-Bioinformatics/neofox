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

# import modules
import math
from typing import List
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.neoantigen import Annotation, PredictedEpitope, Neoantigen
from neofox.model.factories import AnnotationFactory
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)


class PriorityScore:
    def calc_logistic_function(self, mhc_score):
        """
        Calculates negative logistic function given mhc score
        """
        try:
            log_score = 1.0 / (1.0 + math.exp(5.0 * (mhc_score - 2.0)))
        except (OverflowError, ValueError):
            log_score = None
        return log_score

    def calc_priority_score(
        self,
        vaf,
        transcript_gene_expr,
        no_mismatch,
        score_mut,
        score_wt,
        mut_not_in_prot,
    ):
        """
        This function calculates the Priority Score using parameters for mhc I.
        Bjerregard, 2017, Cancer Immunol Immunother
        https://doi.org/10.1007/s00262-017-2001-3
        """
        priority_score = None
        try:
            if vaf is not None and vaf != -1:
                l_mut = self.calc_logistic_function(score_mut)
                l_wt = self.calc_logistic_function(score_wt)
                priority_score = self.mupexi(
                    l_mut=l_mut,
                    l_wt=l_wt,
                    mut_not_in_prot=mut_not_in_prot,
                    no_mismatch=no_mismatch,
                    transcript_gene_expr=transcript_gene_expr,
                    vaf_tumor=vaf
                )
        except (TypeError, ValueError):
            pass
        return priority_score

    def mupexi(
        self, l_mut, l_wt, mut_not_in_prot, no_mismatch, transcript_gene_expr, vaf_tumor
    ):
        priority_score = (l_mut * vaf_tumor * math.tanh(transcript_gene_expr)) * (
            float(mut_not_in_prot) * (1 - 2 ** (-no_mismatch) * l_wt)
        )
        return priority_score

    def get_annotations(
        self,
        netmhcpan: BestAndMultipleBinder,
        mut_not_in_prot,
        neoantigen: Neoantigen
    ) -> List[Annotation]:
        """
        returns number of mismatches between best MHCI / MHC II epitopes (rank) and their corresponding WTs
        """
        num_mismatches_mhc1 = None
        priority_score_dna = None
        priority_score_rna = None
        priority_score_imputed_dna = None
        priority_score_imputed_rna = None
        if netmhcpan.best_epitope_by_rank.wild_type_peptide and netmhcpan.best_epitope_by_rank.mutated_peptide:
            num_mismatches_mhc1 = EpitopeHelper.number_of_mismatches(
                epitope_wild_type=netmhcpan.best_epitope_by_rank.wild_type_peptide,
                epitope_mutation=netmhcpan.best_epitope_by_rank.mutated_peptide,
            )

            priority_score_dna = self.calc_priority_score(
                        vaf=neoantigen.dna_variant_allele_frequency,
                        transcript_gene_expr=neoantigen.rna_expression,
                        no_mismatch=num_mismatches_mhc1,
                        score_mut=netmhcpan.best_epitope_by_rank.rank_mutated,
                        score_wt=netmhcpan.best_epitope_by_rank.rank_wild_type,
                        mut_not_in_prot=mut_not_in_prot,
                    )
            priority_score_rna = self.calc_priority_score(
                        vaf=neoantigen.rna_variant_allele_frequency,
                        transcript_gene_expr=neoantigen.rna_expression,
                        no_mismatch=num_mismatches_mhc1,
                        score_mut=netmhcpan.best_epitope_by_rank.rank_mutated,
                        score_wt=netmhcpan.best_epitope_by_rank.rank_wild_type,
                        mut_not_in_prot=mut_not_in_prot,
                    )
            priority_score_imputed_dna = self.calc_priority_score(
                        vaf=neoantigen.dna_variant_allele_frequency,
                        transcript_gene_expr=neoantigen.imputed_gene_expression,
                        no_mismatch=num_mismatches_mhc1,
                        score_mut=netmhcpan.best_epitope_by_rank.rank_mutated,
                        score_wt=netmhcpan.best_epitope_by_rank.rank_wild_type,
                        mut_not_in_prot=mut_not_in_prot,
                    )
            priority_score_imputed_rna = self.calc_priority_score(
                        vaf=neoantigen.rna_variant_allele_frequency,
                        transcript_gene_expr=neoantigen.imputed_gene_expression,
                        no_mismatch=num_mismatches_mhc1,
                        score_mut=netmhcpan.best_epitope_by_rank.rank_mutated,
                        score_wt=netmhcpan.best_epitope_by_rank.rank_wild_type,
                        mut_not_in_prot=mut_not_in_prot,
                    )
        annotations = [
            AnnotationFactory.build_annotation(
                value=num_mismatches_mhc1, name="Number_of_mismatches_MCHI"
            ),
            # priority score with rank score
            AnnotationFactory.build_annotation(
                value=priority_score_dna,
                name="Priority_score_fromDNA",
            ),
            # imputed priority score with rank score
            AnnotationFactory.build_annotation(
                value=priority_score_imputed_rna,
                name="Priority_score_imputed_fromRNA"
            ),
            # priority score with rank score f
            AnnotationFactory.build_annotation(
                value=priority_score_rna,
                name="Priority_score_fromRNA",
            ),
            # imputed priority score with rank score
            AnnotationFactory.build_annotation(
                value=priority_score_imputed_dna,
                name="Priority_score_imputed_fromDNA"
            )
        ]
        return annotations

    def get_annotations_epitope_mhci(self, epitope: PredictedEpitope, vaf_tumor, transcript_exp, vaf_rna, gene_exp) -> \
            List[Annotation]:

        priority_score_dna = None
        priority_score_rna = None
        priority_score_imputed_dna = None
        priority_score_imputed_rna = None

        if epitope.mutated_peptide and epitope.rank_wild_type: 
            priority_score_dna = self.calc_priority_score(
                    vaf=vaf_tumor,
                    transcript_gene_expr=transcript_exp,
                    no_mismatch=int(EpitopeHelper.get_annotation_by_name(
                        epitope.neofox_annotations.annotations, name='number_of_mismatches')),
                    score_mut=epitope.rank_mutated,
                    score_wt=epitope.rank_wild_type,
                    mut_not_in_prot=bool(EpitopeHelper.get_annotation_by_name(
                        epitope.neofox_annotations.annotations, name='mutation_not_found_in_proteome')))

        if epitope.mutated_peptide and epitope.rank_wild_type: 
            priority_score_imputed_dna = priority_score_rna = self.calc_priority_score(
                    vaf=vaf_tumor,
                    transcript_gene_expr=gene_exp,
                    no_mismatch=int(EpitopeHelper.get_annotation_by_name(
                        epitope.neofox_annotations.annotations, name='number_of_mismatches')),
                    score_mut=epitope.rank_mutated,
                    score_wt=epitope.rank_wild_type,
                    mut_not_in_prot=bool(EpitopeHelper.get_annotation_by_name(
                        epitope.neofox_annotations.annotations, name='mutation_not_found_in_proteome')))

        if epitope.mutated_peptide and epitope.rank_wild_type: 
            priority_score_rna = self.calc_priority_score(
                    vaf=vaf_rna,
                    transcript_gene_expr=transcript_exp,
                    no_mismatch=int(EpitopeHelper.get_annotation_by_name(
                        epitope.neofox_annotations.annotations, name='number_of_mismatches')),
                    score_mut=epitope.rank_mutated,
                    score_wt=epitope.rank_wild_type,
                    mut_not_in_prot=bool(EpitopeHelper.get_annotation_by_name(
                        epitope.neofox_annotations.annotations, name='mutation_not_found_in_proteome')))

        if epitope.mutated_peptide and epitope.rank_wild_type: 
            priority_score_imputed_rna = self.calc_priority_score(
                    vaf=vaf_rna,
                    transcript_gene_expr=gene_exp,
                    no_mismatch=int(EpitopeHelper.get_annotation_by_name(
                        epitope.neofox_annotations.annotations, name='number_of_mismatches')),
                    score_mut=epitope.rank_mutated,
                    score_wt=epitope.rank_wild_type,
                    mut_not_in_prot=bool(EpitopeHelper.get_annotation_by_name(
                        epitope.neofox_annotations.annotations, name='mutation_not_found_in_proteome')))

        return [
            AnnotationFactory.build_annotation(
                value=priority_score_dna,
                name='Priority_score_fromDNA'),
            AnnotationFactory.build_annotation(
                value=priority_score_imputed_dna,
                name='Priority_score_imputed_fromDNA'),
            AnnotationFactory.build_annotation(
                value=priority_score_rna,
                name='Priority_score_fromRNA'),
            AnnotationFactory.build_annotation(
                value=priority_score_imputed_rna,
                name='Priority_score_imputed_fromRNA'),
            ]
