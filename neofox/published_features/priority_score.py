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
from neofox.model.neoantigen import Annotation, PredictedEpitope, Neoantigen, Patient
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
        vaf_dna,
        vaf_rna,
        transcript_expr,
        no_mismatch,
        score_mut,
        score_wt,
        mut_not_in_prot,
    ):
        """
        This function calculates the Priority Score using parameters for mhc I.
        """
        priority_score = None
        vaf = None
        try:
            if vaf_dna is not None and vaf_dna != -1:
                vaf = vaf_dna
            elif vaf_rna is not None and vaf_rna != -1:
                vaf = vaf_rna
            if vaf:
                l_mut = self.calc_logistic_function(score_mut)
                l_wt = self.calc_logistic_function(score_wt)
                priority_score = self.mupexi(
                    l_mut=l_mut,
                    l_wt=l_wt,
                    mut_not_in_prot=mut_not_in_prot,
                    no_mismatch=no_mismatch,
                    transcript_expr=transcript_expr,
                    vaf_tumor=vaf
                )
        except (TypeError, ValueError):
            pass
        return priority_score

    def mupexi(
        self, l_mut, l_wt, mut_not_in_prot, no_mismatch, transcript_expr, vaf_tumor
    ):
        priority_score = (l_mut * vaf_tumor * math.tanh(transcript_expr)) * (
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
        priority_score = None
        if netmhcpan.best_epitope_by_rank.wild_type_peptide and netmhcpan.best_epitope_by_rank.mutated_peptide:
            num_mismatches_mhc1 = EpitopeHelper.number_of_mismatches(
                epitope_wild_type=netmhcpan.best_epitope_by_rank.wild_type_peptide,
                epitope_mutation=netmhcpan.best_epitope_by_rank.mutated_peptide,
            )
            vaf_rna = neoantigen.dna_variant_allele_frequency
            if vaf_rna is None:
                vaf_rna = neoantigen.rna_variant_allele_frequency

            priority_score = self.calc_priority_score(
                        vaf_dna=neoantigen.dna_variant_allele_frequency,
                        vaf_rna=vaf_rna,
                        transcript_expr=neoantigen.rna_expression,
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
                value=priority_score,
                name="Priority_score",
            ),
        ]
        return annotations

    def get_annotations_epitope_mhci(self, epitope: PredictedEpitope, vaf_tumor, transcript_exp, vaf_rna) -> \
            List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.calc_priority_score(
                    vaf_dna=vaf_tumor,
                    vaf_rna=vaf_rna,
                    transcript_expr=transcript_exp,
                    no_mismatch=int(EpitopeHelper.get_annotation_by_name(
                        epitope.neofox_annotations.annotations, name='number_of_mismatches')),
                    score_mut=epitope.rank_mutated,
                    score_wt=epitope.rank_wild_type,
                    mut_not_in_prot=bool(EpitopeHelper.get_annotation_by_name(
                        epitope.neofox_annotations.annotations, name='mutation_not_found_in_proteome'))),
                name='Priority_score')
            ]
