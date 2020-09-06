#!/usr/bin/env python

# import modules
import math
from typing import List

from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.predictors.netmhcpan4.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.predictors.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder


class PriorityScore:

    def calc_logistic_function(self, mhc_score):
        """
        Calculates negative logistic function given mhc score
        """
        try:
            log_score = 1.0 / (1.0 + math.exp(5.0 * (mhc_score - 2.0)))
        except (OverflowError, ValueError) as e:
            log_score = None
        return log_score

    def calc_priority_score(self, vaf_tumor, vaf_rna, transcript_expr, no_mismatch, score_mut, score_wt,
                            mut_not_in_prot):
        """
        This function calculates the Priority Score using parameters for mhc I.
        """
        l_mut = self.calc_logistic_function(score_mut)
        l_wt = self.calc_logistic_function(score_wt)
        priority_score = None
        try:
            if vaf_tumor is not None and vaf_tumor != -1:
                priority_score = self.mupexi(l_mut, l_wt, mut_not_in_prot, no_mismatch, transcript_expr, vaf_tumor)
            elif vaf_rna is not None and vaf_rna != -1:
                priority_score = self.mupexi(l_mut, l_wt, mut_not_in_prot, no_mismatch, transcript_expr, vaf_rna)
        except (TypeError, ValueError) as e:
            pass
        return priority_score

    def mupexi(self, l_mut, l_wt, mut_not_in_prot, no_mismatch, transcript_expr, vaf_tumor):
        priority_score = (l_mut * vaf_tumor * math.tanh(transcript_expr)) * (
                float(mut_not_in_prot) * (1 - 2 ** (-no_mismatch) * l_wt))
        return priority_score

    def get_annotations(
            self, netmhcpan: BestAndMultipleBinder, netmhcpan2: BestAndMultipleBinderMhcII,
            mut_not_in_prot, expr, vaf_tum, vaf_transcr) -> List[Annotation]:
        """
        returns number of mismatches between best MHCI / MHC II epitopes (rank) and their corresponding WTs
        """
        num_mismatches_mhc1 = EpitopeHelper.number_of_mismatches(
            epitope_wild_type=netmhcpan.best4_mhc_epitope_WT, epitope_mutation=netmhcpan.best4_mhc_epitope)
        return [
            AnnotationFactory.build_annotation(
                value=num_mismatches_mhc1,
                name="Number_of_mismatches_MCHI"),
            AnnotationFactory.build_annotation(
                value=EpitopeHelper.number_of_mismatches(
                    epitope_wild_type=netmhcpan2.best_mhcII_pan_epitope_WT,
                    epitope_mutation=netmhcpan2.best_mhcII_pan_epitope),
                name="Number_of_mismatches_MHCII"),
            # priority score with rank score
            AnnotationFactory.build_annotation(value=self.calc_priority_score(
                vaf_tumor=vaf_tum, vaf_rna=vaf_transcr, transcript_expr=expr, no_mismatch=num_mismatches_mhc1,
                score_mut=netmhcpan.best4_mhc_score, score_wt=netmhcpan.best4_mhc_score_WT,
                mut_not_in_prot=mut_not_in_prot),
                name="Priority_score"),
            # priority score using multiplexed representation score
            AnnotationFactory.build_annotation(value=self.calc_priority_score(
                vaf_tumor=vaf_tum, vaf_rna=vaf_transcr, transcript_expr=expr, no_mismatch=num_mismatches_mhc1,
                score_mut=netmhcpan.MHC_score_top10[1], score_wt=netmhcpan.MHC_score_top10_WT[1],
                mut_not_in_prot=mut_not_in_prot),
                name="Priority_score_multiple_binding")
            ]
