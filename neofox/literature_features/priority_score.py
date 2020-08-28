#!/usr/bin/env python

# import modules
import math
from typing import List

from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


class PriorityScore:

    def calc_logistic_function(self, mhc_score):
        """
        Calculates negative logistic function given mhc score
        """
        try:
            log_score = 1.0 / (1.0 + math.exp(5.0 * (mhc_score - 2.0)))
        except (OverflowError, ValueError) as e:
            log_score = "NA"
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

    def get_annotations(self, epi_wt_mhci, epi_mut_mhci, epi_wt_mhcii, epi_mut_mhcii, mut_not_in_prot,
                        rank_mut, rank_wt, mb_mut, mb_wt, expr, vaf_tum, vaf_transcr) -> List[Annotation]:
        """
        returns number of mismatches between best MHCI / MHC II epitopes (rank) and their corresponding WTs
        """
        num_mismatches_mhc1 = EpitopeHelper.number_of_mismatches(
            epitope_wild_type=epi_wt_mhci, epitope_mutation=epi_mut_mhci)
        return [
            AnnotationFactory.build_annotation(
                value=num_mismatches_mhc1,
                name="Number_of_mismatches_mhcI"),
            AnnotationFactory.build_annotation(
                value=EpitopeHelper.number_of_mismatches(
                    epitope_wild_type=epi_wt_mhcii, epitope_mutation=epi_mut_mhcii),
                name="Number_of_mismatches_mhcII"),
            # priority score with rank score
            AnnotationFactory.build_annotation(value=self.calc_priority_score(
                vaf_tumor=vaf_tum, vaf_rna=vaf_transcr, transcript_expr=expr, no_mismatch=num_mismatches_mhc1,
                score_mut=rank_mut, score_wt=rank_wt, mut_not_in_prot=mut_not_in_prot),
                name="Priority_score"),
            # priority score using multiplexed representation score
            AnnotationFactory.build_annotation(value=self.calc_priority_score(
                vaf_tumor=vaf_tum, vaf_rna=vaf_transcr, transcript_expr=expr, no_mismatch=num_mismatches_mhc1,
                score_mut=mb_mut, score_wt=mb_wt, mut_not_in_prot=mut_not_in_prot),
                name="Priority_score_MB")
            ]
