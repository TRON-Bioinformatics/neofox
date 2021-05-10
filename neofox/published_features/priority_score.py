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
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
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
        vaf_tumor,
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
        l_mut = self.calc_logistic_function(score_mut)
        l_wt = self.calc_logistic_function(score_wt)
        priority_score = None
        try:
            if vaf_tumor is not None and vaf_tumor != -1:
                priority_score = self.mupexi(
                    l_mut,
                    l_wt,
                    mut_not_in_prot,
                    no_mismatch,
                    transcript_expr,
                    vaf_tumor,
                )
            elif vaf_rna is not None and vaf_rna != -1:
                priority_score = self.mupexi(
                    l_mut, l_wt, mut_not_in_prot, no_mismatch, transcript_expr, vaf_rna
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
        expr,
        vaf_tum,
        vaf_transcr,
    ) -> List[Annotation]:
        """
        returns number of mismatches between best MHCI / MHC II epitopes (rank) and their corresponding WTs
        """
        num_mismatches_mhc1 = None
        priority_score = None
        if netmhcpan.best_wt_epitope_by_rank.peptide and netmhcpan.best_epitope_by_rank.peptide:
            num_mismatches_mhc1 = EpitopeHelper.number_of_mismatches(
                epitope_wild_type=netmhcpan.best_wt_epitope_by_rank.peptide,
                epitope_mutation=netmhcpan.best_epitope_by_rank.peptide,
            )
            priority_score = self.calc_priority_score(
                        vaf_tumor=vaf_tum,
                        vaf_rna=vaf_transcr,
                        transcript_expr=expr,
                        no_mismatch=num_mismatches_mhc1,
                        score_mut=netmhcpan.best_epitope_by_rank.rank,
                        score_wt=netmhcpan.best_wt_epitope_by_rank.rank,
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
