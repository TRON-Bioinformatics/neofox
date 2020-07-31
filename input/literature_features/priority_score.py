#!/usr/bin/env python

# import modules
import math


class PriorityScore:

    def number_of_mismatches(self, epitope_wild_type, epitope_mutation):
        """
        This function calculates the number of mismatches between the wt and the mutated epitope
        """
        p1 = 0
        for aa_mut, aa_wt in zip(epitope_mutation, epitope_wild_type):
            if aa_mut != aa_wt:
                p1 += 1
        return p1

    def calc_logistic_function(self, mhc_score):
        """
        Calculates negative logistic function given mhc score
        """
        try:
            log_score = 1.0 / (1.0 + math.exp(5.0 * (float(mhc_score) - 2.0)))
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
        priority_score = "NA"
        try:
            if vaf_tumor not in ["-1", "NA"]:
                priority_score = self.mupexi(l_mut, l_wt, mut_not_in_prot, no_mismatch, transcript_expr,
                                             vaf_tumor)
            elif vaf_rna not in ["-1", "NA"]:
                priority_score = self.mupexi(l_mut, l_wt, mut_not_in_prot, no_mismatch, transcript_expr,
                                   vaf_rna)
        except (TypeError, ValueError) as e:
            pass
        return priority_score

    def mupexi(self, l_mut, l_wt, mut_not_in_prot, no_mismatch, transcript_expr, vaf_tumor):
        priority_score = (l_mut * float(vaf_tumor) * math.tanh(float(transcript_expr))) * (
                float(mut_not_in_prot) * (1 - 2 ** (-float(no_mismatch)) * l_wt))
        return priority_score




