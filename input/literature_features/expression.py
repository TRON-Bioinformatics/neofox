#!/usr/bin/env python


class Expression:

    def rna_expression_mutation(self, transcript_expression, vaf_rna):
        """
        This function calculates the product of VAF in RNA and transcript expression
        to reflect the expression of the mutated transcript
        """
        try:
            expression_mut = transcript_expression * float(vaf_rna) if float(vaf_rna) >= 0.0 else "NA"
        except ValueError:
            expression_mut = "NA"
        return expression_mut

    def rna_expression_mutation_tc(self, expression_mutation, tumor_content):
        """
        calculated expression of mutation corrected by tumour content
        """
        try:
            expression_mut_tc = expression_mutation / tumor_content
        except (TypeError, ZeroDivisionError) as e:
            expression_mut_tc = "NA"
        return expression_mut_tc

