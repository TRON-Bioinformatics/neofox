
from unittest import TestCase

from input.literature_features.expression import Expression



class TestExpression(TestCase):

    def setUp(self):
        self.expression_calculator = Expression()

    def test_calculate_expression_mutation(self):
        result = self.expression_calculator.rna_expression_mutation(transcript_expression=12.0,
                                                                    vaf_rna=0.2)
        self.assertGreater(result, 0.0)
        # no reads for mut
        result = self.expression_calculator.rna_expression_mutation(transcript_expression=12.0,
                                                                    vaf_rna=0.0)
        self.assertEqual(result, 0.0)
        # no reads for mut/wt
        result = self.expression_calculator.rna_expression_mutation(transcript_expression=12.0,
                                                                    vaf_rna=-1)
        self.assertEqual(result, "NA")
        result = self.expression_calculator.rna_expression_mutation(transcript_expression="NA",
                                                                    vaf_rna=-1)
        self.assertEqual(result, "NA")


    def test_expression_mutation_tumor_content(self):
        result = self.expression_calculator.rna_expression_mutation_tc(expression_mutation=12.0,
                                                                       tumor_content=0.7)
        self.assertGreater(result, 0.0)
        # expression mutation = "NA"
        result = self.expression_calculator.rna_expression_mutation_tc(expression_mutation="NA",
                                                                       tumor_content=0.7)
        self.assertEqual(result, "NA")
        # tumor content = "NA"
        result = self.expression_calculator.rna_expression_mutation_tc(expression_mutation=12.0,
                                                                       tumor_content="NA")
        self.assertEqual(result, "NA")