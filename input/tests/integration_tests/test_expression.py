
from unittest import TestCase

from input.literature_features.expression import Expression



class TestExpression(TestCase):

    def setUp(self):
        self.expression_calculator = Expression()

    def test_calculate_expression(self):
        self.expression_calculator.calculate_expression(transcript_expression=12.0,
                                                        vaf_rna=0.2, tumor_content=0.7)
        self.assertGreater(float(self.expression_calculator.expression_mutation), 0.0)
        self.assertGreater(float(self.expression_calculator.expression_mutation_tc), 0.0)
        # no reads for mut
        self.expression_calculator.calculate_expression(transcript_expression=12.0,
                                                        vaf_rna=0, tumor_content=0.7)
        self.assertEqual(float(self.expression_calculator.expression_mutation), 0.0)
        self.assertEqual(float(self.expression_calculator.expression_mutation_tc), 0.0)
        # no reads for mut/wt
        self.expression_calculator.calculate_expression( 14.0, -1, 0.7)
        self.assertEqual(self.expression_calculator.expression_mutation, "NA")
        self.assertEqual(self.expression_calculator.expression_mutation_tc, "NA")
        # tumor content = "NA"
        self.expression_calculator.calculate_expression(14.0, 0.2, "NA")
        self.assertGreater(float(self.expression_calculator.expression_mutation), 0.0)
        self.assertEqual(self.expression_calculator.expression_mutation_tc, "NA")