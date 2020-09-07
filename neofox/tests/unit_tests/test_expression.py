
from unittest import TestCase

from neofox.published_features.expression import Expression



class TestExpression(TestCase):

    def test_calculate_expression_mutation(self):
        result = Expression(transcript_expression=12.0, vaf_rna=0.2, tumor_content=None).get_annotations()[0]
        self.assertGreater(result.value, "0.0")
        # no reads for mut
        result = Expression(transcript_expression=12.0, vaf_rna=0.0, tumor_content=None).get_annotations()[0]
        self.assertEqual(result.value, "0.0")
        # no reads for mut/wt
        result = Expression(transcript_expression=12.0, vaf_rna=-1, tumor_content=None).get_annotations()[0]
        self.assertEqual(result.value, "NA")
        result = Expression(transcript_expression=None, vaf_rna=-1, tumor_content=None).get_annotations()[0]
        self.assertEqual(result.value, "NA")

    def test_expression_mutation_tumor_content(self):
        result = Expression(transcript_expression=12.0, tumor_content=0.7, vaf_rna=1.0).get_annotations()[1]
        self.assertGreater(float(result.value), 0.0)
        result = Expression(transcript_expression=12.0, tumor_content=0.7, vaf_rna=0.0).get_annotations()[1]
        self.assertEqual(float(result.value), 0.0)
        result = Expression(transcript_expression=None, tumor_content=0.7, vaf_rna=None).get_annotations()[1]
        self.assertEqual(result.value, "NA")
        result = Expression(transcript_expression=12.0, tumor_content=None, vaf_rna=None).get_annotations()[1]
        self.assertEqual(result.value, "NA")
