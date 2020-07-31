from unittest import TestCase

from input.annotation_resources.gtex.gtex import GTEx


class TestGtex(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.gtex = GTEx()

    def test_get_metrics(self):
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics("BRCA2", "skin")
        self.assertEqual(0.146363636363636, mean_expression)
        self.assertEqual(27.37, sum_expression)
        self.assertEqual(0.0950531054943768, sd_expression)
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics("CFTR", "bladder")
        self.assertEqual(0.35454545454545505, mean_expression)
        self.assertEqual(3.9, sum_expression)
        self.assertEqual(1.11291835851839, sd_expression)

    def test_get_metrics_from_non_existing_tissue(self):
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics("BRCA2", "chinichinchin")
        self.assertEqual("NA", mean_expression)
        self.assertEqual("NA", sum_expression)
        self.assertEqual("NA", sd_expression)
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics("BRCA2", None)
        self.assertEqual("NA", mean_expression)
        self.assertEqual("NA", sum_expression)
        self.assertEqual("NA", sd_expression)

    def test_get_metrics_from_non_existing_gene(self):
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics("NOOOOOOOPE", "skin")
        self.assertEqual("NA", mean_expression)
        self.assertEqual("NA", sum_expression)
        self.assertEqual("NA", sd_expression)
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics(None, "skin")
        self.assertEqual("NA", mean_expression)
        self.assertEqual("NA", sum_expression)
        self.assertEqual("NA", sd_expression)


