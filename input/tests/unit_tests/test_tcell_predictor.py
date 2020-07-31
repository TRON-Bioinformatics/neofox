from unittest import TestCase
from input.predictors.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction


class TestTCellPredictor(TestCase):

    def setUp(self) -> None:
        self.tcell_predictor = TcellPrediction()

    # TODO: Franzis maybe you can add some more sensible tests here?

    def test_non_existing_gene(self):
        result = self.tcell_predictor.calculate_tcell_predictor_score(
            gene="BLAH", substitution='blaaaah', epitope="BLAHBLAH", score=5, threshold=10)
        self.assertEqual("NA", result)

    def test_existing_gene_with_too_short_epitope(self):
        result = self.tcell_predictor.calculate_tcell_predictor_score(
            gene="BRCA2", substitution='C', epitope="CCCCCC", score=5, threshold=10)
        self.assertEqual("NA", result)

    def test_existing_gene_with_too_long_epitope(self):
        result = self.tcell_predictor.calculate_tcell_predictor_score(
            gene="BRCA2", substitution='C', epitope="CCCCCCCCCC", score=5, threshold=10)
        self.assertEqual("NA", result)

    def test_existing_gene(self):
        result = self.tcell_predictor.calculate_tcell_predictor_score(
            gene="BRCA2", substitution='CCCCVCCCC', epitope="CCCCCCCCC", score=5, threshold=10)
        self.assertAlmostEqual(0.2453409331088489, float(result))
