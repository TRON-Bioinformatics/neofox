from unittest import TestCase
from neofox.literature_features.priority_score import PriorityScore


class TestPriorityScore(TestCase):

    def setUp(self):
        self.priority_calculator = PriorityScore()

    def test_priority(self):
        result = self.priority_calculator.calc_priority_score(vaf_tumor=0.35, vaf_rna=0.33, transcript_expr=12,
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=True)
        self.assertGreater(result, 0)
        result = self.priority_calculator.calc_priority_score(vaf_tumor=None, vaf_rna=0.33, transcript_expr=12,
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=True)
        self.assertGreater(result, 0)
        result = self.priority_calculator.calc_priority_score(vaf_tumor=0.35, vaf_rna=None, transcript_expr=12,
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=True)
        self.assertGreater(result, 0)
        result = self.priority_calculator.calc_priority_score(vaf_tumor=None, vaf_rna=-1, transcript_expr=12,
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=True)
        self.assertEqual(result, None)
        result = self.priority_calculator.calc_priority_score(vaf_tumor=0.35, vaf_rna=0.33, transcript_expr=None,
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=True)
        self.assertEqual(result, None)
        result = self.priority_calculator.calc_priority_score(vaf_tumor=0.35, vaf_rna=0.33, transcript_expr=None,
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=True)
        self.assertEqual(result, None)
