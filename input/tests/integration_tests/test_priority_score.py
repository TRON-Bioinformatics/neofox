from unittest import TestCase

from input.literature_features.priority_score import PriorityScore
from input.predict_all_epitopes import BunchEpitopes
from input.references import ReferenceFolder


class TestPriorityScore(TestCase):

    def setUp(self):
        self.priority_calculator = PriorityScore()
        self.references = ReferenceFolder()
        self.db = BunchEpitopes.load_proteome(self.references.uniprot)


    def test_proteome_match(self):
        # sequence in proteome
        result = self.priority_calculator.match_not_in_proteome(sequence="MAELCPLA", db=self.db)
        self.assertEqual(result, "0")
        # sequence not int proteome
        result = self.priority_calculator.match_not_in_proteome(sequence="FIAGLIAIV", db=self.db)
        self.assertEqual(result, "1")

    def test_number_of_mismatches(self):
        # one mismatch
        result = self.priority_calculator.number_of_mismatches("FIAGLIAIV", "FIAGDIAIV")
        self.assertEqual(result, 1)
        # three mismatch
        result = self.priority_calculator.number_of_mismatches("FIAGLIAIV", "FLAGDIAIN")
        self.assertEqual(result, 3)
        # no mismatch
        result = self.priority_calculator.number_of_mismatches("FIAGLIAIV", "FIAGLIAIV")
        self.assertEqual(result, 0)
        # wt sequence shorter than mut sequence
        result = self.priority_calculator.number_of_mismatches("FIAGI", "FIAGLIAIV")
        self.assertEqual(result, 1)
        # mut sequence shorter than wt sequence
        result = self.priority_calculator.number_of_mismatches("FIAGLIAIV", "FIAGI" )
        self.assertEqual(result, 1)

    def test_priority(self):
        result = self.priority_calculator.calc_priority_score(vaf_tumor=0.35, vaf_rna=0.33, transcript_expr=12,
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=1)
        self.assertGreater(result, 0)
        result = self.priority_calculator.calc_priority_score(vaf_tumor="NA", vaf_rna=0.33, transcript_expr=12,
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=1)
        self.assertGreater(result, 0)
        result = self.priority_calculator.calc_priority_score(vaf_tumor=0.35, vaf_rna="NA", transcript_expr=12,
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=1)
        self.assertGreater(result, 0)
        result = self.priority_calculator.calc_priority_score(vaf_tumor="NA", vaf_rna="-1", transcript_expr=12,
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=1)
        self.assertEqual(result, "NA")
        result = self.priority_calculator.calc_priority_score(vaf_tumor=0.35, vaf_rna=0.33, transcript_expr="NA",
                                                              no_mismatch=1, score_mut=1.1, score_wt=10,
                                                              mut_not_in_prot=1)
        self.assertEqual(result, "NA")
        result = self.priority_calculator.calc_priority_score(vaf_tumor=0.35, vaf_rna=0.33, transcript_expr="NA",
                                                              no_mismatch=1, score_mut="NA", score_wt=10,
                                                              mut_not_in_prot=1)
        self.assertEqual(result, "NA")