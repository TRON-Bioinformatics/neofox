from unittest import TestCase

import input.self_similarity.self_similarity as self_similarity
from input import MHC_I, MHC_II


class TestSelfSimilarity(TestCase):

    def test_get_self_similarity(self):
        result = self_similarity.get_self_similarity(wild_type="DDD", mutation="DDD")
        self.assertEqual('1.0', result)

    def test_is_improved_binder_mhci(self):
        result = self_similarity.is_improved_binder(
            props={'best%Rank_netmhcpan4': '1.0', 'best%Rank_netmhcpan4_WT': '1.3'},
            mhc=MHC_I)
        self.assertEqual('1', result)

    def test_is_improved_binder_mhcii(self):
        result = self_similarity.is_improved_binder(
            props={'MHC_II_score_.best_prediction.': '1,0', 'MHC_II_score_.WT.': '1,3'},
            mhc=MHC_II)
        self.assertEqual('1', result)

    def test_position_mutation(self):
        position = self_similarity.position_of_mutation_epitope(wild_type="AAAAAA", mutation="AAANAA")
        self.assertEqual(position, "4")
        position = self_similarity.position_of_mutation_epitope(wild_type="AAAAAA", mutation="AAAAAA")
        self.assertEqual(position, "-1")
        position = self_similarity.position_of_mutation_epitope(wild_type="AAAAAA", mutation="AANNNN")
        self.assertEqual(position, "6")

    def test_compute_self_similarity_calculator(self):

        s = self_similarity.SelfSimilarityCalculator()
        self.assertEqual(s.compute_k_hat_3("AAAAA", "AAAAA"), 1.0)
        for i in range(5):
            self.assertTrue(s.compute_k_hat_3("AAAAA", "WWWWW" * (i + 1)) < 1.0)
        for i in list(s.k1.keys()):
            if i == "A":
                self.assertTrue(s.compute_k_hat_3("AAAAA", "AA" + i + "AA") == 1.0)
            else:
                self.assertTrue(s.compute_k_hat_3("AAAAA", "AA" + i + "AA") < 1.0)


