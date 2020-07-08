from unittest import TestCase

import input.self_similarity.self_similarity as self_similarity
from input import MHC_I, MHC_II


class TestSelfSimilarity(TestCase):

    def test_get_self_similarity(self):
        result = self_similarity.get_self_similarity(wild_type="DDD", mutation="DDD")
        self.assertEqual('1.0', result)

    def test_is_improved_binder(self):
        result = self_similarity.is_improved_binder(
            score_mutation='1.0', score_wild_type='1.3')
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


