from unittest import TestCase

from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator


class TestSelfSimilarity(TestCase):

    def test_get_self_similarity(self):
        result = SelfSimilarityCalculator().get_self_similarity(wild_type="DDD", mutation="DDD")
        self.assertEqual('1.0', result)

    def test_is_improved_binder(self):
        result = SelfSimilarityCalculator().is_improved_binder(score_mutation=1.0, score_wild_type=1.3)
        self.assertTrue(result)

    def test_compute_self_similarity_calculator(self):

        s = SelfSimilarityCalculator()
        self.assertEqual(s.compute_k_hat_3("AAAAA", "AAAAA"), 1.0)
        for i in range(5):
            self.assertTrue(s.compute_k_hat_3("AAAAA", "WWWWW" * (i + 1)) < 1.0)
        for i in list(s.k1.keys()):
            if i == "A":
                self.assertTrue(s.compute_k_hat_3("AAAAA", "AA" + i + "AA") == 1.0)
            else:
                self.assertTrue(s.compute_k_hat_3("AAAAA", "AA" + i + "AA") < 1.0)


