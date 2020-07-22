from unittest import TestCase

from input.literature_features.differential_binding import DifferentialBinding


class TestDifferentialBinding(TestCase):

    def setUp(self):
        self.diffbdg_calculator = DifferentialBinding()

    def test_dai(self):
        result = self.diffbdg_calculator.dai(score_mutation=50, score_wild_type=500, affin_filtering=False)
        self.assertEqual(result, 450)
        result = self.diffbdg_calculator.dai(score_mutation=550, score_wild_type=2000, affin_filtering=True)
        self.assertEqual(result, "NA")
        result = self.diffbdg_calculator.dai(score_mutation=50, score_wild_type=10, affin_filtering=True)
        self.assertLess(result, 0.0)

    def test_diff_number_binders(self):
        result = self.diffbdg_calculator.diff_number_binders(num_mutation=5, num_wild_type=2)
        self.assertEqual(result, 3)
        result = self.diffbdg_calculator.diff_number_binders(num_mutation=2, num_wild_type=4)
        self.assertLess(result, 0)
        result = self.diffbdg_calculator.diff_number_binders(num_mutation="NA", num_wild_type=4)
        self.assertEqual(result, "NA")

    def test_ratio_number_binders(self):
        result = self.diffbdg_calculator.ratio_number_binders(num_mutation=10, num_wild_type=5)
        self.assertEqual(result, 2)
        result = self.diffbdg_calculator.ratio_number_binders(num_mutation=0, num_wild_type=5)
        self.assertEqual(result, 0)
        result = self.diffbdg_calculator.ratio_number_binders(num_mutation=5, num_wild_type=0)
        self.assertEqual(result, "NA")
        result = self.diffbdg_calculator.ratio_number_binders(num_mutation="NA", num_wild_type=3)
        self.assertEqual(result, "NA")

    def test_classify_adn_cdn(self):
        result = self.diffbdg_calculator.classify_adn_cdn(score_mutation=2, amplitude=11, bdg_cutoff_classical=50,
                                                          bdg_cutoff_alternative=5000, amplitude_cutoff=10,
                                                          category="CDN")
        self.assertEqual(result, "True")
        result = self.diffbdg_calculator.classify_adn_cdn(score_mutation=70, amplitude=11, bdg_cutoff_classical=50,
                                                          bdg_cutoff_alternative=5000, amplitude_cutoff=10,
                                                          category="CDN")
        self.assertEqual(result, "False")
        result = self.diffbdg_calculator.classify_adn_cdn(score_mutation=70, amplitude=11, bdg_cutoff_classical=50,
                                                          bdg_cutoff_alternative=5000, amplitude_cutoff=10,
                                                          category="ADN")
        self.assertEqual(result, "True")
        result = self.diffbdg_calculator.classify_adn_cdn(score_mutation=70, amplitude=5, bdg_cutoff_classical=50,
                                                          bdg_cutoff_alternative=5000, amplitude_cutoff=10,
                                                          category="ADN")
        self.assertEqual(result, "False")
        result = self.diffbdg_calculator.classify_adn_cdn(score_mutation="NA", amplitude=5, bdg_cutoff_classical=50,
                                                          bdg_cutoff_alternative=5000, amplitude_cutoff=10,
                                                          category="ADN")
        self.assertEqual(result, "NA")
