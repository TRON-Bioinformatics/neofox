from unittest import TestCase

from neofox.potential_features.nmer_frequency.nmer_frequency import AminoacidFrequency, FourmerFrequency


class TestNmerFrequency(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.aa_frequencies = AminoacidFrequency()
        cls.fourmer_frequencies = FourmerFrequency()

    def test_get_aminoacid_frequency(self):
        self.assertEqual(7.27339103101, self.aa_frequencies._get_frequency("A"))

    def test_non_existing_aminoacid_frequency(self):
        self.assertIsNone(self.aa_frequencies._get_frequency("X"))

    def test_get_4mer_frequency(self):
        self.assertEqual(0.0148452054899, self.fourmer_frequencies._get_frequency_4mer("AAAAAAAA"))
        self.assertEqual(0.000233496235874, self.fourmer_frequencies._get_frequency_4mer("AAAACCAA"))

    def test_non_existing_4mer_frequency(self):
        self.assertIsNone(self.fourmer_frequencies._get_frequency_4mer("AA"))
        self.assertIsNone(self.fourmer_frequencies._get_frequency_4mer(""))
        self.assertIsNone(self.fourmer_frequencies._get_frequency_4mer("XXXXXXXX"))

    def test_get_4mer_product_frequency(self):
        self.assertEqual(2798.6445730350238, self.aa_frequencies._get_product_4mer_frequencies("AAAAAAAA"))
        self.assertEqual(215.14851896633044, self.aa_frequencies._get_product_4mer_frequencies("AAAACCAA"))

    def test_non_existing_4mer_product_frequency(self):
        self.assertIsNone(self.aa_frequencies._get_product_4mer_frequencies("AA"))
        self.assertIsNone(self.aa_frequencies._get_product_4mer_frequencies(""))
        self.assertIsNone(self.aa_frequencies._get_product_4mer_frequencies("XXXXXXXX"))
