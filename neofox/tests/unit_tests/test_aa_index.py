from unittest import TestCase

from neofox.potential_features.aa_index.aa_index import AminoacidIndex


class TestAaIndex(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.aaindex = AminoacidIndex()

    def test_get_aaindex1_frequency(self):
        self.assertEqual(0.946, self.aaindex.aaindex1["KARP850102"]['A'])

    def test_get_aaindex2_frequency(self):
        self.assertEqual(5.7, self.aaindex.aaindex2["VOGG950101"]["A"]["C"])
