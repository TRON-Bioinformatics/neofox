from unittest import TestCase

from input.aa_index.aa_index import AaIndex


class TestAaIndex(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.aaindex = AaIndex()

    def test_get_aaindex1_frequency(self):
        self.assertEqual(0.946, self.aaindex.get_aaindex1()["KARP850102"]['A'])

    def test_get_aaindex2_frequency(self):
        self.assertEqual(5.7, self.aaindex.get_aaindex2()["VOGG950101"]["A"]["C"])
