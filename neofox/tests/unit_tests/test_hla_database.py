import unittest
from neofox.tests.fake_classes import FakeHlaDatabase


class TestHlaDatabase(unittest.TestCase):

    def test_hla_database(self):
        hla_database = FakeHlaDatabase()
        self.assertIsNotNone(hla_database)
        self.assertIsNotNone(hla_database.alleles)
        self.assertTrue(len(hla_database.alleles) > 1000)
        self.assertTrue("DPB1*104:01" in hla_database.alleles)
        self.assertFalse("DPB1*10:401" in hla_database.alleles)
        self.assertTrue("B*15:228" in hla_database.alleles)
        self.assertFalse("B*152:28" in hla_database.alleles)
