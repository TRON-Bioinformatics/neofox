import unittest
from neofox.model.neoantigen import MhcAllele
from neofox.tests.fake_classes import FakeHlaDatabase, FakeH2Database


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

    def test_exists(self):
        hla_database = FakeHlaDatabase()
        self.assertTrue(hla_database.exists(MhcAllele(gene="DPB1", group="104", protein="01")))
        self.assertFalse(hla_database.exists(MhcAllele(gene="DPB1", group="10", protein="401")))
        self.assertTrue(hla_database.exists(MhcAllele(gene="B", group="15", protein="228")))
        self.assertFalse(hla_database.exists(MhcAllele(gene="B", group="152", protein="28")))
        # badly formed HLA alleles
        self.assertFalse(hla_database.exists(MhcAllele(gene="B", group="15", protein=None)))
        self.assertFalse(hla_database.exists(MhcAllele(gene="B", group=None, protein="228")))
        self.assertFalse(hla_database.exists(MhcAllele(gene=None, group="15", protein="228")))
        self.assertFalse(hla_database.exists(MhcAllele(gene="Z", group="15", protein="228")))


class TestH2Database(unittest.TestCase):

    def test_h2_database(self):
        h2_database = FakeH2Database()
        self.assertIsNotNone(h2_database)
        self.assertIsNotNone(h2_database.alleles)
        self.assertTrue(len(h2_database.alleles) == 69)
        self.assertTrue("H2Kp" in h2_database.alleles)
        self.assertFalse("H2Kx" in h2_database.alleles)
        self.assertTrue("H2Lf" in h2_database.alleles)
        self.assertFalse("Nope" in h2_database.alleles)

    def test_exists(self):
        h2_database = FakeH2Database()
        self.assertTrue(h2_database.exists(MhcAllele(gene="H2K", protein="p")))
        self.assertFalse(h2_database.exists(MhcAllele(gene="H2K", protein="x")))
        self.assertTrue(h2_database.exists(MhcAllele(gene="H2L", protein="f")))
        # badly formed H2 alleles
        self.assertFalse(h2_database.exists(MhcAllele(gene="H2K", group="p", protein=None)))
        self.assertFalse(h2_database.exists(MhcAllele(gene="H2K", group=None, protein=None)))
        self.assertFalse(h2_database.exists(MhcAllele(gene=None, protein="p")))
        self.assertFalse(h2_database.exists(MhcAllele(gene="Z", group="15", protein="228")))
