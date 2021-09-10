import unittest
from neofox.model.mhc_parser import MhcParser
from neofox.tests.fake_classes import FakeH2Database


class TestH2Parser(unittest.TestCase):

    def setUp(self) -> None:
        self.mhc_parser = MhcParser.get_mhc_parser(FakeH2Database())

    def test_parse_h2_allele(self):
        mhc = self.mhc_parser.parse_mhc_allele("H2Kb")
        self.assertEqual("H2K", mhc.gene)
        self.assertEqual("", mhc.group)
        self.assertEqual("b", mhc.protein)

    def test_parse_h2_allele_with_number(self):
        mhc = self.mhc_parser.parse_mhc_allele("H2Kq2")
        self.assertEqual("H2K", mhc.gene)
        self.assertEqual("", mhc.group)
        self.assertEqual("q2", mhc.protein)

    def test_parse_h2_non_existing_in_h2_database_does_not_fail(self):
        self.mhc_parser.parse_mhc_allele("H2Kx")

    def test_parse_h2_mhc2(self):
        mhc = self.mhc_parser.parse_mhc_allele("H2Eb")
        self.assertEqual("H2E", mhc.gene)
        self.assertEqual("", mhc.group)
        self.assertEqual("b", mhc.protein)

    def test_parse_mhc2_non_existing_in_h2_database_does_not_fail(self):
        self.mhc_parser.parse_mhc_allele("H2Ex")
