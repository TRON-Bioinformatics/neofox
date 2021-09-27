import unittest

from neofox.exceptions import NeofoxDataValidationException
from neofox.model.mhc_parser import MhcParser
from neofox.model.validation import ModelValidator
from neofox.references.references import ORGANISM_MUS_MUSCULUS
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

    def test_mhc_i_allele_parsing(self):
        # adds the star
        self._assert_allele_parsing(
            expected="H2Kd", allele=self.mhc_parser.parse_mhc_allele("H2Kd")
        )
        self._assert_allele_parsing(
            expected="H2Kd1", allele=self.mhc_parser.parse_mhc_allele("H2Kd1")
        )

    def _assert_allele_parsing(self, allele, expected, organism=ORGANISM_MUS_MUSCULUS):
        ModelValidator.validate_mhc_allele_representation(allele, organism)
        self.assertEqual(expected, allele.name)

    def test_mhc_ii_allele_parsing(self):
        self._assert_allele_parsing(
            expected="H2Aa", allele=self.mhc_parser.parse_mhc_allele("H2Aa")
        )
        self._assert_allele_parsing(
            expected="H2Aa1", allele=self.mhc_parser.parse_mhc_allele("H2Aa1")
        )

    def test_invalid_mhcalleles(self):
        # Z gene is not valid
        self._assert_invalid_allele("H2Zd")
        # protein information not valid
        self._assert_invalid_allele("H2D1")
        # no H2 prefix
        self._assert_invalid_allele("Dd")
        # wrong organism, only mice
        self._assert_invalid_allele("GOGO-A01:01")
        # no gene
        self._assert_invalid_allele("H2d")
        # nonsense
        self._assert_invalid_allele("nonsense")

    def _assert_invalid_allele(self, allele):
        self.assertRaises(
            NeofoxDataValidationException,
            self.mhc_parser.parse_mhc_allele,
            allele
        )
