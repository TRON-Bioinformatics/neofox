import unittest

from neofox.exceptions import NeofoxDataValidationException
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import MhcAllele
from neofox.model.validation import ModelValidator
from neofox.references.references import ORGANISM_HOMO_SAPIENS
from neofox.tests.fake_classes import FakeHlaDatabase


class TestHlaParser(unittest.TestCase):

    def setUp(self) -> None:
        self.mhc_parser = MhcParser.get_mhc_parser(FakeHlaDatabase())

    def test_parse_mhc_with_3_digits_in_second_place(self):
        mhc = self.mhc_parser.parse_mhc_allele("B15:228")
        self.assertEqual("B", mhc.gene)
        self.assertEqual("15", mhc.group)
        self.assertEqual("228", mhc.protein)

    def test_parse_mhc_non_existing_in_hla_database_does_not_fail(self):
        self.mhc_parser.parse_mhc_allele("B152:28")

    def test_parse_mhc_no_separator_with_3_digits_in_second_place(self):
        mhc = self.mhc_parser.parse_mhc_allele("B15228")
        self.assertEqual("B", mhc.gene)
        self.assertEqual("15", mhc.group)
        self.assertEqual("228", mhc.protein)

    def test_parse_mhc2_with_3_digits_in_first_place(self):
        mhc = self.mhc_parser.parse_mhc_allele("DPB1104:01")
        self.assertEqual("DPB1", mhc.gene)
        self.assertEqual("104", mhc.group)
        self.assertEqual("01", mhc.protein)

    def test_parse_mhc2_non_existing_in_hla_database_does_not_fail(self):
        self.mhc_parser.parse_mhc_allele("DPB110:401")

    def test_parse_mhc2_no_separator_with_3_digits_in_first_place(self):
        mhc = self.mhc_parser.parse_mhc_allele("DPB110401")
        self.assertEqual("DPB1", mhc.gene)
        self.assertEqual("104", mhc.group)
        self.assertEqual("01", mhc.protein)

    def test_mhc_i_allele_parsing(self):
        # adds the star
        self._assert_allele_parsing(
            expected="HLA-A*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-A01:01")
        )
        # adds the HLA-
        self._assert_allele_parsing(
            expected="HLA-A*01:01", allele=self.mhc_parser.parse_mhc_allele("A01:01")
        )
        # adds the colon to homogenise representation
        self._assert_allele_parsing(
            expected="HLA-A*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-A01:01")
        )
        # does not modify an originally good representation
        self._assert_allele_parsing(
            expected="HLA-A*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-A*01:01")
        )
        # removes further information
        self._assert_allele_parsing(
            expected="HLA-A*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-A01:01:02:03N")
        )
        self._assert_allele_parsing(
            expected="HLA-A*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-A01:01:02N")
        )
        self._assert_allele_parsing(
            expected="HLA-A*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-A01:01N")
        )
        self._assert_allele_parsing(
            expected="HLA-A*01:01", allele=MhcAllele(gene="A", group="01", protein="01", name="HLA-A*01:01")
        )

    def _assert_allele_parsing(self, allele, expected, organism=ORGANISM_HOMO_SAPIENS):
        ModelValidator.validate_mhc_allele_representation(allele, organism)
        self.assertEqual(expected, allele.name)

    def test_mhc_ii_allele_parsing(self):
        # add the star
        self._assert_allele_parsing(
            expected="HLA-DPB1*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-DPB101:01")
        )
        # adds the HLA-
        self._assert_allele_parsing(
            expected="HLA-DPB1*01:01", allele=self.mhc_parser.parse_mhc_allele("DPB1*01:01")
        )
        # adds the colon to homogenise representation
        self._assert_allele_parsing(
            expected="HLA-DPA1*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-DPA101:01")
        )
        # does not reove the star
        self._assert_allele_parsing(
            expected="HLA-DPA1*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-DPA1*01:01")
        )
        # removes further information
        self._assert_allele_parsing(
            expected="HLA-DPA1*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-DPA101:01:02:03N")
        )
        self._assert_allele_parsing(
            expected="HLA-DPA1*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-DPA101:01:02N")
        )
        self._assert_allele_parsing(
            expected="HLA-DPB1*01:01", allele=self.mhc_parser.parse_mhc_allele("HLA-DPB101:01")
        )
        self._assert_allele_parsing(
            expected="HLA-DPA1*01:01",
            allele=MhcAllele(gene="DPA1", group="01", protein="01", name="HLA-DPA1*01:01"),
        )

    def test_invalid_mhc_i_alleles(self):
        # P gene is not valid
        self._assert_invalid_allele("HLA-P01:01")
        # serotype 1 is not valid
        self._assert_invalid_allele("HLA-A1:01")
        # no protein information
        self._assert_invalid_allele("HLA-A01")
        # bad protein format
        self._assert_invalid_allele("HLA-A01:ABC")
        # wrong organism, only human
        self._assert_invalid_allele("GOGO-A01:01")
        # no gene
        self._assert_invalid_allele("HLA-0123456")
        # nonsense
        self._assert_invalid_allele("nonsense")

    def _assert_invalid_allele(self, allele):
        self.assertRaises(
            NeofoxDataValidationException,
            self.mhc_parser.parse_mhc_allele,
            allele
        )

    def test_invalid_mhc_ii_alleles(self):
        # P gene is not valid
        self._assert_invalid_allele("HLA-DPR01:01")
        # serotype 1 is not valid
        self._assert_invalid_allele("HLA-DPA11:01")
        # no protein information
        self._assert_invalid_allele("HLA-DPA101")
        # bad protein format
        self._assert_invalid_allele("HLA-DPA101:ABC")
        # wrong organism, only human
        self._assert_invalid_allele("GOGO-DPA101:01")
        # no gene
        self._assert_invalid_allele("HLA-0123456")
        # nonsense
        self._assert_invalid_allele("nonsense")
