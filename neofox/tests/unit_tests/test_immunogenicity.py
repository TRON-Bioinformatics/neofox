from unittest import TestCase

from neofox.published_features.iedb_immunogenicity.iedb import IEDBimmunogenicity


class TestImmunogenicity(TestCase):

    def setUp(self):
        self.immunogenicity_calculator = IEDBimmunogenicity()

    def test_immunogenicity(self):
        result = self.immunogenicity_calculator.calculate_iedb_immunogenicity(epitope="ENPVVHFF",
                                                                              mhc_allele="HLA-A*68:01", mhc_score=600)
        self.assertGreater(result, 0)
        result = self.immunogenicity_calculator.calculate_iedb_immunogenicity(epitope="ENPVVHFF",
                                                                              mhc_allele="HLA-A*68:01", mhc_score=600,
                                                                              affin_filtering=True)
        self.assertIsNone(result)
