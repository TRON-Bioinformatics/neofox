from unittest import TestCase

from input.IEDB_Immunogenicity.predict_immunogenicity_simple import IEDBimmunogenicity


class TestImmunogenicity:

    def setUp(self):
        self.immunogenicity_calculator = IEDBimmunogenicity()

    def test_immunogenicity(self):
        result = self.immunogenicity_calculator.calc_IEDB_immunogenicity(epitope="NLVPMVATV",
                                                                         mhc_allele="HLA-A*68:01", mhc_score=600)
        self.assertGreater(result, 0)
        result = self.immunogenicity_calculator.calc_IEDB_immunogenicity(epitope="NLVPMVATV",
                                                                         mhc_allele="HLA-A*68:01", mhc_score=600,
                                                                         affin_filtering=True)
        self.assertGreater(result, "NA")
