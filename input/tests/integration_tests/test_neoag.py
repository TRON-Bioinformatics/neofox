from unittest import TestCase, SkipTest

from input.predictors.neoag.neoag_gbm_model import NeoagCalculator
from input.helpers.runner import Runner
import input.tests.integration_tests.integration_test_tools as integration_test_tools


class TestNeoantigenFitness(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        self.runner = Runner()

    def test_neoag(self):
        result = NeoagCalculator(runner=self.runner, configuration=self.configuration).wrapper_neoag(
            sample_id="12345",
            mut_peptide="DDDDDV",
            score_mut="0",
            ref_peptide="DDDDDD",
            peptide_variant_position="123")
        self.assertTrue(isinstance(result, str))
        self.assertTrue(float(result) > 0)
