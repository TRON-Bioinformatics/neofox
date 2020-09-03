from unittest import TestCase, SkipTest

from neofox.model.neoantigen import Annotation
from neofox.predictors.neoag.neoag_gbm_model import NeoagCalculator
from neofox.helpers.runner import Runner
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.tests.fake_classes import FakeBestAndMultipleBinder


class TestNeoantigenFitness(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        self.runner = Runner()

    def test_neoag(self):
        result = NeoagCalculator(runner=self.runner, configuration=self.configuration).get_annotation(
            sample_id="12345",
            netmhcpan=FakeBestAndMultipleBinder(mutated_epitope="DDDDDV", wild_type_epitope="DDDDDD", affinity=0),
            peptide_variant_position="123")
        self.assertTrue(isinstance(result, Annotation))
        self.assertTrue(float(result.value) > 0)
