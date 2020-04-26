from unittest import TestCase

import input.tests.integration_tests.integration_test_tools as integration_test_tools
from input.helpers.runner import Runner
from input.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator


class TestDissimilarity(TestCase):

    def setUp(self):
        self.references = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        self.runner = Runner()

    def test_dissimilarity(self):
        result = DissimilarityCalculator(runner=self.runner).calculate_dissimilarity(
            props={'best_affinity_epitope_netmhcpan4': 'hey', 'best_affinity_netmhcpan4': 'ho'},
            fastafile=self.fastafile.name,
            references=self.references)
        self.assertEqual('0', result)
