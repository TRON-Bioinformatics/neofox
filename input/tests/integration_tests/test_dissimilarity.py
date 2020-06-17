from unittest import TestCase

import input.tests.integration_tests.integration_test_tools as integration_test_tools
from input.helpers.runner import Runner
from input.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator


class TestDissimilarity(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        self.runner = Runner()

    def test_dissimilar_sequences(self):
        result = DissimilarityCalculator(
            runner=self.runner, configuration=self.configuration)\
            .calculate_dissimilarity(
            mhc_mutation='tocino', mhc_affinity='velocidad', references=self.references)
        self.assertEqual(1, result)

    def test_similar_sequences(self):
        result = DissimilarityCalculator(
            runner=self.runner, configuration=self.configuration)\
            .calculate_dissimilarity(
            mhc_mutation='DDDDDD', mhc_affinity='DDDDDD', references=self.references)
        self.assertTrue(result < 0.000001)
