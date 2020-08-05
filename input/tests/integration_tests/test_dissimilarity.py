from unittest import TestCase

import input.tests.integration_tests.integration_test_tools as integration_test_tools
from input.helpers.runner import Runner
from input.predictors.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator


class TestDissimilarity(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        self.runner = Runner()

    def test_dissimilar_sequences(self):
        result = DissimilarityCalculator(
            runner=self.runner, configuration=self.configuration, proteome_db=self.references.proteome_db)\
            .calculate_dissimilarity(
            mhc_mutation='tocino', mhc_affinity='velocidad')
        self.assertEqual(1, result)

    def test_similar_sequences(self):
        result = DissimilarityCalculator(
            runner=self.runner, configuration=self.configuration, proteome_db=self.references.proteome_db)\
            .calculate_dissimilarity(
            mhc_mutation='DDDDDD', mhc_affinity='DDDDDD')
        self.assertTrue(result < 0.000001)

    def test_dissimilarity_mhcii(self):
        # peptide with point mutation
        result = DissimilarityCalculator(
            runner=self.runner, configuration=self.configuration, proteome_db=self.references.proteome_db) \
            .calculate_dissimilarity(
            mhc_mutation='LGLSDSQFLQTFLFM', mhc_affinity='430')
        self.assertEqual(result, 0)
        # unsimmilar peptide
        result = DissimilarityCalculator(
            runner=self.runner, configuration=self.configuration, proteome_db=self.references.proteome_db) \
            .calculate_dissimilarity(
            mhc_mutation='LFTSPIMTKSAEMIV', mhc_affinity='430')
        print(result)
        self.assertGreater(result, 0.0)
