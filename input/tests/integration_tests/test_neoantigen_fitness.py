from collections import defaultdict
from unittest import TestCase

from input.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from input.helpers.runner import Runner
import input.tests.integration_tests.integration_test_tools as integration_test_tools
from input import MHC_II, MHC_I


class TestNeoantigenFitness(TestCase):

    def setUp(self):
        self.references, self.configuration, self.fastafile = self._load_references()
        self.neoantigen_fitness_calculator = NeoantigenFitnessCalculator(
            runner=Runner(), configuration=self.configuration)

    def _load_references(self):
        references, configuration = integration_test_tools.load_references()
        fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        return references, configuration, fastafile

    def test_pathogen_similarity(self):
        result = self.neoantigen_fitness_calculator.wrap_pathogensimilarity(
            mutation='hey',
            fastafile=self.fastafile.name,
            iedb=self.references.iedb)
        self.assertEqual('0', result)

    def test_amplitude_mhc(self):
        self.assertEqual('1.0', self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            score_mutation="1.0", score_wild_type="1.0"))
        self.assertEqual('0.9997000899730081', self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            score_mutation="1.0", score_wild_type="1.0", apply_correction=True))

    def test_recognition_potential(self):
        props = defaultdict(lambda: "1.0")
        props['Mutation_in_anchor_netmhcpan'] = '0'
        props['Mutation_in_anchor_netmhcpan_9mer'] = '0'
        self.assertEqual('1.0', self.neoantigen_fitness_calculator.calculate_recognition_potential(
            amplitude="1.0", pathogen_similarity="1.0", mutation_in_anchor="0"))
