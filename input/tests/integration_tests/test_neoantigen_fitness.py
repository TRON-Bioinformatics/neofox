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
        props = defaultdict(lambda: "1,0")
        self.assertEqual('1.0', self.neoantigen_fitness_calculator.calculate_amplitude_mhc(props=props, mhc=MHC_I))
        self.assertEqual('1.0', self.neoantigen_fitness_calculator.calculate_amplitude_mhc(props=props, mhc=MHC_II))
        self.assertEqual('1.0',
                         self.neoantigen_fitness_calculator.calculate_amplitude_mhc(props=props, mhc=MHC_I, multiple_binding=True))
        self.assertEqual('1.0',
                         self.neoantigen_fitness_calculator.calculate_amplitude_mhc(props=props, mhc=MHC_II, multiple_binding=True))
        self.assertEqual('0.9997000899730081',
                         self.neoantigen_fitness_calculator.calculate_amplitude_mhc(props=props, mhc=MHC_I, affinity=True))
        self.assertEqual('0.9997000899730081',
                         self.neoantigen_fitness_calculator.calculate_amplitude_mhc(props=props, mhc=MHC_II, affinity=True))
        self.assertEqual('1.0', self.neoantigen_fitness_calculator.calculate_amplitude_mhc(props=props, mhc=MHC_I, netmhcscore=True))
        self.assertEqual('1.0', self.neoantigen_fitness_calculator.calculate_amplitude_mhc(props=props, mhc=MHC_II, netmhcscore=True))
        self.assertEqual('0.9997000899730081',
                         self.neoantigen_fitness_calculator.calculate_amplitude_mhc(props=props, mhc=MHC_I, nine_mer=True))
        self.assertEqual('0.9997000899730081',
                         self.neoantigen_fitness_calculator.calculate_amplitude_mhc(props=props, mhc=MHC_II, nine_mer=True))

    def test_recognition_potential(self):
        props = defaultdict(lambda: "1.0")
        props['Mutation_in_anchor_netmhcpan'] = '0'
        props['Mutation_in_anchor_netmhcpan_9mer'] = '0'
        self.assertEqual('1.0', self.neoantigen_fitness_calculator.calculate_recognition_potential(props=props, mhc=MHC_I))
        self.assertEqual('1.0', self.neoantigen_fitness_calculator.calculate_recognition_potential(props=props, mhc=MHC_II))
        self.assertEqual('1.0',
                         self.neoantigen_fitness_calculator.calculate_recognition_potential(props=props, mhc=MHC_I, affinity=True))
        self.assertEqual('1.0',
                         self.neoantigen_fitness_calculator.calculate_recognition_potential(props=props, mhc=MHC_II, affinity=True))
        self.assertEqual('1.0',
                         self.neoantigen_fitness_calculator.calculate_recognition_potential(props=props, mhc=MHC_I, netmhcscore=True))
        self.assertEqual('1.0',
                         self.neoantigen_fitness_calculator.calculate_recognition_potential(props=props, mhc=MHC_II, netmhcscore=True))
        self.assertEqual('1.0', self.neoantigen_fitness_calculator.calculate_recognition_potential(
            props=props, mhc=MHC_I, nine_mer=True))
        self.assertEqual('1.0', self.neoantigen_fitness_calculator.calculate_recognition_potential(
            props=props, mhc=MHC_II, nine_mer=True))
