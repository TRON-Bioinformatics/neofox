from unittest import TestCase
from collections import defaultdict
import input.neoantigen_fitness.neoantigen_fitness as neoantigen_fitness
import input.tests.integration_tests.integration_test_tools as integration_test_tools
from input import MHC_II, MHC_I


class TestNeoantigenFitness(TestCase):

    def setUp(self):
        self.references = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()

    def test_pathogen_similarity_mhcI(self):
        result = neoantigen_fitness.wrap_pathogensimilarity(
            props={"MHC_I_epitope_.best_prediction.": "hey"},
            mhc=MHC_I,
            fastafile=self.fastafile.name,
            iedb=self.references.iedb)
        self.assertEqual('0', result)

    def test_pathogen_similarity_mhcII(self):
        result = neoantigen_fitness.wrap_pathogensimilarity(
            props={"MHC_II_epitope_.best_prediction.": "hey"},
            mhc=MHC_II,
            fastafile=self.fastafile.name,
            iedb=self.references.iedb)
        self.assertEqual('0', result)

    def test_amplitude_mhc(self):
        props = defaultdict(lambda: "1,0")
        self.assertEqual('1.0', neoantigen_fitness.calculate_amplitude_mhc(props=props, mhc=MHC_I))
        self.assertEqual('1.0', neoantigen_fitness.calculate_amplitude_mhc(props=props, mhc=MHC_II))
        self.assertEqual('1.0', neoantigen_fitness.calculate_amplitude_mhc(props=props, mhc=MHC_I, multiple_binding=True))
        self.assertEqual('1.0', neoantigen_fitness.calculate_amplitude_mhc(props=props, mhc=MHC_II, multiple_binding=True))
        self.assertEqual('0.999700089973', neoantigen_fitness.calculate_amplitude_mhc(props=props, mhc=MHC_I, affinity=True))
        self.assertEqual('0.999700089973', neoantigen_fitness.calculate_amplitude_mhc(props=props, mhc=MHC_II, affinity=True))
        self.assertEqual('1.0', neoantigen_fitness.calculate_amplitude_mhc(props=props, mhc=MHC_I, netmhcscore=True))
        self.assertEqual('1.0', neoantigen_fitness.calculate_amplitude_mhc(props=props, mhc=MHC_II, netmhcscore=True))
        self.assertEqual('0.999700089973', neoantigen_fitness.calculate_amplitude_mhc(props=props, mhc=MHC_I, nine_mer=True))
        self.assertEqual('0.999700089973', neoantigen_fitness.calculate_amplitude_mhc(props=props, mhc=MHC_II, nine_mer=True))

    def test_recognition_potential(self):
        props = defaultdict(lambda: "1.0")
        props['Mutation_in_anchor'] = '0'
        self.assertEqual('1.0', neoantigen_fitness.calculate_recognition_potential(props=props, mhc=MHC_I))
        self.assertEqual('1.0', neoantigen_fitness.calculate_recognition_potential(props=props, mhc=MHC_II))
        self.assertEqual('1.0', neoantigen_fitness.calculate_recognition_potential(props=props, mhc=MHC_I, affinity=True))
        self.assertEqual('1.0', neoantigen_fitness.calculate_recognition_potential(props=props, mhc=MHC_II, affinity=True))
        self.assertEqual('1.0', neoantigen_fitness.calculate_recognition_potential(props=props, mhc=MHC_I, netmhcscore=True))
        self.assertEqual('1.0', neoantigen_fitness.calculate_recognition_potential(props=props, mhc=MHC_II, netmhcscore=True))
        self.assertEqual('1.0', neoantigen_fitness.calculate_recognition_potential(props=props, mhc=MHC_I, nine_mer=True))
        self.assertEqual('1.0', neoantigen_fitness.calculate_recognition_potential(props=props, mhc=MHC_II, nine_mer=True))


