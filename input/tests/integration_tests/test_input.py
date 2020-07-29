from unittest import TestCase
import pkg_resources
import input.tests
from input.immunogenicity_neoantigen_prediction import ImmunogenicityNeoantigenPredictionToolbox
from input.tests.integration_tests import integration_test_tools


class TestInput(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        # self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        # self.runner = Runner()

    def test_input(self):
        """
        This test is equivalent to the command line call:
        input --icam-file /projects/SUMMIT/WP1.2/input/development/Pt29.sequences4testing.txt --patient-id Pt29
        --patients-data ../resources/patient.pt29.csv

        NOTE: we will need to check the output when the calculation of resuls and printing to stdout have been decoupled
        :return:
        """
        patient_id = 'Pt29'
        input_file = pkg_resources.resource_filename(input.tests.__name__, "resources/test_data.txt")
        patients_file = pkg_resources.resource_filename(input.tests.__name__, "resources/patient.Pt29.csv")
        ImmunogenicityNeoantigenPredictionToolbox(
            icam_file=input_file,
            patient_id=patient_id,
            patients_file=patients_file).get_annotations()
