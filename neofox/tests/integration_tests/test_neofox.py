from unittest import TestCase
from datetime import datetime
import pkg_resources
import neofox.tests
from neofox.neofox import NeoFox
from neofox.tests.integration_tests import integration_test_tools


class TestNeofox(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        # self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        # self.runner = Runner()

    def test_neofox(self):
        """
        This test is equivalent to the command line call:
        neofox --icam-file /projects/SUMMIT/WP1.2/neofox/development/Pt29.sequences4testing.txt --patient-id Pt29
        --patients-data ../resources/patient.pt29.csv

        NOTE: we will need to check the output when the calculation of resuls and printing to stdout have been decoupled
        :return:
        """
        patient_id = 'Pt29'
        input_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/test_data.txt")
        output_file = pkg_resources.resource_filename(neofox.tests.__name__,
                                                      "resources/output_{:%Y%m%d%H%M%S}.txt".format(datetime.now()))
        patients_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/patient.Pt29.csv")
        annotations, header = NeoFox(
            icam_file=input_file,
            patient_id=patient_id,
            patients_file=patients_file).get_annotations()
        NeoFox.write_to_file_sorted(annotations, header, output_file=output_file)
