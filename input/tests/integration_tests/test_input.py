from unittest import TestCase

from input.predict_all_epitopes import BunchEpitopes
from input.tests.integration_tests import integration_test_tools


class TestInput(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        # self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        # self.runner = Runner()

    def test_input(self):
        """
        This test is equivalent to the command line call:
        input -i /projects/SUMMIT/WP1.2/input/development/Pt29.sequences4testing.txt -p Pt29
        -a /projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/vanallen/output_tables/20200106_alleles_extended.csv
        -tc /projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/vanallen/output_tables/vanallen_patient_overview.csv

        NOTE: we will need to check the output when the calculation of resuls and printing to stdout have been decoupled
        :return:
        """
        patient_id = 'Pt29'
        input_file = '/projects/SUMMIT/WP1.2/input/development/Pt29.sequences4testing.txt'
        alleles_file = '/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/vanallen/output_tables/20200106_alleles_extended.csv'
        tumor_content_file = '/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/vanallen/output_tables/vanallen_patient_overview.csv'
        BunchEpitopes().wrapper_table_add_feature_annotation(
            icam_file=input_file,
            patient_id=patient_id,
            indel=False,
            hla_file=alleles_file,
            tissue='skin',
            tumour_content_file=tumor_content_file)
