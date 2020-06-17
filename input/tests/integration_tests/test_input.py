from unittest import TestCase

from input.predict_all_epitopes import BunchEpitopes


class TestInput(TestCase):

    def test_input(self):
        """
        This test is equivalent to the command line call:
        input -i /projects/SUMMIT/WP1.2/input/development/Pt29.sequences4testing.txt
        -a /projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/vanallen/output_tables/20200106_alleles_extended.csv
        -tc /projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/vanallen/output_tables/vanallen_patient_overview.csv

        NOTE: we will need to check the output when the calculation of resuls and printing to stdout have been decoupled
        :return:
        """
        input_file = '/projects/SUMMIT/WP1.2/input/development/Pt29.sequences4testing.txt'
        alleles_file = '/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/vanallen/output_tables/20200106_alleles_extended.csv'
        tumor_content_file = '/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/vanallen/output_tables/vanallen_patient_overview.csv'
        BunchEpitopes().wrapper_table_add_feature_annotation(
            file=input_file,
            indel=False,
            path_to_hla_file=alleles_file,
            tissue='skin',
            tumour_content_file=tumor_content_file)
