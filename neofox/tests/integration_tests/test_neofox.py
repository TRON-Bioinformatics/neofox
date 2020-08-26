from unittest import TestCase
from datetime import datetime
import pkg_resources
import neofox.tests
from neofox.neofox import NeoFox
from neofox.tests.integration_tests import integration_test_tools
import pandas as pd
from logzero import logger
import os
import shutil


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
        self._regression_test_on_output_file(new_file=output_file)

    def _regression_test_on_output_file(self, new_file):
        previous_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/output_previous.txt")
        if os.path.exists(previous_file):
            previous_df = pd.read_csv(previous_file, sep="\t")
            new_df = pd.read_csv(new_file, sep="\t")
            shared_columns = self._check_columns(new_df, previous_df)
            has_error = False
            for c in shared_columns:
                has_error |= self._check_values(c, new_df, previous_df)
            self.assertFalse(has_error)
        else:
            logger.warning("No previous file to compare output with")
        # copies the new file to the previous file for the next execution only if no values were different
        shutil.copyfile(new_file, previous_file)

    def _check_values(self, column_name, new_df, previous_df):
        error = False
        value_counts = previous_df[column_name].isin(new_df[column_name]).value_counts(dropna=False)

        try:
            ok_values_count = value_counts[True]
        except KeyError:
            ok_values_count = 0
        if ok_values_count == 0:
            logger.error("There no equal values at all for column {}".format(column_name))

        try:
            ko_values_count = value_counts[False]
        except KeyError:
            ko_values_count = 0
        if ko_values_count > 0:
            logger.error("There are {} different values for column {}".format(ko_values_count, column_name))
            logger.error("Previous version: {}".format(previous_df[column_name].get_values()))
            logger.error("New version: {}".format(new_df[column_name].get_values()))
            error = True

        return error

    def _check_columns(self, new_df, previous_df):
        shared_columns = set(previous_df.columns).intersection(set(new_df.columns))
        lost_columns = set(previous_df.columns).difference(set(new_df.columns))
        if len(lost_columns) > 0:
            logger.warning("There are {} lost columns: {}".format(len(lost_columns), lost_columns))
        gained_columns = set(new_df.columns).difference(set(previous_df.columns))
        if len(gained_columns) > 0:
            logger.warning("There are {} gained columns: {}".format(len(gained_columns), gained_columns))
        # fails the test if there are no shared columns
        self.assertTrue(len(shared_columns) > 0)
        return shared_columns

