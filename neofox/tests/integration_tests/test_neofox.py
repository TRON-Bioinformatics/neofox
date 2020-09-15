#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#
from unittest import TestCase
from datetime import datetime
import pkg_resources
import neofox.tests
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import NeoantigenAnnotations
from neofox.neofox import NeoFox
from neofox.tests.integration_tests import integration_test_tools
import pandas as pd
from logzero import logger
import os
import shutil
import math
import numpy as np


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
        output_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/output_{:%Y%m%d%H%M%S}.txt".format(datetime.now()))
        output_file_tall_skinny = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/output_{:%Y%m%d%H%M%S}.annotations.txt".format(datetime.now()))
        output_file_neoantigens = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/output_{:%Y%m%d%H%M%S}.neoantigens.txt".format(datetime.now()))
        output_json_neoantigens = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/output_{:%Y%m%d%H%M%S}.neoantigens.json".format(datetime.now()))
        output_json_annotations = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/output_{:%Y%m%d%H%M%S}.annotations.json".format(datetime.now()))
        patients_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/patient.Pt29.csv")
        neoantigens = ModelConverter.parse_icam_file(input_file)
        patients = ModelConverter.parse_patients_file(patients_file)
        annotations = NeoFox(
            neoantigens=neoantigens, patient_id=patient_id, patients=patients, num_cpus=2).get_annotations()

        # writes output
        ModelConverter.annotations2short_wide_table(
            neoantigen_annotations=annotations, neoantigens=neoantigens).to_csv(output_file, sep='\t', index=False)
        ModelConverter.annotations2tall_skinny_table(annotations).to_csv(output_file_tall_skinny, sep='\t', index=False)
        ModelConverter.objects2dataframe(neoantigens).to_csv(output_file_neoantigens, sep='\t', index=False)
        ModelConverter.objects2json(annotations, output_json_annotations)
        ModelConverter.objects2json(neoantigens, output_json_neoantigens)

        # regression test
        self._regression_test_on_output_file(new_file=output_file)

    def test_neofox_only_one_neoantigen(self):
        """
        """
        patient_id = 'Pt29'
        input_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/test_data_only_one.txt")
        patients_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/patient.Pt29.csv")
        neoantigens = ModelConverter.parse_icam_file(input_file)
        patients = ModelConverter.parse_patients_file(patients_file)
        annotations = NeoFox(
            neoantigens=neoantigens, patient_id=patient_id, patients=patients, num_cpus=2).get_annotations()
        self.assertEqual(1, len(annotations))
        self.assertIsInstance(annotations[0], NeoantigenAnnotations)
        self.assertTrue(len(annotations[0].annotations) > 10)

    def test_neofox_model_input(self):
        """
        """
        patient_id = 'Pt29'
        input_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/test_data_model.txt")
        patients_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/patient.Pt29.csv")
        neoantigens = ModelConverter.parse_neoantigens_file(input_file)
        patients = ModelConverter.parse_patients_file(patients_file)
        annotations = NeoFox(
            neoantigens=neoantigens, patient_id=patient_id, patients=patients, num_cpus=2).get_annotations()
        self.assertEqual(5, len(annotations))
        self.assertIsInstance(annotations[0], NeoantigenAnnotations)
        self.assertTrue(len(annotations[0].annotations) > 10)

    def _regression_test_on_output_file(self, new_file):
        previous_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/output_previous.txt")
        if os.path.exists(previous_file):
            previous_df = pd.read_csv(previous_file, sep="\t")
            new_df = pd.read_csv(new_file, sep="\t")
            self._check_rows(new_df, previous_df)
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
        ok_values_count = 0
        ko_values_count = 0
        for s1, s2 in zip(previous_df[column_name], new_df[column_name]):
            if self._check_single_value(s1, s2):
                ok_values_count += 1
            else:
                ko_values_count += 1

        if ok_values_count == 0:
            logger.error("There no equal values at all for column {}".format(column_name))

        if ko_values_count > 0:
            logger.error("There are {} different values for column {}".format(ko_values_count, column_name))
            logger.error("Previous version: {}".format(previous_df[column_name].transform(
                lambda x: x[0:20] + "..." if len(x) > 20 else x).get_values()))
            logger.error("New version: {}".format(new_df[column_name].transform(
                lambda x: x[0:20] + "..." if len(x) > 20 else x).get_values()))
            error = True

        return error

    def _check_single_value(self, s1, s2):
        if isinstance(s1, float) or isinstance(s1, np.float):
            # equality of NaN is never true so we force it
            # relative tolerance set to consider equal very close floats
            is_equal = True if np.isnan(s1) and np.isnan(s2) else math.isclose(s1, s2, rel_tol=0.0001)
        else:
            is_equal = s1 == s2
        return is_equal

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

    def _check_rows(self, new_df, previous_df):
        # fails the test if the number of rows differ
        self.assertEqual(previous_df.shape[0], new_df.shape[0], "Mismatching number of rows")
