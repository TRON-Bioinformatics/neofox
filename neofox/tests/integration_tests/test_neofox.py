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
import unittest
from unittest import TestCase
from datetime import datetime
import pkg_resources
from neofox.exceptions import NeofoxConfigurationException

from neofox import NEOFOX_MIXMHCPRED_ENV, NEOFOX_MIXMHC2PRED_ENV

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
import timeit


class TestNeofox(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        # self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        # self.runner = Runner()
        self.patient_id = "Pt29"
        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_candidate_file.txt"
        )
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_patient_file.txt"
        )
        self.patients = ModelConverter.parse_patients_file(patients_file)
        self.neoantigens, external_annotations = ModelConverter.parse_candidate_file(
            input_file
        )

    def test_neofox(self):
        """
        This test is equivalent to the command line call:
        neofox --candidate-file /projects/SUMMIT/WP1.2/neofox/development/Pt29.sequences4testing.txt --patient-id Pt29
        --patients-data ../resources/patient.pt29.csv

        NOTE: we will need to check the output when the calculation of resuls and printing to stdout have been decoupled
        """
        output_file = pkg_resources.resource_filename(
            neofox.tests.__name__,
            "resources/output_{:%Y%m%d%H%M%S}_neoantigen_candidates_annotated.tsv".format(
                datetime.now()
            ),
        )
        output_file_tall_skinny = pkg_resources.resource_filename(
            neofox.tests.__name__,
            "resources/output_{:%Y%m%d%H%M%S}.neoantigen_features.tsv".format(
                datetime.now()
            ),
        )
        output_file_neoantigens = pkg_resources.resource_filename(
            neofox.tests.__name__,
            "resources/output_{:%Y%m%d%H%M%S}.neoantigens.tsv".format(datetime.now()),
        )
        output_json_neoantigens = pkg_resources.resource_filename(
            neofox.tests.__name__,
            "resources/output_{:%Y%m%d%H%M%S}.neoantigen_candidates.json".format(
                datetime.now()
            ),
        )
        output_json_annotations = pkg_resources.resource_filename(
            neofox.tests.__name__,
            "resources/output_{:%Y%m%d%H%M%S}.neoantigen_features.json".format(
                datetime.now()
            ),
        )
        annotations = NeoFox(
            neoantigens=self.neoantigens,
            patient_id=self.patient_id,
            patients=self.patients,
            num_cpus=4,
        ).get_annotations()
        annotation_names = [a.name for n in annotations for a in n.annotations]

        # check it does contain any of the MixMHCpred annotations
        self.assertIn("MixMHC2pred_best_peptide", annotation_names)
        self.assertIn("MixMHC2pred_best_rank", annotation_names)
        self.assertIn("MixMHC2pred_best_allele", annotation_names)
        self.assertIn("MixMHCpred_best_peptide", annotation_names)
        self.assertIn("MixMHCpred_best_score", annotation_names)
        self.assertIn("MixMHCpred_best_rank", annotation_names)
        self.assertIn("MixMHCpred_best_allele", annotation_names)
        # checks it does have some of the NetMHCpan annotations
        self.assertIn("Best_affinity_MHCI_9mer_position_mutation", annotation_names)
        self.assertIn("Best_rank_MHCII_score", annotation_names)

        # writes output
        ModelConverter.annotations2short_wide_table(
            neoantigen_annotations=annotations, neoantigens=self.neoantigens
        ).to_csv(output_file, sep="\t", index=False)
        ModelConverter.annotations2tall_skinny_table(annotations).to_csv(
            output_file_tall_skinny, sep="\t", index=False
        )
        ModelConverter.objects2dataframe(self.neoantigens).to_csv(
            output_file_neoantigens, sep="\t", index=False
        )
        ModelConverter.objects2json(annotations, output_json_annotations)
        ModelConverter.objects2json(self.neoantigens, output_json_neoantigens)

        # regression test
        self._regression_test_on_output_file(new_file=output_file)

    def test_neofox_only_one_neoantigen(self):
        """"""
        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_only_one.txt"
        )
        neoantigens, external_annotations = ModelConverter.parse_candidate_file(
            input_file
        )
        annotations = NeoFox(
            neoantigens=neoantigens,
            patient_id=self.patient_id,
            patients=self.patients,
            num_cpus=4,
        ).get_annotations()
        self.assertEqual(1, len(annotations))
        self.assertIsInstance(annotations[0], NeoantigenAnnotations)
        self.assertTrue(len(annotations[0].annotations) > 10)

    def test_neofox_model_input(self):
        """"""
        neoantigens, patients, patient_id = self._get_test_data()
        annotations = NeoFox(
            neoantigens=neoantigens,
            patient_id=patient_id,
            patients=patients,
            num_cpus=2,
        ).get_annotations()
        self.assertEqual(5, len(annotations))
        self.assertIsInstance(annotations[0], NeoantigenAnnotations)
        self.assertTrue(len(annotations[0].annotations) > 10)

    def test_neofox_without_mixmhcpreds(self):
        """
        This test aims at testing neofox when MixMHCpred and MixMHC2pred are not configured. As these are optional it
        shoudl just run, but without these annotations in the output
        """
        del os.environ[NEOFOX_MIXMHCPRED_ENV]
        del os.environ[NEOFOX_MIXMHC2PRED_ENV]
        annotations = NeoFox(
            neoantigens=self.neoantigens,
            patient_id=self.patient_id,
            patients=self.patients,
            num_cpus=1,
        ).get_annotations()
        annotation_names = [a.name for n in annotations for a in n.annotations]
        # check it does not contain any of the MixMHCpred annotations
        self.assertNotIn("MixMHC2pred_best_peptide", annotation_names)
        self.assertNotIn("MixMHC2pred_best_rank", annotation_names)
        self.assertNotIn("MixMHC2pred_best_allele", annotation_names)
        self.assertNotIn("MixMHCpred_best_peptide", annotation_names)
        self.assertNotIn("MixMHCpred_best_score", annotation_names)
        self.assertNotIn("MixMHCpred_best_rank", annotation_names)
        self.assertNotIn("MixMHCpred_best_allele", annotation_names)
        # checks it does have some of the NetMHCpan annotations
        self.assertIn("Best_affinity_MHCI_9mer_position_mutation", annotation_names)
        self.assertIn("Best_rank_MHCII_score", annotation_names)

    @unittest.skip
    def test_neofox_performance(self):
        def compute_annotations():
            return NeoFox(
                neoantigens=self.neoantigens,
                patient_id=self.patient_id,
                patients=self.patients,
                num_cpus=4,
            ).get_annotations()

        print("Average time: {}".format(timeit.timeit(compute_annotations, number=5)))

    @unittest.skip
    def test_neofox_performance_single_neoantigen(self):

        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_only_one.txt"
        )
        neoantigens, _ = ModelConverter.parse_candidate_file(input_file)

        def compute_annotations():
            return NeoFox(
                neoantigens=neoantigens,
                patient_id=self.patient_id,
                patients=self.patients,
                num_cpus=4,
            ).get_annotations()

        print("Average time: {}".format(timeit.timeit(compute_annotations, number=10)))

    def test_neofox_with_config(self):
        neoantigens, patients, patient_id = self._get_test_data()
        config_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/neofox_config.txt"
        )
        try:
            NeoFox(
                neoantigens=neoantigens,
                patient_id=patient_id,
                patients=patients,
                num_cpus=1,
                configuration_file=config_file,
            )
        except NeofoxConfigurationException as e:
            assert "/neofox/testing/reference_data" in str(e)
            return
        assert False

    def test_neofox_without_mhc2(self):
        """"""
        neoantigens, patients, patient_id = self._get_test_data()
        for p in patients:
            p.mhc2 = None
        annotations = NeoFox(
            neoantigens=neoantigens,
            patient_id=self.patient_id,
            patients=patients,
            num_cpus=1,
        ).get_annotations()
        self.assertEqual(5, len(annotations))
        self.assertIsInstance(annotations[0], NeoantigenAnnotations)
        self.assertTrue(len(annotations[0].annotations) > 10)

    def test_neofox_without_mhc1(self):
        """"""
        neoantigens, patients, patient_id = self._get_test_data()
        for p in patients:
            p.mhc1 = None
        annotations = NeoFox(
            neoantigens=neoantigens,
            patient_id=patient_id,
            patients=patients,
            num_cpus=1,
        ).get_annotations()
        self.assertEqual(5, len(annotations))
        self.assertIsInstance(annotations[0], NeoantigenAnnotations)
        self.assertTrue(len(annotations[0].annotations) > 10)

    def _get_test_data(self):
        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_model_file.txt"
        )
        neoantigens, external_annotations = ModelConverter.parse_neoantigens_file(
            input_file
        )
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_patient_file.txt"
        )
        patients = ModelConverter.parse_patients_file(patients_file)
        patient_id = "Pt29"
        return neoantigens, patients, patient_id

    def _regression_test_on_output_file(self, new_file):
        previous_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/output_previous.txt"
        )
        if os.path.exists(previous_file):
            NeofoxChecker(previous_file, new_file)
        else:
            logger.warning("No previous file to compare output with")
        # copies the new file to the previous file for the next execution only if no values were different
        shutil.copyfile(new_file, previous_file)


class NeofoxChecker:
    def __init__(self, previous_file, new_file):
        previous_df = pd.read_csv(previous_file, sep="\t")
        new_df = pd.read_csv(new_file, sep="\t")
        self._check_rows(new_df, previous_df)
        shared_columns = self._check_columns(new_df, previous_df)
        differing_columns = []
        for c in shared_columns:
            if self._check_values(c, new_df, previous_df):
                differing_columns.append(c)
        if len(differing_columns) > 0:
            differing_columns.sort()
            logger.error(
                "There are {} columns with differing values".format(
                    len(differing_columns)
                )
            )
            logger.error("Differing columns {}".format(differing_columns))
        assert len(differing_columns) == 0, "The regression test contains errors"

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
            logger.error(
                "There are no equal values at all for column {}".format(column_name)
            )

        if ko_values_count > 0:
            logger.error(
                "There are {} different values for column {}".format(
                    ko_values_count, column_name
                )
            )
            logger.error(
                "Previous version: {}".format(
                    previous_df[column_name]
                    .transform(lambda x: x[0:20] + "..." if len(x) > 20 else x)
                    .get_values()
                )
            )
            logger.error(
                "New version: {}".format(
                    new_df[column_name]
                    .transform(lambda x: x[0:20] + "..." if len(x) > 20 else x)
                    .get_values()
                )
            )
            error = True

        return error

    def _check_single_value(self, s1, s2):
        is_float_s1 = self._is_float(s1)
        is_float_s2 = self._is_float(s2)
        if is_float_s1 and is_float_s2:
            # equality of NaN is never true so we force it
            # relative tolerance set to consider equal very close floats
            is_equal = (
                True
                if np.isnan(s1) and np.isnan(s2)
                else math.isclose(s1, s2, rel_tol=0.0001)
            )
        else:
            is_equal = s1 == s2
        return is_equal

    def _is_float(self, s1):
        return isinstance(s1, float) or isinstance(s1, np.float)

    def _check_columns(self, new_df, previous_df):
        shared_columns = set(previous_df.columns).intersection(set(new_df.columns))
        lost_columns = set(previous_df.columns).difference(set(new_df.columns))
        if len(lost_columns) > 0:
            logger.warning(
                "There are {} lost columns: {}".format(len(lost_columns), lost_columns)
            )
        gained_columns = set(new_df.columns).difference(set(previous_df.columns))
        if len(gained_columns) > 0:
            logger.warning(
                "There are {} gained columns: {}".format(
                    len(gained_columns), gained_columns
                )
            )
        # fails the test if there are no shared columns
        assert len(shared_columns) > 0, "No shared columns"
        return shared_columns

    def _check_rows(self, new_df, previous_df):
        # fails the test if the number of rows differ
        assert previous_df.shape[0] == new_df.shape[0], "Mismatching number of rows"
