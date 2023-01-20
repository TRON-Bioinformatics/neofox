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
import os
import unittest
from unittest import TestCase

import pkg_resources

from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import Neoantigen, Patient

import neofox
from neofox.exceptions import (
    NeofoxConfigurationException,
    NeofoxDataValidationException,
)
from neofox.neofox import NeoFox
from neofox.tests.fake_classes import FakeReferenceFolder, FakeDependenciesConfiguration, FakeHlaDatabase


class TestNeofox(TestCase):

    def setUp(self) -> None:
        self.hla_database = FakeHlaDatabase()

    def test_missing_input_raises_exception(self):
        with self.assertRaises(NeofoxConfigurationException):
            NeoFox(
                neoantigens=None,
                patients=None,
                num_cpus=1,
                reference_folder=FakeReferenceFolder(),
            )
        with self.assertRaises(NeofoxConfigurationException):
            NeoFox(
                neoantigens=[],
                patients=[],
                num_cpus=1,
                reference_folder=FakeReferenceFolder(),
            )

    def test_not_set_reference_folder_fails(self):
        with self.assertRaises(NeofoxConfigurationException):
            NeoFox(
                neoantigens=[self._get_test_neoantigen()],
                patients=[self._get_test_patient()],
                num_cpus=1,
                reference_folder=FakeReferenceFolder(),
            ).get_annotations()

    def test_empty_reference_folder_fails(self):
        os.environ[neofox.REFERENCE_FOLDER_ENV] = "dummy"
        with self.assertRaises(NeofoxConfigurationException):
            NeoFox(
                neoantigens=[self._get_test_neoantigen()],
                patients=[self._get_test_patient()],
                num_cpus=1,
                reference_folder=FakeReferenceFolder(),
            ).get_annotations()

    def test_validation_captures_bad_wild_type_xmer(self):
        neoantigen = self._get_test_neoantigen()
        neoantigen.wild_type_xmer = "123"  # should be a valid aminoacid
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(
                neoantigens=[neoantigen],
                patients=[self._get_test_patient()],
                num_cpus=1,
                reference_folder=FakeReferenceFolder(),
                configuration=FakeDependenciesConfiguration(),
            )

    def test_validation_captures_bad_mutated_xmer(self):
        neoantigen = self._get_test_neoantigen()
        neoantigen.mutated_xmer = "123"  # should be a valid aminoacid
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(
                neoantigens=[neoantigen],
                patients=[self._get_test_patient()],
                num_cpus=1,
                reference_folder=FakeReferenceFolder(),
                configuration=FakeDependenciesConfiguration(),
            )

    def test_validation_captures_bad_patient(self):
        patient = self._get_test_patient()
        patient.identifier = 12345  # should be a string
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(
                neoantigens=[self._get_test_neoantigen()],
                patients=[patient],
                num_cpus=1,
                reference_folder=FakeReferenceFolder(),
                configuration=FakeDependenciesConfiguration(),
            )

    def test_valid_data_does_not_raise_exceptions(self):
        NeoFox(
            neoantigens=[self._get_test_neoantigen()],
            patients=[self._get_test_patient()],
            num_cpus=1,
            reference_folder=FakeReferenceFolder(),
            configuration=FakeDependenciesConfiguration(),
        )

    def test_neoantigens_referring_to_non_existing_patients(self):
        neoantigen = self._get_test_neoantigen()
        neoantigen.patient_identifier = (
            "I am not patient"  # should be a valid aminoacid
        )
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(
                neoantigens=[neoantigen],
                patients=[self._get_test_patient()],
                num_cpus=1,
                reference_folder=FakeReferenceFolder(),
                configuration=FakeDependenciesConfiguration(),
            )
        neoantigen.patient_identifier = None
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(
                neoantigens=[neoantigen],
                patients=[self._get_test_patient()],
                num_cpus=1,
                reference_folder=FakeReferenceFolder(),
                configuration=FakeDependenciesConfiguration(),
            )
        neoantigen.patient_identifier = ""
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(
                neoantigens=[neoantigen],
                patients=[self._get_test_patient()],
                num_cpus=1,
                reference_folder=FakeReferenceFolder(),
                configuration=FakeDependenciesConfiguration(),
            )

    def test_no_expression_imputation(self):
        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_model_realistic.txt"
        )
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_patient_file.txt"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        neoantigens = ModelConverter.parse_candidate_file(input_file)
        neofox_runner = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            reference_folder=FakeReferenceFolder(),
            configuration=FakeDependenciesConfiguration(),
        )
        for neoantigen in neoantigens:
            for neoantigen_imputed in neofox_runner.neoantigens:
                if neoantigen.mutated_xmer == neoantigen_imputed.mutated_xmer:
                    self.assertEqual(
                        neoantigen.rna_expression, neoantigen_imputed.rna_expression
                    )

    def test_with_expression_imputation(self):
        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_model_realistic_Pty.txt"
        )
        neoantigens= ModelConverter.parse_candidate_file(input_file)
        import copy

        original_neoantigens = copy.deepcopy(neoantigens)
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_patient_file.txt"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        neofox_runner = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            reference_folder=FakeReferenceFolder(),
            configuration=FakeDependenciesConfiguration(),
        )
        for neoantigen, neoantigen_imputed in zip(original_neoantigens, neofox_runner.neoantigens):
            self.assertIsNotNone(neoantigen_imputed.imputed_gene_expression)
            if neoantigen.rna_expression is None:
                self.assertNotEqual(neoantigen.rna_expression, neoantigen_imputed.rna_expression)
            else:
                self.assertEqual(neoantigen.rna_expression, neoantigen_imputed.rna_expression)

    def _get_test_neoantigen(self):
        return Neoantigen(
            gene="GENE",
            mutated_xmer="AAAAAAAIAAAAAAAA",
            wild_type_xmer="AAAAAAALAAAAAAAA",
            patient_identifier="12345",
            rna_expression=0.12345,
        )

    def _get_test_patient(self):
        return Patient(identifier="12345")


if __name__ == "__main__":
    unittest.main()
