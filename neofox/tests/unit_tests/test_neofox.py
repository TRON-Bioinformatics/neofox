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

from neofox.model.neoantigen import Neoantigen, Transcript, Mutation, Patient

import neofox
from neofox.exceptions import NeofoxConfigurationException, NeofoxDataValidationException
from neofox.neofox import NeoFox
from neofox.tests.fake_classes import FakeReferenceFolder, FakeDependenciesConfiguration


class TestNeofox(TestCase):

    def test_missing_input_raises_exception(self):
        with self.assertRaises(NeofoxConfigurationException):
            NeoFox(neoantigens=None, patient_id=None, patients=None, num_cpus=1, reference_folder=FakeReferenceFolder())
        with self.assertRaises(NeofoxConfigurationException):
            NeoFox(neoantigens=[], patient_id=None, patients=[], num_cpus=1, reference_folder=FakeReferenceFolder())

    def test_not_set_reference_folder_fails(self):
        with self.assertRaises(NeofoxConfigurationException):
            NeoFox(neoantigens=[self._get_test_neoantigen()], patient_id=None,
                   patients=[self._get_test_patient()], num_cpus=1, reference_folder=FakeReferenceFolder()).get_annotations()

    def test_empty_reference_folder_fails(self):
        os.environ[neofox.REFERENCE_FOLDER_ENV] = 'dummy'
        with self.assertRaises(NeofoxConfigurationException):
            NeoFox(neoantigens=[self._get_test_neoantigen()], patient_id=None,
                   patients=[self._get_test_patient()], num_cpus=1, reference_folder=FakeReferenceFolder()).get_annotations()

    def test_validation_captures_bad_neoantigen(self):
        neoantigen = self._get_test_neoantigen()
        neoantigen.mutation.wild_type_aminoacid = "XXX"     # should be a valid aminoacid
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(neoantigens=[neoantigen], patient_id=None, patients=[self._get_test_patient()], num_cpus=1,
                   reference_folder=FakeReferenceFolder(), configuration=FakeDependenciesConfiguration())

    def test_validation_captures_bad_patient(self):
        patient = self._get_test_patient()
        patient.identifier = 12345      # should be a string
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(neoantigens=[self._get_test_neoantigen()], patient_id=None, patients=[patient], num_cpus=1,
                   reference_folder=FakeReferenceFolder(), configuration=FakeDependenciesConfiguration())

    def test_valid_data_does_not_raise_exceptions(self):
        NeoFox(neoantigens=[self._get_test_neoantigen()], patient_id=None, patients=[self._get_test_patient()],
               num_cpus=1, reference_folder=FakeReferenceFolder(), configuration=FakeDependenciesConfiguration())

    def test_neoantigens_referring_to_non_existing_patients(self):
        neoantigen = self._get_test_neoantigen()
        neoantigen.patient_identifier = "I am not patient"     # should be a valid aminoacid
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(neoantigens=[neoantigen], patient_id=None, patients=[self._get_test_patient()], num_cpus=1,
                   reference_folder=FakeReferenceFolder(), configuration=FakeDependenciesConfiguration())
        neoantigen.patient_identifier = None
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(neoantigens=[neoantigen], patient_id=None, patients=[self._get_test_patient()], num_cpus=1,
                   reference_folder=FakeReferenceFolder(), configuration=FakeDependenciesConfiguration())
        neoantigen.patient_identifier = ""
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(neoantigens=[neoantigen], patient_id=None, patients=[self._get_test_patient()], num_cpus=1,
                   reference_folder=FakeReferenceFolder(), configuration=FakeDependenciesConfiguration())

    def test_repeated_neoantigens(self):
        neoantigen = self._get_test_neoantigen()
        with self.assertRaises(NeofoxDataValidationException):
            NeoFox(neoantigens=[neoantigen, neoantigen], patient_id=None, patients=[self._get_test_patient()],
                   num_cpus=1, reference_folder=FakeReferenceFolder(), configuration=FakeDependenciesConfiguration())

    def _get_test_neoantigen(self):
        return Neoantigen(
            transcript=Transcript(assembly="hg19", identifier="ENST12345", gene="GENE"),
            mutation=Mutation(position=123, mutated_aminoacid="I", wild_type_aminoacid="L",
                              left_flanking_region="AAAAAAA", right_flanking_region="AAAAAAAA"),
            patient_identifier="12345",
            rna_expression=0.12345
        )

    def _get_test_patient(self):
        return Patient(
            identifier="12345",
            is_rna_available=True
        )


if __name__ == "__main__":
    unittest.main()
