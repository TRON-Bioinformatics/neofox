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
import pkg_resources

import neofox.tests
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import PredictedEpitope, MhcAllele, Mhc2Isoform
from neofox.neofox_epitope import NeoFoxEpitope
from neofox.references.references import ORGANISM_MUS_MUSCULUS
from neofox.tests.integration_tests import integration_test_tools


class TestNeofoxEpitope(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.references_mouse, self.configuration_mouse = integration_test_tools.load_references(
            organism=ORGANISM_MUS_MUSCULUS)

        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_candidate_file.txt"
        )
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_patient_file.txt"
        )
        patients_file_mouse = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_patient_file_mouse.txt"
        )
        self.hla_database = self.references.get_mhc_database()
        self.h2_database = self.references_mouse.get_mhc_database()
        self.patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        self.patients_mouse = ModelConverter.parse_patients_file(patients_file_mouse, self.h2_database)
        self.neoantigens = ModelConverter.parse_candidate_file(input_file)
        self.neoantigens_mouse = ModelConverter.parse_candidate_file(input_file)

    def test_neofox_epitope(self):

        neoepitopes = [
            PredictedEpitope(
                mutated_peptide="AAAAAAAAAAAAA",
                wild_type_peptide="AAAAAAADAAAAA",
                allele_mhc_i=MhcAllele(name='HLA-A*01:01')
            ),
            PredictedEpitope(
                mutated_peptide="AAAAAAAAAAAAAAAAAA",
                wild_type_peptide="AAAAAAAAADAAAAAAAA",
                isoform_mhc_i_i=Mhc2Isoform(name='DRB1*01:01')
            )
        ]

        annotated_neoepitopes = NeoFoxEpitope(
            neoepitopes=neoepitopes,
            num_cpus=4,
        ).get_annotations()

        self.assertEqual(len(annotated_neoepitopes), 2)
