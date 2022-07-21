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
import pkg_resources
import neofox.tests
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import PredictedEpitope, Patient
from neofox.model.validation import ModelValidator
from neofox.neofox_epitope import NeoFoxEpitope
from neofox.tests.integration_tests.integration_test_tools import get_hla_one_test, get_hla_two_test, \
    BaseIntegrationTest


class TestNeofoxEpitope(BaseIntegrationTest):

    def setUp(self):
        super().setUp()

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
                mutated_peptide="DILVTDQTR",
                wild_type_peptide="DILVIDQTR",
                allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            ),
            PredictedEpitope(
                mutated_peptide="DEVLGEPSQDILVTDQTR",
                wild_type_peptide="DEVLGEPSQDILVIDQTR",
                isoform_mhc_i_i=self._get_test_mhcii_isoform("HLA-DRB1*01:01")
            )
        ]

        annotated_neoepitopes = NeoFoxEpitope(
            neoepitopes=neoepitopes,
            num_cpus=4,
        ).get_annotations()

        self.assertEqual(len(annotated_neoepitopes), 2)    # 16 from the patient neoepitope
        for e, ae in zip(neoepitopes, annotated_neoepitopes):
            if ModelValidator.is_mhci_epitope(e):
                self.assert_neoepitope_mhci(original_neoepitope=e, annotated_neoepitope=ae)
            elif ModelValidator.is_mhcii_epitope(e):
                self.assert_neoepitope_mhcii(original_neoepitope=e, annotated_neoepitope=ae)
            else:
                self.assertTrue(False)

    def test_neofox_epitope_with_patients(self):

        neoepitopes = [
            PredictedEpitope(
                mutated_peptide="DILVTDQTR",
                wild_type_peptide="DILVIDQTR",
                allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            ),
            PredictedEpitope(
                mutated_peptide="DEVLGEPSQDILVTDQTR",
                wild_type_peptide="DEVLGEPSQDILVIDQTR",
                isoform_mhc_i_i=self._get_test_mhcii_isoform("HLA-DRB1*01:01")
            ),
            PredictedEpitope(
                mutated_peptide="DILVTDQTR",
                wild_type_peptide="DILVIDQTR",
                patient_identifier="123"
            ),
        ]

        patients = [
            Patient(
                identifier="123",
                mhc1=get_hla_one_test(self.references.get_mhc_database()),
                mhc2=get_hla_two_test(self.references.get_mhc_database())
            )
        ]

        annotated_neoepitopes = NeoFoxEpitope(
            neoepitopes=neoepitopes,
            patients=patients,
            num_cpus=4,
        ).get_annotations()

        self.assertEqual(len(annotated_neoepitopes), 18)    # 16 from the patient neoepitope

    def _assert_neeoepitope(self, neoepitope: PredictedEpitope):
        # netMHCpan or netMHC2pan annotations
        self.assertIsInstance(neoepitope.rank_mutated, float)
        self.assertIsInstance(neoepitope.rank_wild_type, float)
        self.assertIsInstance(neoepitope.affinity_mutated, float)
        self.assertIsInstance(neoepitope.affinity_wild_type, float)

        # MixMHCpred annotations
        self._assert_float_annotation(neoepitope, annotation_name="MixMHCpred_affinity_score")
        self._assert_float_annotation(neoepitope, annotation_name="MixMHCpred_rank")
        self._assert_float_annotation(neoepitope, annotation_name="MixMHCpred_WT_affinity_score")
        self._assert_float_annotation(neoepitope, annotation_name="MixMHCpred_WT_rank")

        # PRIME annotations
        self._assert_float_annotation(neoepitope, annotation_name="PRIME_affinity_score")
        self._assert_float_annotation(neoepitope, annotation_name="PRIME_rank")
        self._assert_float_annotation(neoepitope, annotation_name="PRIME_WT_affinity_score")
        self._assert_float_annotation(neoepitope, annotation_name="PRIME_WT_rank")

        # additional annotations
        self._assert_annotation(neoepitope, annotation_name="position_mutation")
        self._assert_annotation(neoepitope, annotation_name="anchor_mutated")
        self._assert_annotation(neoepitope, annotation_name="amplitude")
        self._assert_annotation(neoepitope, annotation_name="pathogen_similarity")
        self._assert_annotation(neoepitope, annotation_name="recognition_potential")
        self._assert_annotation(neoepitope, annotation_name="DAI")
        self._assert_annotation(neoepitope, annotation_name="Improved_Binder_MHCI")
        self._assert_annotation(neoepitope, annotation_name="Selfsimilarity")
        self._assert_annotation(neoepitope, annotation_name="Selfsimilarity_conserved_binder")
        self._assert_annotation(neoepitope, annotation_name="mutation_not_found_in_proteome")
        self._assert_annotation(neoepitope, annotation_name="dissimilarity_score")
        self._assert_annotation(neoepitope, annotation_name="number_of_mismatches")
        self._assert_annotation(neoepitope, annotation_name="IEDB_Immunogenicity")
        self._assert_annotation(neoepitope, annotation_name="hex_alignment_score")

        # others to comes
        self._assert_annotation(neoepitope, annotation_name="Priority_score")
        self._assert_annotation(neoepitope, annotation_name="Tcell_predictor_score")