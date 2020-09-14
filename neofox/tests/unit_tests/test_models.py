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
import struct
from unittest import TestCase
import pkg_resources
import pandas as pd

import neofox.tests
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import Neoantigen, Gene, Mutation, Patient, Annotation, NeoantigenAnnotations
from neofox.model.validation import ModelValidator
from neofox.tests.unit_tests.tools import get_random_neoantigen, get_random_patient


class ModelConverterTest(TestCase):

    def test_model2json(self):
        neoantigens = [get_random_neoantigen() for _ in range(5)]
        json_data = [n.to_json() for n in neoantigens]
        self.assertIsInstance(json_data, list)
        self.assertEqual(5, len(json_data))
        neoantigens2 = [Neoantigen().from_json(j) for j in json_data]
        self._assert_lists_equal(neoantigens, neoantigens2)

    def test_model2dict(self):
        neoantigens = [get_random_neoantigen() for _ in range(5)]
        json_data = [n.to_dict() for n in neoantigens]
        self.assertIsInstance(json_data, list)
        self.assertEqual(5, len(json_data))
        neoantigens2 = [Neoantigen().from_dict(j) for j in json_data]
        self._assert_lists_equal(neoantigens, neoantigens2)

    def test_model2csv(self):
        neoantigen = get_random_neoantigen()
        csv_data = ModelConverter.object2series(neoantigen)
        self.assertIsNotNone(csv_data)
        self.assertIsInstance(csv_data, pd.Series)
        self.assertEqual(neoantigen.dna_variant_allele_frequency, csv_data.dna_variant_allele_frequency)

    def test_model2flat_dict(self):
        neoantigen = get_random_neoantigen()
        flat_dict = ModelConverter.object2flat_dict(neoantigen)
        self.assertIsNotNone(flat_dict)
        self.assertEqual(neoantigen.dna_variant_allele_frequency, flat_dict['dna_variant_allele_frequency'])
        self.assertEqual(neoantigen.gene.transcript_identifier, flat_dict['gene.transcript_identifier'])

    def test_model2csv2model(self):
        neoantigen = get_random_neoantigen()
        csv_data = ModelConverter.object2series(neoantigen)
        neoantigen2 = ModelConverter.neoantigens_csv2object(csv_data)
        self.assertEqual(neoantigen, neoantigen2)

    def test_patient_csv2model(self):
        patients = [get_random_patient() for _ in range(5)]
        csv_data = ModelConverter.objects2dataframe(patients)
        patients2 = ModelConverter.patient_metadata_csv2objects(csv_data)
        self._assert_lists_equal(patients, patients2)

    def test_neoantigen_annotations(self):
        neoantigen = get_random_neoantigen()
        annotations = NeoantigenAnnotations()
        annotations.neoantigen_identifier = ModelValidator.generate_neoantigen_identifier(neoantigen)
        annotations.annotations = [Annotation(name='string_annotation', value='blabla'),
                                      Annotation(name='integer_annotation', value=1),
                                      Annotation(name='float_annotation', value=1.1)]
        annotations_dict = annotations.to_dict()
        self.assertTrue(len(annotations_dict.get('annotations')) == 3)
        self.assertEqual(annotations_dict.get('annotations')[0].get('value'), 'blabla')
        # this does not fail, but it will fail validation
        self.assertEqual(annotations_dict.get('annotations')[1].get('value'), 1)
        self.assertEqual(annotations_dict.get('annotations')[2].get('value'), 1.1)
        self.assertTrue(annotations_dict.get('neoantigenIdentifier'), ModelValidator.generate_neoantigen_identifier(neoantigen))

    def test_icam2model(self):
        icam_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/test_data.txt")
        with open(icam_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens = ModelConverter().parse_icam_file(icam_file)
        self.assertIsNotNone(neoantigens)
        # NOTE: the file contains 2 indels that are filtered out
        self.assertEqual(self.count_lines - 1 - 2, len(neoantigens))
        for n in neoantigens:
            self.assertIsInstance(n, Neoantigen)
            self.assertIsInstance(n.gene, Gene)
            self.assertIsInstance(n.mutation, Mutation)
            self.assertTrue(n.gene.transcript_identifier is not None and len(n.gene.transcript_identifier) > 0)
            self.assertTrue(n.mutation.mutated_aminoacid is not None and len(n.mutation.mutated_aminoacid) == 1)
            self.assertTrue(n.rna_variant_allele_frequency is None or
                            (0 <= n.rna_variant_allele_frequency <= 1))
            self.assertTrue(n.rna_expression is None or n.rna_expression >= 0)
            self.assertTrue(0 <= n.dna_variant_allele_frequency <= 1)

    def test_csv2model(self):
        neoantigens_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/test_data_model.txt")
        neoantigens = ModelConverter.parse_neoantigens_file(neoantigens_file)
        self.assertEqual(5, len(neoantigens))
        for n in neoantigens:
            self.assertTrue(isinstance(n, Neoantigen))
            self.assertTrue(n.gene.gene is not None)
            self.assertTrue(n.gene.transcript_identifier is not None)
            self.assertTrue(n.gene.assembly is not None)

    def test_overriding_patient_id(self):
        icam_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/test_data.txt")
        with open(icam_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens = ModelConverter().parse_icam_file(icam_file, patient_id='patientX')
        for n in neoantigens:
            self.assertEqual(n.patient_identifier, 'patientX')
        neoantigens = ModelConverter().parse_icam_file(icam_file)
        for n in neoantigens:
            self.assertEqual(n.patient_identifier, None)

    def test_patients_csv_file2model(self):
        patients_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/alleles.Pt29.csv")
        patients = ModelConverter.parse_patients_file(patients_file)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Pt29")
        self.assertEqual(len(patients[0].mhc_i_alleles), 6)
        self.assertEqual(len(patients[0].mhc_i_i_alleles), 10)
        self.assertEqual(patients[0].is_rna_available, False)

    def test_patients_csv_file2model2(self):
        patients_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/patient.Pt29.csv")
        patients = ModelConverter.parse_patients_file(patients_file)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Pt29")
        self.assertEqual(len(patients[0].mhc_i_alleles), 6)
        self.assertEqual(len(patients[0].mhc_i_i_alleles), 10)
        self.assertEqual(patients[0].is_rna_available, True)

    def test_patients_csv_file2model3(self):
        patients_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/test_patient_info.txt")
        patients = ModelConverter.parse_patients_file(patients_file)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Ptx")
        self.assertEqual(len(patients[0].mhc_i_alleles), 6)
        self.assertEqual(len(patients[0].mhc_i_i_alleles), 10)
        self.assertEqual(patients[0].is_rna_available, True)

    def test_annotations2short_wide_df(self):
        annotations = [
            NeoantigenAnnotations(
                neoantigen_identifier='12345',annotations=[
                    Annotation(name='this_name', value='this_value'), Annotation(name='that_name', value='that_value'),
                    Annotation(name='diese_name', value='diese_value'), Annotation(name='das_name', value='das_value')]
            ),
            NeoantigenAnnotations(
                neoantigen_identifier='6789', annotations=[
                    Annotation(name='this_name', value='0'), Annotation(name='that_name', value='1'),
                    Annotation(name='diese_name', value='2'), Annotation(name='das_name', value='3')]
            )
        ]
        neoantigens = [
            Neoantigen(
                identifier="12345", mutation=Mutation(position=10, wild_type_aminoacid="A", mutated_aminoacid="C")),
            Neoantigen(
                identifier="6789", mutation=Mutation(position=20, wild_type_aminoacid="G", mutated_aminoacid="Z"))
        ]
        df = ModelConverter.annotations2short_wide_table(neoantigen_annotations=annotations, neoantigens=neoantigens)
        self.assertEqual(df.shape[0], 2)
        self.assertEqual(df.shape[1], 22)

        df_annotations = ModelConverter.annotations2tall_skinny_table(annotations)
        self.assertEqual(df_annotations.shape[0], 8)
        self.assertEqual(df_annotations.shape[1], 3)

    def _assert_lists_equal(self, neoantigens, neoantigens2):
        self.assertEqual(len(neoantigens), len(neoantigens2))
        for n1, n2, in zip(neoantigens, neoantigens2):
            self.assertEqual(n1, n2)
