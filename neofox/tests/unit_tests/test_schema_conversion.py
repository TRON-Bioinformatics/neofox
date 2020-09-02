import random
import struct
from unittest import TestCase
import pkg_resources
import pandas as pd

from Bio.Data import IUPACData
import numpy as np

import neofox.tests
from neofox.model.schema_conversion import SchemaConverter
from neofox.model.neoantigen import Neoantigen, Gene, Mutation, Patient, Annotation


class SchemaConverterTest(TestCase):

    def test_model2json(self):
        neoantigens = [_get_random_neoantigen() for _ in range(5)]
        json_data = [n.to_json() for n in neoantigens]
        self.assertIsInstance(json_data, list)
        self.assertEqual(5, len(json_data))
        neoantigens2 = [Neoantigen().from_json(j) for j in json_data]
        self._assert_lists_equal(neoantigens, neoantigens2)

    def test_model2dict(self):
        neoantigens = [_get_random_neoantigen() for _ in range(5)]
        json_data = [n.to_dict() for n in neoantigens]
        self.assertIsInstance(json_data, list)
        self.assertEqual(5, len(json_data))
        neoantigens2 = [Neoantigen().from_dict(j) for j in json_data]
        self._assert_lists_equal(neoantigens, neoantigens2)

    def test_model2csv(self):
        neoantigen = _get_random_neoantigen()
        csv_data = SchemaConverter.object2series(neoantigen)
        self.assertIsNotNone(csv_data)
        self.assertIsInstance(csv_data, pd.Series)
        self.assertEqual(neoantigen.dna_variant_allele_frequency, csv_data.dna_variant_allele_frequency)

    def test_model2flat_dict(self):
        neoantigen = _get_random_neoantigen()
        flat_dict = SchemaConverter.object2flat_dict(neoantigen)
        self.assertIsNotNone(flat_dict)
        self.assertEqual(neoantigen.dna_variant_allele_frequency, flat_dict['dna_variant_allele_frequency'])
        self.assertEqual(neoantigen.gene.transcript_identifier, flat_dict['gene.transcript_identifier'])

    def test_csv2model(self):
        neoantigen = _get_random_neoantigen()
        csv_data = SchemaConverter.object2series(neoantigen)
        neoantigen2 = SchemaConverter.neoantigens_csv2object(csv_data)
        self.assertEqual(neoantigen, neoantigen2)

    def test_patient_csv2model(self):
        patients = [_get_random_patient() for _ in range(5)]
        csv_data = SchemaConverter.objects2dataframe(patients)
        patients2 = SchemaConverter.patient_metadata_csv2objects(csv_data)
        self._assert_lists_equal(patients, patients2)

    def test_annotate_neoantigen(self):
        neoantigen = _get_random_neoantigen()
        neoantigen.annotations = [Annotation(name='string_annotation', value='blabla'),
                                      Annotation(name='integer_annotation', value=1),
                                      Annotation(name='float_annotation', value=1.1)]
        neoantigen_dict = neoantigen.to_dict()
        self.assertTrue(len(neoantigen_dict.get('annotations')) == 3)
        self.assertEqual(neoantigen_dict.get('annotations')[0].get('value'), 'blabla')
        self.assertEqual(neoantigen_dict.get('annotations')[1].get('value'), 1)
        self.assertEqual(neoantigen_dict.get('annotations')[2].get('value'), 1.1)

    def test_icam2model(self):
        icam_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/test_data.txt")
        with open(icam_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens = SchemaConverter().parse_icam_file(icam_file)
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

    def test_overriding_patient_id(self):
        icam_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/test_data.txt")
        with open(icam_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens = SchemaConverter().parse_icam_file(icam_file, patient_id='patientX')
        for n in neoantigens:
            self.assertEqual(n.patient_identifier, 'patientX')
        neoantigens = SchemaConverter().parse_icam_file(icam_file)
        for n in neoantigens:
            self.assertEqual(n.patient_identifier, None)

    def test_patients_csv_file2model(self):
        patients_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/alleles.Pt29.csv")
        patients = SchemaConverter.parse_patients_file(patients_file)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Pt29")
        self.assertEqual(len(patients[0].mhc_i_alleles), 6)
        self.assertEqual(len(patients[0].mhc_i_i_alleles), 10)
        self.assertEqual(patients[0].estimated_tumor_content, 0.0)
        self.assertEqual(patients[0].is_rna_available, False)
        self.assertEqual(patients[0].tissue, '')

    def test_patients_csv_file2model2(self):
        patients_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/patient.Pt29.csv")
        patients = SchemaConverter.parse_patients_file(patients_file)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Pt29")
        self.assertEqual(len(patients[0].mhc_i_alleles), 6)
        self.assertEqual(len(patients[0].mhc_i_i_alleles), 10)
        self.assertEqual(patients[0].estimated_tumor_content, 0.69)
        self.assertEqual(patients[0].is_rna_available, True)
        self.assertEqual(patients[0].tissue, 'skin')

    def test_patients_csv_file2model3(self):
        patients_file = pkg_resources.resource_filename(neofox.tests.__name__, "resources/test_patient_info.txt")
        patients = SchemaConverter.parse_patients_file(patients_file)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Ptx")
        self.assertEqual(len(patients[0].mhc_i_alleles), 6)
        self.assertEqual(len(patients[0].mhc_i_i_alleles), 10)
        self.assertEqual(patients[0].estimated_tumor_content, 0.84)
        self.assertEqual(patients[0].is_rna_available, True)
        self.assertEqual(patients[0].tissue, 'skin')

    def _assert_lists_equal(self, neoantigens, neoantigens2):
        self.assertEqual(len(neoantigens), len(neoantigens2))
        for n1, n2, in zip(neoantigens, neoantigens2):
            self.assertEqual(n1, n2)


class SchemaValidationTest(TestCase):
    
    def test_validation(self):
        neoantigens = [_get_random_neoantigen() for _ in range(5)]
        for n in neoantigens:
            SchemaConverter.validate(n)

    def test_field_invalid_type(self):
        neoantigen = _get_random_neoantigen()
        neoantigen.rna_expression = "5.7"  # should be a float
        with self.assertRaises(struct.error):
            SchemaConverter.validate(neoantigen)

    def test_annnotation_invalid_type(self):
        annotation = Annotation(name='invalid_annotation', value=123)
        with self.assertRaises(Exception):  # NOTE: when  the offending value is a literal exception is not captured
            SchemaConverter.validate(annotation)


def _get_random_neoantigen():
    neoantigen = Neoantigen()
    neoantigen.variant_allele_frequency = np.random.uniform(0, 1)
    neoantigen.expression_value = np.random.uniform(0, 50)
    mutation = Mutation()
    mutation.mutated_aminoacid = random.choices(list(IUPACData.protein_letters), k=1)[0]
    mutation.wild_type_aminoacid = random.choices(list(IUPACData.protein_letters), k=1)[0]
    mutation.position = np.random.randint(0, 1000)
    neoantigen.mutation = mutation
    gene = Gene()
    gene.gene = "BRCA2"
    gene.transcript_identifier = "ENST1234567"
    gene.assembly = "hg19"
    neoantigen.gene = gene
    return neoantigen


def _get_random_patient():
    patient = Patient()
    patient.estimated_tumor_content = np.random.uniform(0, 1)
    patient.is_rna_available = np.random.choice([True, False], 1)[0]
    patient.identifier = 'Pt12345'
    patient.mhc_i_alleles = ['A', 'B', 'C']
    patient.mhc_i_i_alleles = ['X', 'Y']
    patient.tissue = 'skin'
    return patient