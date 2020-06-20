import random
import struct
from unittest import TestCase
import pkg_resources

from Bio.Data import IUPACData
import numpy as np

import input.tests
from input.model.schema_conversion import SchemaConverter
from input.model.neoantigen import Neoantigen, Gene, Mutation, Patient


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
        neoantigens = [_get_random_neoantigen() for _ in range(5)]
        csv_data = SchemaConverter.model2csv(neoantigens)
        self.assertIsNotNone(csv_data)
        self.assertEqual(csv_data.shape[0], len(neoantigens))
        for n in neoantigens:
            self.assertEqual(n.variant_allele_frequency,
                             csv_data[
                                 (csv_data['mutation.position'] == n.mutation.position) &
                                 (csv_data['mutation.mutatedAminoacid'] == n.mutation.mutated_aminoacid)
                             ].variantAlleleFrequency.iloc[0])

    def test_csv2model(self):
        neoantigens = [_get_random_neoantigen() for _ in range(5)]
        csv_data = SchemaConverter.model2csv(neoantigens)
        neoantigens2 = SchemaConverter.neoantigens_csv2model(csv_data)
        self._assert_lists_equal(neoantigens, neoantigens2)

    def test_patient_csv2model(self):
        patients = [_get_random_patient() for _ in range(5)]
        csv_data = SchemaConverter.model2csv(patients)
        patients2 = SchemaConverter.patient_metadata_csv2model(csv_data)
        self._assert_lists_equal(patients, patients2)

    def test_patient_csv_file2model(self):
        patients_file = pkg_resources.resource_filename(input.tests.__name__, "resources/Alleles.Pt29.csv")
        patients = [_get_random_patient() for _ in range(5)]
        csv_data = SchemaConverter.model2csv(patients)
        patients2 = SchemaConverter.patient_metadata_csv2model(csv_data)
        self._assert_lists_equal(patients, patients2)
        
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
        neoantigen.expression_value = "5.7"  # should be a float
        with self.assertRaises(struct.error):
            SchemaConverter.validate(neoantigen)


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