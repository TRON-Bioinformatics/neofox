from unittest import TestCase
import pkg_resources
import input.tests
from input.model.neoantigen import Patient
import input.helpers.data_import as data_import


class DataImportTest(TestCase):

    def test_patient_alleles_csv_file2model(self):
        patients_file = pkg_resources.resource_filename(input.tests.__name__, "resources/alleles.Pt29.csv")
        patients = data_import.import_patients_data(patients_file)
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

    def test_patient_csv_file2model(self):
        patients_file = pkg_resources.resource_filename(input.tests.__name__, "resources/patient.Pt29.csv")
        patients = data_import.import_patients_data(patients_file)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Pt29")
        self.assertEqual(len(patients[0].mhc_i_alleles), 6)
        self.assertEqual(len(patients[0].mhc_i_i_alleles), 10)
        self.assertEqual(patients[0].estimated_tumor_content, 69.0)
        self.assertEqual(patients[0].is_rna_available, True)
        self.assertEqual(patients[0].tissue, 'skin')

