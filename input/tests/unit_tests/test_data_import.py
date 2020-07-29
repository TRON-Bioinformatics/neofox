from unittest import TestCase
import pkg_resources
import input.tests
from input.model.neoantigen import Patient
import input.helpers.data_import as data_import


class DataImportTest(TestCase):

    def test_patients_csv_file2model(self):
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

    def test_patients_csv_file2model2(self):
        patients_file = pkg_resources.resource_filename(input.tests.__name__, "resources/patient.Pt29.csv")
        patients = data_import.import_patients_data(patients_file)
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
        patients_file = pkg_resources.resource_filename(input.tests.__name__, "resources/test_patient_info.txt")
        patients = data_import.import_patients_data(patients_file)
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

    def test_data_import(self):
        icam_file = pkg_resources.resource_filename(input.tests.__name__, "resources/test_data.txt")
        header, rows = data_import.import_dat_icam(icam_file)
        self.assertIsNotNone(header)
        self.assertEqual(len(header), 44)
        self.assertIsNotNone(rows)
        self.assertEqual(len(rows), 9)  # 2 indels are excluded from the input file
        self.assertIsInstance(rows, list)
        for r in rows:
            self.assertEqual(len(r), 44)

