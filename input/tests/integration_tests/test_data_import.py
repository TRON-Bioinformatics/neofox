from unittest import TestCase

from input.helpers import data_import


class TestDataImport(TestCase):

    def test_data_import(self):
        icam_file = '/projects/SUMMIT/WP1.2/input/development/Pt29.sequences4testing.txt'
        # icam_file = '\\\\192.168.171.199\\projects$\\SUMMIT\\WP1.2\\input\\development\\Pt29.sequences4testing.txt'
        header, rows = data_import.import_dat_icam(icam_file)
        self.assertIsNotNone(header)
        self.assertEqual(len(header), 44)
        self.assertIsNotNone(rows)
        self.assertEqual(len(rows), 9)
        self.assertIsInstance(rows, list)
        for r in rows:
            self.assertEqual(len(r), 44)
